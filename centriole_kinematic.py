"""Kinematic centriole model: rigid bodies joined by strength-weighted bonds.

All lengths are in **nanometres**, all angles in degrees.

The centriole is modelled as a **mechanical network**: rigid bodies (each
with a free pose) joined at connection points by springs whose stiffness
is the bond's strength. Perturbing the system therefore has a well-posed
answer -- bodies re-arrange, and the model reports how far each joint had
to rotate, what had to buckle, which bonds are closest to rupture, and
what collides with what.

Bodies and bonds of one repeating unit (dimensions measured off the
cryo-ET-derived schematic by `svg_calibration.py`):

    SAS-6 head (on the hub ring) -> coiled-coil spoke
        |                                   |
        |                             [pinhead-spoke bond]
        |                                   v
        |                                pinhead ---[base-pinhead bond]--> triplet base
        |                                   |                                    |
        |                          [pinhead-triplet bond]              [base-linker bond]
        |                                   v                                    v
        `--------------------------> MT triplet (A+B+C) ---[linker-C bond]--> A-C linker
                                            ^                                    |
                                            `----[linker-A bond, next unit]------'

**The triplet ring is the primary scaffold.** Triplets are held at even
angular spacing by strong springs and are joined head-to-tail by the A-C
linkers; the cartwheel attaches to that ring wherever it can reach. So
when the two symmetries disagree the triplets keep their wild-type
spacing and the *cartwheel* absorbs the mismatch -- spokes strain to
reach, and any surplus triplets are simply left unattached. Spokes and
triplets are paired one-to-one by optimal angular assignment, so with
`N_cw < N_mt` the excess triplets carry no pinhead at all.

**Mechanics.** Segments never stretch; they may *buckle* (bow, shortening
the end-to-end chord at fixed contour length), bounded to shorten only.
Joint rotations are penalised in **per-joint bands** (see
:data:`JOINT_BANDS`): within the OK limit rotation is nearly free, between
OK and HARD it is costly, beyond HARD it is heavily penalised. Contacts on
microtubules are the tightest, because they grip a rigid ordered lattice;
the triplet axis and base are the most permissive. Band values are
reasoned heuristics rather than measurements, so
:func:`band_sensitivity` is provided to check whether any conclusion
actually depends on them.

**Bonds and rupture.** Every connection is a spring whose stiffness is its
strength (see :data:`BOND_STRENGTH`), so weak bonds yield first and the
load distribution is an output rather than an assumption. Bond force is
reported per connection, along with which bond is closest to rupture.

**Sterics.** Tubule-tubule overlap is penalised, and the clearance of the
A-C linker, triplet base and spoke against every microtubule they do not
bind is measured and reported.

**Protofilament register shift** (off by default) lets the linker and
pinhead contacts slide to neighbouring protofilaments, testing whether a
mutant could relieve strain by re-registering rather than deforming.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Optional

import numpy as np
from scipy.optimize import least_squares, linear_sum_assignment

# --------------------------------------------------------------------------
COLOR_SPOKE = "#90b2fc"
COLOR_PINHEAD = "#9e42f1"
COLOR_BASE = "#688950"
COLOR_LINKER = "#9eeadf"
COLOR_MT = "#ababab"
COLOR_HUB = "#5b9bd5"

TUBULES = ("A", "B", "C")

#: Relative bond strengths, strongest to weakest. Doubles as spring
#: stiffness: a stronger bond both resists opening and ruptures later, so
#: load naturally concentrates in the stiff bonds while the weak ones give.
BOND_STRENGTH = {
    "linker-A": 1.00,     # A-C linker to A-tubule  (strongest)
    "linker-C": 0.85,     # A-C linker to C-tubule
    "pinhead-triplet": 0.70,
    "pinhead-spoke": 0.55,
    "base-pinhead": 0.40,
    "base-linker": 0.25,  # (weakest)
}

#: Per-joint rotation bands, as (OK-limit, HARD-limit) in degrees. Beyond
#: the HARD limit a rotation is treated as forbidden.
#:
#: The *ordering* here is defensible; the absolute values are reasoned
#: heuristics, not measurements (no sub-tomogram angular variance is
#: available to calibrate against). Use :func:`band_sensitivity` to check
#: whether a conclusion actually depends on them.
#:
#: Rationale:
#:  - contacts on microtubules are tightest: they grip a rigid, ordered
#:    lattice at defined protofilaments;
#:  - the SAS-6 spoke is NOT treated as unusually soft -- in 3D many SAS-6
#:    rings stack, which stiffens what looks floppy in a single 2D slice;
#:  - the triplet axis and base are the most permissive: neither is a
#:    single interface, and both are implicated in iris-like motion.
JOINT_BANDS = {
    "linker-A": (8.0, 20.0),
    "linker-C": (8.0, 20.0),
    "pinhead-A": (10.0, 22.0),
    "spoke": (15.0, 35.0),
    "pinhead": (15.0, 30.0),
    "triplet": (20.0, 40.0),
    "base": (20.0, 40.0),
}

#: Rest orientation of each strand within the body frame of the tubule it
#: grips (degrees), taken from the wild-type solution of the measured
#: geometry. Deviation from these is what the contact bands penalise.
CONTACT_REST = {"linker-A": -13.199, "linker-C": -53.199, "pinhead-A": 77.394}

W_FREE, W_HARD, W_FORBID = 1.0, 5.0, 30.0

PF1_ANCHOR = {"A": -43.131, "B": 123.036, "C": 110.525}


def _u(deg):
    t = np.radians(deg)
    return np.array([np.cos(t), np.sin(t)])


def _R(deg):
    t = np.radians(deg)
    return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])


def _wrap(a):
    return (a + 180.0) % 360.0 - 180.0


def angle_penalty(dev_deg, joint="spoke", bands=None):
    """Piecewise-linear banded penalty on a joint deviation (degrees).

    Continuous, so the optimiser stays well behaved, but its slope steps
    up at that joint's OK limit and again at its HARD limit.
    """
    free_lim, hard_lim = (bands or JOINT_BANDS)[joint]
    a = np.abs(dev_deg)
    free = np.minimum(a, free_lim)
    hard = np.clip(a - free_lim, 0.0, hard_lim - free_lim)
    forbid = np.maximum(a - hard_lim, 0.0)
    return W_FREE * free + W_HARD * hard + W_FORBID * forbid


def grade(dev_deg, joint, bands=None):
    """Grade the worst deviation of a joint as OK / HARD / FORBIDDEN."""
    if len(np.atleast_1d(dev_deg)) == 0:
        return "-"
    free_lim, hard_lim = (bands or JOINT_BANDS)[joint]
    m = float(np.max(np.abs(dev_deg)))
    return "OK" if m <= free_lim else ("HARD" if m <= hard_lim else "FORBIDDEN")


def _seg_point_dist(a, b, p):
    """Distance from point(s) p to segment ab."""
    ab = b - a
    L2 = float(ab @ ab)
    if L2 < 1e-12:
        return np.linalg.norm(p - a, axis=-1)
    t = np.clip(((p - a) @ ab) / L2, 0.0, 1.0)
    proj = a + np.outer(t, ab) if p.ndim > 1 else a + t * ab
    return np.linalg.norm(p - proj, axis=-1)


# --------------------------------------------------------------------------
@dataclass
class Geometry:
    """Rigid-body dimensions of one unit, in nanometres / degrees."""

    N_cw: int = 9                  # cartwheel: SAS-6 dimers / spokes
    N_mt: int = 9                  # microtubule triplets
    MTn: int = 3                   # 3 = triplet, 2 = doublet, 1 = singlet

    head_contact: float = 10.40    # head-head spacing on the hub ring
    head_length: float = 4.98
    spoke_rod: float = 45.03

    pinhead_span: float = 21.30
    pinhead_base_frac: tuple = (0.340, -0.411)

    n_pf: dict = field(default_factory=lambda: {"A": 13, "B": 9, "C": 9})
    pf_width: float = 5.747
    ab_spacing: float = 18.75
    bc_spacing: float = 19.63

    base_length: float = 34.68
    base_thickness: float = 2.23

    linker_arm_C: float = 14.64
    linker_arm_A: float = 11.64
    linker_vertex_deg: float = 138.7

    rest_spoke: float = 0.0
    rest_pinhead: float = 21.5
    rest_triplet: float = -56.9
    rest_base: float = -54.3

    pin_pf: tuple = (3, 4)
    linkA_pf: tuple = (8,)
    linkC_pf: tuple = (8, 9)

    # --- derived ----------------------------------------------------------
    @property
    def tubule_radius(self) -> float:
        return self.n_pf["A"] * self.pf_width / (2 * np.pi)

    @property
    def _radius_scale(self) -> float:
        return self.tubule_radius / (13 * 5.747 / (2 * np.pi))

    @property
    def spacing_ab(self) -> float:
        return self.ab_spacing * self._radius_scale

    @property
    def spacing_bc(self) -> float:
        return self.bc_spacing * self._radius_scale

    @property
    def pf_pitch(self) -> float:
        return 360.0 / self.n_pf["A"]

    @property
    def outer_tubule(self) -> str:
        return TUBULES[self.MTn - 1]

    @property
    def hub_radius(self) -> float:
        return self.head_contact / (2.0 * np.sin(np.pi / self.N_cw))

    def pf_angle(self, tubule: str, idx, shift: float = 0.0) -> float:
        """Mean protofilament angle relative to the triplet axis.

        `shift` moves the contact by that many protofilament steps, used
        by the register-shift option.
        """
        a = [
            np.radians(PF1_ANCHOR[tubule] - self.pf_pitch * (i - 1 + shift))
            for i in np.atleast_1d(idx)
        ]
        return float(np.degrees(np.arctan2(np.mean(np.sin(a)), np.mean(np.cos(a)))))

    def pf_angles(self, tubule: str) -> dict:
        return {i: _wrap(PF1_ANCHOR[tubule] - self.pf_pitch * (i - 1))
                for i in range(1, self.n_pf[tubule] + 1)}

    def pairing(self):
        """One-to-one spoke<->triplet assignment minimising angular offset.

        Leaves surplus triplets (or spokes) unattached rather than forcing
        one spoke to serve two, which is what lets the triplet ring keep
        its spacing when the symmetries disagree.
        """
        t_mt = np.arange(self.N_mt) * 360.0 / self.N_mt
        t_cw = np.arange(self.N_cw) * 360.0 / self.N_cw
        cost = np.abs(_wrap(t_cw[:, None] - t_mt[None, :]))
        rows, cols = linear_sum_assignment(cost)
        return list(zip(rows.tolist(), cols.tolist()))   # (spoke, triplet)


# --------------------------------------------------------------------------
# Body-local geometry
# --------------------------------------------------------------------------
def _triplet_local(g: Geometry, reg=(0.0, 0.0)):
    """Points on a triplet, in its own frame (origin A centre, +x = axis)."""
    Rt = g.tubule_radius
    pin_shift, linkC_shift = reg
    pts = {
        "A": np.zeros(2),
        "B": np.array([g.spacing_ab, 0.0]),
        "C": np.array([g.spacing_ab + g.spacing_bc, 0.0]),
    }
    ot = g.outer_tubule
    return dict(
        centres=pts,
        A_pf34=Rt * _u(g.pf_angle("A", g.pin_pf, pin_shift)),
        A_pf8=Rt * _u(g.pf_angle("A", g.linkA_pf)),
        C_pf89=pts[ot] + Rt * _u(g.pf_angle(ot, g.linkC_pf, linkC_shift)),
    )


def _place(local, pose):
    """Map body-local points by pose = (x, y, theta_deg)."""
    R = _R(pose[2])
    o = pose[:2]
    return {k: (o + R @ v if isinstance(v, np.ndarray) else
                {kk: o + R @ vv for kk, vv in v.items()})
            for k, v in local.items()}


# --------------------------------------------------------------------------
@dataclass
class Layout:
    geom: Geometry

    def __post_init__(self):
        g = self.geom
        self.pairs = g.pairing()
        self.nm, self.ncw, self.npair = g.N_mt, g.N_cw, len(self.pairs)
        i = 0

        def take(n):
            nonlocal i
            s = slice(i, i + n)
            i += n
            return s

        self.i_trip = take(3 * self.nm)      # triplet poses
        self.i_link = take(3 * self.nm)      # linker poses
        self.i_pin = take(3 * self.npair)    # pinhead poses
        self.i_base = take(3 * self.npair)   # base poses
        self.i_spoke = take(self.ncw)        # spoke direction offsets
        self.n_pose = i
        self.b_spoke = take(self.ncw)
        self.b_pin = take(self.npair)
        self.b_base = take(self.npair)
        self.b_linkC = take(self.nm)
        self.b_linkA = take(self.nm)
        self.n_total = i


def assemble(z, g: Geometry, lay: Layout, reg=(0.0, 0.0)) -> dict:
    trip = z[lay.i_trip].reshape(-1, 3)
    link = z[lay.i_link].reshape(-1, 3)
    pin = z[lay.i_pin].reshape(-1, 3)
    base = z[lay.i_base].reshape(-1, 3)
    a_spoke = z[lay.i_spoke]
    bs, bp, bb = z[lay.b_spoke], z[lay.b_pin], z[lay.b_base]
    bC, bA = z[lay.b_linkC], z[lay.b_linkA]

    loc = _triplet_local(g, reg)
    T = [_place(loc, trip[i]) for i in range(lay.nm)]

    phi_cw = np.arange(g.N_cw) * 360.0 / g.N_cw
    spoke_dir = phi_cw + a_spoke
    head = np.stack([g.hub_radius * _u(p) for p in phi_cw])
    spoke_tip = np.stack([head[j] + g.spoke_rod * (1 - bs[j]) * _u(spoke_dir[j])
                          for j in range(g.N_cw)])

    L = []
    for i in range(lay.nm):
        o, th = link[i, :2], link[i, 2]
        L.append(dict(vertex=o,
                      end_A=o + g.linker_arm_A * (1 - bA[i]) * _u(th),
                      end_C=o + g.linker_arm_C * (1 - bC[i]) * _u(th + g.linker_vertex_deg),
                      theta=th))

    P, B = [], []
    for p, (j, t) in enumerate(lay.pairs):
        o, th = pin[p, :2], pin[p, 2]
        span = g.pinhead_span * (1 - bp[p])
        ex, ey = _u(th), _u(th + 90)
        P.append(dict(spoke_end=o, A_end=o + span * ex,
                      base_pt=o + g.pinhead_base_frac[0] * span * ex
                              + g.pinhead_base_frac[1] * span * ey, theta=th))
        o2, th2 = base[p, :2], base[p, 2]
        B.append(dict(pin_end=o2,
                      link_end=o2 + g.base_length * (1 - bb[p]) * _u(th2), theta=th2))

    return dict(T=T, L=L, P=P, B=B, head=head, spoke_tip=spoke_tip,
                spoke_dir=spoke_dir, trip=trip, link=link, pin=pin, base=base)


# --------------------------------------------------------------------------
def bond_gaps(st, g: Geometry, lay: Layout):
    """Vector gap at every bond, keyed by bond type."""
    nm = lay.nm
    prev = (np.arange(nm) - 1) % nm
    out = {
        "linker-A": [st["L"][i]["end_A"] - st["T"][prev[i]]["A_pf8"] for i in range(nm)],
        "linker-C": [st["L"][i]["end_C"] - st["T"][i]["C_pf89"] for i in range(nm)],
        "pinhead-triplet": [], "pinhead-spoke": [], "base-pinhead": [], "base-linker": [],
    }
    for p, (j, t) in enumerate(lay.pairs):
        out["pinhead-triplet"].append(st["P"][p]["A_end"] - st["T"][t]["A_pf34"])
        out["pinhead-spoke"].append(st["P"][p]["spoke_end"] - st["spoke_tip"][j])
        out["base-pinhead"].append(st["B"][p]["pin_end"] - st["P"][p]["base_pt"])
        out["base-linker"].append(st["B"][p]["link_end"] - st["L"][t]["vertex"])
    return {k: (np.array(v).reshape(-1, 2) if len(v) else np.zeros((0, 2)))
            for k, v in out.items()}


def joint_deviations(st, g: Geometry, lay: Layout):
    """Deviation of each joint from its wild-type rest angle (degrees)."""
    dev = {"spoke": _wrap(st["spoke_dir"] - np.arange(g.N_cw) * 360.0 / g.N_cw - g.rest_spoke)}
    pin, trip, base = [], [], []
    for p, (j, t) in enumerate(lay.pairs):
        sd = st["spoke_dir"][j]
        pin.append(_wrap(st["P"][p]["theta"] - sd - g.rest_pinhead))
        trip.append(_wrap(st["trip"][t][2] - sd - g.rest_triplet))
        base.append(_wrap(st["B"][p]["theta"] - sd - g.rest_base))
    dev["pinhead"] = np.array(pin) if pin else np.zeros(0)
    dev["triplet"] = np.array(trip) if trip else np.zeros(0)
    dev["base"] = np.array(base) if base else np.zeros(0)
    dev.update(contact_deviations(st, g, lay))
    return dev


def contact_deviations(st, g: Geometry, lay: Layout):
    """Rotation of each strand within the frame of the tubule it grips.

    Microtubules are effectively rigid, so what matters at these contacts
    is how far the strand has twisted relative to the lattice it binds --
    measured in the tubule's own body frame and compared to
    :data:`CONTACT_REST`.
    """
    nm = lay.nm
    prev = (np.arange(nm) - 1) % nm
    lA = np.array([_wrap(st["L"][i]["theta"] - st["trip"][prev[i]][2] - CONTACT_REST["linker-A"])
                   for i in range(nm)])
    lC = np.array([_wrap(st["L"][i]["theta"] - st["trip"][i][2] - CONTACT_REST["linker-C"])
                   for i in range(nm)])
    pA = np.array([_wrap(st["P"][p]["theta"] - st["trip"][t][2] - CONTACT_REST["pinhead-A"])
                   for p, (j, t) in enumerate(lay.pairs)]) if lay.pairs else np.zeros(0)
    return {"linker-A": lA, "linker-C": lC, "pinhead-A": pA}


def tubule_positions(st, g: Geometry, lay: Layout):
    names = TUBULES[: g.MTn]
    P = np.array([st["T"][i]["centres"][t] for t in names for i in range(lay.nm)])
    unit = np.tile(np.arange(lay.nm), len(names))
    tub = np.repeat(np.arange(len(names)), lay.nm)
    return P, unit, tub


def tubule_overlaps(st, g, lay):
    P, unit, _ = tubule_positions(st, g, lay)
    d = np.linalg.norm(P[:, None] - P[None, :], axis=2)
    n = len(P)
    iu = np.triu_indices(n, 1)
    diff = unit[iu[0]] != unit[iu[1]]
    return np.where(diff, np.maximum(0.0, 2 * g.tubule_radius - d[iu]), 0.0)


def strand_clearances(st, g: Geometry, lay: Layout):
    """Clearance of linker / base / spoke strands against microtubules.

    Returns a dict of (min clearance, overlap depth) per strand type.
    Negative clearance means the strand is inside a tubule wall. The
    tubules a strand legitimately binds are excluded.
    """
    P, unit, tub = tubule_positions(st, g, lay)
    Rt = g.tubule_radius
    nm = lay.nm
    prev = (np.arange(nm) - 1) % nm
    ot_idx = g.MTn - 1          # index of the outer tubule the linker binds
    res = {}

    def scan(segments, exclude_per_seg):
        """Exclusions are (unit, tubule_index) pairs -- only the tubules a
        strand genuinely binds are ignored, so a strand cutting across a
        *different* tubule of the same triplet is still caught."""
        best, over, worst = np.inf, 0.0, None
        for (a, b), excl in zip(segments, exclude_per_seg):
            keep = np.ones(len(P), bool)
            for (uu, tt) in excl:
                keep &= ~((unit == uu) & (tub == tt))
            if not keep.any():
                continue
            d = _seg_point_dist(a, b, P[keep]) - Rt
            k = int(np.argmin(d))
            if float(d[k]) < best:
                best = float(d[k])
                idx = np.where(keep)[0][k]
                worst = f"{TUBULES[tub[idx]]}{unit[idx]}"
            over = max(over, float(max(0.0, -d.min())))
        return dict(min_clearance=None if best is np.inf else best, overlap=over, nearest=worst)

    # A-C linker: binds its own triplet's outer tubule and the previous A
    segs, excl = [], []
    for i in range(nm):
        e = [(i, ot_idx), (prev[i], 0)]
        segs += [(st["L"][i]["end_C"], st["L"][i]["vertex"]),
                 (st["L"][i]["vertex"], st["L"][i]["end_A"])]
        excl += [e, e]
    res["linker"] = scan(segs, excl)

    segs, excl = [], []
    for p, (j, t) in enumerate(lay.pairs):
        segs.append((st["B"][p]["pin_end"], st["B"][p]["link_end"]))
        excl.append([])
    res["base"] = scan(segs, excl) if segs else dict(min_clearance=None, overlap=0.0, nearest=None)

    segs, excl = [], []
    for p, (j, t) in enumerate(lay.pairs):
        segs.append((st["head"][j], st["spoke_tip"][j]))
        excl.append([])
    res["spoke"] = scan(segs, excl) if segs else dict(min_clearance=None, overlap=0.0, nearest=None)
    return res


# --------------------------------------------------------------------------
def residuals(z, g, lay, k_bond, k_angle, k_buckle, k_steric, k_uniform, reg, bands=None):
    st = assemble(z, g, lay, reg)
    parts = []

    gaps = bond_gaps(st, g, lay)
    for name, gv in gaps.items():
        if len(gv):
            parts.append(k_bond * BOND_STRENGTH[name] * gv.ravel())

    dev = joint_deviations(st, g, lay)
    for jname, v in dev.items():
        if len(v):
            parts.append(k_angle * angle_penalty(v, jname, bands))

    # keep the triplet ring evenly spaced -- the primary scaffold
    pos = st["trip"][:, :2]
    phi = np.degrees(np.arctan2(pos[:, 1], pos[:, 0]))
    target = phi[0] + np.arange(lay.nm) * 360.0 / lay.nm
    parts.append(k_uniform * _wrap(phi - target))
    r = np.linalg.norm(pos, axis=1)
    parts.append(0.3 * k_uniform * (r - r.mean()))

    parts.append(k_buckle * z[lay.n_pose:] * 100.0)
    parts.append(k_steric * tubule_overlaps(st, g, lay))
    return np.concatenate(parts)


def _initial_guess(g: Geometry, lay: Layout, reg) -> np.ndarray:
    """Seed from the ideal wild-type-like arrangement."""
    z = np.zeros(lay.n_total)
    Rt = g.tubule_radius
    # place triplets on a ring whose radius follows the spoke reach
    r0 = g.hub_radius + g.spoke_rod + g.pinhead_span * np.cos(np.radians(g.rest_pinhead))
    for i in range(lay.nm):
        phi = i * 360.0 / lay.nm
        z[lay.i_trip][3 * i: 3 * i + 3] = [r0 * np.cos(np.radians(phi)),
                                           r0 * np.sin(np.radians(phi)),
                                           phi + g.rest_triplet]
    z[lay.i_trip] = z[lay.i_trip]  # (slice assignment above already applied)

    tri = z[lay.i_trip].reshape(-1, 3)
    loc = _triplet_local(g, reg)
    for i in range(lay.nm):
        T = _place(loc, tri[i])
        nxt = (i + 1) % lay.nm
        v = T["C_pf89"]
        z[lay.i_link][3 * i: 3 * i + 3] = [v[0], v[1], tri[i][2]]
    lk = z[lay.i_link].reshape(-1, 3)
    for i in range(lay.nm):
        lk[i, :2] = lk[i, :2] - g.linker_arm_C * _u(lk[i, 2] + g.linker_vertex_deg)
    z[lay.i_link] = lk.ravel()

    phi_cw = np.arange(g.N_cw) * 360.0 / g.N_cw
    for p, (j, t) in enumerate(lay.pairs):
        tip = (g.hub_radius + g.spoke_rod) * _u(phi_cw[j])
        z[lay.i_pin][3 * p: 3 * p + 3] = [tip[0], tip[1], phi_cw[j] + g.rest_pinhead]
        ex = _u(phi_cw[j] + g.rest_pinhead)
        ey = _u(phi_cw[j] + g.rest_pinhead + 90)
        bp = tip + g.pinhead_base_frac[0] * g.pinhead_span * ex \
             + g.pinhead_base_frac[1] * g.pinhead_span * ey
        z[lay.i_base][3 * p: 3 * p + 3] = [bp[0], bp[1], phi_cw[j] + g.rest_base]
    return z


# --------------------------------------------------------------------------
@dataclass
class Solution:
    geom: Geometry
    z: np.ndarray
    lay: Layout
    state: dict
    reg: tuple
    success: bool
    outer_diameter: float
    a_ring_diameter: float
    lumen_diameter: float        # central aperture (inner tubule wall)
    triplet_tilt: float          # mean triplet axis vs radial, degrees
    joint_strain: dict
    joint_bands: dict
    buckling: dict
    bond_force: dict
    worst_bond: str
    n_clashes: int
    max_overlap: float
    strand: dict
    unattached_triplets: list

    def report(self) -> str:
        g = self.geom
        L = [f"N_cw={g.N_cw}  N_mt={g.N_mt}  MTn={g.MTn}  "
             f"{'converged' if self.success else 'DID NOT CONVERGE'}",
             f"  hub radius     : {g.hub_radius:6.2f} nm   tubule radius {g.tubule_radius:5.2f} nm",
             f"  centriole diam : {self.outer_diameter:6.2f} nm  "
             f"(A-tubule ring {self.a_ring_diameter:.2f} nm, lumen {self.lumen_diameter:.2f} nm)",
             f"  triplet tilt   : {self.triplet_tilt:6.2f} deg from radial"]
        if self.reg != (0.0, 0.0):
            L.append(f"  register shift : pinhead {self.reg[0]:+.0f} pf, linker-C {self.reg[1]:+.0f} pf")
        if self.unattached_triplets:
            L.append(f"  unattached triplets (no pinhead): {self.unattached_triplets}")
        L.append("  joint rotation (deg from wild-type rest)   [band]:")
        for k, v in self.joint_strain.items():
            L.append(f"      {k:<9}: rms {v['rms']:6.2f}  max {v['max']:+7.2f}   "
                     f"{self.joint_bands[k]}")
        buck = {k: v for k, v in self.buckling.items() if abs(v) > 0.01}
        L.append("  buckling (% contour lost): "
                 + (", ".join(f"{k} {v:.2f}%" for k, v in buck.items()) if buck else "none"))
        L.append("  bond load (nm gap x strength; higher = closer to rupture):")
        for k, v in sorted(self.bond_force.items(), key=lambda kv: -kv[1]["force"]):
            L.append(f"      {k:<16}: force {v['force']:7.3f}  gap {v['gap']:6.3f} nm")
        L.append(f"    -> closest to rupture: {self.worst_bond}")
        L.append(f"  MT-MT clashes  : {self.n_clashes} pairs, worst overlap {self.max_overlap:.3f} nm")
        L.append("  strand clearance vs microtubules (negative = clashing):")
        for k, v in self.strand.items():
            c = v["min_clearance"]
            L.append(f"      {k:<8}: min {c:+7.3f} nm  (nearest {v.get('nearest')})"
                     if c is not None else f"      {k:<8}: n/a")
        return "\n".join(L)


def solve(
    geom: Geometry,
    k_bond: float = 12.0,
    k_angle: float = 0.25,
    k_buckle: float = 3.0,
    k_steric: float = 30.0,
    k_uniform: float = 6.0,
    max_buckle: float = 0.35,
    register_shift: bool = False,
    shift_range=(-1, 0, 1),
    bands: Optional[dict] = None,
    base_buckle: Optional[float] = None,
) -> Solution:
    """Relax the network. Optionally search protofilament register shifts.

    With `register_shift=True` the pinhead and linker-C contacts are
    allowed to slide to neighbouring protofilaments; every combination in
    `shift_range` is solved and the lowest-cost one returned.

    `bands` overrides :data:`JOINT_BANDS` -- used by
    :func:`band_sensitivity` to test whether a result depends on them.

    `base_buckle` *drives* the triplet base rather than letting the solver
    choose it: 0.0 holds the base fully extended (its contour length, the
    state the schematic depicts), and larger values bend it so its
    end-to-end chord shortens to `base_length * (1 - base_buckle)`. This is
    the degree of freedom proposed for iris-like blooming, so driving it
    turns blooming into a controlled experiment -- see
    :func:`blooming_scan`.
    """
    combos = [(a, b) for a in shift_range for b in shift_range] if register_shift else [(0.0, 0.0)]
    best = None
    for reg in combos:
        s = _solve_one(geom, reg, k_bond, k_angle, k_buckle, k_steric,
                       k_uniform, max_buckle, bands, base_buckle)
        if best is None or s[1] < best[1]:
            best = (s[0], s[1])
    return best[0]


def _solve_one(geom, reg, k_bond, k_angle, k_buckle, k_steric, k_uniform,
               max_buckle, jbands=None, base_buckle=None):
    g = geom
    lay = Layout(g)
    z0 = _initial_guess(g, lay, reg)
    lo = np.full(lay.n_total, -np.inf)
    hi = np.full(lay.n_total, np.inf)
    lo[lay.n_pose:] = 0.0
    hi[lay.n_pose:] = max_buckle
    if base_buckle is not None:
        # pin the base bend instead of letting the solver pick it
        lo[lay.b_base] = base_buckle
        hi[lay.b_base] = base_buckle + 1e-9
        z0[lay.b_base] = base_buckle

    out = least_squares(residuals, z0, bounds=(lo, hi), method="trf",
                        x_scale="jac", ftol=1e-8, xtol=1e-8, gtol=1e-8, max_nfev=8000,
                        args=(g, lay, k_bond, k_angle, k_buckle, k_steric, k_uniform, reg, jbands))
    z = out.x
    st = assemble(z, g, lay, reg)

    gaps = bond_gaps(st, g, lay)
    bond = {}
    for k, v in gaps.items():
        if len(v) == 0:
            continue
        d = np.linalg.norm(v, axis=1)
        bond[k] = dict(gap=float(d.max()), force=float(d.max() * BOND_STRENGTH[k]))
    worst = max(bond, key=lambda k: bond[k]["force"]) if bond else "-"

    dev = joint_deviations(st, g, lay)
    strain, bands = {}, {}
    for k, v in dev.items():
        strain[k] = dict(rms=float(np.sqrt(np.mean(v**2))) if len(v) else 0.0,
                         max=float(v[np.argmax(np.abs(v))]) if len(v) else 0.0)
        bands[k] = grade(v, k, jbands)

    bl = dict(spoke_rod=lay.b_spoke, pinhead=lay.b_pin, base=lay.b_base,
              linker_armC=lay.b_linkC, linker_armA=lay.b_linkA)
    buckling = {k: float(100 * np.max(z[s])) if z[s].size else 0.0 for k, s in bl.items()}

    ov = tubule_overlaps(st, g, lay)
    P, _, _ = tubule_positions(st, g, lay)
    outer = 2.0 * (np.linalg.norm(P, axis=1).max() + g.tubule_radius)
    a_ring = 2.0 * float(np.mean([np.linalg.norm(st["T"][i]["centres"]["A"])
                                  for i in range(lay.nm)]))
    lumen = 2.0 * max(0.0, float(np.linalg.norm(P, axis=1).min() - g.tubule_radius))
    tilt = float(np.mean([_wrap(st["trip"][i][2] - np.degrees(np.arctan2(
        st["T"][i]["centres"]["A"][1], st["T"][i]["centres"]["A"][0])))
        for i in range(lay.nm)]))
    attached = {t for _, t in lay.pairs}
    sol = Solution(
        geom=g, z=z, lay=lay, state=st, reg=reg, success=bool(out.success),
        outer_diameter=float(outer), a_ring_diameter=a_ring,
        lumen_diameter=lumen, triplet_tilt=tilt,
        joint_strain=strain, joint_bands=bands, buckling=buckling,
        bond_force=bond, worst_bond=worst,
        n_clashes=int((ov > 1e-6).sum()), max_overlap=float(ov.max()) if ov.size else 0.0,
        strand=strand_clearances(st, g, lay),
        unattached_triplets=sorted(set(range(lay.nm)) - attached),
    )
    return sol, float(0.5 * np.sum(out.fun**2))


# --------------------------------------------------------------------------
def draw(sol: Solution, ax=None, show_pf_labels=False, title=None):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, Polygon

    g, st, lay = sol.geom, sol.state, sol.lay
    if ax is None:
        _, ax = plt.subplots(figsize=(7.5, 7.5))

    ax.add_patch(Polygon(st["head"], closed=True, facecolor="none",
                         edgecolor=COLOR_HUB, lw=3, zorder=1))
    for j in range(g.N_cw):
        ax.plot(*np.c_[st["head"][j], st["spoke_tip"][j]], color=COLOR_SPOKE, lw=5,
                solid_capstyle="round", zorder=2)
        ax.add_patch(Circle(st["head"][j], g.head_length / 2, facecolor=COLOR_SPOKE,
                            edgecolor="#231f20", lw=0.6, zorder=3))
    for p in range(lay.npair):
        P, B = st["P"][p], st["B"][p]
        ax.add_patch(Polygon([P["spoke_end"], P["base_pt"], P["A_end"]], closed=True,
                             facecolor=COLOR_PINHEAD, edgecolor="#231f20", lw=0.6, zorder=3))
        ax.plot(*np.c_[B["pin_end"], B["link_end"]], color=COLOR_BASE, lw=3,
                solid_capstyle="round", zorder=2)
    for i in range(lay.nm):
        L = st["L"][i]
        ax.plot(*np.c_[L["end_C"], L["vertex"], L["end_A"]], color=COLOR_LINKER, lw=5,
                solid_capstyle="round", solid_joinstyle="round", zorder=2)

    ov = tubule_overlaps(st, g, lay)
    P, unit, _ = tubule_positions(st, g, lay)
    bad = set()
    if ov.size:
        iu = np.triu_indices(len(P), 1)
        for k in np.where(ov > 1e-6)[0]:
            bad.add(int(unit[iu[0][k]])); bad.add(int(unit[iu[1][k]]))
    for i in range(lay.nm):
        face = "#e08a8a" if i in bad else COLOR_MT
        for t in TUBULES[: g.MTn]:
            c = st["T"][i]["centres"][t]
            for idx, rel in g.pf_angles(t).items():
                pos = c + g.tubule_radius * _u(st["trip"][i][2] + rel)
                ax.add_patch(Circle(pos, g.tubule_radius * 0.23, facecolor=face,
                                    edgecolor="#231f20", lw=0.5, zorder=4))
                if show_pf_labels and i == 0:
                    ax.annotate(str(idx), pos, ha="center", va="center", fontsize=4, zorder=5)

    lim = 1.1 * (np.abs(P).max() + g.tubule_radius)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_aspect("equal"); ax.set_xlabel("nm")
    if title is None:
        rms = np.sqrt(np.mean([v["rms"] ** 2 for v in sol.joint_strain.values()]))
        title = (f"N_cw={g.N_cw} N_mt={g.N_mt}  diam {sol.outer_diameter:.0f} nm  "
                 f"joint rms {rms:.1f}$\\degree$  clashes {sol.n_clashes}")
    ax.set_title(title, fontsize=10)
    return ax


# --------------------------------------------------------------------------
DERIVED_PARAMS = {"linker_length": "end-to-end A-C linker span; scales both arms",
                  "n_pf_A": "A-tubule protofilament count; also sets tubule radius",
                  "n_pf_B": "B-tubule protofilament count",
                  "n_pf_C": "C-tubule protofilament count"}


def _linker_span(g: Geometry) -> float:
    a, b = g.linker_arm_A, g.linker_arm_C
    return float(np.sqrt(a * a + b * b - 2 * a * b * np.cos(np.radians(g.linker_vertex_deg))))


def set_param(geom: Geometry, name: str, value) -> Geometry:
    if name == "linker_length":
        f = float(value) / _linker_span(geom)
        return replace(geom, linker_arm_C=geom.linker_arm_C * f,
                       linker_arm_A=geom.linker_arm_A * f)
    if name in ("n_pf_A", "n_pf_B", "n_pf_C"):
        n = dict(geom.n_pf); n[name[-1]] = int(value)
        return replace(geom, n_pf=n)
    if not hasattr(geom, name):
        raise ValueError(f"unknown parameter {name!r}; expected a Geometry field "
                         f"or one of {sorted(DERIVED_PARAMS)}")
    return replace(geom, **{name: value})


def summarise(sol: Solution) -> dict:
    g, js = sol.geom, sol.joint_strain
    rms = float(np.sqrt(np.mean([v["rms"] ** 2 for v in js.values()])))
    rec = {"N_cw": g.N_cw, "N_mt": g.N_mt, "MTn": g.MTn, "spoke_rod": g.spoke_rod,
           "base_length": g.base_length, "pinhead_span": g.pinhead_span,
           "linker_length": round(_linker_span(g), 3), "n_pf_A": g.n_pf["A"],
           "diameter_nm": round(sol.outer_diameter, 2),
           "A_ring_nm": round(sol.a_ring_diameter, 2),
           "joint_rms_deg": round(rms, 2)}
    for k, v in js.items():
        rec[f"{k}_deg"] = round(v["max"], 2)
        rec[f"{k}_band"] = sol.joint_bands[k]
    rec["lumen_nm"] = round(sol.lumen_diameter, 2)
    rec["triplet_tilt_deg"] = round(sol.triplet_tilt, 2)
    rec["max_buckle_pct"] = round(max(sol.buckling.values()), 2)
    for k, v in sol.buckling.items():
        rec[f"buckle_{k}_pct"] = round(v, 3)
    for k, v in sol.bond_force.items():
        rec[f"bond_{k}"] = round(v["force"], 4)
    rec["worst_bond"] = sol.worst_bond
    rec["worst_bond_force"] = round(sol.bond_force.get(sol.worst_bond, {"force": 0})["force"], 3)
    rec["n_clashes"] = sol.n_clashes
    rec["max_overlap_nm"] = round(sol.max_overlap, 3)
    for k, v in sol.strand.items():
        rec[f"{k}_clear_nm"] = None if v["min_clearance"] is None else round(v["min_clearance"], 3)
    rec["n_unattached"] = len(sol.unattached_triplets)
    rec["converged"] = sol.success
    return rec


@dataclass
class ModeAnalysis:
    """Result of :func:`mode_analysis` -- the network's softest motions."""

    geom: Geometry
    lay: "Layout"
    z0: np.ndarray
    reg: tuple
    eigenvalues: np.ndarray          # ascending
    eigenvectors: np.ndarray         # columns, pose-space
    table: "object"                  # pandas DataFrame summary

    def __repr__(self):
        return f"ModeAnalysis(N_mt={self.geom.N_mt}, {len(self.eigenvalues)} modes)"


def _rotation_generator(z0, lay) -> np.ndarray:
    """Unit vector for 'rotate the whole assembly about the centriole axis'.

    Not an exact symmetry of the model -- the SAS-6 head positions are
    pinned to fixed radial directions -- so this appears as a very soft
    'pinwheel' rather than a true zero mode. Reported per mode so
    rotation-like modes are obvious rather than mistaken for iris motion.
    """
    v = np.zeros(lay.n_pose)
    for sl, n in ((lay.i_trip, lay.nm), (lay.i_link, lay.nm),
                  (lay.i_pin, lay.npair), (lay.i_base, lay.npair)):
        arr = z0[sl].reshape(n, 3)
        out = np.zeros((n, 3))
        for i in range(n):
            out[i] = [-arr[i, 1], arr[i, 0], 1.0]
        v[sl] = out.ravel()
    v[lay.i_spoke] = 1.0
    nrm = np.linalg.norm(v)
    return v / nrm if nrm > 0 else v


def _ring_wavenumber(disp, pos):
    """Dominant angular wavenumber m of a per-unit displacement pattern.

    m = 0 means every unit moves alike (breathing or uniform twist);
    m = 1 the whole ring shifts sideways; m = 2 an elliptical distortion;
    higher m are ripples that alternate around the ring. Far more
    informative than a single 'is it collective?' number.
    """
    n = len(pos)
    c = np.zeros(n, dtype=complex)
    for i in range(n):
        r = pos[i] / max(np.linalg.norm(pos[i]), 1e-12)
        t = np.array([-r[1], r[0]])
        c[i] = (disp[i] @ r) + 1j * (disp[i] @ t)
    amp = np.abs(np.fft.fft(c))
    m = int(np.argmax(amp))
    return min(m, n - m)          # fold to 0..N/2


def _describe_mode(m, d_diam, d_tilt, rot_overlap):
    if rot_overlap > 0.55:
        return "global pinwheel (whole assembly rotates against the hub)"
    if m == 0:
        if abs(d_diam) > 0.02:
            return "breathing / IRIS-LIKE (in-phase, diameter changes)"
        if abs(d_tilt) > 0.02:
            return "uniform twist (in-phase, diameter fixed)"
        return "in-phase, but neither diameter nor tilt changes"
    if m == 1:
        return "whole ring shifts sideways"
    if m == 2:
        return "elliptical distortion"
    return f"ripple around the ring (m={m})"


def mode_analysis(geom: Optional[Geometry] = None, n_modes: int = 6,
                  eps: float = 1e-4, amp: float = 1.0, **kw) -> ModeAnalysis:
    """Find and rank the softest collective motions of the relaxed network.

    Poke the wild-type solution in every direction and measure how hard it
    pushes back. Directions that barely resist are motions the structure
    can actually perform; stiff ones it effectively cannot. This is normal
    mode / elastic network analysis applied to the linkage: the Hessian is
    approximated as ``J^T J`` at the solution and eigendecomposed.

    Each mode is reported with:

    - ``stiffness_rel``   eigenvalue relative to the softest mode
    - ``wavenumber_m``    how the motion varies around the ring (0 = every
                          unit alike, 1 = ring shifts, 2 = ellipse, ...)
    - ``rotation_overlap`` how much of the mode is just global rotation
    - ``d_diameter_nm`` / ``d_tilt_deg`` what it actually changes
    - ``description``     the above in words

    An iris mode is ``wavenumber_m = 0`` with low rotation overlap and
    non-zero ``d_diameter_nm``.

    Restricted to body-pose coordinates: buckling variables sit against
    their lower bound at wild type (nothing may stretch), so perturbing
    them symmetrically is unphysical -- drive those explicitly instead.
    """
    import pandas as pd

    g = geom if geom is not None else Geometry()
    lay = Layout(g)
    sol = solve(g, **kw)
    z0, reg = sol.z.copy(), sol.reg
    npose = lay.n_pose

    def resid(zz):
        return residuals(zz, g, lay, 12.0, 0.25, 3.0, 30.0, 6.0, reg, None)

    r0 = resid(z0)
    J = np.zeros((len(r0), npose))
    for i in range(npose):
        dz = np.zeros(len(z0))
        dz[i] = eps
        J[:, i] = (resid(z0 + dz) - resid(z0 - dz)) / (2 * eps)
    w, V = np.linalg.eigh(J.T @ J)
    rot = _rotation_generator(z0, lay)

    def measure(zz):
        st = assemble(zz, g, lay, reg)
        P, _, _ = tubule_positions(st, g, lay)
        d = 2.0 * (np.linalg.norm(P, axis=1).max() + g.tubule_radius)
        tilt = float(np.mean([_wrap(st["trip"][i][2] - np.degrees(np.arctan2(
            st["T"][i]["centres"]["A"][1], st["T"][i]["centres"]["A"][0])))
            for i in range(lay.nm)]))
        A = np.array([st["T"][i]["centres"]["A"] for i in range(lay.nm)])
        return d, tilt, A

    _, _, A0 = measure(z0)
    rows = []
    for k in range(min(n_modes, npose)):
        v = np.zeros(len(z0))
        v[:npose] = V[:, k]
        dp, tp, Ap = measure(z0 + amp * v)
        dm, tm, Am = measure(z0 - amp * v)
        disp = (Ap - Am) / 2.0
        m = _ring_wavenumber(disp, A0)
        ro = abs(float(rot @ V[:, k]))
        dd, dt = (dp - dm) / 2, (tp - tm) / 2
        rows.append({"mode": k,
                     "stiffness": float(w[k]),
                     "stiffness_rel": round(float(w[k] / max(w[0], 1e-30)), 1),
                     "wavenumber_m": m,
                     "rotation_overlap": round(ro, 2),
                     "d_diameter_nm": round(dd, 3),
                     "d_tilt_deg": round(dt, 3),
                     "description": _describe_mode(m, dd, dt, ro)})
    return ModeAnalysis(geom=g, lay=lay, z0=z0, reg=reg,
                        eigenvalues=w, eigenvectors=V, table=pd.DataFrame(rows))


def soft_modes(geom: Optional[Geometry] = None, n_modes: int = 6, **kw):
    """Convenience wrapper: just the summary table from :func:`mode_analysis`."""
    return mode_analysis(geom, n_modes=n_modes, **kw).table


def draw_mode(ma: ModeAnalysis, k: int, ax=None, target_nm: float = 9.0,
              arrow_scale: float = 1.0, title: Optional[str] = None):
    """Draw one mode: wild type in grey, the deformed shape over it, arrows.

    The eigenvector is rescaled so the largest tubule displacement equals
    `target_nm`, purely so the motion is visible -- modes have no intrinsic
    amplitude.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle

    g, lay = ma.geom, ma.lay
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))

    v = np.zeros(len(ma.z0))
    v[:lay.n_pose] = ma.eigenvectors[:, k]

    def tub(zz):
        st = assemble(zz, g, lay, ma.reg)
        P, _, _ = tubule_positions(st, g, lay)
        return st, P

    st0, P0 = tub(ma.z0)
    _, P1 = tub(ma.z0 + v)
    step = target_nm / max(np.linalg.norm(P1 - P0, axis=1).max(), 1e-12)
    stD, PD = tub(ma.z0 + step * v)

    for st, P, col, lw, alpha, z in ((st0, P0, "#b8bcc4", 1.0, 1.0, 1),
                                     (stD, PD, "#6b3fa0", 1.6, 0.95, 3)):
        for j in range(g.N_cw):
            ax.plot(*np.c_[st["head"][j], st["spoke_tip"][j]],
                    color=col, lw=lw, alpha=alpha, zorder=z)
        for i in range(lay.nm):
            L = st["L"][i]
            ax.plot(*np.c_[L["end_C"], L["vertex"], L["end_A"]],
                    color=col, lw=lw, alpha=alpha, zorder=z)
        for p in range(lay.npair):
            ax.plot(*np.c_[st["P"][p]["spoke_end"], st["P"][p]["A_end"]],
                    color=col, lw=lw, alpha=alpha, zorder=z)
            ax.plot(*np.c_[st["B"][p]["pin_end"], st["B"][p]["link_end"]],
                    color=col, lw=lw, alpha=alpha, zorder=z)
        for c in P:
            ax.add_patch(Circle(c, g.tubule_radius, facecolor="none",
                                edgecolor=col, lw=lw, alpha=alpha, zorder=z))

    for a, b in zip(P0, PD):
        d = (b - a) * arrow_scale
        if np.linalg.norm(d) > 0.35:
            ax.arrow(a[0], a[1], d[0], d[1], head_width=4.0, head_length=5.0,
                     fc="#c0392b", ec="#c0392b", length_includes_head=True, zorder=5)

    row = ma.table.iloc[k]
    lim = 1.12 * (np.abs(PD).max() + g.tubule_radius)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title or f"Mode {k} - {row['description']}\n"
                          f"stiffness {row['stiffness_rel']}x softest   "
                          f"m={row['wavenumber_m']}   "
                          f"$\\Delta$diam {row['d_diameter_nm']:+.2f} nm", fontsize=9)
    return ax


def plot_modes(ma: ModeAnalysis, n: int = 6, path: Optional[str] = None,
               ncols: int = 3):
    """Grid of the softest modes plus the stiffness spectrum."""
    import matplotlib.pyplot as plt

    n = min(n, len(ma.table))
    nrows = -(-n // ncols)
    fig = plt.figure(figsize=(5.2 * ncols, 5.2 * nrows + 3.2))
    gs = fig.add_gridspec(nrows + 1, ncols, height_ratios=[1] * nrows + [0.62])
    for k in range(n):
        draw_mode(ma, k, ax=fig.add_subplot(gs[k // ncols, k % ncols]))

    ax = fig.add_subplot(gs[nrows, :])
    t = ma.table
    cols = ["#c0392b" if r["wavenumber_m"] == 0 and r["rotation_overlap"] <= 0.55
            else ("#95a5a6" if r["rotation_overlap"] > 0.55 else "#6b3fa0")
            for _, r in t.iterrows()]
    ax.bar(t["mode"], t["stiffness_rel"], color=cols)
    ax.set_xlabel("mode")
    ax.set_ylabel("stiffness\n(x softest)")
    ax.set_title("Stiffness spectrum   "
                 "grey = global rotation,  red = in-phase (m=0),  purple = m>0",
                 fontsize=9)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=130, bbox_inches="tight")
    return fig


def blooming_scan(geom: Optional[Geometry] = None, n: int = 9,
                  max_bend: float = 0.32, **kw):
    """Drive the triplet base from bent to fully extended and record the result.

    Tests the proposed iris/blooming mechanism, in which the triplet base
    is bent in the closed state and straightens to its full contour length
    as the centriole opens. The schematic the model is calibrated on shows
    the base **fully extended**, so extension corresponds to
    ``base_buckle = 0`` and bending to positive values.

    Returns a DataFrame ordered from most bent to fully extended, with the
    base chord, the resulting diameters and triplet tilt, and whether the
    cartwheel was disturbed -- the proposal is that it should not be.
    """
    import pandas as pd

    base = geom if geom is not None else Geometry()
    rows = []
    for b in np.linspace(max_bend, 0.0, n):
        sol = solve(base, base_buckle=float(b), **kw)
        rows.append({
            "base_bend_frac": round(float(b), 4),
            "base_chord_nm": round(base.base_length * (1 - b), 2),
            "outer_diam_nm": round(sol.outer_diameter, 2),
            "lumen_diam_nm": round(sol.lumen_diameter, 2),
            "A_ring_nm": round(sol.a_ring_diameter, 2),
            "triplet_tilt_deg": round(sol.triplet_tilt, 2),
            "spoke_dev_deg": round(sol.joint_strain["spoke"]["max"], 2),
            "spoke_band": sol.joint_bands["spoke"],
            "pinhead_dev_deg": round(sol.joint_strain["pinhead"]["max"], 2),
            "linkerA_dev_deg": round(sol.joint_strain["linker-A"]["max"], 2),
            "linkerA_band": sol.joint_bands["linker-A"],
            "base_linker_force": round(sol.bond_force.get("base-linker", {"force": 0})["force"], 3),
            "worst_bond": sol.worst_bond,
            "n_clashes": sol.n_clashes,
        })
    return pd.DataFrame(rows)


def band_sensitivity(geom: Optional[Geometry] = None,
                     factors=(0.5, 0.75, 1.0, 1.5, 2.0), **kw):
    """Re-solve with every joint band scaled, to test robustness.

    The band values are heuristics, so a conclusion is only reportable if
    it survives them moving. Returns a DataFrame with one row per scale
    factor: the key metrics plus each joint's grade.

    >>> band_sensitivity(Geometry(N_cw=8, N_mt=9))
    """
    import pandas as pd

    base = geom if geom is not None else Geometry()
    rows = []
    for f in factors:
        bands = {k: (a * f, b * f) for k, (a, b) in JOINT_BANDS.items()}
        sol = solve(base, bands=bands, **kw)
        rec = {"band_scale": f, "diameter_nm": round(sol.outer_diameter, 2),
               "joint_rms_deg": round(float(np.sqrt(np.mean(
                   [v["rms"] ** 2 for v in sol.joint_strain.values()]))), 2),
               "n_clashes": sol.n_clashes, "worst_bond": sol.worst_bond}
        rec.update({f"{k}_band": v for k, v in sol.joint_bands.items()})
        rows.append(rec)
    return pd.DataFrame(rows)


def sweep(param: str, values, geom: Optional[Geometry] = None, **kw):
    import pandas as pd
    base = geom if geom is not None else Geometry()
    rows = []
    for v in values:
        rec = {param: v}
        rec.update(summarise(solve(set_param(base, param, v), **kw)))
        rows.append(rec)
    return pd.DataFrame(rows)


def sweep2(pa, va, pb, vb, geom: Optional[Geometry] = None, **kw):
    import pandas as pd
    base = geom if geom is not None else Geometry()
    rows = []
    for a in va:
        for b in vb:
            rec = {pa: a, pb: b}
            rec.update(summarise(solve(set_param(set_param(base, pa, a), pb, b), **kw)))
            rows.append(rec)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    print(solve(Geometry()).report())
