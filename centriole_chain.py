"""Chain formulation: connection points hold by construction, not by penalty.

`centriole_kinematic.solve()` gives every rigid body a free position and ties
the bodies together with springs. Springs can open, and under a large
perturbation they do -- shortening the SAS-6 spoke by 17 nm opens the
spoke-to-pinhead connection by ~2.8 nm, which is visible as a gap in the
drawing. A connection point separating is wrong at any load.

The cause is redundant coordinates, so the fix is to remove them rather than
to penalise them harder. (Penalising harder was tried: a banded bond penalty
made the objective non-smooth and the solve either crawled or diverged.)

Most of the structure is a *tree*, and every tree edge can be satisfied
exactly by construction:

    spoke tip            -> pinhead origin        EXACT
    pinhead A-end        -> A-tubule contact      EXACT
    pinhead base-point   -> triplet base origin   EXACT
    triplet C pf8/9      -> linker C-end          EXACT

Only two connections close a *loop*, and those genuinely cannot be satisfied
by construction -- they are where a real structure has to strain:

    triplet base tip     -> linker vertex         loop closure
    linker A-end         -> next triplet's A pf8  loop closure

So four of the six connections become unbreakable, including the one that was
visibly separating, and the residual strain is confined to the two places
where it belongs. It is also much smaller: 45 unknowns for a 9-fold
centriole instead of 162.

Each body carries one angle (its rotation about the joint that feeds it)
instead of a full pose. Unattached triplets -- which occur when the cartwheel
and triplet symmetries differ, so some triplets have no pinhead -- have no
parent in the tree and keep a free pose.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.optimize import least_squares

from centriole_kinematic import (
    BOND_STRENGTH, JOINT_BANDS, TUBULES, CONTACT_REST, Geometry, _R, _u, _wrap,
    angle_penalty, grade, tubule_overlaps, tubule_positions, strand_clearances,
    _linker_span,
)


def _reg3(reg):
    """Accept a 2-tuple (pinhead, linker-C) or a 3-tuple, return the 3-tuple."""
    r = tuple(float(v) for v in reg)
    return r if len(r) == 3 else (r[0], 0.0, r[1])


@dataclass
class ChainLayout:
    geom: Geometry

    def __post_init__(self):
        g = self.geom
        self.pairs = g.pairing()
        self.nm, self.ncw, self.npair = g.N_mt, g.N_cw, len(self.pairs)
        self.trip_of_spoke = {t: j for j, t in self.pairs}
        self.free_trips = [i for i in range(self.nm) if i not in self.trip_of_spoke]
        i = 0

        def take(n):
            nonlocal i
            s = slice(i, i + n)
            i += n
            return s

        self.i_spoke = take(self.ncw)                    # spoke vs radial
        self.i_pin = take(self.npair)                    # pinhead vs spoke
        self.i_trip = take(self.npair)                   # triplet vs spoke
        self.i_base = take(self.npair)                   # base vs spoke
        self.i_link = take(self.nm)                      # linker orientation
        self.i_free = take(3 * len(self.free_trips))     # unparented triplets
        self.n_pose = i
        self.b_spoke = take(self.ncw)
        self.b_base = take(self.npair)
        self.n_total = i


def assemble_chain(z, g: Geometry, lay: ChainLayout, reg=(0.0, 0.0, 0.0)) -> dict:
    """Forward kinematics. Every tree edge is exact by construction.

    `reg` shifts each contact by whole protofilaments:
    (pinhead on A, linker on A, linker on C).
    """
    Rt = g.tubule_radius
    pin_shift, linkA_shift, linkC_shift = _reg3(reg)
    a_spoke, a_pin = z[lay.i_spoke], z[lay.i_pin]
    a_trip, a_base = z[lay.i_trip], z[lay.i_base]
    th_link = z[lay.i_link]
    free = z[lay.i_free].reshape(-1, 3) if len(lay.free_trips) else np.zeros((0, 3))
    bs, bb = z[lay.b_spoke], z[lay.b_base]

    phi_cw = np.arange(g.N_cw) * 360.0 / g.N_cw
    spoke_dir = phi_cw + a_spoke
    head = np.stack([g.hub_radius * _u(p) for p in phi_cw])
    spoke_tip = np.stack([head[j] + g.spoke_rod * (1 - bs[j]) * _u(spoke_dir[j])
                          for j in range(g.N_cw)])

    ot = g.outer_tubule
    aA_pin = g.pf_angle("A", g.pin_pf, pin_shift)
    aA_link = g.pf_angle("A", g.linkA_pf, linkA_shift)
    aC_link = g.pf_angle(ot, g.linkC_pf, linkC_shift)

    T = [None] * lay.nm
    axis = np.zeros(lay.nm)
    P, B = [], []

    for p, (j, t) in enumerate(lay.pairs):
        sd = spoke_dir[j]
        th_p = sd + a_pin[p]
        ex, ey = _u(th_p), _u(th_p + 90)
        # --- pinhead: its origin IS the spoke tip. Cannot separate.
        origin = spoke_tip[j]
        A_end = origin + g.pinhead_span * ex
        base_pt = (origin + g.pinhead_base_frac[0] * g.pinhead_span * ex
                   + g.pinhead_base_frac[1] * g.pinhead_span * ey)
        P.append(dict(spoke_end=origin, A_end=A_end, base_pt=base_pt, theta=th_p))

        # --- triplet: placed so its A-contact IS the pinhead A-end. Cannot separate.
        ax = sd + a_trip[p]
        axis[t] = ax
        A_c = A_end - Rt * _u(ax + aA_pin)
        T[t] = _triplet_from(g, A_c, ax, Rt, aA_link, aC_link, ot)

        # --- base: its origin IS the pinhead base-point. Cannot separate.
        bdir = sd + a_base[p]
        tip = base_pt + g.base_length * (1 - bb[p]) * _u(bdir)
        B.append(dict(pin_end=base_pt, link_end=tip, theta=bdir))

    for k, t in enumerate(lay.free_trips):          # no pinhead to hang from
        ax = free[k, 2]
        axis[t] = ax
        T[t] = _triplet_from(g, free[k, :2], ax, Rt, aA_link, aC_link, ot)

    # --- linker: its C-end IS the triplet's C pf8/9 point. Cannot separate.
    L = []
    for i in range(lay.nm):
        th = th_link[i]
        end_C = T[i]["C_pf89"]
        vertex = end_C - g.linker_arm_C * _u(th + g.linker_vertex_deg)
        end_A = vertex + g.linker_arm_A * _u(th)
        L.append(dict(end_C=end_C, vertex=vertex, end_A=end_A, theta=th))

    trip = np.stack([np.array([T[i]["centres"]["A"][0], T[i]["centres"]["A"][1],
                               axis[i]]) for i in range(lay.nm)])
    return dict(T=T, L=L, P=P, B=B, head=head, spoke_tip=spoke_tip,
                spoke_dir=spoke_dir, trip=trip, axis=axis, src=g.pairing())


def _triplet_from(g, A_c, ax, Rt, aA_link, aC_link, ot):
    cen = {"A": np.asarray(A_c, dtype=float)}
    cen["B"] = cen["A"] + g.spacing_ab * _u(ax)
    cen["C"] = cen["B"] + g.spacing_bc * _u(ax)
    return dict(centres=cen,
                A_pf8=cen["A"] + Rt * _u(ax + aA_link),
                C_pf89=cen[ot] + Rt * _u(ax + aC_link))


def loop_gaps(st, g: Geometry, lay: ChainLayout):
    """Only the two connections that close a loop can open at all."""
    prev = (np.arange(lay.nm) - 1) % lay.nm
    out = {"linker-A": np.array([st["L"][i]["end_A"] - st["T"][prev[i]]["A_pf8"]
                                 for i in range(lay.nm)]).reshape(-1, 2)}
    bl = [st["B"][p]["link_end"] - st["L"][t]["vertex"]
          for p, (j, t) in enumerate(lay.pairs)]
    out["base-linker"] = (np.array(bl).reshape(-1, 2) if bl else np.zeros((0, 2)))
    return out


def chain_deviations(st, g: Geometry, lay: ChainLayout):
    dev = {"spoke": _wrap(st["spoke_dir"] - np.arange(g.N_cw) * 360.0 / g.N_cw
                          - g.rest_spoke)}
    pin, trip, base = [], [], []
    for p, (j, t) in enumerate(lay.pairs):
        sd = st["spoke_dir"][j]
        pin.append(_wrap(st["P"][p]["theta"] - sd - g.rest_pinhead))
        trip.append(_wrap(st["axis"][t] - sd - g.rest_triplet))
        base.append(_wrap(st["B"][p]["theta"] - sd - g.rest_base))
    dev["pinhead"] = np.array(pin) if pin else np.zeros(0)
    dev["triplet"] = np.array(trip) if trip else np.zeros(0)
    dev["base"] = np.array(base) if base else np.zeros(0)
    nm, prev = lay.nm, (np.arange(lay.nm) - 1) % lay.nm
    dev["linker-A"] = np.array([_wrap(st["L"][i]["theta"] - st["axis"][prev[i]]
                                      - CONTACT_REST["linker-A"]) for i in range(nm)])
    dev["linker-C"] = np.array([_wrap(st["L"][i]["theta"] - st["axis"][i]
                                      - CONTACT_REST["linker-C"]) for i in range(nm)])
    dev["pinhead-A"] = np.array(
        [_wrap(st["P"][p]["theta"] - st["axis"][t] - CONTACT_REST["pinhead-A"])
         for p, (j, t) in enumerate(lay.pairs)]) if lay.pairs else np.zeros(0)
    return dev


def chain_residuals(z, g, lay, k_bond, k_angle, k_buckle, k_steric, k_uniform,
                    reg, bands=None):
    st = assemble_chain(z, g, lay, reg)
    parts = []
    for name, gv in loop_gaps(st, g, lay).items():
        if len(gv):
            parts.append(k_bond * BOND_STRENGTH[name] * gv.ravel())
    for jname, v in chain_deviations(st, g, lay).items():
        if len(v):
            parts.append(k_angle * angle_penalty(v, jname, bands))
    pos = st["trip"][:, :2]
    phi = np.degrees(np.arctan2(pos[:, 1], pos[:, 0]))
    parts.append(k_uniform * _wrap(phi - (phi[0] + np.arange(lay.nm) * 360.0 / lay.nm)))
    parts.append(k_buckle * z[lay.n_pose:] * 100.0)
    parts.append(k_steric * tubule_overlaps(st, g, lay))
    return np.concatenate(parts)


@dataclass
class ChainSolution:
    geom: Geometry
    z: np.ndarray
    lay: ChainLayout
    state: dict
    success: bool
    outer_diameter: float
    a_ring_diameter: float
    lumen_diameter: float
    triplet_tilt: float
    joint_strain: dict
    joint_bands: dict
    buckling: dict
    bond_force: dict
    worst_bond: str
    n_clashes: int
    max_overlap: float
    strand: dict
    unattached_triplets: list
    reg: tuple = (0.0, 0.0, 0.0)
    register_scan: list = field(default_factory=list)
    cost: float = 0.0
    ruptured: list = field(default_factory=list)
    exact_bonds: tuple = ("pinhead-spoke", "pinhead-triplet", "base-pinhead", "linker-C")

    def report(self) -> str:
        g = self.geom
        L = [f"N_cw={g.N_cw}  N_mt={g.N_mt}  MTn={g.MTn}  "
             f"{'converged' if self.success else 'DID NOT CONVERGE'}   [chain form]",
             f"  centriole diam : {self.outer_diameter:6.2f} nm  "
             f"(A-tubule ring {self.a_ring_diameter:.2f} nm, lumen {self.lumen_diameter:.2f} nm)",
             f"  triplet tilt   : {self.triplet_tilt:6.2f} deg from radial",
             "  joint rotation (deg from wild-type rest)   [band]:"]
        for k, v in self.joint_strain.items():
            L.append(f"      {k:<9}: rms {v['rms']:6.2f}  max {v['max']:+7.2f}   "
                     f"{self.joint_bands[k]}")
        L.append("  connections held EXACTLY (cannot separate): "
                 + ", ".join(self.exact_bonds))
        L.append("  loop closures (the only ones that can open):")
        for k, v in sorted(self.bond_force.items(), key=lambda kv: -kv[1]["gap"]):
            L.append(f"      {k:<14}: gap {v['gap']:6.3f} nm")
        if any(self.reg):
            L.append(f"  register        : pinhead {self.reg[0]:+.0f} pf, "
                     f"linker-A {self.reg[1]:+.0f} pf, linker-C {self.reg[2]:+.0f} pf")
        if self.register_scan:
            L.append("  register options, best model cost first "
                     "(NOTE: lowest cost is not evidence -- compare against data):")
            for r in self.register_scan[:5]:
                L.append(f"      A{r['linkA']:+.0f}/C{r['linkC']:+.0f}  "
                         f"cost {r['cost']:9.2f}  outer {r['outer']:7.2f} nm  "
                         f"tilt {r['tilt']:+6.2f} deg")
        buck = {k: v for k, v in self.buckling.items() if abs(v) > 0.01}
        L.append("  buckling: " + (", ".join(f"{k} {v:.2f}%" for k, v in buck.items())
                                   if buck else "none"))
        L.append(f"  MT-MT clashes  : {self.n_clashes} pairs, "
                 f"worst overlap {self.max_overlap:.3f} nm")
        return "\n".join(L)


def solve_chain(geom: Optional[Geometry] = None, k_bond: float = 600.0,
                k_angle: float = 0.25, k_buckle: float = 3.0,
                k_steric: float = 30.0, k_uniform: float = 6.0,
                max_buckle: float = 0.35, bands: Optional[dict] = None,
                reg=(0.0, 0.0, 0.0), max_nfev: int = 8000,
                register_shift: bool = False,
                shift_range=(-2, -1, 0, 1, 2)) -> ChainSolution:
    """Solve, optionally searching protofilament registers.

    With `register_shift=True` both A-C linker contacts are slid over
    `shift_range` whole protofilaments (25 combinations by default, ~1 s
    each) and the lowest-cost solution is returned. The pinhead contact
    stays at whatever `reg[0]` specifies.

    **The ranking is by model cost, which is not evidence.** Every candidate
    is kept in `register_scan` precisely because the cheapest is often not
    the interesting one: for a spoke shortened to 28 nm the model prefers
    wild-type register, yet the shifted register is what reproduces the
    measured diameter. Read the whole scan and compare against data rather
    than trusting the winner.
    """
    if not register_shift:
        return _solve_chain_once(geom, k_bond, k_angle, k_buckle, k_steric,
                                 k_uniform, max_buckle, bands, reg, max_nfev)[0]

    pin = _reg3(reg)[0]
    best, scan = None, []
    for dA in shift_range:
        for dC in shift_range:
            sol, cost = _solve_chain_once(geom, k_bond, k_angle, k_buckle,
                                          k_steric, k_uniform, max_buckle, bands,
                                          (pin, float(dA), float(dC)), max_nfev)
            scan.append(dict(linkA=float(dA), linkC=float(dC), cost=cost,
                             outer=sol.outer_diameter, tilt=sol.triplet_tilt,
                             worst_gap=max(v["gap"] for v in sol.bond_force.values())))
            if best is None or cost < best[1]:
                best = (sol, cost)
    scan.sort(key=lambda r: r["cost"])
    best[0].register_scan = scan
    return best[0]


def _solve_chain_once(geom, k_bond, k_angle, k_buckle, k_steric, k_uniform,
                      max_buckle, bands, reg, max_nfev):
    """One solve at a fixed register. Returns (solution, least-squares cost).

    `k_bond` defaults to 600, which makes loop closure strictly dominate the
    angular penalties. That matters: the chain form has 45 pose unknowns
    against only 36 loop constraints, so 9 spare degrees of freedom exist and
    exact closure IS attainable -- but only if closing is preferred to holding
    an angle. At the old k_bond=12 the solver traded a 6.4 nm gap against a
    rotation instead, which is exactly the visible separation this formulation
    was meant to remove. Raising it drives every gap below 0.03 nm, pushes the
    strain into the joint angles where it can be read and graded, and happens
    to be ~150x faster because the problem becomes well determined.
    """
    g = geom if geom is not None else Geometry()
    lay = ChainLayout(g)

    z0 = np.zeros(lay.n_total)
    z0[lay.i_spoke] = g.rest_spoke
    z0[lay.i_pin] = g.rest_pinhead
    z0[lay.i_trip] = g.rest_triplet
    z0[lay.i_base] = g.rest_base
    step = 360.0 / lay.nm
    z0[lay.i_link] = [g.rest_triplet + CONTACT_REST["linker-C"] + i * step
                      for i in range(lay.nm)]
    if len(lay.free_trips):
        r0 = g.hub_radius + g.spoke_rod + g.pinhead_span
        fr = np.zeros((len(lay.free_trips), 3))
        for k, t in enumerate(lay.free_trips):
            phi = t * step
            fr[k] = [r0 * np.cos(np.radians(phi)), r0 * np.sin(np.radians(phi)),
                     phi + g.rest_triplet]
        z0[lay.i_free] = fr.ravel()

    lo = np.full(lay.n_total, -np.inf)
    hi = np.full(lay.n_total, np.inf)
    lo[lay.n_pose:] = 0.0
    hi[lay.n_pose:] = max_buckle

    out = least_squares(chain_residuals, z0, bounds=(lo, hi), method="trf",
                        x_scale="jac", ftol=1e-9, xtol=1e-9, gtol=1e-9,
                        max_nfev=max_nfev,
                        args=(g, lay, k_bond, k_angle, k_buckle, k_steric,
                              k_uniform, reg, bands))
    z = out.x
    st = assemble_chain(z, g, lay, reg)

    bond = {}
    for k, v in loop_gaps(st, g, lay).items():
        if len(v):
            d = np.linalg.norm(v, axis=1)
            bond[k] = dict(gap=float(d.max()), force=float(d.max() * BOND_STRENGTH[k]))
    for k in BOND_STRENGTH:                     # exact connections: zero by construction
        bond.setdefault(k, dict(gap=0.0, force=0.0))
    worst = max(bond, key=lambda k: bond[k]["force"]) if bond else "-"

    dev = chain_deviations(st, g, lay)
    strain, bnds = {}, {}
    for k, v in dev.items():
        strain[k] = dict(rms=float(np.sqrt(np.mean(v**2))) if len(v) else 0.0,
                         max=float(v[np.argmax(np.abs(v))]) if len(v) else 0.0)
        bnds[k] = grade(v, k, bands)

    # keep the same five keys as the network solver so summarise() and the app
    # need no special-casing; pinhead and the linker arms are rigid here
    buckling = {"spoke_rod": float(100 * np.max(z[lay.b_spoke])) if lay.ncw else 0.0,
                "pinhead": 0.0,
                "base": float(100 * np.max(z[lay.b_base])) if lay.npair else 0.0,
                "linker_armC": 0.0, "linker_armA": 0.0}

    ov = tubule_overlaps(st, g, lay)
    P, _, _ = tubule_positions(st, g, lay)
    outer = 2.0 * (np.linalg.norm(P, axis=1).max() + g.tubule_radius)
    a_ring = 2.0 * float(np.mean([np.linalg.norm(st["T"][i]["centres"]["A"])
                                  for i in range(lay.nm)]))
    lumen = 2.0 * max(0.0, float(np.linalg.norm(P, axis=1).min() - g.tubule_radius))
    tilt = float(np.mean([_wrap(st["axis"][i] - np.degrees(np.arctan2(
        st["T"][i]["centres"]["A"][1], st["T"][i]["centres"]["A"][0])))
        for i in range(lay.nm)]))

    return (ChainSolution(
        geom=g, z=z, lay=lay, state=st, success=bool(out.success),
        outer_diameter=float(outer), a_ring_diameter=a_ring, lumen_diameter=lumen,
        triplet_tilt=tilt, joint_strain=strain, joint_bands=bnds, buckling=buckling,
        bond_force=bond, worst_bond=worst,
        n_clashes=int((ov > 1e-6).sum()), max_overlap=float(ov.max()) if ov.size else 0.0,
        strand=strand_clearances(st, g, lay),
        unattached_triplets=sorted(lay.free_trips), reg=_reg3(reg),
        cost=float(out.cost)),
        float(out.cost))
