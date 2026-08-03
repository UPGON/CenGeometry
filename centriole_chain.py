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
from types import SimpleNamespace
from typing import Optional

import numpy as np
from scipy.optimize import least_squares

from centriole_kinematic import (
    BOND_STRENGTH, CLEAR_TOL, Geometry, _approach_report,
    _reg3, _strand_report, _u, _wrap, angle_penalty, contact_approach,
    contact_rest, grade, strand_clearances, strand_overlaps, tubule_overlaps,
    tubule_positions,
)


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


def chain_deviations(st, g: Geometry, lay: ChainLayout, reg=(0.0, 0.0, 0.0)):
    """Deviation of each joint from its rest angle (degrees).

    The three microtubule contacts are measured against the rest angle of the
    protofilament they are actually on, not against wild type -- see
    :func:`centriole_kinematic.contact_rest`. Without that, sliding a contact
    by one protofilament looked like 27.7 deg of strain that was never there.

    Note that in this formulation the linker is rigid with a single degree of
    freedom, so `linker-A` and `linker-C` are two readings of one angle and are
    numerically identical whenever the ring is symmetric and both contacts
    carry the same shift. That is not a bug -- the linker really is one body --
    but it does mean its orientation is charged twice, so the effective weight
    on the linker is double that on any other joint.
    """
    pin_s, dA, dC = _reg3(reg)
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
    rA, rC = contact_rest("linker-A", g, dA), contact_rest("linker-C", g, dC)
    rP = contact_rest("pinhead-A", g, pin_s)
    dev["linker-A"] = np.array([_wrap(st["L"][i]["theta"] - st["axis"][prev[i]] - rA)
                                for i in range(nm)])
    dev["linker-C"] = np.array([_wrap(st["L"][i]["theta"] - st["axis"][i] - rC)
                                for i in range(nm)])
    dev["pinhead-A"] = np.array(
        [_wrap(st["P"][p]["theta"] - st["axis"][t] - rP)
         for p, (j, t) in enumerate(lay.pairs)]) if lay.pairs else np.zeros(0)
    return dev


def chain_residuals(z, g, lay, k_bond, k_angle, k_buckle, k_steric, k_uniform,
                    reg, bands=None, k_strand=30.0):
    st = assemble_chain(z, g, lay, reg)
    parts = []
    for name, gv in loop_gaps(st, g, lay).items():
        if len(gv):
            parts.append(k_bond * BOND_STRENGTH[name] * gv.ravel())
    for jname, v in chain_deviations(st, g, lay, reg).items():
        if len(v):
            parts.append(k_angle * angle_penalty(v, jname, bands))
    pos = st["trip"][:, :2]
    phi = np.degrees(np.arctan2(pos[:, 1], pos[:, 0]))
    parts.append(k_uniform * _wrap(phi - (phi[0] + np.arange(lay.nm) * 360.0 / lay.nm)))
    parts.append(k_buckle * z[lay.n_pose:] * 100.0)
    parts.append(k_steric * tubule_overlaps(st, g, lay))
    # Strands are solid too. Until this term existed, penetration of a
    # microtubule by the A-C linker was measured and reported but never
    # penalised, so nothing stopped a solve from laying an arm along -- or
    # through -- a tubule wall to reach a shifted protofilament.
    parts.append(k_strand * strand_overlaps(st, g, lay))
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
    approach: dict = field(default_factory=dict)
    spoke_pivot: bool = True
    register_scan: list = field(default_factory=list)
    cost: float = 0.0
    ruptured: list = field(default_factory=list)
    exact_bonds: tuple = ("pinhead-spoke", "pinhead-triplet", "base-pinhead", "linker-C")

    @property
    def reachable(self) -> bool:
        """Does every microtubule contact get reached from outside its tubule?"""
        vals = [float(np.min(v)) for v in self.approach.values() if len(v)]
        return bool(vals) and min(vals) > 0.0

    @property
    def worst_strand_clearance(self) -> float:
        c = [v["min_clearance"] for v in self.strand.values()
             if v["min_clearance"] is not None]
        return float(min(c)) if c else 0.0

    @property
    def feasible(self) -> bool:
        """Reachable *and* no strand buried in a microtubule wall.

        Reachability alone is too weak. A contact approached at 17 deg above
        tangential scores a positive cosine while the arm still lies along the
        wall for several nm -- which is what the drawing shows as a clash. Both
        conditions have to hold.

        The two halves do not rest on the same footing, and the difference
        matters when reporting a result. **Reachability is assumption-free**:
        it asks only whether a strand arrives from outside its tubule or from
        within, and comes out the same whatever strand thickness is assumed and
        whatever weight the steric term carries. **Clearance is not**: it
        counts the strand's own half-width, which is assumed rather than
        measured (:attr:`Geometry.strand_half_width`), so registers sitting
        within a few tenths of a nanometre of zero flip either way if that
        number moves. Treat those as marginal, and prefer to rest a conclusion
        on reachability.
        """
        return self.reachable and self.worst_strand_clearance > -CLEAR_TOL

    def report(self) -> str:
        g = self.geom
        L = [f"N_cw={g.N_cw}  N_mt={g.N_mt}  MTn={g.MTn}  "
             f"{'converged' if self.success else 'DID NOT CONVERGE'}   [chain form]",
             f"  centriole diam : {self.outer_diameter:6.2f} nm  "
             f"(A-tubule ring {self.a_ring_diameter:.2f} nm, lumen {self.lumen_diameter:.2f} nm)",
             f"  triplet tilt   : {self.triplet_tilt:6.2f} deg from radial",
             "  joint rotation (deg from the rest angle of the bound protofilament)   [band]:"]
        for k, v in self.joint_strain.items():
            note = "   (LOCKED: spoke held radial)" if (
                k == "spoke" and not self.spoke_pivot) else ""
            L.append(f"      {k:<9}: rms {v['rms']:6.2f}  max {v['max']:+7.2f}   "
                     f"{self.joint_bands[k]}{note}")
        L.append("  connections held EXACTLY (cannot separate): "
                 + ", ".join(self.exact_bonds))
        L.append("  loop closures (the only ones that can open):")
        for k, v in sorted(self.bond_force.items(), key=lambda kv: -kv[1]["gap"]):
            L.append(f"      {k:<14}: gap {v['gap']:6.3f} nm")
        if any(self.reg):
            L.append(f"  register        : pinhead {self.reg[0]:+.0f} pf, "
                     f"linker-A {self.reg[1]:+.0f} pf, linker-C {self.reg[2]:+.0f} pf")
        buck = {k: v for k, v in self.buckling.items() if abs(v) > 0.01}
        L.append("  buckling: " + (", ".join(f"{k} {v:.2f}%" for k, v in buck.items())
                                   if buck else "none"))
        L.append(f"  MT-MT clashes  : {self.n_clashes} pairs, "
                 f"worst overlap {self.max_overlap:.3f} nm")
        L += _strand_report(self.strand)
        L += _approach_report(self.approach)
        if self.register_scan:
            L += _register_scan_report(self.register_scan)
        return "\n".join(L)


def _register_scan_report(scan: list, n: int = 5) -> list:
    """The register table, reachable candidates first.

    Ordering is by model cost *within* the reachable set, because an
    unreachable register is not a cheaper answer -- it is not an answer. Cost
    itself remains a model quantity and not evidence; see :func:`solve_chain`.
    """
    ok = [r for r in scan if r["feasible"]]
    bad = [r for r in scan if not r["feasible"]]
    L = [f"  register options ({len(ok)} of {len(scan)} feasible), "
         "cheapest feasible first:",
         "      A / C     cost      outer    tilt   worst-gap  clearance  approach"]

    def row(r, mark=" "):
        return (f"    {mark} {r['linkA']:+.0f}/{r['linkC']:+.0f}  {r['cost']:9.1f}  "
                f"{r['outer']:7.2f}  {r['tilt']:+6.2f}  {r['worst_gap']:7.3f}  "
                f"{r['strand_clear']:+8.3f}  {r['approach']:+7.3f}")
    L += [row(r) for r in ok[:n]] or ["      (none -- every register clashes)"]
    if bad:
        L.append(f"      -- {len(bad)} rejected: a strand would arrive through a "
                 "tubule wall, or ends up buried in one. Worst first:")
        L += [row(r, "x") for r in sorted(bad, key=lambda r: r["strand_clear"])[:3]]
    L.append("      Lowest cost is NOT evidence. Compare the outer and tilt "
             "columns against your own measurement.")
    return L


def solve_chain(geom: Optional[Geometry] = None, k_bond: float = 600.0,
                k_angle: float = 0.25, k_buckle: float = 3.0,
                k_steric: float = 30.0, k_uniform: float = 6.0,
                max_buckle: float = 0.35, bands: Optional[dict] = None,
                reg=(0.0, 0.0, 0.0), max_nfev: int = 8000,
                register_shift: bool = False, k_strand: float = 30.0,
                spoke_pivot: bool = True,
                shift_range=(-2, -1, 0, 1, 2)) -> ChainSolution:
    """Solve, optionally searching protofilament registers.

    `spoke_pivot=False` forbids the SAS-6 coiled coil from turning at its
    head, so each spoke points strictly along its own radius -- see
    :func:`_solve_chain_once`.

    With `register_shift=True` both A-C linker contacts are slid over
    `shift_range` whole protofilaments (25 combinations by default) and the
    best solution is returned. The pinhead contact stays at whatever `reg[0]`
    specifies. Every candidate is kept in `register_scan`, and
    :func:`best_registers` extracts the leading few with their metrics.

    **Read this before using the ranking.** Candidates are ordered by model
    cost *within the reachable set*: a register that requires a strand to
    arrive at its protofilament from inside the tubule is not a cheap answer,
    it is not an answer, so it is listed separately rather than ranked.
    Among the reachable ones, cost is still a model quantity built from
    reasoned-not-measured joint bands and bond strengths, so **the cheapest
    register is not the most likely one**. Compare the outer-diameter and tilt
    columns against your own data and treat the ranking as a shortlist.

    Two earlier defects made this search actively misleading, and both are
    fixed:

    - contact rest angles were not re-referenced when the register shifted,
      so every shifted candidate was charged a spurious 27.7 deg of strain
      (:func:`centriole_kinematic.contact_rest`);
    - nothing penalised a strand lying along or through a microtubule, so the
      search could and did return registers whose A-C linker reached its site
      tangentially through the tubule wall.
    """
    if not register_shift:
        return _solve_chain_once(geom, k_bond, k_angle, k_buckle, k_steric,
                                 k_uniform, max_buckle, bands, reg, max_nfev,
                                 k_strand, spoke_pivot)[0]

    pin = _reg3(reg)[0]
    scan = []
    for dA in shift_range:
        for dC in shift_range:
            sol, cost = _solve_chain_once(geom, k_bond, k_angle, k_buckle,
                                          k_steric, k_uniform, max_buckle, bands,
                                          (pin, float(dA), float(dC)), max_nfev,
                                          k_strand, spoke_pivot)
            worst_app = min((float(np.min(v)) for v in sol.approach.values() if len(v)),
                            default=1.0)
            clear = min((v["min_clearance"] for v in sol.strand.values()
                         if v["min_clearance"] is not None), default=0.0)
            rms = float(np.sqrt(np.mean([v["rms"] ** 2
                                         for v in sol.joint_strain.values()])))
            scan.append(dict(
                linkA=float(dA), linkC=float(dC), cost=cost, sol=sol,
                outer=sol.outer_diameter, a_ring=sol.a_ring_diameter,
                lumen=sol.lumen_diameter, tilt=sol.triplet_tilt, joint_rms=rms,
                worst_gap=max(v["gap"] for v in sol.bond_force.values()),
                strand_clear=float(clear), approach=float(worst_app),
                reachable=bool(worst_app > 0.0),
                feasible=bool(worst_app > 0.0 and clear > -CLEAR_TOL),
                n_clashes=sol.n_clashes, converged=sol.success,
                worst_band=max(sol.joint_bands.values(),
                               key=lambda b: ["-", "OK", "HARD", "SEVERE"].index(b))))
    # feasible first, then by cost -- an infeasible register is not a rival
    scan.sort(key=lambda r: (not r["feasible"], r["cost"]))
    best = scan[0]["sol"]
    best.register_scan = scan
    return best


def best_registers(sol: ChainSolution, n: int = 3, feasible_only: bool = True):
    """The leading `n` candidates from a register search, with their metrics.

    Returns a list of plain dicts (no solution objects), ready to tabulate.
    `feasible_only` drops registers that would need a strand to reach its
    protofilament through a tubule wall, or that leave one buried in a wall;
    pass False to see them anyway and read the `reachable` and `strand_clear`
    columns yourself.
    """
    scan = sol.register_scan or []
    rows = [r for r in scan if r["feasible"]] if feasible_only else list(scan)
    return [{k: v for k, v in r.items() if k != "sol"} for r in rows[:n]]


def register_solution(sol: ChainSolution, linkA: float, linkC: float):
    """Pull one candidate's full ChainSolution back out of a register scan."""
    for r in sol.register_scan or []:
        if r["linkA"] == linkA and r["linkC"] == linkC:
            return r["sol"]
    return None


def _solve_chain_once(geom, k_bond, k_angle, k_buckle, k_steric, k_uniform,
                      max_buckle, bands, reg, max_nfev, k_strand=30.0,
                      spoke_pivot=True):
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

    The solve runs in two stages when the strand-vs-tubule term is on. That
    term is a barrier, and the rest-angle initial guess starts *inside* it --
    at a 28 nm spoke the guess buries a linker arm 7.0 nm into a tubule wall.
    Handed the full problem from there the solver crawls out over ~1500
    iterations. Solving first without the barrier costs about 30 iterations
    and lands somewhere already clear, so the second stage starts feasible and
    converges immediately. Both stages minimise the same objective apart from
    that one term, and the reported cost is the second stage's.

    `spoke_pivot=False` holds every spoke at `rest_spoke`, i.e. exactly along
    the radius through its own SAS-6 head. The coiled coil may then still
    shorten by buckling, but it cannot turn at the head -- the assumption most
    treatments of the cartwheel make, in which the spoke strains radially and
    nothing hinges where it meets the hub. It removes `N_cw` degrees of
    freedom, so for a 9-fold centriole the chain form drops from 45 unknowns
    to 36 against 36 loop constraints: exactly determined, with no slack left
    to absorb a mismatch. Expect strain to appear elsewhere, and expect the
    loop closures to open where they previously did not -- which is the point
    of being able to compare the two.
    """
    g = geom if geom is not None else Geometry()
    reg = _reg3(reg)
    lay = ChainLayout(g)

    z0 = np.zeros(lay.n_total)
    z0[lay.i_spoke] = g.rest_spoke
    z0[lay.i_pin] = g.rest_pinhead
    z0[lay.i_trip] = g.rest_triplet
    z0[lay.i_base] = g.rest_base
    step = 360.0 / lay.nm
    z0[lay.i_link] = [g.rest_triplet + contact_rest("linker-C", g, reg[2]) + i * step
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
    if not spoke_pivot:
        # pin rather than remove: least_squares needs lo strictly below hi
        lo[lay.i_spoke] = g.rest_spoke
        hi[lay.i_spoke] = g.rest_spoke + 1e-9
        z0[lay.i_spoke] = g.rest_spoke

    def run(start, ks, nfev):
        """Solve, falling back to an iterative trust-region step if needed.

        `trf` factorises the augmented Jacobian with a dense SVD, and on badly
        conditioned geometries LAPACK occasionally fails to converge and raises
        -- a singlet centriole (MTn=1) at a shifted register does it. That is a
        linear-algebra failure, not a modelling one, so retry with `lsmr`,
        which takes the same step iteratively and never factorises. If even
        that fails, return the starting point flagged unconverged rather than
        taking down a whole 25-register scan for one bad cell.
        """
        kw = dict(bounds=(lo, hi), method="trf", ftol=1e-9, xtol=1e-9,
                  gtol=1e-9, max_nfev=nfev,
                  args=(g, lay, k_bond, k_angle, k_buckle, k_steric,
                        k_uniform, reg, bands, ks))
        for solver in ("exact", "lsmr"):
            try:
                return least_squares(chain_residuals, start, x_scale="jac",
                                     tr_solver=solver, **kw)
            except (np.linalg.LinAlgError, ValueError):
                continue
        return SimpleNamespace(
            x=np.asarray(start, float), success=False,
            cost=float(0.5 * np.sum(chain_residuals(
                start, g, lay, k_bond, k_angle, k_buckle, k_steric,
                k_uniform, reg, bands, ks) ** 2)))

    out = run(z0, 0.0 if k_strand else k_strand, max_nfev)
    if k_strand and strand_overlaps(assemble_chain(out.x, g, lay, reg),
                                    g, lay).max(initial=0.0) > 0.0:
        # Only pay for the second stage when something is actually inside a
        # wall. Where nothing is, the barrier contributes no residual and no
        # gradient, so re-solving would return the same point after one
        # iteration -- and that is the common case.
        out = run(out.x, k_strand, max_nfev)
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

    dev = chain_deviations(st, g, lay, reg)
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
        unattached_triplets=sorted(lay.free_trips), reg=reg,
        approach=contact_approach(st, g, lay), spoke_pivot=bool(spoke_pivot),
        cost=float(out.cost)),
        float(out.cost))
