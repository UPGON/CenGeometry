"""Geometry model with explicit cartwheel/pinhead/triplet-base sub-structure.

`centriole.py` is a faithful port of the original MATLAB `centriole.m` and
treats the entire path from the cartwheel hub to the microtubule wall as a
single length `r`. Real centrioles (see cryo-ET cross-section
reconstructions) resolve that path into distinct substructures:

    CID (center) -> SAS-6 cartwheel (hub of SAS-6 dimers + spokes)
    -> pinhead -> triplet base -> MT triplet (A/B/C tubules)

The A-tubule attaches directly to the pinhead (pinhead tip = A's near
edge) -- so the MT ring's radius depends only on `HUB+SPOKE+PINHEAD`, not
on `BASE_LENGTH`. The triplet base is a separate branch: it starts at the
pinhead's MIDPOINT (not its tip), bends to run parallel to the local
triplet axis (A-B-C direction), and its far end (the "tip") is the A-C
linker attachment. Unlike `centriole.py`'s A-C linker (which anchors
directly on the tubule wall at a fixed angle GAMMA), this model's linker
connects neighboring triplets' base tips directly -- i.e. the achieved
linker length is an OUTPUT of the base/pinhead/spoke geometry, not an
independently anchored quantity. GAMMA is retired; there's no longer a
separate wall-anchor angle to tune.

Two independent symmetries are supported: `SYM_MT` (microtubule triplet
count) and `SYM_CW` (cartwheel SAS-6 dimer / spoke count), matching a real
failure mode where a SAS-6 mutation changes cartwheel symmetry without
necessarily changing the MT wall symmetry (or vice versa). When they
differ, each triplet's base is fed by its nearest-neighbor spoke (angular
nearest-neighbor matching) -- some spokes may feed more than one triplet,
others none. In that mismatched case the "base parallel to triplet axis"
property is only approximate (it is exact when SYM_CW == SYM_MT), since a
spoke and its assigned triplet are then generally at different angles.

CAVEAT ON DEFAULT NUMBERS: the default segment lengths are rough,
EYEBALLED proportions read off a reference cryo-ET cross-section image
relative to its scale bar, not precise pixel measurements. Treat them as
a starting point to calibrate against real data, not verified ground
truth.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge

from centriole import CentrioleResult, _num2str, _round2xdigit

COLOR_CID = "#1b2a6b"  # dark blue: central inner density
COLOR_CARTWHEEL = "#bcd9f0"  # light blue: cartwheel hub disc
COLOR_DIMER = "#5b9bd5"  # slightly darker blue: individual SAS-6 dimers
COLOR_SPOKE = "#7fb3e0"  # SAS-6 coiled-coil spoke
COLOR_PINHEAD = "#9b59d0"  # purple: pinhead
COLOR_BASE = "#7c9a3a"  # green: triplet base
COLOR_LINKER = "#3fc9c9"  # cyan: A-C linker
COLOR_MT = "#a6a6a6"  # grey: MT triplet wall
COLOR_LUMEN = "white"  # tubule lumen ("hole")


def _nearest_neighbor_map(angles_from: np.ndarray, angles_to: np.ndarray) -> np.ndarray:
    """For each angle in `angles_from`, index of the nearest angle in `angles_to`."""

    diff = np.angle(np.exp(1j * (angles_from[:, None] - angles_to[None, :])))
    return np.argmin(np.abs(diff), axis=1)


def centriole_v2(
    SYM_MT: int = 9,
    SYM_CW: Optional[int] = None,
    LK: float = 25,
    CID_RADIUS: float = 8,
    HUB_RADIUS: float = 15,
    SPOKE_LENGTH: float = 25,
    PINHEAD_LENGTH: float = 8,
    BASE_LENGTH: float = 20,
    MT_RADIUS: float = 11,
    L: float = 50,
    MTn: int = 3,
    MT_LUMEN_RATIO: float = 0.55,
    FUSION_ARC_DEG: Optional[float] = None,
    CORRECT_OVERLAP: bool = False,
    show_result: bool = True,
    ax: "Optional[plt.Axes]" = None,
    return_details: bool = False,
):
    """Model (and optionally plot) a centriole cross-section with explicit
    CID / cartwheel / pinhead / triplet-base / A-C-linker sub-structure.

    Parameters
    ----------
    SYM_MT : microtubule triplet count (rotational symmetry of the MT wall).
    SYM_CW : cartwheel SAS-6 dimer / spoke count. Defaults to `SYM_MT` (the
        WT case). Set independently to model a symmetry-mismatch mutant;
        see module docstring for how mismatched counts are paired.
    LK : target A-C linker length -- the tilt angle `b` is searched so that
        the distance between neighboring triplets' base-rod TIPS matches
        this value (the linker is no longer independently wall-anchored).
    CID_RADIUS, HUB_RADIUS, SPOKE_LENGTH, PINHEAD_LENGTH, BASE_LENGTH,
    MT_RADIUS : lengths of each named layer, center outward.
    L : triplet width, end of A to end of C -- sets how far B (at L/2) and
        C (at L) sit from A along the local triplet axis.
    MTn : 1 (singlet), 2 (doublet), or 3 (triplet -- the biological default).
    MT_LUMEN_RATIO : inner "hole" radius as a fraction of MT_RADIUS
        (cosmetic).
    FUSION_ARC_DEG : degrees of B/C tubule circumference to omit on the
        side fused to its inward neighbor (B fuses to A, C fuses to B).
        If None (default), derived analytically from the actual
        circle-circle overlap geometry given `L` and `MT_RADIUS` -- pass
        an explicit value to override.
    CORRECT_OVERLAP, show_result, ax, return_details : as in `centriole()`.

    Returns
    -------
    CentrioleResult (or (CentrioleResult, details) if return_details=True)
        (success, overlap, b_deg). `success` means a tilt angle was found
        matching `LK` within tolerance; `overlap` means adjacent MT
        triplets' tubules collide at that angle.
    """

    if MTn not in (1, 2, 3):
        raise ValueError(f"MTn must be 1 (singlet), 2 (doublet), or 3 (triplet); got {MTn}")
    if SYM_CW is None:
        SYM_CW = SYM_MT

    # Per-triplet tubule offset fractions of (L - 2*MT_RADIUS) from the
    # A-tubule (frac 0.0): singlet has only A; doublet has A + a second
    # tubule at frac 1.0 (labeled "B", playing C's linker-connecting role);
    # triplet has A, B (frac 0.5), C (frac 1.0). The LAST entry is always
    # the tubule the A-C linker connects to the next triplet's A-tubule.
    offsets = [0.0] if MTn == 1 else ([0.0, 1.0] if MTn == 2 else [0.0, 0.5, 1.0])
    far_frac = offsets[-1]

    # Rc: radius at which A-tubules sit (pinhead tip -- A attaches directly
    # to the pinhead, so this does NOT include BASE_LENGTH).
    # Rc_mid: radius of the pinhead's MIDPOINT, where the triplet base
    # branches off (base runs from here, parallel to the triplet axis, to
    # the A-C linker attachment -- it does not end at the tubule wall).
    Rc = HUB_RADIUS + SPOKE_LENGTH + PINHEAD_LENGTH
    Rc_mid = HUB_RADIUS + SPOKE_LENGTH + PINHEAD_LENGTH / 2

    print()
    print(f"-------------- desired LK: {LK}")

    t_mt = np.linspace(0, 2 * np.pi, SYM_MT, endpoint=False)  # triplet angles
    t_cw = np.linspace(0, 2 * np.pi, SYM_CW, endpoint=False)  # spoke/dimer angles
    bb = np.linspace(0, np.pi / 2, 362)

    # --- Angle search: representative triplet pair 0 & 1. The A-C linker
    # connects triplet 0's "far" tubule (C, or B for a doublet, or A
    # itself for a singlet) to triplet 1's A-tubule, anchored on each
    # tubule's facing perimeter point -- so the achieved linker length is
    # simply the tubule CENTER-to-center distance minus 2*MT_RADIUS (the
    # two anchor points sit MT_RADIUS in from each center, along the same
    # connecting line). A-tubule positions don't depend on bb (frac=0
    # term has no bb dependence), only the far tubule's does.
    theta0, theta1 = t_mt[0], t_mt[1]
    dir0 = np.pi / 2 - theta0 + bb
    A0_x = (Rc + MT_RADIUS) * np.cos(theta0)
    A0_y = (Rc + MT_RADIUS) * np.sin(theta0)
    A1_x = (Rc + MT_RADIUS) * np.cos(theta1)
    A1_y = (Rc + MT_RADIUS) * np.sin(theta1)
    far0_span = far_frac * (L - 2 * MT_RADIUS)
    far0_x = A0_x - far0_span * np.cos(dir0)
    far0_y = A0_y + far0_span * np.sin(dir0)

    LKmV = np.sqrt((A1_x - far0_x) ** 2 + (A1_y - far0_y) ** 2) - 2 * MT_RADIUS

    idx_candidates = np.where(np.abs(LKmV - LK) < 0.15)[0]
    if idx_candidates.size > 0:
        gaps = np.where(np.diff(idx_candidates) > 1)[0]
        if gaps.size > 0:
            idx_candidates = idx_candidates[gaps[-1] + 1 :]
        best = np.argmin(np.abs(LKmV[idx_candidates] - LK))
        idx = idx_candidates[best]
        success = True
    else:
        idx = int(np.argmin(np.abs(LKmV - LK)))
        success = False
        print("It is not possible to reach the desired A-C LINKER.")

    LKm = _round2xdigit(LKmV[idx])
    b = bb[idx]
    b_deg = _round2xdigit(np.degrees(b))

    # --- Overlap check: unchanged from centriole.py's algorithm -- still
    # about MT tubule circles between adjacent triplets colliding, which
    # is independent of how the base/linker are drawn.
    overlap = False
    if success:
        r1 = np.array([(Rc + MT_RADIUS) * np.cos(t_mt[0]), (Rc + MT_RADIUS) * np.sin(t_mt[0])])
        r2 = np.array([(Rc + MT_RADIUS) * np.cos(t_mt[1]), (Rc + MT_RADIUS) * np.sin(t_mt[1])])
        d2 = np.linalg.norm(r1 - r2)

        if d2 < 2 * MT_RADIUS:
            print("It is not possible to avoid MT overlap.. Increase the SAS-6/pinhead/base reach.")
            overlap = True
        else:

            def _tubule_centers(angle0, angle1, beta):
                x0V = np.array(
                    [
                        (Rc + MT_RADIUS) * np.cos(angle0),
                        (Rc + MT_RADIUS) * np.cos(angle0) - (L - 2 * MT_RADIUS) * np.cos(np.pi / 2 - angle0 + beta),
                        (Rc + MT_RADIUS) * np.cos(angle0) - 0.5 * (L - 2 * MT_RADIUS) * np.cos(np.pi / 2 - angle0 + beta),
                    ]
                )
                y0V = np.array(
                    [
                        (Rc + MT_RADIUS) * np.sin(angle0),
                        (Rc + MT_RADIUS) * np.sin(angle0) + (L - 2 * MT_RADIUS) * np.sin(np.pi / 2 - angle0 + beta),
                        (Rc + MT_RADIUS) * np.sin(angle0) + 0.5 * (L - 2 * MT_RADIUS) * np.sin(np.pi / 2 - angle0 + beta),
                    ]
                )
                x0pV = np.array(
                    [
                        (Rc + MT_RADIUS) * np.cos(angle1),
                        (Rc + MT_RADIUS) * np.cos(angle1) - (L - 2 * MT_RADIUS) * np.cos(np.pi / 2 - angle1 + beta),
                        (Rc + MT_RADIUS) * np.cos(angle1) - 0.5 * (L - 2 * MT_RADIUS) * np.cos(np.pi / 2 - angle1 + beta),
                    ]
                )
                y0pV = np.array(
                    [
                        (Rc + MT_RADIUS) * np.sin(angle1),
                        (Rc + MT_RADIUS) * np.sin(angle1) + (L - 2 * MT_RADIUS) * np.sin(np.pi / 2 - angle1 + beta),
                        (Rc + MT_RADIUS) * np.sin(angle1) + 0.5 * (L - 2 * MT_RADIUS) * np.sin(np.pi / 2 - angle1 + beta),
                    ]
                )
                return x0V, y0V, x0pV, y0pV

            def _any_overlap(x0V, y0V, x0pV, y0pV):
                for jj in range(MTn):
                    for kk in range(MTn):
                        dist = np.hypot(x0V[jj] - x0pV[kk], y0V[jj] - y0pV[kk])
                        if dist <= 2 * MT_RADIUS:
                            return True
                return False

            x0V, y0V, x0pV, y0pV = _tubule_centers(t_mt[0], t_mt[1], b)
            overlap = _any_overlap(x0V, y0V, x0pV, y0pV)

            if overlap and CORRECT_OVERLAP:
                print(
                    "It is not possible to reach the specified LK without overlap "
                    "(showing closest solution)... Calculating (wait)..."
                )
                for jjj in range(idx + 1, len(bb)):
                    x0V, y0V, x0pV, y0pV = _tubule_centers(t_mt[0], t_mt[1], bb[jjj])
                    if not _any_overlap(x0V, y0V, x0pV, y0pV):
                        overlap = False
                        idx = jjj
                        b = bb[idx]
                        b_deg = _round2xdigit(np.degrees(b))
                        LKm = _round2xdigit(LKmV[idx])
                        break

            if overlap and CORRECT_OVERLAP:
                print("No better solution was found... Showing with overlap.")
            elif overlap and not CORRECT_OVERLAP:
                print("Showing with overlap!")

    # --- Derive the B/C fusion cut angle from actual circle-circle overlap
    # geometry (consecutive A-B / B-C center spacing is always
    # 0.5*(L - 2*MT_RADIUS) given the fixed 0/0.5/1.0 offset fractions).
    consecutive_spacing = 0.5 * (L - 2 * MT_RADIUS)
    if FUSION_ARC_DEG is None:
        if 0 < consecutive_spacing < 2 * MT_RADIUS:
            fusion_arc_deg = 2 * np.degrees(np.arccos(consecutive_spacing / (2 * MT_RADIUS)))
        else:
            fusion_arc_deg = 0.0  # tubules don't actually touch at these L/MT_RADIUS -- draw complete rings
            print(
                "Note: L and MT_RADIUS put consecutive tubule centers "
                f"{consecutive_spacing:.1f} apart, >= 2*MT_RADIUS ({2*MT_RADIUS}) -- "
                "B/C would not actually be fused to their neighbor at these values; "
                "drawing them as complete rings."
            )
    else:
        fusion_arc_deg = FUSION_ARC_DEG

    # --- Plotting ---
    if show_result:
        owns_fig = ax is None
        if owns_fig:
            fig1, ax1 = plt.subplots(figsize=(6, 6))
        else:
            ax1 = ax

        # CID + cartwheel hub, drawn once.
        ax1.add_patch(Circle((0, 0), HUB_RADIUS, facecolor=COLOR_CARTWHEEL, edgecolor="black", linewidth=0.5, zorder=1))
        ax1.add_patch(Circle((0, 0), CID_RADIUS, facecolor=COLOR_CID, edgecolor="black", linewidth=0.5, zorder=2))
        ax1.add_patch(Circle((0, 0), CID_RADIUS * 0.15, facecolor="black", zorder=3))

        # SAS-6 dimers: SYM_CW discrete subunits around the hub ring.
        dimer_r = 0.3 * (2 * np.pi * 0.72 * HUB_RADIUS / SYM_CW)
        for m in range(SYM_CW):
            dx, dy = 0.72 * HUB_RADIUS * np.cos(t_cw[m]), 0.72 * HUB_RADIUS * np.sin(t_cw[m])
            ax1.add_patch(Circle((dx, dy), dimer_r, facecolor=COLOR_DIMER, edgecolor="black", linewidth=0.4, zorder=2))

        # Spokes + pinheads: SYM_CW of them, each collinear (hub edge ->
        # pinhead tip, where the A-tubule attaches). The pinhead's own
        # MIDPOINT (marked separately) is where the triplet base branches
        # off -- a distinct point from the tip.
        pinhead_tip = np.zeros((SYM_CW, 2))
        pinhead_mid = np.zeros((SYM_CW, 2))
        for m in range(SYM_CW):
            angle = t_cw[m]
            ux, uy = np.cos(angle), np.sin(angle)
            p0 = (HUB_RADIUS * ux, HUB_RADIUS * uy)
            p1 = ((HUB_RADIUS + SPOKE_LENGTH) * ux, (HUB_RADIUS + SPOKE_LENGTH) * uy)
            p2 = (Rc * ux, Rc * uy)
            ax1.plot([p0[0], p1[0]], [p0[1], p1[1]], color=COLOR_SPOKE, linewidth=4, solid_capstyle="round", zorder=1)
            ax1.plot([p1[0], p2[0]], [p1[1], p2[1]], color=COLOR_PINHEAD, linewidth=7, solid_capstyle="round", zorder=1)
            pinhead_tip[m] = p2
            pinhead_mid[m] = (Rc_mid * ux, Rc_mid * uy)

        # Nearest-neighbor map: which spoke feeds each triplet's base.
        spoke_for_triplet = _nearest_neighbor_map(t_mt, t_cw)

        # Pre-pass: every triplet's tilt direction and tubule centers, so
        # the A-C linker (computed below) can reference a NEIGHBORING
        # triplet's A-tubule center before that triplet's own turn in the
        # per-triplet drawing loop.
        dir_rad_of = np.zeros(SYM_MT)
        centers_of = np.zeros((SYM_MT, len(offsets), 2))  # [triplet, offset-index, xy]
        for k in range(SYM_MT):
            theta = t_mt[k]
            dir_rad_of[k] = np.pi / 2 - theta + b
            ux, uy = np.cos(theta), np.sin(theta)
            for oi, frac in enumerate(offsets):
                span = frac * (L - 2 * MT_RADIUS)
                centers_of[k, oi] = (
                    (Rc + MT_RADIUS) * ux - span * np.cos(dir_rad_of[k]),
                    (Rc + MT_RADIUS) * uy + span * np.sin(dir_rad_of[k]),
                )

        # MT triplets (attached directly to their spoke's pinhead tip) +
        # triplet bases (branching from the pinhead's MIDPOINT).
        base_tip = np.zeros((SYM_MT, 2))
        for k in range(SYM_MT):
            dir_rad = dir_rad_of[k]

            src = pinhead_mid[spoke_for_triplet[k]]
            tip_x = src[0] - BASE_LENGTH * np.cos(dir_rad)
            tip_y = src[1] + BASE_LENGTH * np.sin(dir_rad)
            base_tip[k] = (tip_x, tip_y)
            ax1.plot([src[0], tip_x], [src[1], tip_y], color=COLOR_BASE, linewidth=4, solid_capstyle="round", zorder=1)

            for oi, frac in enumerate(offsets):
                cx, cy = centers_of[k, oi]
                if frac == 0.0:
                    # A tubule: complete ring.
                    ax1.add_patch(Circle((cx, cy), MT_RADIUS, facecolor=COLOR_MT, edgecolor="black", linewidth=0.5, zorder=2))
                    ax1.add_patch(Circle((cx, cy), MT_RADIUS * MT_LUMEN_RATIO, facecolor=COLOR_LUMEN, zorder=3))
                elif fusion_arc_deg <= 0:
                    # Not actually touching at these L/MT_RADIUS -- complete ring.
                    ax1.add_patch(Circle((cx, cy), MT_RADIUS, facecolor=COLOR_MT, edgecolor="black", linewidth=0.5, zorder=2))
                    ax1.add_patch(Circle((cx, cy), MT_RADIUS * MT_LUMEN_RATIO, facecolor=COLOR_LUMEN, zorder=3))
                else:
                    # B/C tubule: incomplete ring. The A->B (and B->C)
                    # center displacement vector is (-cos(dir_rad),
                    # sin(dir_rad)) -- note the sign asymmetry inherited
                    # from centriole.py's own anchor formulas, which means
                    # its true bearing is (180 - dir_rad_deg), NOT
                    # dir_rad_deg. The fusion cut faces back toward the
                    # lower-offset (already-placed) neighbor, i.e. at the
                    # opposite bearing: (180 - dir_rad_deg) + 180 = -dir_rad_deg.
                    fusion_bearing = -np.degrees(dir_rad)
                    theta1 = fusion_bearing + fusion_arc_deg / 2
                    theta2 = theta1 + (360 - fusion_arc_deg)
                    ax1.add_patch(
                        Wedge(
                            (cx, cy),
                            MT_RADIUS,
                            theta1,
                            theta2,
                            width=MT_RADIUS * (1 - MT_LUMEN_RATIO),
                            facecolor=COLOR_MT,
                            edgecolor="black",
                            linewidth=0.5,
                            zorder=2,
                        )
                    )

        # A-C linkers: connect this triplet's C-tubule (its "far" tubule --
        # B for a doublet, itself for a singlet) to the NEXT triplet's
        # A-tubule, anchored on each tubule's facing perimeter point (the
        # point on the C/A tubule surface closest to the other), not the
        # tubule center and not the (separate) triplet-base tip.
        for k in range(SYM_MT):
            k2 = (k + 1) % SYM_MT
            C_k = centers_of[k, -1]
            A_k2 = centers_of[k2, 0]
            vec = A_k2 - C_k
            dist = np.hypot(vec[0], vec[1])
            direction = vec / dist if dist > 1e-9 else np.array([1.0, 0.0])
            anchor_far = C_k + MT_RADIUS * direction
            anchor_near = A_k2 - MT_RADIUS * direction
            ax1.plot(
                [anchor_far[0], anchor_near[0]],
                [anchor_far[1], anchor_near[1]],
                color=COLOR_LINKER,
                linewidth=3,
                solid_capstyle="round",
                zorder=1,
            )

        lim = 1.15 * (Rc + MT_RADIUS + L + BASE_LENGTH)
        ax1.set_xlim(-lim, lim)
        ax1.set_ylim(-lim, lim)
        ax1.set_aspect("equal")
        ax1.set_box_aspect(1)
        for spine in ax1.spines.values():
            spine.set_visible(True)
        ax1.set_title(
            f"SYM_MT: {SYM_MT}  SYM_CW: {SYM_CW}  LK: {_num2str(LK)}  Rc: {_num2str(Rc)}"
            f"  L: {_num2str(L)}  MTr: {_num2str(MT_RADIUS)}",
            fontsize=9,
        )

        if owns_fig:
            base_name = (
                f"v2_SYMMT{SYM_MT}_SYMCW{SYM_CW}_LK{_num2str(LK)}_Rc{_num2str(Rc)}"
                f"_L{_num2str(L)}_MTr{_num2str(MT_RADIUS)}_MTn{MTn}"
            )
            fig1.savefig(f"{base_name}.pdf", bbox_inches="tight")
            fig1.savefig(f"{base_name}.png", bbox_inches="tight")
            plt.close(fig1)

    print(f"Obtained LK: {LKm} with an angle beta = {b_deg + 90} degrees.")

    result = CentrioleResult(success=success, overlap=overlap, b_deg=b_deg)
    if not return_details:
        return result

    details = {
        "LKm": LKm,
        "LKmV": LKmV,
        "bb_deg": np.degrees(bb),
        "t_mt_deg": np.degrees(t_mt),
        "t_cw_deg": np.degrees(t_cw),
        "Rc": Rc,
        "fusion_arc_deg": fusion_arc_deg,
    }
    return result, details


if __name__ == "__main__":
    import matplotlib

    matplotlib.use("Agg")
    centriole_v2()
