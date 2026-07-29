"""Python port of centriole.m.

Models a 2D cross-section of a centriole: a central "cartwheel" hub of
radius CW_RADIUS at the origin, SYM-fold rotationally symmetric spokes
("SAS-6 coiled-coils", length r) extending to SYM microtubule blades
(singlet/doublet/triplet, MTn=1/2/3, tubule radius MT_RADIUS), each blade
of length L (end of A to end of C) connected circumferentially to its
neighbor by an "A-C linker" of length LK anchored at angle GAMMA (degrees)
from the tubule's local radial direction.

The function searches over a candidate tilt angle b in [0, 90] degrees for
the value that makes the modeled linker distance match the target LK, then
optionally checks/corrects for MT-MT circle overlap between adjacent
blades, then plots the result.

This is a faithful port of the original MATLAB implementation, including
several of its original quirks (documented inline below). A few small,
clearly-flagged improvements were made where the spec explicitly allowed
it (vectorized angle sweep, closed-form overlap test instead of symbolic
solve, MTn > 3 raises a clear error, output filename includes MTn/GAMMA).
"""

from __future__ import annotations

from typing import NamedTuple, Optional

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


class CentrioleResult(NamedTuple):
    """Result of a centriole() call: (success, overlap, b_deg)."""

    success: bool
    overlap: bool
    b_deg: float


def _round2xdigit(x: float) -> float:
    """Python stand-in for the original MATLAB `round2xdigit`.

    NOTE: `round2xdigit` is called in the original .m file but is not
    defined anywhere in the codebase -- it appears to be an external
    (MATLAB File Exchange) helper whose exact rounding precision could not
    be verified. `round(x, 2)` (2 decimal places) is used here as a
    reasonable equivalent; this is an assumption, not a verified port.
    """

    return round(x, 2)


def _num2str(x: float) -> str:
    """Python stand-in for MATLAB's `num2str` default formatting.

    MATLAB's `num2str(x)` prints ~4-5 significant digits by default. The
    original title/filename strings were built with `num2str`, so e.g.
    `CW_RADIUS = 11.029437556268347` prints as `"11.0294"`, not the full
    double. Using Python's raw `f"{x}"` instead (as an earlier version of
    this port did) produces far longer strings than the original ever
    did, which visibly overflows the plot title. `%.4g`-style formatting
    approximates `num2str`'s default closely enough for this script's
    purposes (SYM/LK/CW/r/L/MT_RADIUS/GAMMA are all small-to-moderate
    magnitude values).
    """

    return f"{x:.4g}"


def centriole(
    SYM: int = 9,
    LK: float = 28.5,
    CW_RADIUS: Optional[float] = None,
    r: float = 27,
    L: Optional[float] = None,
    MTn: int = 1,
    MT_RADIUS: float = 10,
    GAMMA: float = 60,
    CORRECT_OVERLAP: bool = False,
    show_result: bool = True,
    ax: "Optional[plt.Axes]" = None,
    return_details: bool = False,
):
    """Model and (optionally) plot a centriole cross-section.

    Parameters mirror the original MATLAB signature. ``CW_RADIUS`` and
    ``L`` replicate the MATLAB ``nargin < 1`` default block: when left as
    ``None`` they are computed from the other parameters exactly as the
    original script did (``CW_RADIUS = SYM*7.7/(2*pi)``, and
    ``L = 20 * [1, 1.75, 2.65][MTn-1]``).

    Two additions beyond the original MATLAB script, both opt-in and
    backward compatible (the default call signature/return value is
    unchanged):

    ax : matplotlib.axes.Axes, optional
        Draw the geometry onto a caller-supplied Axes instead of creating
        (and saving) a new figure. Useful for building side-by-side
        comparison panels (e.g. wild-type vs. a perturbed parameter set)
        with :func:`perturbation.plot_grid`. When ``ax`` is given, the
        function does NOT save a PDF/PNG and does NOT create the second
        "LK vs angle" diagnostic figure -- the caller owns that figure's
        lifecycle. When ``ax`` is ``None`` (default), behavior is
        unchanged from the original: a new figure is created, drawn, and
        saved, plus the diagnostic figure.
    return_details : bool, optional
        If True, also return a dict with the full angle-sweep result
        (``LKm``, ``LKmV``, ``bb_deg``) and the blade angle positions
        (``t_deg``). Useful for quantifying how a perturbation shifts the
        *achievable* LK curve, not just the single best-fit point.

    Returns
    -------
    CentrioleResult
        A (success, overlap, b_deg) named tuple. If ``return_details`` is
        True, returns ``(CentrioleResult, details_dict)`` instead.
    """

    if MTn not in (1, 2, 3):
        # The original MATLAB code would silently index out of bounds
        # (L_ratio(MTn) / the plotting if/elif chain) for MTn > 3 or < 1.
        # Raising a clear error here is a deliberate, explicitly-allowed
        # improvement over that raw crash.
        raise ValueError(f"MTn must be 1 (singlet), 2 (doublet), or 3 (triplet); got {MTn}")

    if CW_RADIUS is None:
        CW_RADIUS = SYM * 7.7 / (2 * np.pi)  # central hub radius
    if L is None:
        L_single = 20  # triplet length (from the end of A to the end of C)
        L_ratio = [1, 1.75, 2.65]
        L = L_single * L_ratio[MTn - 1]  # MTn is 1-indexed in the original

    print()
    print(f"-------------- desired LK: {LK}")

    if show_result:
        # MT colors: sample ceil(sqrt(SYM)) evenly-spaced colors from a
        # grayscale colormap. The plotting loop below cycles through these
        # colors, but the cycle *resets* every floor(sqrt(SYM)) blades (see
        # the `j > sqrt(SYM)` comparison there, taken verbatim from the
        # MATLAB `if j > sqrt(SYM)`) -- for a non-perfect-square SYM this
        # means the last sampled color is allocated but never actually
        # used. That mismatch between "colors sampled" and "colors used"
        # is an original-MATLAB quirk, faithfully preserved here.
        n_colors = int(np.ceil(np.sqrt(SYM)))
        col = cm.gray(np.linspace(0, 1, n_colors))
        owns_fig = ax is None
        if owns_fig:
            fig1, ax1 = plt.subplots()
        else:
            ax1 = ax

    t = np.linspace(0, 2 * np.pi, SYM + 1)  # SYM+1 points; t[0]=0, t[-1]=2*pi
    tt = np.linspace(0, 2 * np.pi, 100)  # for drawing circles
    g = np.radians(GAMMA)

    # NOTE: `a = 2*pi/SYM` and `R = (r+CW_RADIUS-MT_RADIUS)/cos(a/2)` are
    # computed in the original MATLAB source but never used again --
    # dead code, intentionally omitted here.

    bb = np.linspace(0, np.pi / 2, 362)

    # --- Calculate LK and angle b (with respect to tangent) ----------------
    # Vectorized over bb with numpy (the original used a MATLAB for-loop);
    # this is a pure elementwise-trig computation so vectorizing changes
    # nothing about the result, only the implementation.
    i = 0  # representative blade pair: "blade 1" (t[0]) and "blade 2" (t[1])
    if MTn > 1:  # doublet / triplet
        anchC_x = (
            (r + CW_RADIUS + MT_RADIUS) * np.cos(t[i])
            - (L - 2 * MT_RADIUS) * np.cos(np.pi / 2 - t[i] + bb)
            - MT_RADIUS * np.cos(np.pi / 2 - t[i] + bb - g)
        )
        anchC_y = (
            (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i])
            + (L - 2 * MT_RADIUS) * np.sin(np.pi / 2 - t[i] + bb)
            + MT_RADIUS * np.sin(np.pi / 2 - t[i] + bb - g)
        )
        anchA_x = (r + CW_RADIUS + MT_RADIUS) * np.cos(t[i + 1]) + MT_RADIUS * np.cos(
            np.pi / 2 - t[i + 1] + bb - g
        )
        anchA_y = (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i + 1]) - MT_RADIUS * np.sin(
            np.pi / 2 - t[i + 1] + bb - g
        )
    else:  # singlet
        # NOTE: for MTn == 1, anchC/anchA do not depend on bb at all, so
        # LKmV is CONSTANT across the entire sweep. This means the angle
        # search is effectively a no-op for singlets -- this is the
        # original model's actual behavior and is faithfully preserved
        # here, not "fixed".
        anchC_x = np.full_like(bb, (r + CW_RADIUS + MT_RADIUS) * np.cos(t[i]) - MT_RADIUS * np.cos(np.pi / 2 - t[i] - g))
        anchC_y = np.full_like(bb, (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i]) + MT_RADIUS * np.sin(np.pi / 2 - t[i] - g))
        anchA_x = np.full_like(bb, (r + CW_RADIUS + MT_RADIUS) * np.cos(t[i + 1]) + MT_RADIUS * np.cos(np.pi / 2 - t[i + 1] - g))
        anchA_y = np.full_like(bb, (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i + 1]) - MT_RADIUS * np.sin(np.pi / 2 - t[i + 1] - g))

    LKmV = np.sqrt((anchC_x - anchA_x) ** 2 + (anchC_y - anchA_y) ** 2)

    # --- Select the best-matching angle index -------------------------------
    idx_candidates = np.where(np.abs(LKmV - LK) < 0.15)[0]
    if idx_candidates.size > 0:
        gaps = np.where(np.diff(idx_candidates) > 1)[0]
        if gaps.size > 0:
            # keep only the last contiguous run (discard earlier, lower-angle runs)
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

    # --- Overlap check -------------------------------------------------------
    # IMPORTANT GOTCHA (preserved from the original): overlap is only
    # evaluated when success is True. If success is False, `overlap`
    # simply stays False -- this does NOT mean "no overlap was found", it
    # means overlap was never checked at all.
    overlap = False
    if success:
        r1 = np.array([(r + CW_RADIUS + MT_RADIUS) * np.cos(t[i]), (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i])])
        r2 = np.array(
            [(r + CW_RADIUS + MT_RADIUS) * np.cos(t[i + 1]), (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i + 1])]
        )
        d2 = np.linalg.norm(r1 - r2)

        if d2 < 2 * MT_RADIUS:
            print("It is not possible to avoid MT overlap.. Increase SAS-6 length (r)")
            overlap = True
        else:

            def _tubule_centers(angle0: float, angle1: float, beta: float):
                """Tubule CENTER positions (A, C, B) for a blade pair.

                Distinct from the anchor (perimeter) points used in the LK
                search above. Returns (x0V, y0V) for the blade at angle0
                and (x0pV, y0pV) for the blade at angle1, each a 3-vector
                indexed [A, C(full offset), B(half offset)].
                """

                x0V = np.array(
                    [
                        (r + CW_RADIUS + MT_RADIUS) * np.cos(angle0),
                        (r + CW_RADIUS + MT_RADIUS) * np.cos(angle0) - (L - 2 * MT_RADIUS) * np.cos(np.pi / 2 - angle0 + beta),
                        (r + CW_RADIUS + MT_RADIUS) * np.cos(angle0) - 0.5 * (L - 2 * MT_RADIUS) * np.cos(np.pi / 2 - angle0 + beta),
                    ]
                )
                y0V = np.array(
                    [
                        (r + CW_RADIUS + MT_RADIUS) * np.sin(angle0),
                        (r + CW_RADIUS + MT_RADIUS) * np.sin(angle0) + (L - 2 * MT_RADIUS) * np.sin(np.pi / 2 - angle0 + beta),
                        (r + CW_RADIUS + MT_RADIUS) * np.sin(angle0) + 0.5 * (L - 2 * MT_RADIUS) * np.sin(np.pi / 2 - angle0 + beta),
                    ]
                )
                x0pV = np.array(
                    [
                        (r + CW_RADIUS + MT_RADIUS) * np.cos(angle1),
                        (r + CW_RADIUS + MT_RADIUS) * np.cos(angle1) - (L - 2 * MT_RADIUS) * np.cos(np.pi / 2 - angle1 + beta),
                        (r + CW_RADIUS + MT_RADIUS) * np.cos(angle1) - 0.5 * (L - 2 * MT_RADIUS) * np.cos(np.pi / 2 - angle1 + beta),
                    ]
                )
                y0pV = np.array(
                    [
                        (r + CW_RADIUS + MT_RADIUS) * np.sin(angle1),
                        (r + CW_RADIUS + MT_RADIUS) * np.sin(angle1) + (L - 2 * MT_RADIUS) * np.sin(np.pi / 2 - angle1 + beta),
                        (r + CW_RADIUS + MT_RADIUS) * np.sin(angle1) + 0.5 * (L - 2 * MT_RADIUS) * np.sin(np.pi / 2 - angle1 + beta),
                    ]
                )
                return x0V, y0V, x0pV, y0pV

            def _any_overlap(x0V, y0V, x0pV, y0pV) -> bool:
                """Closed-form circle-circle intersection test.

                Replaces the original MATLAB `solve()` symbolic approach:
                two circles of radius MT_RADIUS intersect iff the distance
                between their centers is <= 2*MT_RADIUS.
                """

                for j in range(MTn):
                    for k in range(MTn):
                        dist = np.hypot(x0V[j] - x0pV[k], y0V[j] - y0pV[k])
                        if dist <= 2 * MT_RADIUS:
                            return True
                return False

            x0V, y0V, x0pV, y0pV = _tubule_centers(t[0], t[1], b)
            overlap = _any_overlap(x0V, y0V, x0pV, y0pV)

            if overlap and CORRECT_OVERLAP:
                print(
                    "It is not possible to reach the specified LK without overlap "
                    "(showing closest solution)... Calculating (wait)..."
                )
                # Forward-only greedy search: scan bb[idx+1:], never backward.
                # If no angle without overlap is found, `overlap` simply
                # stays True and idx/b/b_deg/LKm keep their pre-search values.
                for j in range(idx + 1, len(bb)):
                    x0V, y0V, x0pV, y0pV = _tubule_centers(t[0], t[1], bb[j])
                    if not _any_overlap(x0V, y0V, x0pV, y0pV):
                        overlap = False
                        idx = j
                        b = bb[idx]
                        b_deg = _round2xdigit(np.degrees(b))
                        LKm = _round2xdigit(LKmV[idx])
                        break

            if overlap and CORRECT_OVERLAP:
                print("No better solution was found... Showing with overlap.")
            elif overlap and not CORRECT_OVERLAP:
                print("Showing with overlap!")

    # --- Plotting --------------------------------------------------------------
    if show_result:
        # 1-indexed MT color counter (kept 1-indexed, matching the MATLAB
        # loop exactly): using `col[j - 1]` and resetting on `j > sqrt(SYM)`
        # (not `n_colors = ceil(sqrt(SYM))`) reproduces the original's
        # color-cycle length precisely, including for non-perfect-square
        # SYM where the two thresholds differ.
        j = 1
        for i in range(SYM):
            ax1.plot([0, (r + CW_RADIUS) * np.cos(t[i])], [0, (r + CW_RADIUS) * np.sin(t[i])], color="black")

            if MTn > 2:  # triplet
                ax1.plot(
                    [
                        (r + CW_RADIUS + MT_RADIUS) * np.cos(t[i]) - (L - 2 * MT_RADIUS) * np.cos(np.pi / 2 - t[i] + b) - MT_RADIUS * np.cos(np.pi / 2 - t[i] + b - g),
                        (r + CW_RADIUS + MT_RADIUS) * np.cos(t[i + 1]) + MT_RADIUS * np.cos(np.pi / 2 - t[i + 1] + b - g),
                    ],
                    [
                        (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i]) + (L - 2 * MT_RADIUS) * np.sin(np.pi / 2 - t[i] + b) + MT_RADIUS * np.sin(np.pi / 2 - t[i] + b - g),
                        (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i + 1]) - MT_RADIUS * np.sin(np.pi / 2 - t[i + 1] + b - g),
                    ],
                    color="black",
                )
                # C (full offset)
                ax1.fill(
                    (r + CW_RADIUS + MT_RADIUS) * np.cos(t[i]) - (L - 2 * MT_RADIUS) * np.cos(np.pi / 2 - t[i] + b) + MT_RADIUS * np.cos(tt),
                    (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i]) + (L - 2 * MT_RADIUS) * np.sin(np.pi / 2 - t[i] + b) + MT_RADIUS * np.sin(tt),
                    color=col[j - 1],
                    edgecolor="black",
                    linewidth=0.5,
                )
                # B (half offset)
                ax1.fill(
                    (r + CW_RADIUS + MT_RADIUS) * np.cos(t[i]) - 0.5 * (L - 2 * MT_RADIUS) * np.cos(np.pi / 2 - t[i] + b) + MT_RADIUS * np.cos(tt),
                    (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i]) + 0.5 * (L - 2 * MT_RADIUS) * np.sin(np.pi / 2 - t[i] + b) + MT_RADIUS * np.sin(tt),
                    color=col[j - 1],
                    edgecolor="black",
                    linewidth=0.5,
                )
                # A (zero offset)
                ax1.fill(
                    (r + CW_RADIUS + MT_RADIUS) * np.cos(t[i]) + MT_RADIUS * np.cos(tt),
                    (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i]) + MT_RADIUS * np.sin(tt),
                    color=col[j - 1],
                    edgecolor="black",
                    linewidth=0.5,
                )
            elif MTn > 1:  # doublet
                ax1.plot(
                    [
                        (r + CW_RADIUS + MT_RADIUS) * np.cos(t[i]) - (L - 2 * MT_RADIUS) * np.cos(np.pi / 2 - t[i] + b) - MT_RADIUS * np.cos(np.pi / 2 - t[i] + b - g),
                        (r + CW_RADIUS + MT_RADIUS) * np.cos(t[i + 1]) + MT_RADIUS * np.cos(np.pi / 2 - t[i + 1] + b - g),
                    ],
                    [
                        (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i]) + (L - 2 * MT_RADIUS) * np.sin(np.pi / 2 - t[i] + b) + MT_RADIUS * np.sin(np.pi / 2 - t[i] + b - g),
                        (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i + 1]) - MT_RADIUS * np.sin(np.pi / 2 - t[i + 1] + b - g),
                    ],
                    color="black",
                )
                # B (full offset -- same geometry as triplet's "C" position;
                # relabeled in the original comments only, preserved as-is)
                ax1.fill(
                    (r + CW_RADIUS + MT_RADIUS) * np.cos(t[i]) - (L - 2 * MT_RADIUS) * np.cos(np.pi / 2 - t[i] + b) + MT_RADIUS * np.cos(tt),
                    (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i]) + (L - 2 * MT_RADIUS) * np.sin(np.pi / 2 - t[i] + b) + MT_RADIUS * np.sin(tt),
                    color=col[j - 1],
                    edgecolor="black",
                    linewidth=0.5,
                )
                # A (zero offset)
                ax1.fill(
                    (r + CW_RADIUS + MT_RADIUS) * np.cos(t[i]) + MT_RADIUS * np.cos(tt),
                    (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i]) + MT_RADIUS * np.sin(tt),
                    color=col[j - 1],
                    edgecolor="black",
                    linewidth=0.5,
                )
            elif MTn == 1:  # singlet
                ax1.plot(
                    [
                        (r + CW_RADIUS + MT_RADIUS) * np.cos(t[i]) - MT_RADIUS * np.cos(np.pi / 2 - t[i] + b - g),
                        (r + CW_RADIUS + MT_RADIUS) * np.cos(t[i + 1]) + MT_RADIUS * np.cos(np.pi / 2 - t[i + 1] + b - g),
                    ],
                    [
                        (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i]) + MT_RADIUS * np.sin(np.pi / 2 - t[i] + b - g),
                        (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i + 1]) - MT_RADIUS * np.sin(np.pi / 2 - t[i + 1] + b - g),
                    ],
                    color="black",
                )
                # A (zero offset)
                ax1.fill(
                    (r + CW_RADIUS + MT_RADIUS) * np.cos(t[i]) + MT_RADIUS * np.cos(tt),
                    (r + CW_RADIUS + MT_RADIUS) * np.sin(t[i]) + MT_RADIUS * np.sin(tt),
                    color=col[j - 1],
                    edgecolor="black",
                    linewidth=0.5,
                )

            j += 1  # MT color index
            if j > np.sqrt(SYM):
                j = 1

        # White cartwheel hub disk. NOTE: original uses
        # (CW_RADIUS*sin(tt), CW_RADIUS*cos(tt)) -- sin for x, cos for y,
        # reversed relative to the tubule patches (which use cos/sin).
        # Harmless for a full circle; preserved as in the original.
        ax1.fill(CW_RADIUS * np.sin(tt), CW_RADIUS * np.cos(tt), color="white", edgecolor="black", linewidth=0.5)

        lim = 1.1 * (r + CW_RADIUS + 6 * MT_RADIUS)
        ax1.set_xlim(-lim, lim)
        ax1.set_ylim(-lim, lim)
        ax1.set_aspect("equal")
        ax1.set_box_aspect(1)
        for spine in ax1.spines.values():
            spine.set_visible(True)
        ax1.set_title(
            f"SYM: {SYM}   LK: {_num2str(LK)}   CW: {_num2str(CW_RADIUS)}   "
            f"SAS-6: {_num2str(r)}  MTl: {_num2str(L)}   MTr: {_num2str(MT_RADIUS)}",
            fontsize=10,
        )

        # Saving a PDF/PNG and creating the second diagnostic figure only
        # make sense when this call owns its figure. When `ax` was
        # supplied by the caller (e.g. building a multi-panel comparison
        # grid), the caller owns the figure's lifecycle -- don't save
        # files or spawn an extra standalone diagnostic figure per panel.
        if owns_fig:
            # Output filename encodes the same parameters as the original,
            # PLUS MTn and GAMMA. The original filename omitted MTn/GAMMA/
            # CORRECT_OVERLAP entirely, so two runs differing only in those
            # parameters would silently overwrite each other's output file.
            # Including MTn and GAMMA here is a deliberate fix for that
            # collision; CORRECT_OVERLAP is left out since it does not affect
            # the final geometry when a solution is found. Values are
            # formatted with `_num2str` (not raw Python float repr) so e.g. a
            # computed CW_RADIUS doesn't turn into a 17-digit filename.
            base_name = (
                f"SYM{SYM}_LK{_num2str(LK)}_CW{_num2str(CW_RADIUS)}_SAS{_num2str(r)}"
                f"_MTl{_num2str(L)}_MTr{_num2str(MT_RADIUS)}_MTn{MTn}_GAMMA{_num2str(GAMMA)}"
            )
            # bbox_inches="tight" guards against title/label clipping (the bug
            # that made earlier renders show a cut-off/overflowing title).
            fig1.savefig(f"{base_name}.pdf", bbox_inches="tight")
            # .fig has no direct Python equivalent; PNG is saved as a reasonable substitute.
            fig1.savefig(f"{base_name}.png", bbox_inches="tight")

            # Second, diagnostic figure: LK(angle)
            fig2, ax2 = plt.subplots()
            ax2.plot(np.degrees(bb), LKmV)
            ax2.set_xlabel("Angle (degrees)")
            ax2.set_ylabel("LK")

            # Both figures have been saved (fig1) or are purely diagnostic
            # (fig2); close them so repeated calls (e.g. sweeping parameters
            # in a loop) don't accumulate open figures in memory.
            plt.close(fig1)
            plt.close(fig2)

    # --- Final report ------------------------------------------------------
    # NOTE: the printed angle is b_deg + 90, but the RETURNED b_deg is the
    # plain value (without +90). These are deliberately different values;
    # both are preserved distinctly, exactly as in the original.
    print(f"Obtained LK: {LKm} with an angle beta = {b_deg + 90} degrees.")

    result = CentrioleResult(success=success, overlap=overlap, b_deg=b_deg)
    if not return_details:
        return result

    details = {
        "LKm": LKm,
        "LKmV": LKmV,
        "bb_deg": np.degrees(bb),
        "t_deg": np.degrees(t),
    }
    return result, details


if __name__ == "__main__":
    # Force a non-interactive backend for CLI/headless runs (e.g. no
    # display available). Done here rather than at module import time so
    # importing this module (e.g. from a Jupyter notebook that already has
    # an interactive backend configured) doesn't silently override it.
    import matplotlib

    matplotlib.use("Agg")
    centriole()
