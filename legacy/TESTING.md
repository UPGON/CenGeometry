# Manual test log for `centriole.py`

No automated test suite exists for this project. This log records manual
smoke tests run against `centriole.py`, with the exact parameters and
observed results, so the coverage is reproducible.

Environment: Python 3.11.9, numpy 2.4.6, matplotlib 3.11.1, `Agg` backend.

## Basic configurations

| Case | Params | success | overlap | b_deg |
|---|---|---|---|---|
| Default (singlet) | `centriole()` | True | False | 0.0 |
| Doublet | `MTn=2` | True | False | 81.77 |
| Triplet | `MTn=3` | True | False | 51.11 |
| Non-perfect-square SYM, triplet | `SYM=13, MTn=3` | True | False | 52.35 |

All ran without exceptions and produced `LKm` within the 0.15 tolerance of
the requested `LK`.

## Overlap paths

| Case | Params | success | overlap | Notes |
|---|---|---|---|---|
| Hard-fail overlap | `SYM=4, r=1, MT_RADIUS=15, MTn=1, LK≈28.55` | False | False | `LK` unreachable at these dimensions, so overlap is never evaluated (matches the documented `success=False` ⇒ overlap-not-checked gotcha) |
| Soft overlap, uncorrected | `r=10, MTn=2, LK=8, CORRECT_OVERLAP=False` | True | **True** | Prints "Showing with overlap!"; `b_deg=41.14` |
| Soft overlap, corrected | `r=10, MTn=2, LK=8, CORRECT_OVERLAP=True` | True | **False** | Forward scan found `b_deg=44.38` (LKm relaxed to 9.11 to escape overlap) — confirms the forward-only greedy correction search works |
| Unreachable LK | `LK=1000` | False | False | Prints "It is not possible to reach the desired A-C LINKER." |
| Invalid MTn | `MTn=4` | — | — | Raises `ValueError` as designed, instead of MATLAB's raw index-out-of-bounds crash |

To reliably trigger the *hard-fail* overlap path (`d2 < 2*MT_RADIUS`,
independent of `b`), reduce `r` until adjacent A-tubule centers are closer
than `2*MT_RADIUS` for the given `SYM`/`CW_RADIUS`. To trigger the
*correctable/soft* overlap path, keep `r` above that hard-fail threshold but
pick an `LK` that forces a low tilt angle `b`, where the offset tubules of
neighboring blades (e.g. doublet's "B" position) come within `2*MT_RADIUS`
of each other.

## Visual checks

Rendered and visually inspected the saved `.png` for:
- Default singlet (`SYM=9`) — 9 single circles around the hub, connected by
  linker lines, spokes converge at the origin. Correct.
- Triplet (`SYM=9, MTn=3`) — 9 three-circle blades (A/B/C), correct relative
  spacing and linker connections. Correct.

**Note:** the blade color palette is sampled from a grayscale colormap
including pure white (`cm.gray(np.linspace(0, 1, n_colors))`), a direct
carry-over from the original MATLAB `colormap(gray)` call. On the white
plot background, any blade assigned the white end of the palette renders
as invisible. This is not a porting bug — it reproduces the original
script's behavior — but it can make the figure look like it's missing
blades at a glance.

## Resource management

Confirmed `plt.close(fig1)` / `plt.close(fig2)` prevent figure
accumulation across repeated calls (checked via
`len(plt.get_fignums())` before/after 5 successive calls — returns to 0
each time).

## Regression found via user-provided reference image (2026-07-28)

The user supplied a screenshot of an actual `centriole.m` run
(`SYM=9, LK=25, CW_RADIUS=11.0294, r=77, MTl=50, MTr=10`, triplet) to
compare against. Two real bugs were found and fixed as a result:

1. **Title/filename overflow.** The title and output filename embedded
   `CW_RADIUS` (and other params) using Python's raw float repr (e.g.
   `11.029437556268347`), producing far longer strings than MATLAB's
   `num2str` (~4-5 significant digits) ever would. This overflowed the
   figure title off the canvas. Fixed by adding a `_num2str()` helper
   (`%.4g` formatting) used for all title/filename numeric fields, plus
   `bbox_inches="tight"` on `savefig` as a safeguard.
2. **Missing patch outlines (the actual "geometry looks wrong" report).**
   MATLAB's `patch()` draws a black edge around every filled shape by
   default. The Python port's `ax1.fill(...)` calls didn't set an edge
   color, so any blade whose color landed at the white end of the
   grayscale palette rendered as fully invisible (no fill, no border)
   instead of an outlined white circle. This made renders look like
   they were missing 1/3 of the blades. Fixed by adding
   `edgecolor="black", linewidth=0.5` to every `ax1.fill(...)` call
   (blade tubule circles and the hub disc).

Verification: reproducing the user's exact parameters
(`SYM=9, LK=25, CW_RADIUS=11.0294, r=77, L=50, MTn=3, MT_RADIUS=10,
GAMMA=60`) produced `angle: 98.98` — numerically identical to the
reference screenshot's title — confirming the underlying geometry/search
math was already correct and only the rendering was at fault. The
corrected PDF render now visually matches the reference: same 9-blade
triplet pinwheel, same black/gray/white color cycling, same connector and
spoke layout, hub circle now visible with its outline.

## Not tested

- No MATLAB/Octave installation was available in this environment, so
  outputs were verified by formula-by-formula comparison against
  `centriole.m` and by manual/geometric reasoning rather than direct
  numeric cross-checking against a live MATLAB run.
- `round2xdigit`'s exact original precision remains unverified (see
  `README.md`); `round(x, 2)` is used as a documented stand-in.
