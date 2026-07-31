# Open question: what does the A-C linker do under radial compression?

**Status: largely resolved, needs experimental confirmation.** This began as a
disagreement between model and measurement and ended with a specific,
testable prediction. Not blocked.

---

## The observation

A mutant with the SAS-6 spoke shortened by **17 nm** (45.03 → 28.03 nm),
measured by U-ExM on the **outer wall**:

| | Wild type | Mutant | Change |
|---|---|---|---|
| Outer diameter | 250 nm | 216 nm | **−34 nm** |

## The initial disagreement

Solved with `solve_chain` at wild-type protofilament register:

| Ring | WT | Mutant | Change |
|---|---|---|---|
| A-tubule centres | 183.9 | 150.8 | −33.11 |
| Tubulin centroid | 205.2 | 178.5 | −26.62 |
| **Outer wall** | **255.0** | **233.8** | **−21.16** |
| *measured* | *250* | *216* | *−34* |

The mutant solve was otherwise clean: zero clashes, all connections closed,
a buildable geometry.

**The contradiction that made this worth chasing:** the measured *absolute*
wild-type diameter (250) matches the outer wall (255.0) and nothing else,
while the measured *change* (−34) matches the A-tubule ring (−33.11) and
nothing else. No single ring reproduced both, so it could not be waved away
as "U-ExM is really sampling the A-ring" — that would put wild type at
~184 nm, not 250.

## What was ruled out

The outer wall is buffered because the blade **rotates**: tilt goes −60.9° →
−48.1°, a change of +12.9°, swinging the C-tubule outward as the A-tubule
comes inward.

The obvious suspect was the `triplet` band, the most permissive in the model
at 20°/40°. **This was tested and is not the cause:**

| Triplet band | Outer change | Tilt change |
|---|---|---|
| 20°/40° (default) | −21.16 | +12.85° |
| 8°/16° | −21.15 | +12.84° |
| 3°/6° | −21.15 | +12.84° |
| 1°/2° (near-rigid) | −21.18 | +12.81° |

Tightening the band by a factor of twenty changes nothing. The rotation is
**kinematically forced**, not energetically chosen: the A-C linker is a
fixed-length tie between neighbouring triplets, so as the A-ring shrinks the
chord between adjacent A-tubules shortens, and the only way the linker can
still span is for the blade to swing outward. No band setting can prevent it,
because no alternative configuration closes.

## The resolution

If the rotation is forced by the linker at wild-type register, then the
linker must be doing something different in the mutant. Shifting the
C-tubule contact by whole protofilaments:

| linkC shift | Outer Ø change | Tilt change | Worst loop gap |
|---|---|---|---|
| −2 | **−32.82 nm** | **+2.67°** | 0.203 nm |
| −1 | −25.09 | +8.64° | 0.119 nm |
| 0 (wild type) | −21.16 | +12.85° | 0.023 nm |
| +1 | −21.30 | +12.80° | 0.014 nm |
| +2 | −27.04 | +3.45° | 0.007 nm |
| *measured* | *−34* | — | — |

**A −2 protofilament shift reproduces the measurement to within 1.2 nm**, and
does so by letting the blade very nearly keep its wild-type tilt (+2.7°
instead of +12.9°). That is a coherent mechanism rather than a fitted
coincidence: re-registering removes the kinematic forcing, the blade stops
rotating, and the outer wall then follows the A-ring inward.

Worth noting this is the mechanism the user proposed independently before it
was tested.

## The prediction to test

**In the mutant, the A-C linker should bind roughly 2 protofilaments away
from its wild-type position on the C-tubule, and the triplet tilt should be
essentially unchanged.**

Both are checkable by cryo-ET, and they are sharply distinguishable from the
wild-type-register alternative, which requires a +12.9° tilt change. Either
observation alone settles it:

- tilt unchanged + linker shifted → mechanism confirmed
- tilt rotated ~13° → mechanism wrong, and the outer-wall discrepancy returns

Also cheap and worth doing: confirm what the U-ExM signal is and how the
diameter is fitted from the ring. Peak-to-peak of a signal ring is not the
outer wall envelope, and the ~24 nm wall may be unresolved at typical
expansion factors.

## Register search (added after the manual scan)

`solve_chain(..., register_shift=True)` now searches both A-C linker contacts
over +-2 protofilaments (25 combinations, ~240 s) and keeps every candidate in
`register_scan`. Results for the 17 nm-shorter spoke:

| Register (linker-A / linker-C) | Model cost | Outer change | Tilt change |
|---|---|---|---|
| +2 / +0 — *model's pick* | 127 | −28.26 nm | +6.19° |
| +1 / +1 | 1,780 | −23.36 nm | +9.23° |
| +0 / +2 | 3,372 | −27.04 nm | +3.45° |
| +2 / +1 | 4,247 | −38.56 nm | −1.84° |
| **+1 / −2 — *matches the data*** | **> 15,000** | **−34.20 nm** | **+4.47°** |
| *measured* | — | *−34 nm* | — |

**The register that reproduces the measurement is one the model strongly
disfavours.** It reproduces −34 nm to within 0.2 nm and lets the blade nearly
keep its tilt, but ranks outside the top eight on cost — roughly 100x the
cheapest option.

That is a genuine tension, and it has only two readings:

1. **The cost function is wrong.** The joint bands and bond strengths are
   reasoned heuristics, never measured. If they misprice a linker-contact
   rotation, the ranking is meaningless and the +1/−2 register may be
   perfectly ordinary. This is quite likely, and it is the same weakness the
   critical-evaluation deck identifies as fundamental.
2. **Re-registering is not the mechanism.** The register reproduces the
   number, but if the model's energetics are even roughly right, the mutant
   would not adopt this configuration, and something else explains the
   measurement.

The experiment below distinguishes these without needing to fix the cost
function first, which is why it is worth doing before any further modelling.

## Caveats

- The −2 register gives a slightly worse loop closure (0.203 nm against
  0.023 nm at wild-type register), so the model does not *prefer* it on
  energy grounds — it is selected here because it matches the data. Treat the
  register value as a hypothesis to test, not a result.
- **Do not tune anything to fit this measurement.** One number agreeing does
  not validate a mechanism, and fitting to it would reintroduce exactly the
  circularity documented in the critical-evaluation deck.
- The register search ranks by model cost, and **cost ranking is not
  evidence** — it does not select the data-matching register. Always read the
  whole `register_scan`, never just the returned winner.
- The search takes ~240 s for 25 registers, far slower than a single solve,
  because shifted registers converge less easily.
- Never use `solve()` for this perturbation: it opens its bond springs 1–3 nm,
  which visibly separates the spoke from the pinhead and corrupts the answer.

## Reproducing

```python
from centriole_kinematic import Geometry, set_param
from centriole_chain import solve_chain

wt  = solve_chain(Geometry())
mut = solve_chain(set_param(Geometry(), "spoke_rod", 28.03), reg=(0, 1, -2))
print(wt.outer_diameter, mut.outer_diameter)          # 255.0, 220.8
print(mut.triplet_tilt - wt.triplet_tilt)             # +4.47 deg

# or search, then read the whole ranking rather than the winner
m = solve_chain(g, register_shift=True)
for r in m.register_scan: print(r)
```

Under a second each. Use `solve_chain`, never `solve`, for this perturbation.

## Related

- `docs/parked_soft_modes.md` — abandoned normal-mode analysis
- `docs/blooming_handoff.md` — iris/blooming, blocked on the 2D limitation
- Grades are **assumed tolerances, not feasibility verdicts**: this mutant
  grades SEVERE at the linker contacts yet builds real centrioles.
