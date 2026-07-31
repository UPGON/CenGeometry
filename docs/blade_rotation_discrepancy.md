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

## Caveats

- The −2 register gives a slightly worse loop closure (0.203 nm against
  0.023 nm at wild-type register), so the model does not *prefer* it on
  energy grounds — it is selected here because it matches the data. Treat the
  register value as a hypothesis to test, not a result.
- **Do not tune anything to fit this measurement.** One number agreeing does
  not validate a mechanism, and fitting to it would reintroduce exactly the
  circularity documented in the critical-evaluation deck.
- `solve_chain` takes a fixed `reg=(pinhead, linkC)`; it does not search
  registers. A search exists only in the older spring solver, which cannot be
  trusted here because it opens its bonds 1–3 nm under this perturbation.
  Adding a proper register search to the chain solver is the obvious next
  piece of work.

## Reproducing

```python
from centriole_kinematic import Geometry, set_param
from centriole_chain import solve_chain

wt  = solve_chain(Geometry())
mut = solve_chain(set_param(Geometry(), "spoke_rod", 28.03), reg=(0.0, -2.0))
print(wt.outer_diameter, mut.outer_diameter)          # 255.0, 222.1
print(mut.triplet_tilt - wt.triplet_tilt)             # +2.67 deg
```

Under a second each. Use `solve_chain`, never `solve`, for this perturbation.

## Related

- `docs/parked_soft_modes.md` — abandoned normal-mode analysis
- `docs/blooming_handoff.md` — iris/blooming, blocked on the 2D limitation
- Grades are **assumed tolerances, not feasibility verdicts**: this mutant
  grades SEVERE at the linker contacts yet builds real centrioles.
