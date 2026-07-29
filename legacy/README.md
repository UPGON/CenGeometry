# Legacy code

Superseded by [`../centriole_kinematic.py`](../centriole_kinematic.py).
Kept for provenance — none of it is needed to run the current model, and
none of it is maintained.

| File | What it is |
|---|---|
| `centriole.m` | The original MATLAB script this project grew from. Draws a parametric centriole cross-section and searches for the blade tilt angle that matches a target A-C linker length. |
| `centriole.py` | A faithful Python port of `centriole.m`, preserving its behaviour (including several of its quirks) so results could be compared directly. |
| `centriole_v2.py` | First attempt at an anatomically explicit model — splits the single lumped `r` into CID, cartwheel, pinhead and triplet base. Still a parametric drawing rather than a mechanical model. |
| `perturbation.py` | Sweep and comparison helpers written against `centriole.py` / `centriole_v2.py`. |
| `perturbation_analysis.ipynb` | Notebook demonstrating those helpers. |
| `TESTING.md` | Manual test log for the MATLAB port. |

## Why they were replaced

These versions *draw* a shape from parameters. The current model instead
treats the centriole as rigid bodies joined by rotatable, strength-weighted
bonds and solves for the least-strained arrangement, which is what makes
perturbations answerable rather than merely drawable.

Two specific problems drove the rewrite:

1. **The cartwheel was primary.** Triplets hung off their nearest spoke, so
   a symmetry mismatch dragged them out of position. Biologically the
   triplet ring is the scaffold and the cartwheel should adapt to it.
2. **All dimensions were guessed.** Segment lengths were eyeballed
   proportions. They are now measured from a cryo-ET-derived schematic by
   [`../svg_calibration.py`](../svg_calibration.py).
