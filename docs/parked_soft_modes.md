# Soft-mode (normal-mode) analysis — parked notes

**Status: parked.** The code works and is committed, but the output did not
prove useful or interpretable in practice and is not surfaced in the app or
the standard report. Revisit only if there is a specific question it would
answer.

---

## What it does

`mode_analysis()` approximates the Hessian as `J^T J` at the relaxed
wild-type solution and eigendecomposes it. Low eigenvalues are directions
in which the structure barely resists deformation — motions it could
perform at biological energy scales. This is normal-mode / elastic-network
analysis applied to the linkage.

```python
from centriole_kinematic import Geometry, mode_analysis, plot_modes
ma = mode_analysis(Geometry(), n_modes=6)
print(ma.table)
plot_modes(ma, n=6, path="modes.png")
```

Each mode is labelled by its **angular wavenumber `m`** around the ring
(0 = every unit alike, 1 = ring shifts sideways, 2 = ellipse, ≥3 =
ripples), plus its overlap with global rotation and the diameter/tilt
change it produces. `plot_modes()` draws wild type in grey, the deformed
shape over it, and per-tubule arrows.

## What it found

| Mode | Stiffness | m | Character |
|---|---|---|---|
| 0 | 1.0× | 0 | global pinwheel — whole assembly rotates against the hub (rotation overlap 0.77) |
| 1–2 | 28× | 4 | ripples around the ring |
| 3–4 | 34× | 3 | ripples |
| 5 | 44× | 2 | elliptical distortion |

So: one nearly-free rotational sloppiness against the hub, then a stiff,
rumple-dominated spectrum, and **no soft mode that changes the diameter**.

One internal check did pass — modes come out in degenerate pairs for every
`m > 0` (a sine and a cosine of the same ring wave), as ring symmetry
requires. So the decomposition itself was behaving.

## Why it was parked

The results were not interpretable in biological terms. Concretely:

- The only genuinely soft mode is a trivial artefact of the model having no
  absolute angular reference — informative about the parameterisation, not
  the biology.
- Everything else is a high-wavenumber ripple with no obvious
  correspondence to anything observable by cryo-ET or mutation.
- Mode amplitudes are arbitrary (eigenvectors have no intrinsic scale), so
  "how much" a mode moves is set by an arbitrary display factor rather
  than by any energy scale. Without a real energy calibration — which
  would need measured bond energies, not the current relative strengths —
  the stiffness numbers cannot be converted into anything testable like a
  thermal amplitude.

That last point is the substantive blocker. **If this is revisited, the
prerequisite is an absolute energy scale for the bonds** (kT units), which
would turn eigenvalues into predicted fluctuation amplitudes comparable
against cryo-ET particle-to-particle variance. Without that, the analysis
can only ever rank modes relative to one another.

## Related parked work

`docs/blooming_handoff.md` — the iris/blooming question, which is where
this analysis was originally aimed.
