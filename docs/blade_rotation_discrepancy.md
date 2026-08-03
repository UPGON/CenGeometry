# Open question: what does the A-C linker do under radial compression?

**Status: reopened.** This began as a disagreement between model and
measurement and appeared to end with a specific, testable prediction. Two
defects in the register machinery have since been found and fixed, and the
fix removes the ground the prediction stood on. The experiment is still the
right one to do; what it would distinguish has changed. Read
[the 2026 revision](#2026-revision-the-register-search-was-measuring-the-wrong-thing)
before acting on anything above it.

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

> **Superseded.** The register this rests on turns out to need the linker's
> C-arm to reach its site through the C-tubule wall. See
> [the 2026 revision](#2026-revision-the-register-search-was-measuring-the-wrong-thing).

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

> **Superseded — the numbers in this section came from a search with two
> defects** (rest angles not re-referenced when the register shifted, and no
> check that a contact could be reached from outside its tubule). Kept for the
> record. The corrected table is in
> [the 2026 revision](#2026-revision-the-register-search-was-measuring-the-wrong-thing).

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

## 2026 revision: the register search was measuring the wrong thing

Everything above this heading was produced by a register search with two
defects. Both are now fixed, and the conclusions above do not survive them
unchanged.

### Defect 1 — shifted registers were charged a strain they never had

`CONTACT_REST` gives each strand's rest orientation *in the body frame of the
tubule it grips*. `pf_angle` puts protofilament `i + k` at
`PF1_ANCHOR - pitch*(i - 1 + k)`, so sliding a contact by one filament rotates
the binding site by `-pitch` in that same frame. The rest angle was not slid
with it. A strand making the identical interaction one protofilament round was
therefore recorded as having twisted by a whole `pitch` — **27.7°**, already
past the 20° HARD limit on both linker contacts — before the structure had
moved at all.

Every shifted register carried that penalty; the wild-type register did not.
The ranking was mostly bookkeeping. `contact_rest()` now re-references the rest
angle by `-pf_pitch * shift`.

### Defect 2 — nothing stopped a strand from lying in a tubule wall

Strand-versus-tubule penetration was measured but never entered any residual,
and the measurement itself trimmed 3 nm off *both* ends of *every* segment
against *every* tubule — more than half of an 11.6 nm linker arm, including
exactly the stretch next to the binding site where a bad register shows up.
Nothing, therefore, distinguished a contact reached head-on from outside from
one reached tangentially along the wall, or from inside it.

The search exploited this. Its winner for the 28 nm spoke, linker-A `+2`, put
the arm's contact 55° round the A-tubule and brought it in at 1.7° above
tangential — which is what the drawing was showing as the A-C linker clashing
into the A-tubule.

Three things now guard it. Strands have a thickness, so grazing counts as
overlap; the overlap enters the residual; and `contact_approach()` reports the
cosine of the angle at which each strand meets its protofilament, which is
negative when a strand would have to arrive through the wall. The cross-section
draws strands at their modelled width and colours a clashing one red, so what
is drawn matches what is checked.

### What the corrected search says

Wild type 254.96 nm, tilt −60.92°. Spoke shortened to 28.03 nm:

| Register (A / C) | Outer Ø change | Tilt change | Clearance | Approach | |
|---|---|---|---|---|---|
| +0 / +0 (wild type) | −21.16 | +12.85° | +1.06 | +0.79 | reachable |
| +1 / +1 — *model's pick* | −23.86 | +10.24° | +0.39 | +0.49 | reachable |
| +2 / +0 | **−35.11** | **+14.07°** | −0.03 | +0.30 | reachable, grazing |
| +1 / −2 — *the old pick* | **−34.04** | **+4.42°** | −0.88 | **−0.06** | **unreachable** |
| +0 / −2 | −32.62 | +2.63° | −1.29 | **−0.22** | **unreachable** |
| *measured* | *−34* | — | | | |

The registers that reproduced the measurement *by letting the blade keep its
tilt* — the `C−2` family, and the mechanism this document proposed — are the
ones that need the linker's C-arm to arrive through the C-tubule wall.

A register still reproduces −34 nm: `+2 / +0`, at −35.11 nm. But it gets there
by rotating the blade **more**, +14.1° against wild type's +12.9°, not less.
That is the opposite of the proposed mechanism.

### How firm is that?

Firm for `+0 / −2`, borderline for `+1 / −2`. Reachability is the
assumption-free half of the test — it asks only whether a strand arrives from
outside its tubule or from within, and is unchanged by the assumed strand
thickness or the weight on the steric term. But it does depend on the
protofilament angles, which carry about ±2° of calibration noise. Shifting
`PF1_ANCHOR` bodily over ±3°:

| anchor shift | +1 / −2 | +0 / −2 | +2 / +0 |
|---|---|---|---|
| −3° | +0.013 | −0.158 | +0.245 |
| −2° | −0.011 | −0.179 | +0.263 |
| 0° | −0.058 | −0.221 | +0.298 |
| +2° | −0.105 | −0.262 | +0.332 |
| +3° | −0.128 | −0.282 | +0.350 |

`+1 / −2` crosses zero at about −2.5°, i.e. **inside the measurement noise**.
So it is not ruled out; it is pushed to the edge of the possible, where the
C-arm would have to run essentially tangentially along the C-tubule. `+0 / −2`
stays negative throughout and is ruled out. The clearance column is the
*assumed* half — it counts a 1 nm strand half-width that was never measured —
so registers within a few tenths of zero (`+2 / +0` at −0.03) are marginal
either way.

The clashes are geometric, not numerical: raising the steric weight a
hundredfold moves `+2 / +0` from −0.032 nm to −0.014 nm. In the chain form both
ends of the linker are fixed by construction, so a straight arm between them
has no freedom left — clearance follows from where the contacts sit, and the
solver cannot relieve it.

### The revised prediction

The cryo-ET test is unchanged and now discriminates more sharply:

- **tilt essentially unchanged (~+3 to +4°) and the linker shifted ~2 filaments
  on the C-tubule** → the old mechanism, and the model's reachability check is
  wrong or its protofilament angles are off by more than 2°;
- **tilt rotated ~+14° with the linker shifted +2 on the A-tubule** → `+2 / +0`,
  which reproduces the diameter and is reachable;
- **tilt rotated ~+13° at wild-type register** → the original discrepancy
  returns and re-registering is not the explanation.

The middle case is new, and is the one the corrected model favours among the
registers that match the measured diameter.

Still worth doing first, and still cheap: confirm what the U-ExM signal is and
how the diameter is fitted from the ring. Peak-to-peak of a signal ring is not
the outer-wall envelope.

## Caveats

- **Do not tune anything to fit this measurement.** One number agreeing does
  not validate a mechanism, and fitting to it would reintroduce exactly the
  circularity documented in the critical-evaluation deck.
- The register search ranks by model cost, and **cost ranking is not
  evidence**. It is a shortlist built from joint bands and bond strengths that
  were reasoned, never measured. Read the whole `register_scan`, never just
  the returned winner.
- Registers that are *geometrically impossible* are a different matter and are
  excluded rather than ranked. That test is not a cost judgement, but it does
  inherit the ±2° noise on the protofilament angles — see the sensitivity
  table above before resting anything on a borderline case.
- The `feasible` flag combines an assumption-free half (reachability) with an
  assumed one (clearance, which counts a 1 nm strand half-width that was never
  measured). Prefer to rest a conclusion on the former.
- The search takes ~290 s for 25 registers, far slower than a single solve,
  because shifted registers converge less easily.
- Never use `solve()` for this perturbation: it opens its bond springs 1–3 nm,
  which visibly separates the spoke from the pinhead and corrupts the answer.
  `run_analysis.py` and the app both use `solve_chain` for this reason.

## Reproducing

```python
from centriole_kinematic import Geometry, set_param
from centriole_chain import solve_chain, best_registers

wt  = solve_chain(Geometry())
g   = set_param(Geometry(), "spoke_rod", 28.03)

mut = solve_chain(g, reg=(0, 1, -2))                  # the old pick
print(mut.outer_diameter, mut.triplet_tilt - wt.triplet_tilt)   # 220.9, +4.42
print(mut.feasible, mut.reachable, mut.worst_strand_clearance)  # False False -0.88
print(mut.report())                  # lists the approach cosine per contact

# search, then read the whole ranking -- feasible ones first, then by cost
m = solve_chain(g, register_shift=True)               # ~290 s
print(m.report())
for r in best_registers(m, 3): print(r)
```

A single solve is under a second; the search is ~290 s. Use `solve_chain`,
never `solve`, for this perturbation.

## Related

- `docs/parked_soft_modes.md` — abandoned normal-mode analysis
- `docs/blooming_handoff.md` — iris/blooming, blocked on the 2D limitation
- Grades are **assumed tolerances, not feasibility verdicts**: this mutant
  grades SEVERE at the linker contacts yet builds real centrioles.
