# CenGeometry

**Georgios Hatzopoulos** · UPGON

A geometric and mechanical model of the centriole cross-section, built to
predict how structural perturbations — SAS-6 mutations that change
cartwheel symmetry, altered coiled-coil lengths, loss of the C-tubule —
change the architecture of the whole organelle.

**New here?** Double-click the launcher for your platform:

| Platform | File |
|---|---|
| **macOS** | `START HERE (Mac).command` |
| **Windows** | `START HERE (Windows).bat` |

It sets everything up the first time and opens the interactive app in your
browser — no terminal, no code. Use only the file for your own platform:
double-clicking the `.bat` on a Mac just opens it in a text editor.

**First run on macOS** — Gatekeeper blocks unsigned scripts, so the very
first launch needs **right-click → Open → Open** rather than a double-click.
After that, double-clicking works normally.

**If the Mac launcher won't run at all**, the executable permission has
probably been stripped (Dropbox and some file transfers do this). Restore it
once with:

```bash
chmod +x "START HERE (Mac).command"
```

For the command-line route instead, open
[`QUICKSTART.html`](QUICKSTART.html) in a web browser.

---

## What the model does

A centriole is a ring of microtubule triplets held around a central
cartwheel. This model represents one repeating unit as a set of **rigid
bodies joined by connections that can rotate**, then asks what happens to
the whole ring when you change something:

```
SAS-6 head (hub ring) → coiled-coil spoke → pinhead → microtubule triplet
                                    ↓                          ↓
                              triplet base ────→ A-C linker ──→ next triplet
```

Because the units must close into a ring, changing one part forces the
others to adapt. The model finds the least-strained arrangement that still
fits together and reports:

| Output | Meaning |
|---|---|
| **Joint rotation** | How far each connection turned from its wild-type angle, graded `OK` / `HARD` / `SEVERE` against per-joint limits (see below) |
| **Bond load** | Which connection carries the most strain, i.e. which would rupture first |
| **Buckling** | Which segments had to bow because space became too tight (segments never stretch) |
| **Clashes** | Microtubules overlapping in space — physically impossible |
| **Clearance** | How close the linker, base and spoke pass to microtubules, counting the strand's own thickness; negative means it is inside a tubule wall, and it is drawn in red |
| **Approach** | How squarely each strand meets the protofilament it binds: +1 head-on from outside the tubule, 0 tangential, below 0 impossible |
| **Diameter** | Resulting centriole diameter in nm |

### Design principles

- **The triplet ring is the primary scaffold.** Triplets keep even
  spacing; the cartwheel attaches where it can reach. When the two
  symmetries disagree, the cartwheel strains and surplus triplets are left
  unattached — not the other way round.
- **Nothing stretches.** Segments have fixed contour length. Under
  compression they may *buckle* (bow, shortening end-to-end distance).
- **Bonds have strengths.** Each connection is a spring whose stiffness is
  its bond strength, so weak bonds yield first and the load distribution is
  a result rather than an assumption. Strength order (strong → weak):
  linker–A-tubule, linker–C-tubule, pinhead–triplet, pinhead–spoke,
  base–pinhead, base–linker.
- **Joints have individual tolerances.** A rotation means different things
  at different connections, so each is graded against its own limits.
  Contacts on microtubules are tightest — they grip a rigid, ordered
  lattice at defined protofilaments — while the triplet axis and base are
  most permissive, being composites rather than single interfaces. The
  SAS-6 spoke is deliberately *not* treated as unusually soft: many SAS-6
  rings stack along the centriole axis, which stiffens what looks floppy
  in a single 2D slice.

  | Joint / contact | OK | HARD | SEVERE |
  |---|---|---|---|
  | Linker ↔ A-tubule | ≤ 8° | 8–20° | > 20° |
  | Linker ↔ C-tubule | ≤ 8° | 8–20° | > 20° |
  | Pinhead ↔ A-tubule | ≤ 10° | 10–22° | > 22° |
  | Spoke vs radial | ≤ 15° | 15–35° | > 35° |
  | Pinhead vs spoke | ≤ 15° | 15–30° | > 30° |
  | Triplet axis vs spoke | ≤ 20° | 20–40° | > 40° |
  | Base vs spoke | ≤ 20° | 20–40° | > 40° |

  **These are assumed tolerances, not verdicts on feasibility.** Centrioles
  are known to assemble in conditions the model grades SEVERE — for example a
  28 nm SAS-6 spoke builds real centrioles, yet grades SEVERE at the linker
  contacts. Read the grades as *how far from wild type*, and judge from the
  geometry and the numbers.

  These values are reasoned heuristics, not measurements — no
  sub-tomogram angular variance was available to calibrate them. The rank
  ordering is defensible; the absolute numbers are not. Use
  `band_sensitivity()` to confirm a conclusion survives the bands being
  scaled, and only report conclusions that do:

  ```python
  from centriole_kinematic import Geometry, band_sensitivity
  band_sensitivity(Geometry(N_cw=8, N_mt=9))   # scales all bands 0.5x - 2x
  ```

- **The spoke may hinge at its head, or not — your choice.** By default the
  SAS-6 coiled coil can turn where it meets its head on the hub ring, and the
  `spoke` grade reports how far. `spoke_pivot=False` (app: *Spoke cannot bend
  at the SAS-6 head*; CLI: `--lock-spoke`) instead holds each spoke on the
  radius through its own head, so it strains outwards but never hinges — the
  assumption most treatments of the cartwheel make.

  Wild type is nearly indifferent: 255 nm either way, overall strain 1.4° vs
  1.6°, because a wild-type spoke barely pivots. A *mismatched* cartwheel is
  transformed. With 8 spokes on 9 triplets the free spoke swings −20.9° to
  reach; locked, it cannot, so it buckles 27.5% instead, the diameter opens
  from 256 to 279 nm, and the loop closures begin to gape.

  | 8 cw / 9 mt | free spoke | locked spoke |
  |---|---|---|
  | Outer diameter | 256.4 nm | 279.0 nm |
  | Spoke rotation | −20.9° | 0° (by construction) |
  | Spoke buckling | 1.1% | 27.5% |
  | Joint strain, RMS | 7.1° | 19.2° |
  | Worst loop gap | 0.002 nm | 0.186 nm |

  Run both. A conclusion that survives the toggle does not depend on the
  assumption; one that flips is being carried by it. Note that with the spoke
  locked the `spoke` grade reads 0° **by construction, not by relaxation** —
  and that locking removes `N_cw` degrees of freedom, leaving a 9-fold chain
  exactly determined at 36 unknowns against 36 loop constraints, with no
  slack to absorb a mismatch.

- **Hub size follows symmetry.** SAS-6 heads oligomerise at a fixed
  head-to-head spacing, so hub radius = `d / (2·sin(π/N))`. An 8-fold
  cartwheel is automatically tighter than a 9-fold one.

---

## Where the numbers come from

Every dimension is measured from a cryo-ET-derived schematic
(`data/2x_units_wScale.svg`) by [`svg_calibration.py`](svg_calibration.py),
which is re-runnable so the provenance of each constant is traceable:

```bash
python svg_calibration.py
```

Key measured values (nm):

| Quantity | Value |
|---|---|
| SAS-6 dimer (head + coiled coil) | 50.0 (4.98 + 45.0) |
| Pinhead span | 21.3 |
| Triplet base | 34.7 |
| Tubule radius | 11.9 |
| A→B / B→C spacing | 18.8 / 19.6 |
| Triplet axis vs spoke | −56.9° |

Two findings worth noting: the B- and C-tubules carry 9 protofilaments at
the **same 27.7° lattice pitch as the A-tubule's 13** (not 360/9), and the
A-C linker is **bent, not straight**.

**Three values do not come from `svg_calibration.py`.** The A-C linker's two
arm lengths and its vertex angle (14.6 nm to C, 11.6 nm to A, meeting at
138.7°) were read off the schematic by hand; the calibration measures only
where the linker *attaches*. A fourth, the strand half-width used for
clearance and steric checks, is assumed outright at 1 nm — half the diameter
of a two-stranded coiled coil — and nothing kinematic depends on it.

**One value is fitted, not measured.** `head_contact` (the SAS-6 head-head
spacing, which sets the hub radius) cannot be read off the schematic,
because the drawing's centre is ambiguous — its spokes sit at exactly 40°
apart while its triplets sit at 41.4°. It is instead fitted as the value
that lets the measured wild-type unit close with least strain. Re-derive it
with `calibrate_head_contact()`. Every other dimension in the table above
comes straight from `svg_calibration.py`.

**Calibration caveat.** Protofilament positions were read from label
anchors in the SVG, which carries roughly ±2° (~0.4 nm) of noise. That sets
the floor on how exactly wild type can be made to close, and is a
measurement limit rather than a solver limit.

---

## Validation

The model was never told that 9-fold symmetry is special. Sweeping each
parameter independently with `solve_chain`, the measured wild-type value sits
at a strain minimum in every case:

| Parameter | Measured | Strain-minimising |
|---|---|---|
| Symmetry | 9 | **9** |
| SAS-6 coiled coil | 45.0 nm | **45 nm** |
| Triplet base | 34.7 nm | **34.7 nm** |
| A-tubule protofilaments | 13 | **13** |
| Tubules per blade | 3 | **3** |

Wild-type diameter comes out at **255 nm**, against roughly 250 nm measured
in real centrioles.

Re-verified after the 2026 correction to the contact rest angles and the
addition of the steric checks — every minimum is unmoved, and wild type solves
with all seven joints graded OK, zero clashes, and every contact reachable.

---

## Installation

Requires Python 3.9 or newer.

```bash
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Usage

### Interactive app (recommended)

Double-click the launcher for your platform, or run:

```bash
streamlit run app.py
```

Four tabs:

- **Cross-section** — set any parameter, see the geometry redraw with every
  metric: diameters, joint grades, bond loads, buckling, clashes,
  clearances.
- **Parameter scan** — vary one parameter over a range and plot *all*
  metrics against it, grouped into size / joint rotation / buckling /
  steric integrity / bond load / clearance, each with a dotted wild-type
  reference. Cross-sections are rendered above the plots. Downloadable as
  CSV. Symmetry can be scanned coupled (cartwheel and triplets together,
  the biologically coherent case) or either one alone against a fixed
  partner to probe deliberate mismatch. With register shift enabled, the
  protofilament register the solver chose is reported too.
- **2-parameter scan** — cross any two parameters (e.g. SAS-6 coiled-coil
  length against symmetry), choose from 30 metrics, and view them as
  heatmaps, as families of curves with a wild-type reference, or both.
  Choosing the two single-symmetry options reproduces the
  cartwheel-versus-triplet mismatch grid.
- **How to read this** — a plain-English key to every readout.

Symmetry uses sliders, tubules-per-blade a radio button, and every distance
is a typed number box so exact values can be entered.

### Command line

```bash
python run_analysis.py                                  # standard report
python run_analysis.py --cartwheel 8                    # 8 spokes, 9 triplets
python run_analysis.py --triplets 8 --tubules 2         # 8 doublets
python run_analysis.py --sweep spoke_rod --from 35 --to 60
python run_analysis.py --cartwheel 8 --lock-spoke        # spoke cannot hinge
python run_analysis.py --help                           # all options
```

Everything lands in `results/` as PNG figures, CSV tables and a
`summary.txt`. See [`QUICKSTART.html`](QUICKSTART.html) for a fuller guide.

### From Python

```python
from centriole_kinematic import Geometry, draw, sweep, sweep2
from centriole_chain import solve_chain, best_registers

sol = solve_chain(Geometry(N_cw=8, N_mt=9))   # 8-fold cartwheel, 9 triplets
rigid = solve_chain(Geometry(N_cw=8, N_mt=9), spoke_pivot=False)  # no hinge at the head
print(sol.report())
print(sol.outer_diameter, sol.worst_bond, sol.unattached_triplets)
print(sol.feasible, sol.reachable)            # is every contact reachable?

sweep("spoke_rod", [35, 40, 45, 50, 55], solver=solve_chain)
sweep("N_both", range(6, 12), solver=solve_chain)
sweep2("N_both", range(7, 11), "spoke_rod", [40, 45, 50], solver=solve_chain)
```

### Parameters

Commonly varied: `N_cw` (cartwheel symmetry), `N_mt` (triplet symmetry),
`spoke_rod` (SAS-6 coiled-coil length), `base_length` (triplet base).

Less commonly: `pinhead_span`, `linker_length`, `MTn` (3/2/1 for
triplet/doublet/singlet), `n_pf_A` (protofilament count, which also sets
tubule radius).

Optional: `solve_chain(..., register_shift=True)` slides both A-C linker
contacts over ±2 protofilaments (25 combinations, ~290 s), testing whether a
mutant could relieve strain by re-registering instead of deforming.
`best_registers(sol, 3)` returns the leading candidates with their metrics.

Read the result carefully. Registers needing a strand to reach its
protofilament *through* a microtubule are excluded outright — not ranked, as
they are not answers. The rest are ordered by **model cost**, which is built
from joint tolerances and bond strengths that were reasoned rather than
measured, so **the cheapest register is not the most likely one**. Treat the
ranking as a shortlist and compare the diameter and tilt columns against your
own data. See [`docs/blade_rotation_discrepancy.md`](docs/blade_rotation_discrepancy.md)
for a worked case where this distinction changed the conclusion.

---

## Repository layout

```
START HERE (Mac).command / (Windows).bat   double-click launchers
app.py                   interactive Streamlit app (start here)
run_analysis.py          command-line entry point
centriole_chain.py       the solver in use: connection points cannot separate
centriole_kinematic.py   geometry, drawing, metrics, and the older spring solver
svg_calibration.py       extracts every constant from the schematic
QUICKSTART.html          one-page guide for non-programmers
data/                    source schematics (SVG, Illustrator, PNG)
docs/                    self-contained notes on three investigations
presentations/           slide decks
legacy/                  the original MATLAB script and earlier ports
results/                 generated output (not tracked by git)
```

**Two solvers, and only one of them should be used.** `centriole_chain.solve_chain`
makes four of the six connections exact by construction, so parts cannot drift
apart under load. `centriole_kinematic.solve` is the older spring network: its
bonds open by 1–3 nm under a large perturbation, visibly separating the spoke
from the pinhead. The app and `run_analysis.py` both use `solve_chain`; `solve`
is kept for the soft-mode and blooming analyses that were built on it.

`docs/` holds self-contained notes on three investigations, each resumable
without re-deriving anything:

- [**A-C linker register**](docs/blade_rotation_discrepancy.md) — *live, and
  recently reopened.* A measured mutant disagreed with the model until the
  linker was allowed to re-register. Two defects in the register machinery
  have since been fixed, which moved the answer: the register that reproduced
  the measurement turns out to be at or past the edge of what is geometrically
  reachable. Still a cryo-ET-testable prediction, but a different one.
- [blooming / iris motion](docs/blooming_handoff.md) — parked on the 2D
  limitation.
- [soft-mode analysis](docs/parked_soft_modes.md) — parked as
  uninterpretable without an absolute energy scale.

`legacy/` holds the original MATLAB `centriole.m` written by Aitana Neves,
which this project grew from, plus two intermediate Python models. They are
superseded by `centriole_kinematic.py` and kept only for provenance.

---

## Authorship and citation

Written by **Georgios Hatzopoulos** (UPGON). The model, its calibration and
this implementation are his; `legacy/centriole.m`, the parametric MATLAB
script the project grew from, was written by **Aitana Neves**.

If you use the tool or its results, please cite it — GitHub's *Cite this
repository* button reads [`CITATION.cff`](CITATION.cff) and will give you a
formatted reference.

---

## Known limitations

- **Two-dimensional.** A single cross-section; the real centriole twists
  along its length, which is not modelled.
- **Bond strengths are relative,** set by rank order rather than measured
  energies. Rankings between perturbations are more meaningful than
  absolute force values.
- **Joint rotations are unbounded** — they are penalised in bands rather
  than hard-limited, so a `FORBIDDEN` result means "very costly", not
  "geometrically impossible".
- **Steric checks cover** tubule–tubule overlap and strand–tubule
  clearance, not every possible contact. The pinhead is not steric-checked.
- **Strand thickness is assumed,** at 1 nm half-width, so clearances within a
  few tenths of a nanometre of zero are marginal. The companion reachability
  check — does a strand arrive from outside its tubule or from within? — needs
  no such assumption and is the firmer of the two.
- **In the chain formulation the A-C linker has one degree of freedom,** so
  its two contacts report the same angle and its orientation is effectively
  charged twice relative to other joints.
