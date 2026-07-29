# CenGeometry

A geometric and mechanical model of the centriole cross-section, built to
predict how structural perturbations — SAS-6 mutations that change
cartwheel symmetry, altered coiled-coil lengths, loss of the C-tubule —
change the architecture of the whole organelle.

**New here? Open [`QUICKSTART.html`](QUICKSTART.html) in a web browser.**
It is a one-page, no-coding-required guide to installing and running the
model.

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
| **Joint rotation** | How far each connection turned from its wild-type angle, graded `OK` / `HARD` / `FORBIDDEN` against per-joint limits (see below) |
| **Bond load** | Which connection carries the most strain, i.e. which would rupture first |
| **Buckling** | Which segments had to bow because space became too tight (segments never stretch) |
| **Clashes** | Microtubules overlapping in space — physically impossible |
| **Clearance** | How close the linker, base and spoke pass to microtubules; negative means passing through one |
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

  | Joint / contact | OK | HARD | FORBIDDEN |
  |---|---|---|---|
  | Linker ↔ A-tubule | ≤ 8° | 8–20° | > 20° |
  | Linker ↔ C-tubule | ≤ 8° | 8–20° | > 20° |
  | Pinhead ↔ A-tubule | ≤ 10° | 10–22° | > 22° |
  | Spoke vs radial | ≤ 15° | 15–35° | > 35° |
  | Pinhead vs spoke | ≤ 15° | 15–30° | > 30° |
  | Triplet axis vs spoke | ≤ 20° | 20–40° | > 40° |
  | Base vs spoke | ≤ 20° | 20–40° | > 40° |

  These values are reasoned heuristics, not measurements — no
  sub-tomogram angular variance was available to calibrate them. The rank
  ordering is defensible; the absolute numbers are not. Use
  `band_sensitivity()` to confirm a conclusion survives the bands being
  scaled, and only report conclusions that do:

  ```python
  from centriole_kinematic import Geometry, band_sensitivity
  band_sensitivity(Geometry(N_cw=8, N_mt=9))   # scales all bands 0.5x - 2x
  ```

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
| A-C linker arms | 14.6 (to C) and 11.6 (to A), meeting at 138.7° |
| Tubule radius | 11.9 |
| A→B / B→C spacing | 18.8 / 19.6 |
| Triplet axis vs spoke | −56.9° |

Two findings worth noting: the B- and C-tubules carry 9 protofilaments at
the **same 27.7° lattice pitch as the A-tubule's 13** (not 360/9), and the
A-C linker is **bent, not straight**.

**Calibration caveat.** Protofilament positions were read from label
anchors in the SVG, which carries roughly ±2° (~0.4 nm) of noise. That sets
the floor on how exactly wild type can be made to close, and is a
measurement limit rather than a solver limit.

---

## Validation

The model was never told that 9-fold symmetry is special. Sweeping each
parameter independently, the measured wild-type value sits at a strain
minimum in every case:

| Parameter | Measured | Strain-minimising |
|---|---|---|
| Symmetry | 9 | **9** |
| SAS-6 coiled coil | 45.0 nm | **45 nm** |
| Triplet base | 34.7 nm | **34.7 nm** |
| A-tubule protofilaments | 13 | **13** |
| Tubules per blade | 3 | **3** |

Wild-type diameter comes out at **255 nm**, against roughly 250 nm measured
in real centrioles.

---

## Installation

Requires Python 3.9 or newer.

```bash
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Usage

### Without writing code

```bash
python run_analysis.py                                  # standard report
python run_analysis.py --cartwheel 8                    # 8 spokes, 9 triplets
python run_analysis.py --triplets 8 --tubules 2         # 8 doublets
python run_analysis.py --sweep spoke_rod --from 35 --to 60
python run_analysis.py --help                           # all options
```

Everything lands in `results/` as PNG figures, CSV tables and a
`summary.txt`. See [`QUICKSTART.html`](QUICKSTART.html) for a fuller guide.

### From Python

```python
from centriole_kinematic import Geometry, solve, draw, sweep, sweep2

sol = solve(Geometry(N_cw=8, N_mt=9))    # 8-fold cartwheel, 9 triplets
print(sol.report())
print(sol.outer_diameter, sol.worst_bond, sol.unattached_triplets)

sweep("spoke_rod", [35, 40, 45, 50, 55])          # one parameter
sweep2("N_cw", range(7, 12), "N_mt", range(7, 12))  # two-parameter grid
```

### Parameters

Commonly varied: `N_cw` (cartwheel symmetry), `N_mt` (triplet symmetry),
`spoke_rod` (SAS-6 coiled-coil length), `base_length` (triplet base).

Less commonly: `pinhead_span`, `linker_length`, `MTn` (3/2/1 for
triplet/doublet/singlet), `n_pf_A` (protofilament count, which also sets
tubule radius).

Optional: `solve(..., register_shift=True)` lets the pinhead and linker
contacts slide to neighbouring protofilaments, testing whether a mutant
could relieve strain by re-registering instead of deforming.

---

## Repository layout

```
run_analysis.py          command-line entry point (start here)
centriole_kinematic.py   the model
svg_calibration.py       extracts every constant from the schematic
QUICKSTART.html          one-page guide for non-programmers
data/                    source schematics (SVG, Illustrator, PNG)
legacy/                  the original MATLAB script and earlier ports
results/                 generated output (not tracked by git)
```

`legacy/` holds the original MATLAB `centriole.m` this project grew from,
plus two intermediate Python models. They are superseded by
`centriole_kinematic.py` and kept only for provenance.

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
  clearance, not every possible contact.
