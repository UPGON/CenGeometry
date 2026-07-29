# Blooming / iris motion — parked notes

**Status: parked deliberately.** Everything here is self-contained so a
future session can pick it up cold without re-deriving anything. Nothing
in the main model depends on resolving it.

---

## The claim being tested

Laporte et al., *Cell* 2024 — [Time-series reconstruction of the molecular
architecture of human centriole assembly](https://www.cell.com/cell/fulltext/S0092-8674(24)00316-7)
([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11060037/)).

Two features of the paper's actual claim, both of which matter:

1. The **bloom phase is assembly-time growth** of a procentriole — radial
   and longitudinal expansion until tubulin length reaches ~119 nm, after
   which diameter stabilises and only length grows. It is *not* a
   reversible mechanical iris of a mature centriole.
2. The **iris-diaphragm motion is longitudinal**, arising from A-C linker
   *twist*, i.e. a progressive change along the procentriole axis rather
   than an in-plane motion of one cross-section.

Key quantitative anchor: **CEP135 — proposed to be part of the triplet
base — changes diameter from 85 nm to 182 nm** across assembly.

Related: [The A-C linker controls centriole structural integrity and
duplication](https://www.nature.com/articles/s41467-025-62154-6) (Nat
Commun 2025); older precedent in [the iris diaphragm
model](https://pubmed.ncbi.nlm.nih.gov/2268874/) (1990).

---

## What was established

### The scale problem (most important)

85 → 182 nm is a ~48 nm **radial** change. The triplet base in this model
is 34.7 nm of total contour length. Even fully collapsing and re-extending
it cannot span that. Base extension alone therefore cannot account for the
bloom — consistent with the paper's own framing that the base extends *to
accommodate* growth rather than to drive it.

### Driving base extension directly

`blooming_scan()` sweeps the base from bent to fully extended. Note the
schematic the model is calibrated on shows the base **fully extended**, so
extension is `base_buckle = 0` and bending is positive.

| | bent (23.6 nm chord) | extended (34.7 nm) |
|---|---|---|
| Outer diameter | 251.6 nm | 255.0 nm (+3.4) |
| Lumen | 153.8 nm | 160.2 nm (+6.5) |
| Triplet tilt | −57.9° | −61.0° (+3.1°) |
| Spoke deviation | 15.0° (HARD) | 1.2° (OK) |
| Base–linker bond load | 1.047 | 0.019 (55× lower) |

So it produces little expansion, and it does **not** leave the cartwheel
unaffected.

### Reverse causality test

Driving radial expansion instead (spoke 30 → 45 nm, A-ring 163 → 184 nm),
the base stays at **0% bend throughout** — already fully extended, with no
slack to pay out. It cannot be the accommodating element as currently
modelled.

### Soft-mode result

With the cartwheel intact there is **no in-plane iris mode**: the single
very soft mode is the trivial global rotation, and all other soft modes
are non-collective. Weakening the pinhead–spoke bond makes a genuinely
collective diameter+tilt mode appear. **The cartwheel suppresses the
in-plane iris.**

---

## The limitation that blocks progress

**The model is a single 2D cross-section; the paper describes a
longitudinal twist.** The negative in-plane result therefore does *not*
contradict Laporte et al. — it addresses a different motion. This is the
thing to fix before drawing any conclusion about their claim.

---

## Two routes forward

1. **2.5D stack.** Model several cross-sections at fixed axial spacing
   with a prescribed twist increment between them, coupled through the A-C
   linker. Far cheaper than full 3D and directly tests whether A-C linker
   twist can generate the observed axial gradient. This is the route that
   actually engages the paper's claim.
2. **Procentriole states as separate geometries.** Rather than treating
   early states as a bent-base version of the mature one, give them their
   own (smaller, incomplete) architecture — fewer protofilaments, no
   C-tubule, smaller hub — and calibrate each against the U-ExM diameters
   in the paper (PLK4 77, SAS-6 75, STIL 89, CPAP 120, CEP135 85 nm).

---

## Code already in place

Nothing needs writing to resume:

```python
from centriole_kinematic import Geometry, solve, blooming_scan

solve(Geometry(), base_buckle=0.20)   # drive the base; 0 = fully extended
blooming_scan(Geometry(), n=9)        # bent -> extended sweep, returns a DataFrame
```

`Solution` reports `lumen_diameter` and `triplet_tilt` (mean triplet axis
from radial) alongside the usual outputs.

Results from the runs above: `results/blooming_scan.csv`,
`results/bloom_causality.csv`, `results/blooming_no_linker.csv`.
