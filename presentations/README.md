# Presentation decks

Three decks, all generated from source so a slide cannot quietly disagree with
the tool:

| File | For |
|---|---|
| `01_CenGeometry_Introduction.pptx` | Introducing the tool. Biology first, then a signposted six-slide **formalism** section for a mixed audience with physicists in it. |
| `02_CenGeometry_Critical_Evaluation.pptx` | What the model cannot currently be trusted to tell you, ranked by severity. |
| `03_CenGeometry_Expert_Recommendations.pptx` | Where to take it next, and six testable predictions. |

## Rebuilding

Run all three steps, in order, from this directory:

```bash
PYTHONPATH=.. ../.venv/bin/python make_figs.py
../.venv/bin/python make_equations.py
node check_layout.js
```

`check_layout.js` builds the decks *and* checks the geometry; use it instead of
`node build_decks.js` unless you have a reason not to.

Allow about eight minutes — most of it is the 25-register search in
`make_figs.py`.

| Step | Produces |
|---|---|
| `make_figs.py` | every figure (`figs/d_*.png`) and every number the decks quote (`figs/_values.json`) |
| `make_equations.py` | the typeset equations (`figs/eq_*.png`) plus their aspect ratios |
| `check_layout.js` | the three `.pptx` files, and a report of anything off-canvas or overlapping |

## Why it is built this way

**Numbers come from `_values.json`, never typed into a slide.** The decks
originally carried numbers hand-copied from the older spring-network
`solve()`. When the chain solver became the default those numbers silently went
stale, and several were still being presented long after they had stopped being
true — one prediction ("zero clashes from 9 to 18 protofilaments") was a
correction to an earlier claim that had itself gone stale. `build_decks.js` now
reads `figs/_values.json` and interpolates, so re-running the pipeline is the
only way to change a number and every number moves together.

**Equations are images.** pptxgenjs has no maths support, and there is no LaTeX
on this machine, so `make_equations.py` renders each expression with
matplotlib's mathtext. Two of its grammar's gaps shape the result: `\underbrace`
and `\text` do not exist, so the labelling of each residual block lives in the
slide text beside the equation rather than under it. Every expression is parsed
before anything is written, so a typo fails loudly instead of leaving a slide
with a missing image.

**Layout is checked, not eyeballed.** Everything is absolutely positioned, so
the two ways a slide breaks are an element running off the canvas and two
elements overlapping — neither of which makes the build fail, and neither of
which is visible without opening PowerPoint. `check_layout.js` intercepts every
placement call and checks the geometry directly. It caught four real collisions
and one off-canvas equation when the formalism section was first written.

## Editing

Content lives in `build_decks.js`, one function per deck. To change a number,
change what `make_figs.py` computes — not the slide. To add an equation, add it
to the `EQ` dict in `make_equations.py` and place it with `eq()` or `eqCard()`,
which size the image from its measured aspect ratio so nothing is ever
stretched.

`figs/_layout.json` is a by-product of the checker and is only useful for
rendering an approximate preview of a slide during editing.
