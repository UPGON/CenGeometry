#!/usr/bin/env python3
"""Render the model's equations to transparent PNGs for the intro deck.

pptxgenjs has no maths support, and Unicode-in-a-textbox reads badly to an
audience that works with equations. There is no LaTeX on this machine, so
these use matplotlib's mathtext, which covers everything needed here --
with two deliberate omissions. `\\underbrace` and `\\text` are not in its
grammar, so the labelling of each residual block lives in the slide text
beside the equation instead of under it, which reads better in a talk
anyway.

    ../.venv/bin/python make_equations.py

Every expression is parsed before anything is written, so a typo fails
loudly here rather than producing a slide with a missing image.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGS = Path(__file__).resolve().parent / "figs"
FIGS.mkdir(exist_ok=True)

INK = "#1A1A1A"

#: name -> (mathtext body, font size)
EQ = {
    # ---------------------------------------------------------- the problem
    "objective": (
        r"$\mathbf{z}^{\star} = \mathrm{arg\,min}_{\;\mathbf{z}\,\in\,\mathcal{B}}"
        r"\;\; \frac{1}{2}\,\|\,\mathbf{r}(\mathbf{z})\,\|^{2}$", 32),
    "residual": (
        r"$\mathbf{r}(\mathbf{z}) = \left[\; k_b\,s_c\,\mathbf{g}_c \;\;;\;\;"
        r"k_a\,\Phi_j(\delta_j) \;\;;\;\; k_u\,\Delta\varphi_i \;\;;\;\;"
        r"k_\beta\,\beta_s \;\;;\;\; k_t\,o_{ij} \;\;;\;\;"
        r"k_w\,h(-d_{st}) \;\right]$", 21),

    # ------------------------------------------------- forward kinematics
    "hub": (
        r"$R_{hub} = \frac{d_{head}}{2\,\sin(\pi/N_{cw})}"
        r"\qquad \mathbf{x}^{head}_{j} = R_{hub}\,\hat{u}(\varphi_j)"
        r"\qquad \varphi_j = \frac{2\pi j}{N_{cw}}$", 23),
    "chain": (
        r"$\mathbf{x}^{tip}_{j} = \mathbf{x}^{head}_{j} + L_s\,(1-\beta_j)\;"
        r"\hat{u}(\varphi_j + \alpha_j)$" "\n\n"
        r"$\mathbf{x}^{Aend}_{p} = \mathbf{x}^{tip}_{j} + L_p\;\hat{u}(\theta_p)"
        r"\qquad \theta_p = \varphi_j + \alpha_j + \alpha^{pin}_{p}$" "\n\n"
        r"$\mathbf{c}^{A}_{t} = \mathbf{x}^{Aend}_{p} - R_t\;"
        r"\hat{u}\left(\psi_t + a_{pf}(A,\,pin)\right)"
        r"\qquad \psi_t = \varphi_j + \alpha_j + \alpha^{trip}_{p}$", 20),
    "gaps": (
        r"$\mathbf{g}_{linkA} = \mathbf{x}^{endA}_{i} - \mathbf{p}^{A8}_{i-1}"
        r"\qquad\qquad \mathbf{g}_{base} = \mathbf{x}^{tip}_{base} - \mathbf{v}_i$", 22),
    "dof": (
        r"$n = N_{cw} + 3\,N_{pair} + N_{mt} = 45"
        r"\qquad m = 2\,N_{mt} + 2\,N_{pair} = 36$" "\n\n"
        r"$n - m = 9\;\;\mathrm{spare\;DOF}"
        r"\qquad\Rightarrow\qquad \alpha_j \equiv 0"
        r"\;\;\Rightarrow\;\; 36 - 36 = 0$", 21),

    # ------------------------------------------------------- the penalties
    "band": (
        r"$\Phi_j(\delta) = w_1\,\min(|\delta|,\,\delta^{ok}_j)"
        r"\; + \; w_2\,\mathrm{clip}\left(|\delta| - \delta^{ok}_j,\;0,\;"
        r"\delta^{hard}_j - \delta^{ok}_j\right)"
        r"\; + \; w_3\,\max(|\delta| - \delta^{hard}_j,\,0)$" "\n\n"
        r"$(w_1,\,w_2,\,w_3) = (1,\;5,\;30)$", 19),
    "buckle": (
        r"$\|\,\mathbf{b} - \mathbf{a}\,\| = L\,(1 - \beta)"
        r"\qquad \beta \in [\,0,\;\beta_{max}\,]$", 24),
    "steric": (
        r"$o_{ij} = \max\left(0,\;2R_t - \|\,\mathbf{c}_i - \mathbf{c}_j\,\|\right)$"
        "\n\n"
        r"$d_{st} = \min_{\mathbf{q}\,\in\,ab}\|\,\mathbf{q} - \mathbf{c}_t\,\|"
        r"\; - \; R_t \; - \; w_s$", 22),

    # ------------------------------------------- the discrete register axis
    "pfangle": (
        r"$a_{pf}(T,\,i,\,k) = A_T - \frac{2\pi}{n_A}\,(i - 1 + k)"
        r"\qquad R_t = \frac{n_A\,w_{pf}}{2\pi}$", 23),
    "reref": (
        r"$\theta^{rest}_{c}(k) \; = \; \theta^{rest}_{c}(0)"
        r"\; - \; \frac{2\pi}{n_A}\,k"
        r"\qquad\qquad \frac{360^{\circ}}{13} = 27.7^{\circ}$", 24),
    "approach": (
        r"$\cos\chi \; = \; -\,\hat{u}_{arm}\cdot\hat{n}"
        r"\qquad \hat{n} = \frac{\mathbf{x}_{contact} - \mathbf{c}_t}{R_t}$", 24),

    # ------------------------------------------------------------- solver
    "gaussnewton": (
        r"$\left(\mathbf{J}^{\top}\mathbf{J} + \lambda\,\mathbf{D}\right)"
        r"\Delta\mathbf{z} \; = \; -\,\mathbf{J}^{\top}\mathbf{r}"
        r"\qquad \mathbf{J}_{ij} = \frac{\partial r_i}{\partial z_j}$", 23),
}


def check():
    """Parse every expression up front so a typo fails here, not on a slide."""
    from matplotlib.mathtext import MathTextParser
    parser = MathTextParser("agg")
    bad = []
    for name, (body, _) in EQ.items():
        for line in body.split("\n"):
            line = line.strip()
            if not line:
                continue
            try:
                parser.parse(line, dpi=100, prop=None)
            except Exception as exc:                       # noqa: BLE001
                bad.append(f"  {name}: {str(exc).strip().splitlines()[-1]}")
    return bad


def render(name, body, size):
    fig = plt.figure(figsize=(0.01, 0.01))
    txt = fig.text(0, 0, body, fontsize=size, color=INK, ha="left", va="bottom",
                   linespacing=1.9)
    fig.canvas.draw()
    bb = txt.get_window_extent()
    fig.set_size_inches(bb.width / fig.dpi + 0.15, bb.height / fig.dpi + 0.15)
    fig.savefig(FIGS / f"eq_{name}.png", dpi=300, transparent=True,
                bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    return bb.width / max(bb.height, 1e-9)


if __name__ == "__main__":
    problems = check()
    if problems:
        print("mathtext could not parse:", *problems, sep="\n")
        sys.exit(1)
    print(f"rendering {len(EQ)} equations")
    ratios = {}
    for k, (body, size) in EQ.items():
        ratios[k] = round(render(k, body, size), 3)
        print(f"  -> figs/eq_{k}.png   aspect {ratios[k]}")
    (FIGS / "_eq_aspect.json").write_text(
        "{\n" + ",\n".join(f' "{k}": {v}' for k, v in ratios.items()) + "\n}\n")
