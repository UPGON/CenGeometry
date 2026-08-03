#!/usr/bin/env python3
"""Regenerate every figure and number the decks quote.

Run before `node build_decks.js`:

    PYTHONPATH=.. ../.venv/bin/python make_figs.py

Everything here uses `centriole_chain.solve_chain` -- the solver the app and
the CLI actually use. The decks previously carried numbers from the older
spring-network `solve()`, which drifted once the chain form became the
default. Regenerating from one script keeps the slides and the tool telling
the same story, and `figs/_values.json` records exactly what was shown.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from centriole_kinematic import (
    JOINT_BANDS, W_FREE, W_HARD, W_FORBID, Geometry, angle_penalty, draw,
    set_param, summarise,
)
from centriole_chain import solve_chain, best_registers

FIGS = Path(__file__).resolve().parent / "figs"
FIGS.mkdir(exist_ok=True)

PURPLE, CYAN, GREEN, BAD, MUTE = "#6B3FA0", "#2FA8A8", "#5C7A45", "#C62828", "#6B7280"
V: dict = {}


def rms(sol):
    return float(np.sqrt(np.mean([v["rms"] ** 2 for v in sol.joint_strain.values()])))


def save(fig, name, **kw):
    fig.savefig(FIGS / name, dpi=160, bbox_inches="tight", **kw)
    plt.close(fig)
    print(f"  -> figs/{name}")


# ---------------------------------------------------------------- wild type
def wild_type():
    sol = solve_chain(Geometry())
    V["wt"] = summarise(sol)
    V["wt"]["joint_rms_deg"] = round(rms(sol), 2)
    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    draw(sol, ax=ax, title="")
    save(fig, "d_wt.png")
    return sol


# --------------------------------------------------------------- validation
def validation():
    ns = list(range(6, 13))
    strain = []
    for n in ns:
        strain.append(round(rms(solve_chain(Geometry(N_cw=n, N_mt=n))), 2))
    V["symmetry"] = {"N": ns, "strain": strain}

    fig, ax = plt.subplots(figsize=(6.6, 4.1))
    ax.plot(ns, strain, "o-", color=PURPLE, lw=2.4, ms=9)
    ax.axvline(9, color=BAD, ls="--", lw=1.8)
    ax.annotate("measured\nwild type", (9, max(strain) * 0.72), color=BAD,
                ha="center", fontsize=12, fontweight="bold")
    ax.set_xlabel("symmetry (N-fold)", fontsize=12)
    ax.set_ylabel("joint strain (deg, RMS)", fontsize=12)
    ax.grid(alpha=0.3)
    save(fig, "d_validation.png")

    # each parameter's strain minimum, against its measured value
    def minimum(param, vals):
        r = [(v, rms(solve_chain(set_param(Geometry(), param, v)))) for v in vals]
        return {"values": [v for v, _ in r], "strain": [round(s, 2) for _, s in r],
                "best": min(r, key=lambda kv: kv[1])[0]}

    V["minima"] = {
        "N_both": {"measured": 9, "best": ns[int(np.argmin(strain))]},
        "spoke_rod": {"measured": 45.03,
                      **minimum("spoke_rod", [35, 40, 45.03, 50, 55])},
        "base_length": {"measured": 34.68,
                        **minimum("base_length", [25, 30, 34.68, 40, 45])},
        "n_pf_A": {"measured": 13, **minimum("n_pf_A", list(range(11, 16)))},
        "MTn": {"measured": 3, **minimum("MTn", [1, 2, 3])},
    }


# ----------------------------------------------------------------- mismatch
def mismatch():
    cases = [("Wild type", Geometry()),
             ("8-fold cartwheel, 9 triplets", Geometry(N_cw=8, N_mt=9)),
             ("Fully 8-fold", Geometry(N_cw=8, N_mt=8))]
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.9))
    out = {}
    for ax, (tag, g) in zip(axes, cases):
        s = solve_chain(g)
        un = f", {len(s.unattached_triplets)} unattached" if s.unattached_triplets else ""
        ax.set_title(f"{tag}\n{s.outer_diameter:.0f} nm  ·  spoke {s.joint_bands['spoke']}{un}",
                     fontsize=12)
        draw(s, ax=ax, title="")
        out[tag] = dict(diameter=round(s.outer_diameter, 1), rms=round(rms(s), 2),
                        spoke_band=s.joint_bands["spoke"],
                        spoke_deg=round(s.joint_strain["spoke"]["max"], 1),
                        triplet_deg=round(s.joint_strain["triplet"]["max"], 1),
                        unattached=len(s.unattached_triplets))
    fig.tight_layout()
    save(fig, "d_mismatch.png")
    V["mismatch"] = out


# --------------------------------------------------------------- spoke lock
def spoke_lock():
    g = Geometry(N_cw=8, N_mt=9)
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.9))
    out = {}
    for ax, pv in zip(axes, (True, False)):
        s = solve_chain(g, spoke_pivot=pv)
        tag = "Spoke free to hinge at the head" if pv else "Spoke held radial"
        ax.set_title(f"{tag}\n{s.outer_diameter:.0f} nm  ·  spoke buckled "
                     f"{s.buckling['spoke_rod']:.0f}%", fontsize=12)
        draw(s, ax=ax, title="")
        out["free" if pv else "locked"] = dict(
            diameter=round(s.outer_diameter, 1),
            spoke_deg=round(s.joint_strain["spoke"]["max"], 1),
            buckle=round(s.buckling["spoke_rod"], 1), rms=round(rms(s), 1),
            gap=round(max(v["gap"] for v in s.bond_force.values()), 3))
    fig.tight_layout()
    save(fig, "d_spoke_lock.png")
    V["spoke_lock"] = out

    wt = {}
    for pv in (True, False):
        s = solve_chain(Geometry(), spoke_pivot=pv)
        wt["free" if pv else "locked"] = dict(diameter=round(s.outer_diameter, 2),
                                              rms=round(rms(s), 2))
    V["spoke_lock_wt"] = wt


# ------------------------------------------------------------ band penalty
def band_penalty():
    """The piecewise-linear angle penalty, for the maths slide."""
    fig, ax = plt.subplots(figsize=(6.4, 3.7))
    d = np.linspace(0, 30, 600)
    for name, c in (("linker-A", CYAN), ("spoke", PURPLE)):
        ok, hard = JOINT_BANDS[name]
        ax.plot(d, angle_penalty(d, name), lw=2.4, color=c,
                label=f"{name}  (OK {ok:.0f}°, HARD {hard:.0f}°)")
        ax.axvline(ok, color=c, ls=":", lw=1.1, alpha=0.7)
        ax.axvline(hard, color=c, ls="--", lw=1.1, alpha=0.7)
    ax.set_xlabel("joint deviation  |δ|  (deg)", fontsize=12)
    ax.set_ylabel("penalty  Φ(δ)", fontsize=12)
    ax.set_title(f"slopes {W_FREE:.0f} → {W_HARD:.0f} → {W_FORBID:.0f}, "
                 "continuous at every break", fontsize=11, color=MUTE)
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)
    save(fig, "d_bands.png")


# --------------------------------------------------------------- registers
def registers():
    g = set_param(Geometry(), "spoke_rod", 28.03)
    t0 = time.time()
    sol = solve_chain(g, register_shift=True)
    V["register_seconds"] = round(time.time() - t0)
    scan = [{k: v for k, v in r.items() if k != "sol"} for r in sol.register_scan]
    V["register_scan"] = scan
    V["register_best3"] = best_registers(sol, 3)
    V["register_n_feasible"] = sum(1 for r in scan if r["feasible"])

    fig, ax = plt.subplots(figsize=(6.8, 4.3))
    for ok, c, lab, m in ((True, GREEN, "geometrically possible", "o"),
                          (False, BAD, "strand through a tubule wall", "x")):
        pts = [r for r in scan if r["feasible"] == ok]
        ax.scatter([r["outer"] for r in pts], [max(r["cost"], 1) for r in pts],
                   c=c, marker=m, s=64, label=lab, zorder=3)
    from matplotlib.transforms import blended_transform_factory
    ax.axvline(216, color="#111", ls="--", lw=1.8)
    ax.text(216, 0.97, "measured 216 nm", ha="center", va="top", fontsize=10.5,
            fontweight="bold", zorder=5,
            transform=blended_transform_factory(ax.transData, ax.transAxes),
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=2))
    ax.set_yscale("log")
    ax.set_xlabel("outer diameter (nm)", fontsize=12)
    ax.set_ylabel("model cost   (log)", fontsize=12)
    ax.set_title("25 protofilament registers, 28 nm spoke", fontsize=11, color=MUTE)
    ax.legend(fontsize=10, loc="lower left")
    ax.grid(alpha=0.3)
    save(fig, "d_register.png")


# ------------------------------------------------------------- predictions
def predictions():
    V["doublet"] = {}
    for mtn in (3, 2, 1):
        s = solve_chain(Geometry(MTn=mtn))
        V["doublet"][mtn] = dict(diameter=round(s.outer_diameter, 1), rms=round(rms(s), 1),
                                 clashes=s.n_clashes, reachable=bool(s.reachable))

    V["n_pf"] = {}
    for n in (9, 11, 13, 15, 18):
        s = solve_chain(set_param(Geometry(), "n_pf_A", n))
        V["n_pf"][n] = dict(diameter=round(s.outer_diameter, 1), rms=round(rms(s), 1),
                            clashes=s.n_clashes, reachable=bool(s.reachable),
                            clear=round(s.worst_strand_clearance, 2))

    V["spoke_len"] = {}
    for L in (35, 40, 45.03, 50, 55):
        s = solve_chain(set_param(Geometry(), "spoke_rod", L))
        V["spoke_len"][L] = dict(diameter=round(s.outer_diameter, 1), rms=round(rms(s), 1))

    V["bond_order"] = {k: round(v["force"], 4) for k, v in
                       sorted(solve_chain(Geometry(N_cw=8, N_mt=9)).bond_force.items(),
                              key=lambda kv: -kv[1]["force"])}


# ------------------------------------------------------------------ timing
def timing():
    out = {}
    for tag, g in (("wild type", Geometry()), ("8 cw / 9 mt", Geometry(N_cw=8, N_mt=9)),
                   ("n_pf = 9", set_param(Geometry(), "n_pf_A", 9)),
                   ("n_pf = 18", set_param(Geometry(), "n_pf_A", 18)),
                   ("spoke 28 nm", set_param(Geometry(), "spoke_rod", 28.03))):
        t0 = time.time()
        s = solve_chain(g)
        out[tag] = dict(seconds=round(time.time() - t0, 2), converged=bool(s.success))
    V["timing"] = out


if __name__ == "__main__":
    print("regenerating deck figures and values (chain solver)")
    wild_type()
    validation()
    mismatch()
    spoke_lock()
    band_penalty()
    predictions()
    timing()
    registers()
    (FIGS / "_values.json").write_text(json.dumps(V, indent=1, default=str))
    print(f"  -> figs/_values.json  ({len(V)} groups)")
