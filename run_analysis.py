#!/usr/bin/env python3
"""One-command entry point for the centriole geometry model.

Designed to be run without editing any code:

    python run_analysis.py                  # standard report
    python run_analysis.py --cartwheel 8    # 8-fold cartwheel, 9 triplets
    python run_analysis.py --sweep spoke_rod --from 35 --to 60 --steps 6

Everything is written to a `results/` folder: figures as PNG, tables as
CSV, and a plain-text summary.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

RESULTS = Path(__file__).resolve().parent / "results"

FRIENDLY = {
    "N_cw": "cartwheel symmetry (number of SAS-6 spokes)",
    "N_mt": "triplet symmetry (number of microtubule triplets)",
    "MTn": "tubules per blade (3 = triplet, 2 = doublet, 1 = singlet)",
    "spoke_rod": "SAS-6 coiled-coil length (nm)",
    "base_length": "triplet base length (nm)",
    "pinhead_span": "pinhead length (nm)",
    "linker_length": "A-C linker length (nm)",
    "n_pf_A": "protofilaments in the A-tubule",
}


def _need(pkg):
    try:
        __import__(pkg)
    except ImportError:
        sys.exit(
            f"\nMissing package '{pkg}'.\n"
            "Install everything with:\n\n"
            "    pip install -r requirements.txt\n"
        )


for _p in ("numpy", "scipy", "matplotlib", "pandas"):
    _need(_p)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from centriole_kinematic import (  # noqa: E402
    Geometry, draw, set_param, summarise, sweep,
)
# Use the chain formulation, as the app does. The older spring-network solve()
# lets its bonds open by 1-3 nm under a large perturbation -- at a 28 nm SAS-6
# spoke the pinhead visibly separates from the spoke -- which corrupts exactly
# the mutant this CLI exists to model. See docs/blade_rotation_discrepancy.md.
from centriole_chain import solve_chain as solve  # noqa: E402


def _banner(msg):
    print(f"\n{'=' * 62}\n{msg}\n{'=' * 62}", flush=True)


def _rms(sol):
    return float(np.sqrt(np.mean([v["rms"] ** 2 for v in sol.joint_strain.values()])))


def single(cartwheel, triplets, mtn, spoke_pivot=True, **over):
    """Solve and render one configuration."""
    g = Geometry(N_cw=cartwheel, N_mt=triplets, MTn=mtn)
    for k, v in over.items():
        if v is not None:
            g = set_param(g, k, v)
    sol = solve(g, spoke_pivot=spoke_pivot)

    fig, ax = plt.subplots(figsize=(8, 8))
    draw(sol, ax=ax, show_pf_labels=True)
    name = f"cw{cartwheel}_mt{triplets}_MTn{mtn}"
    fig.savefig(RESULTS / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(sol.report())
    print(f"\n  -> figure: results/{name}.png")
    return sol


def standard_report():
    """The default run: wild type, symmetry variants, and two length sweeps."""
    _banner("1/4  Wild type (9-fold cartwheel, 9 triplets, triplets)")
    wt = single(9, 9, 3)

    _banner("2/4  Symmetry variants")
    cases = [(9, 9, "wild type"), (8, 9, "8-fold cartwheel / 9 triplets"),
             (8, 8, "fully 8-fold"), (9, 8, "9-fold cartwheel / 8 triplets")]
    fig, axes = plt.subplots(1, 4, figsize=(26, 7))
    rows = []
    for ax, (ncw, nmt, label) in zip(axes, cases):
        sol = solve(Geometry(N_cw=ncw, N_mt=nmt))
        un = f"  unattached {sol.unattached_triplets}" if sol.unattached_triplets else ""
        draw(sol, ax=ax, title=f"{label}\nN_cw={ncw} N_mt={nmt}  "
                               f"diam {sol.outer_diameter:.0f} nm  "
                               f"spoke {sol.joint_bands['spoke']}{un}")
        rec = {"case": label}
        rec.update(summarise(sol))
        rows.append(rec)
        print(f"  {label:<32} diameter {sol.outer_diameter:6.1f} nm   "
              f"joint rms {_rms(sol):5.2f} deg   spoke {sol.joint_bands['spoke']}", flush=True)
    fig.tight_layout()
    fig.savefig(RESULTS / "symmetry_variants.png", dpi=120, bbox_inches="tight")
    plt.close(fig)

    import pandas as pd
    pd.DataFrame(rows).to_csv(RESULTS / "symmetry_variants.csv", index=False)
    print("\n  -> figure: results/symmetry_variants.png")
    print("  -> table : results/symmetry_variants.csv")

    _banner("3/4  Symmetry sweep (cartwheel and triplets together)")
    rows = []
    for n in range(7, 12):
        sol = solve(Geometry(N_cw=n, N_mt=n))
        rows.append({"symmetry": n, **summarise(sol)})
        print(f"  {n}-fold: diameter {sol.outer_diameter:6.1f} nm   "
              f"joint rms {_rms(sol):5.2f} deg", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "symmetry_sweep.csv", index=False)
    _plot_sweep(df, "symmetry", "symmetry (N-fold)", "symmetry_sweep.png")

    _banner("4/4  SAS-6 coiled-coil length sweep")
    df = sweep("spoke_rod", [35, 40, 45.03, 50, 55, 60], solver=solve)
    df.to_csv(RESULTS / "spoke_length_sweep.csv", index=False)
    _plot_sweep(df, "spoke_rod", "SAS-6 coiled-coil length (nm)", "spoke_length_sweep.png")
    for _, r in df.iterrows():
        print(f"  {r['spoke_rod']:5.1f} nm: diameter {r['diameter_nm']:6.1f} nm   "
              f"joint rms {r['joint_rms_deg']:5.2f} deg", flush=True)

    _summary_file(wt)


def _plot_sweep(df, xcol, xlabel, fname):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    axes[0].plot(df[xcol], df["joint_rms_deg"], "o-", color="#9e42f1")
    axes[0].set_ylabel("joint strain (deg, RMS)")
    axes[1].plot(df[xcol], df["diameter_nm"], "o-", color="#688950")
    axes[1].set_ylabel("centriole diameter (nm)")
    axes[2].plot(df[xcol], df["n_clashes"], "o-", color="#c0392b")
    axes[2].set_ylabel("microtubule clashes")
    for ax in axes:
        ax.set_xlabel(xlabel)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(RESULTS / fname, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  -> figure: results/{fname}")


def _summary_file(wt):
    lines = [
        "CENTRIOLE GEOMETRY MODEL - SUMMARY",
        "=" * 62, "",
        "Wild type (9-fold cartwheel, 9 triplets, microtubule triplets):", "",
        wt.report(), "",
        "-" * 62,
        "Files written to results/:",
        "  cw9_mt9_MTn3.png        wild-type cross-section",
        "  symmetry_variants.png   wild type vs symmetry mutants",
        "  symmetry_variants.csv   the same, as numbers",
        "  symmetry_sweep.png/csv  7- to 11-fold symmetry",
        "  spoke_length_sweep.*    effect of SAS-6 coiled-coil length",
        "",
        "How to read the output:",
        "  joint rotation   how far each connection had to turn from the rest",
        "                   angle of the protofilament it binds, graded OK /",
        "                   HARD / SEVERE against that joint's own limits.",
        "                   These are assumed tolerances, NOT verdicts on",
        "                   whether a structure can form.",
        "  bond load        which connection carries the most strain, i.e.",
        "                   which would break first.",
        "  clashes          microtubules overlapping in space (impossible).",
        "  clearance        how close a strand passes to a microtubule,",
        "                   counting its own thickness; negative means it is",
        "                   inside a tubule wall.",
        "  approach         how squarely a strand meets the protofilament it",
        "                   binds (+1 head-on from outside, below 0 impossible).",
    ]
    (RESULTS / "summary.txt").write_text("\n".join(lines))
    print("\n  -> summary: results/summary.txt")


def main():
    ap = argparse.ArgumentParser(
        description="Model centriole cross-section geometry and test perturbations.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  python run_analysis.py\n"
               "  python run_analysis.py --cartwheel 8\n"
               "  python run_analysis.py --sweep spoke_rod --from 35 --to 60 --steps 6\n",
    )
    ap.add_argument("--cartwheel", type=int, help="number of SAS-6 spokes (default 9)")
    ap.add_argument("--triplets", type=int, help="number of microtubule triplets (default 9)")
    ap.add_argument("--tubules", type=int, choices=(1, 2, 3),
                    help="3=triplet (default), 2=doublet, 1=singlet")
    ap.add_argument("--spoke-length", type=float, help="SAS-6 coiled-coil length in nm")
    ap.add_argument("--base-length", type=float, help="triplet base length in nm")
    ap.add_argument("--pinhead-length", type=float, help="pinhead length in nm")
    ap.add_argument("--linker-length", type=float, help="A-C linker length in nm")
    ap.add_argument("--protofilaments", type=int, help="protofilaments in the A-tubule")
    ap.add_argument("--lock-spoke", action="store_true",
                    help="hold each SAS-6 spoke on its own radius, so it can only "
                         "strain outwards and cannot bend at the head")
    ap.add_argument("--sweep", metavar="PARAM", choices=sorted(FRIENDLY),
                    help="vary one parameter: " + ", ".join(sorted(FRIENDLY)))
    ap.add_argument("--from", dest="lo", type=float, help="sweep start value")
    ap.add_argument("--to", dest="hi", type=float, help="sweep end value")
    ap.add_argument("--steps", type=int, default=6, help="number of sweep values (default 6)")
    args = ap.parse_args()

    RESULTS.mkdir(exist_ok=True)
    print("Centriole geometry model")
    print(f"Results will be written to: {RESULTS}")

    if args.sweep:
        if args.lo is None or args.hi is None:
            sys.exit("\n--sweep also needs --from and --to, e.g.\n"
                     "    python run_analysis.py --sweep spoke_rod --from 35 --to 60\n")
        vals = np.linspace(args.lo, args.hi, args.steps)
        if args.sweep in ("N_cw", "N_mt", "MTn", "n_pf_A"):
            vals = sorted({int(round(v)) for v in vals})
        _banner(f"Sweeping {FRIENDLY[args.sweep]}")
        df = sweep(args.sweep, vals, solver=solve,
                   spoke_pivot=not args.lock_spoke)
        out = RESULTS / f"sweep_{args.sweep}.csv"
        df.to_csv(out, index=False)
        cols = ["diameter_nm", "joint_rms_deg", "max_buckle_pct", "n_clashes", "worst_bond"]
        print(df[[args.sweep] + cols].to_string(index=False))
        _plot_sweep(df, args.sweep, FRIENDLY[args.sweep], f"sweep_{args.sweep}.png")
        print(f"  -> table : results/{out.name}")
        return

    custom = args.lock_spoke or any(v is not None for v in
                 (args.cartwheel, args.triplets, args.tubules, args.spoke_length,
                  args.base_length, args.pinhead_length, args.linker_length,
                  args.protofilaments))
    if custom:
        _banner("Single configuration")
        single(args.cartwheel or 9, args.triplets or 9, args.tubules or 3,
               spoke_pivot=not args.lock_spoke,
               spoke_rod=args.spoke_length, base_length=args.base_length,
               pinhead_span=args.pinhead_length, linker_length=args.linker_length,
               n_pf_A=args.protofilaments)
    else:
        standard_report()

    print("\nDone.")


if __name__ == "__main__":
    main()
