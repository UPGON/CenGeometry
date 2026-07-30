"""CenGeometry — interactive centriole geometry explorer.

Launch by double-clicking the launcher for your platform, or from a
terminal with:

    streamlit run app.py
"""

from __future__ import annotations

import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import streamlit as st

from centriole_kinematic import (
    BOND_STRENGTH, JOINT_BANDS, Geometry, _linker_span, draw, set_param,
    summarise,
)
# The chain solver holds connection points together by construction, so parts
# never visibly separate under load. The older spring-based solve() opened the
# spoke-to-pinhead connection by ~2.8 nm at a 28 nm spoke, which is wrong at
# any load, so the app uses the chain form throughout.
from centriole_chain import solve_chain as solve

st.set_page_config(page_title="CenGeometry", page_icon="🔬", layout="wide")

DEFAULTS = Geometry()
LINKER_DEFAULT = round(_linker_span(DEFAULTS), 2)

JOINTS = list(JOINT_BANDS)
BONDS = list(BOND_STRENGTH)
BUCKLE = ["spoke_rod", "pinhead", "base", "linker_armC", "linker_armA"]

#: parameters offered for scanning: key -> (label, is_integer, lo, hi)
SCANNABLE = {
    "N_both": ("Symmetry — cartwheel and triplets together", True, 6, 11),
    "N_cw": ("Cartwheel symmetry only (triplets stay fixed)", True, 6, 11),
    "N_mt": ("Triplet symmetry only (cartwheel stays fixed)", True, 6, 11),
    "MTn": ("Tubules per blade", True, 1, 3),
    "spoke_rod": ("SAS-6 coiled-coil length (nm)", False, 10.0, 90.0),
    "base_length": ("Triplet base length (nm)", False, 10.0, 70.0),
    "pinhead_span": ("Pinhead length (nm)", False, 5.0, 50.0),
    "linker_length": ("A-C linker length (nm)", False, 5.0, 60.0),
    "head_contact": ("SAS-6 head-head spacing (nm)", False, 3.0, 30.0),
    "n_pf_A": ("A-tubule protofilaments", True, 9, 18),
}

BAND_COLOUR = {"OK": "#2e7d32", "HARD": "#e08a00", "SEVERE": "#c62828", "-": "#888888"}


#: which slot(s) of the PARAMS tuple a scanned parameter writes to
PARAM_SLOTS = {"N_both": (0, 1), "N_cw": (0,), "N_mt": (1,), "MTn": (2,),
               "spoke_rod": (3,), "base_length": (4,), "pinhead_span": (5,),
               "linker_length": (6,), "head_contact": (7,), "n_pf_A": (8,)}


# --------------------------------------------------------------------------
def make_geometry(p: tuple) -> Geometry:
    (N_cw, N_mt, MTn, spoke_rod, base_length, pinhead_span,
     linker_length, head_contact, n_pf_A) = p
    g = Geometry(N_cw=int(N_cw), N_mt=int(N_mt), MTn=int(MTn),
                 spoke_rod=float(spoke_rod), base_length=float(base_length),
                 pinhead_span=float(pinhead_span), head_contact=float(head_contact))
    g = set_param(g, "n_pf_A", int(n_pf_A))
    return set_param(g, "linker_length", float(linker_length))


WT_PARAMS = (9, 9, 3, float(DEFAULTS.spoke_rod), float(DEFAULTS.base_length),
             float(DEFAULTS.pinhead_span), LINKER_DEFAULT,
             float(DEFAULTS.head_contact), int(DEFAULTS.n_pf["A"]))


@st.cache_data(show_spinner=False)
def wt_metrics() -> dict:
    """Wild-type values, drawn as a dotted reference on every scan plot."""
    return summarise(solve(make_geometry(WT_PARAMS)))


@st.cache_data(show_spinner=False, max_entries=128)
def thumbnail(p: tuple, register_shift: bool = False) -> bytes:
    """Small cross-section render, for the strip above the scan plots."""
    sol = solve(make_geometry(p))
    fig, ax = plt.subplots(figsize=(3.1, 3.1))
    draw(sol, ax=ax, title="")
    ax.set_xlabel("")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=95, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


@st.cache_data(show_spinner=False, max_entries=256)
def run_one(p: tuple, register_shift: bool):
    """Solve one configuration. Returns (metrics, report text, figure PNG)."""
    sol = solve(make_geometry(p))
    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    draw(sol, ax=ax)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    return summarise(sol), sol.report(), buf.getvalue()


@st.cache_data(show_spinner=False, max_entries=64)
def run_scan(param: str, values: tuple, p: tuple, register_shift: bool):
    base = make_geometry(p)
    rows = []
    bar = st.progress(0.0, text=f"Solving {len(values)} configurations…")
    for i, v in enumerate(values):
        g = set_param(base, param, v)
        rec = {param: v}
        rec.update(summarise(solve(g)))
        rows.append(rec)
        bar.progress((i + 1) / len(values), text=f"Solved {i+1} of {len(values)}")
    bar.empty()
    return pd.DataFrame(rows)


@st.cache_data(show_spinner=False, max_entries=16)
def run_grid(px: str, xvals: tuple, py: str, yvals: tuple, p: tuple,
             register_shift: bool = False):
    """Cross any two parameters and tabulate every metric for each cell."""
    base = make_geometry(p)
    rows = []
    total = len(xvals) * len(yvals)
    bar = st.progress(0.0, text=f"Solving {total} configurations…")
    k = 0
    for yv in yvals:
        for xv in xvals:
            g = set_param(set_param(base, py, yv), px, xv)
            rec = {px: xv, py: yv}
            rec.update(summarise(solve(g)))
            rows.append(rec)
            k += 1
            bar.progress(k / total, text=f"Solved {k} of {total}")
    bar.empty()
    return pd.DataFrame(rows)


#: metrics offered for 2-parameter display, grouped label -> column
GRID_METRICS = {
    "Outer diameter (nm)": "diameter_nm",
    "A-tubule ring (nm)": "A_ring_nm",
    "Lumen (nm)": "lumen_nm",
    "Triplet tilt (deg)": "triplet_tilt_deg",
    "Joint strain, RMS (deg)": "joint_rms_deg",
    "Max buckling (%)": "max_buckle_pct",
    "MT-MT clashes": "n_clashes",
    "Worst overlap (nm)": "max_overlap_nm",
    "Unattached triplets": "n_unattached",
    "A-C linker clearance (nm)": "linker_clear_nm",
    "Triplet base clearance (nm)": "base_clear_nm",
    "SAS-6 spoke clearance (nm)": "spoke_clear_nm",
    **{f"Joint rotation — {j} (deg)": f"{j}_deg" for j in JOINTS},
    **{f"Bond load — {b}": f"bond_{b}" for b in BONDS},
    **{f"Buckling — {b} (%)": f"buckle_{b}_pct" for b in BUCKLE},
}


def heatmap(df, px, py, col, title, cmap="viridis"):
    piv = df.pivot(index=py, columns=px, values=col)
    fig, ax = plt.subplots(figsize=(5.8, 4.7))
    im = ax.imshow(piv.values, cmap=cmap, origin="lower", aspect="auto")
    ax.set_xticks(range(len(piv.columns)), [f"{v:g}" for v in piv.columns], fontsize=8)
    ax.set_yticks(range(len(piv.index)), [f"{v:g}" for v in piv.index], fontsize=8)
    ax.set_xlabel(SCANNABLE[px][0].split("(")[0].strip(), fontsize=9)
    ax.set_ylabel(SCANNABLE[py][0].split("(")[0].strip(), fontsize=9)
    ax.set_title(title, fontsize=10)
    span = np.nanmax(piv.values) - np.nanmin(piv.values)
    fmt = "{:.0f}" if span > 12 else "{:.2f}"
    if piv.size <= 90:
        for r in range(piv.shape[0]):
            for q in range(piv.shape[1]):
                v = piv.values[r, q]
                if np.isfinite(v):
                    ax.text(q, r, fmt.format(v), ha="center", va="center",
                            color="w", fontsize=6.5)
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    return fig


def family(df, px, py, col, title, ylabel, wt=None):
    """One curve per value of the second parameter -- often easier to read
    than a heatmap when the x axis is continuous."""
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for yv, sub in df.groupby(py):
        sub = sub.sort_values(px)
        ax.plot(sub[px], sub[col], "o-", lw=1.5, ms=3.5,
                label=f"{SCANNABLE[py][0].split('(')[0].strip()} = {yv:g}")
    if wt is not None and col in wt and isinstance(wt[col], (int, float)):
        ax.axhline(wt[col], color="#888", ls=":", lw=1.2, label="wild type")
    ax.set_xlabel(SCANNABLE[px][0])
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    return fig


def axis_values(param, lo, hi, n):
    label, is_int, _, _ = SCANNABLE[param]
    if is_int:
        return list(range(int(round(lo)), int(round(hi)) + 1))
    return [round(float(v), 4) for v in np.linspace(lo, hi, int(n))]


def badge(label: str, grade: str) -> str:
    return (f"<span style='background:{BAND_COLOUR[grade]};color:#fff;padding:2px 9px;"
            f"border-radius:12px;font-size:.78rem;font-weight:700'>{label} {grade}</span>")


def line_panel(df, xcol, xlabel, series, ylabel, title, wt=None):
    """One panel. `wt` adds a faint dotted line at each metric's wild-type
    value, so any departure from wild type is readable at a glance."""
    fig, ax = plt.subplots(figsize=(7, 4))
    for col, lab in series:
        if col not in df:
            continue
        line, = ax.plot(df[xcol], df[col], "o-", label=lab, lw=1.6, ms=4)
        if wt is not None and col in wt and isinstance(wt[col], (int, float)):
            ax.axhline(wt[col], color=line.get_color(), ls=":", lw=1.1, alpha=0.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.3)
    handles, labels = ax.get_legend_handles_labels()
    if wt is not None:
        handles.append(Line2D([0], [0], color="#888", ls=":", lw=1.1))
        labels.append("wild type")
    if len(handles) > 1:
        ax.legend(handles, labels, fontsize=8)
    fig.tight_layout()
    return fig


# --------------------------------------------------------------------------
with st.sidebar:
    st.title("🔬 CenGeometry")
    st.caption("Centriole cross-section geometry and mechanics. All lengths in nm.")

    st.subheader("Symmetry")
    N_cw = st.slider("Cartwheel spokes (SAS-6)", 5, 15, 9, key="N_cw",
                     help="Number of SAS-6 dimers / spokes in the cartwheel. "
                          "Values above ~12 are slow to solve.")
    N_mt = st.slider("Microtubule triplets", 5, 15, 9, key="N_mt",
                     help="Number of microtubule blades. May differ from the "
                          "cartwheel — the triplet ring keeps its spacing and the "
                          "cartwheel adapts. Values above ~12 are slow to solve.")
    MTn = st.radio("Tubules per blade", [3, 2, 1], horizontal=True, key="MTn",
                   format_func=lambda v: {3: "Triplet", 2: "Doublet", 1: "Singlet"}[v])

    st.subheader("Distances (nm)")
    st.caption("Type any value.")
    spoke_rod = st.number_input("SAS-6 coiled coil", key="spoke_rod",
                                min_value=5.0, max_value=200.0,
                                value=float(DEFAULTS.spoke_rod), step=0.5, format="%.2f")
    base_length = st.number_input("Triplet base", key="base_length",
                                  min_value=5.0, max_value=150.0,
                                  value=float(DEFAULTS.base_length), step=0.5, format="%.2f")
    pinhead_span = st.number_input("Pinhead", key="pinhead_span",
                                   min_value=2.0, max_value=100.0,
                                   value=float(DEFAULTS.pinhead_span), step=0.5, format="%.2f")
    linker_length = st.number_input("A-C linker (end to end)", key="linker_length",
                                    min_value=2.0, max_value=100.0,
                                    value=LINKER_DEFAULT, step=0.5, format="%.2f")
    head_contact = st.number_input("SAS-6 head-head spacing", key="head_contact",
                                   min_value=1.0, max_value=60.0,
                                   value=float(DEFAULTS.head_contact), step=0.2, format="%.2f")

    st.subheader("Microtubule")
    n_pf_A = st.number_input("A-tubule protofilaments", key="n_pf_A",
                             min_value=8, max_value=20,
                             value=int(DEFAULTS.n_pf["A"]), step=1,
                             help="Also sets the tubule radius: R = n · 5.747 / 2π nm.")

    register_shift = st.checkbox(
        "Allow protofilament register shift", key="register_shift", value=False,
        disabled=True,
        help="Not available with the current solver, which holds connection "
             "points together by construction. The register search existed only "
             "in the older spring-based solver, which allowed parts to separate.")

    def _reset_to_wt():
        # Must run as an on_click callback: it fires BEFORE the script re-runs,
        # so the widgets are rebuilt from their defaults. Popping the keys
        # inside the button body instead left the old values on screen until
        # the page was manually refreshed.
        for k in ("N_cw", "N_mt", "MTn", "spoke_rod", "base_length", "pinhead_span",
                  "linker_length", "head_contact", "n_pf_A", "register_shift",
                  "scan_done", "grid_done"):
            st.session_state.pop(k, None)

    st.button("Reset to wild type", width='stretch', on_click=_reset_to_wt)

PARAMS = (N_cw, N_mt, MTn, spoke_rod, base_length, pinhead_span,
          linker_length, head_contact, n_pf_A)

tab_x, tab_scan, tab_grid, tab_help = st.tabs(
    ["Cross-section", "Parameter scan", "2-parameter scan", "How to read this"])

# ---------------------------------------------------------------- cross-section
with tab_x:
    with st.spinner("Solving…"):
        m, report, png = run_one(PARAMS, register_shift)

    if not m["converged"]:
        st.warning("The solver did not fully converge at these values — the geometry "
                   "shown is its best attempt, so read the numbers with caution.")
    left, right = st.columns([1, 1.25])
    with left:
        st.image(png, width='stretch')
    with right:
        st.subheader("Size")
        c = st.columns(3)
        c[0].metric("Outer Ø", f"{m['diameter_nm']:.0f} nm",
                    delta=f"{m['diameter_nm'] - 254.99:+.0f} vs WT")
        c[1].metric("A ring", f"{m['A_ring_nm']:.0f} nm")
        c[2].metric("Lumen", f"{m['lumen_nm']:.0f} nm")

        st.subheader("Strain")
        c = st.columns(3)
        c[0].metric("Joint strain", f"{m['joint_rms_deg']:.1f}°")
        c[1].metric("Triplet tilt", f"{m['triplet_tilt_deg']:.0f}°")
        c[2].metric("Buckling", f"{m['max_buckle_pct']:.1f} %")

        st.markdown("**Joint grades**", help="How far each connection has rotated "
                                            "from wild type, against its own limits.")
        st.markdown(" ".join(badge(j, m[f"{j}_band"]) for j in JOINTS),
                    unsafe_allow_html=True)
        st.caption("These are **assumed tolerances, not verdicts on feasibility** — "
                   "the limits were reasoned rather than measured, and structures are "
                   "known to assemble in conditions graded SEVERE. Read them as "
                   "*how far from wild type*, and judge from the geometry and numbers.")

        st.subheader("Integrity")
        c = st.columns(2)
        c[0].metric("MT clashes", m["n_clashes"],
                    delta="impossible" if m["n_clashes"] else "none",
                    delta_color="inverse" if m["n_clashes"] else "off")
        c[1].metric("Unattached", m["n_unattached"])

        st.markdown(f"Closest to rupture: **{m['worst_bond']}** "
                    f"(load {m['worst_bond_force']:.3f})")

        if register_shift:
            rp, rc = m["reg_pinhead_pf"], m["reg_linkerC_pf"]
            if rp == 0 and rc == 0:
                st.info("**Register shift: none chosen.** Staying on the wild-type "
                        "protofilaments was already the lowest-strain option, so "
                        "re-registering would not help here.")
            else:
                st.info(f"**Register shift chosen:** pinhead {rp:+.0f} protofilament(s) "
                        f"on the A-tubule, linker {rc:+.0f} on the C-tubule. "
                        "The structure relieves strain by binding a different "
                        "protofilament rather than by deforming.")
        st.caption(f"Clearance vs microtubules — linker {m['linker_clear_nm']} nm, "
                   f"base {m['base_clear_nm']} nm, spoke {m['spoke_clear_nm']} nm. "
                   f"Negative means passing through one.")

    st.divider()
    a, b = st.columns(2)
    with a:
        st.markdown("**Bond load** (higher = closer to rupture)")
        bl = pd.DataFrame({"bond": BONDS,
                           "load": [m.get(f"bond_{k}", 0.0) for k in BONDS]}
                          ).sort_values("load", ascending=False)
        st.bar_chart(bl, x="bond", y="load", height=240)
    with b:
        st.markdown("**Joint rotation from wild type** (degrees)")
        jd = pd.DataFrame({"joint": JOINTS,
                           "degrees": [abs(m[f"{j}_deg"]) for j in JOINTS]})
        st.bar_chart(jd, x="joint", y="degrees", height=240)

    with st.expander("Full text report"):
        st.code(report, language=None)

# ------------------------------------------------------------------------ scan
with tab_scan:
    st.markdown("Vary one parameter and plot every metric against it. "
                "All other parameters keep their sidebar values.")
    c = st.columns([2.2, 1, 1, 1, 1])
    param = c[0].selectbox("Parameter to vary", list(SCANNABLE),
                           format_func=lambda k: SCANNABLE[k][0])
    label, is_int, lo_lim, hi_lim = SCANNABLE[param]
    fmt = "%d" if is_int else "%.2f"
    step = 1.0 if is_int else 0.5
    lo = c[1].number_input("From", float(lo_lim), float(hi_lim), float(lo_lim),
                           step=step, format=fmt)
    hi = c[2].number_input("To", float(lo_lim), float(hi_lim), float(hi_lim),
                           step=step, format=fmt)
    if is_int:
        # whole-number parameters step one at a time, so the count follows
        # from the range -- for symmetry that is one point per fold
        n = max(int(round(hi)) - int(round(lo)) + 1, 1)
        c[3].number_input("Steps", value=n, disabled=True,
                          help="Set automatically: one point per whole value.")
    else:
        n = int(c[3].number_input("Steps", 2, 40, 8, 1))
    go = c[4].button("Run scan", type="primary", width='stretch')

    if go or st.session_state.get("scan_done"):
        st.session_state["scan_done"] = True
        if is_int:
            vals = list(range(int(round(lo)), int(round(hi)) + 1))
        else:
            vals = [round(float(v), 4) for v in np.linspace(lo, hi, n)]
        df = run_scan(param, tuple(vals), PARAMS, register_shift)
        wt = wt_metrics()

        bad = int((~df["converged"]).sum())
        msg = f"Solved {len(df)} configurations."
        (st.warning if bad else st.success)(
            msg + (f" {bad} did not converge — treat those points with caution."
                   if bad else " All converged."))

        # ---- cross-sections above the plots
        st.markdown("**Cross-sections**")
        shown = vals if len(vals) <= 8 else [vals[i] for i in
                                            np.linspace(0, len(vals) - 1, 8).astype(int)]
        if len(shown) < len(vals):
            st.caption(f"Showing {len(shown)} of {len(vals)} evenly spaced; "
                       f"the plots below use all {len(vals)}.")
        with st.spinner("Rendering cross-sections…"):
            thumbs = []
            for v in shown:
                q = list(PARAMS)
                for idx in PARAM_SLOTS[param]:
                    q[idx] = v
                thumbs.append((v, thumbnail(tuple(q), register_shift)))
        for row_start in range(0, len(thumbs), 4):
            cols = st.columns(4)
            for col, (v, png) in zip(cols, thumbs[row_start:row_start + 4]):
                vtxt = f"{int(v)}" if is_int else f"{v:g}"
                d = df.loc[df[param] == v, "diameter_nm"]
                col.image(png, width='stretch')
                col.caption(f"{label.split('(')[0].strip()} = **{vtxt}**"
                            + (f" · Ø {d.iloc[0]:.0f} nm" if len(d) else ""))
        st.divider()

        g1, g2 = st.columns(2)
        with g1:
            st.pyplot(line_panel(df, param, label,
                                 [("diameter_nm", "outer diameter"),
                                  ("A_ring_nm", "A-tubule ring"),
                                  ("lumen_nm", "lumen")],
                                 "nm", "Size", wt))
            st.pyplot(line_panel(df, param, label,
                                 [(f"{j}_deg", j) for j in JOINTS],
                                 "degrees from wild type", "Joint rotation", wt))
            st.pyplot(line_panel(df, param, label,
                                 [(f"buckle_{b}_pct", b) for b in BUCKLE],
                                 "% of contour lost", "Buckling", wt))
        with g2:
            st.pyplot(line_panel(df, param, label,
                                 [("n_clashes", "MT-MT clashes"),
                                  ("max_overlap_nm", "worst overlap (nm)"),
                                  ("n_unattached", "unattached triplets")],
                                 "count / nm", "Steric integrity", wt))
            st.pyplot(line_panel(df, param, label,
                                 [(f"bond_{k}", k) for k in BONDS],
                                 "load (gap × strength)", "Bond load", wt))
            st.pyplot(line_panel(df, param, label,
                                 [("linker_clear_nm", "A-C linker"),
                                  ("base_clear_nm", "triplet base"),
                                  ("spoke_clear_nm", "SAS-6 spoke")],
                                 "nm (negative = clashing)",
                                 "Strand clearance vs microtubules", wt))
        st.pyplot(line_panel(df, param, label,
                             [("joint_rms_deg", "joint strain RMS"),
                              ("triplet_tilt_deg", "triplet tilt from radial")],
                             "degrees", "Overall strain and tilt", wt))

        if register_shift:
            st.pyplot(line_panel(df, param, label,
                                 [("reg_pinhead_pf", "pinhead on A-tubule"),
                                  ("reg_linkerC_pf", "linker on C-tubule")],
                                 "protofilaments shifted",
                                 "Protofilament register chosen"))
            picked = df[(df["reg_pinhead_pf"] != 0) | (df["reg_linkerC_pf"] != 0)]
            if len(picked):
                st.caption(
                    "Register shifted away from wild type at: "
                    + ", ".join(f"{label.split('(')[0].strip()}={r[param]:g} "
                                f"(pinhead {r['reg_pinhead_pf']:+.0f}, "
                                f"linker {r['reg_linkerC_pf']:+.0f})"
                                for _, r in picked.iterrows()))
            else:
                st.caption("No configuration preferred a shifted register — "
                           "wild-type protofilaments were always lowest-strain.")

        with st.expander("Table of all metrics"):
            st.dataframe(df, width='stretch', height=320)
        st.download_button("Download CSV", df.to_csv(index=False).encode(),
                           file_name=f"scan_{param}.csv", mime="text/csv")

# ------------------------------------------------------------------------ grid
with tab_grid:
    st.markdown("Cross **any two** parameters — e.g. SAS-6 coiled-coil length "
                "against symmetry. Everything else keeps its sidebar values.")

    c = st.columns([2.1, 1, 1, 0.85])
    px = c[0].selectbox("First parameter (x axis)", list(SCANNABLE), index=0,
                        format_func=lambda k: SCANNABLE[k][0], key="gx_param")
    xlab, x_int, xlo_l, xhi_l = SCANNABLE[px]
    xlo = c[1].number_input("From", float(xlo_l), float(xhi_l), float(xlo_l),
                            step=1.0 if x_int else 0.5,
                            format="%d" if x_int else "%.2f", key="gx_lo")
    xhi = c[2].number_input("To", float(xlo_l), float(xhi_l), float(xhi_l),
                            step=1.0 if x_int else 0.5,
                            format="%d" if x_int else "%.2f", key="gx_hi")
    if x_int:
        nx = max(int(round(xhi)) - int(round(xlo)) + 1, 1)
        c[3].number_input("Steps", value=nx, disabled=True, key="gx_n_ro")
    else:
        nx = int(c[3].number_input("Steps", 2, 20, 6, 1, key="gx_n"))

    others = [k for k in SCANNABLE if k != px]
    default_y = "spoke_rod" if px != "spoke_rod" else "N_both"
    c = st.columns([2.1, 1, 1, 0.85])
    py = c[0].selectbox("Second parameter (y axis)", others,
                        index=others.index(default_y),
                        format_func=lambda k: SCANNABLE[k][0], key="gy_param")
    ylab, y_int, ylo_l, yhi_l = SCANNABLE[py]
    ylo = c[1].number_input("From", float(ylo_l), float(yhi_l), float(ylo_l),
                            step=1.0 if y_int else 0.5,
                            format="%d" if y_int else "%.2f", key="gy_lo")
    yhi = c[2].number_input("To", float(ylo_l), float(yhi_l), float(yhi_l),
                            step=1.0 if y_int else 0.5,
                            format="%d" if y_int else "%.2f", key="gy_hi")
    if y_int:
        ny = max(int(round(yhi)) - int(round(ylo)) + 1, 1)
        c[3].number_input("Steps", value=ny, disabled=True, key="gy_n_ro")
    else:
        ny = int(c[3].number_input("Steps", 2, 20, 5, 1, key="gy_n"))

    picked = st.multiselect(
        "Metrics to display", list(GRID_METRICS),
        default=["Outer diameter (nm)", "Joint strain, RMS (deg)",
                 "MT-MT clashes", "Max buckling (%)"],
        key="grid_metrics")
    view = st.radio("View", ["Heatmaps", "Curves", "Both"], horizontal=True,
                    index=2, key="grid_view")

    xv, yv = axis_values(px, xlo, xhi, nx), axis_values(py, ylo, yhi, ny)
    total = len(xv) * len(yv)
    overlap = set(PARAM_SLOTS[px]) & set(PARAM_SLOTS[py])

    if overlap:
        st.error(f"“{xlab}” and “{ylab}” both control the same thing, so they "
                 "cannot be crossed. Pick a different pair — for the "
                 "cartwheel-versus-triplet mismatch grid, choose the two "
                 "single-symmetry options rather than the combined one.")
    elif total > 120:
        st.error(f"{len(xv)} × {len(yv)} = **{total}** configurations is too many "
                 "(the cap is 120, roughly four minutes). Narrow a range or "
                 "reduce the steps.")
    else:
        st.caption(f"{len(xv)} × {len(yv)} = **{total}** configurations, "
                   f"roughly {total * 2 // 60} min {total * 2 % 60} s.")
        if st.button("Run 2-parameter scan", type="primary", key="grid_go") \
                or st.session_state.get("grid_done"):
            st.session_state["grid_done"] = True
            gdf = run_grid(px, tuple(xv), py, tuple(yv), PARAMS, register_shift)
            wt = wt_metrics()

            bad = int((~gdf["converged"]).sum())
            (st.warning if bad else st.success)(
                f"Solved {len(gdf)} configurations."
                + (f" {bad} did not converge — treat those cells with caution."
                   if bad else " All converged."))

            if not picked:
                st.info("Choose at least one metric above to see plots.")
            if view in ("Heatmaps", "Both"):
                st.markdown("#### Heatmaps")
                cols = st.columns(2)
                for i, name in enumerate(picked):
                    col = GRID_METRICS[name]
                    cmap = ("Reds" if "clash" in col or "overlap" in col
                            else "magma_r" if "strain" in col or "_deg" in col
                            else "viridis")
                    cols[i % 2].pyplot(heatmap(gdf, px, py, col, name, cmap))
            if view in ("Curves", "Both"):
                st.markdown("#### Curves — one line per "
                            f"{ylab.split('(')[0].strip().lower()}")
                cols = st.columns(2)
                for i, name in enumerate(picked):
                    cols[i % 2].pyplot(
                        family(gdf, px, py, GRID_METRICS[name], name, name, wt))

            with st.expander("Table of all metrics"):
                st.dataframe(gdf, width='stretch', height=320)
            st.download_button("Download CSV", gdf.to_csv(index=False).encode(),
                               file_name=f"grid_{px}_vs_{py}.csv", mime="text/csv")

# ------------------------------------------------------------------------ help
with tab_help:
    st.markdown(
        """
### What the model does

One repeating centriole unit is a set of **rigid bodies joined by
connections that can rotate**. Because the units must close into a ring,
changing one part forces the others to adapt. The model finds the
least-strained arrangement that still fits together.

### Joint grades

Each connection is graded against **its own** rotation limits — a given
angle means very different things at a 45 nm spoke and a 13 nm linker arm.

| Grade | Meaning |
|---|---|
| **OK** | comfortably accommodated |
| **HARD** | strained but achievable; expect distortion |
| **FORBIDDEN** | very costly; unlikely to form this way |

Contacts on microtubules are tightest (they grip a rigid, ordered
lattice); the triplet axis and base are most permissive.

### The other readouts

| Readout | Meaning |
|---|---|
| **Outer diameter** | Real centrioles are about 250 nm. Wild type here gives 255 nm. |
| **Lumen** | The central aperture. |
| **Bond load** | Which connection carries most strain, i.e. would break first. Bonds have different strengths, so weak ones yield first. |
| **Buckling** | Segments forced to bow because space got tight. Nothing ever stretches. |
| **MT clashes** | Microtubules overlapping in space. Anything above zero is physically impossible. |
| **Strand clearance** | How close the linker, base or spoke passes to a microtubule. Negative means passing through one. |
| **Unattached triplets** | Triplets the cartwheel could not reach — expected when the two symmetries differ. |

### Worth knowing

The triplet ring is the primary scaffold: when the two symmetries
disagree, the triplets keep their spacing and the **cartwheel** absorbs
the mismatch. Try cartwheel 8 with triplets 9 and watch the spoke grade
go HARD while the triplets stay OK.

Every dimension is measured from a cryo-ET-derived schematic. The joint
band values, however, are reasoned heuristics rather than measurements —
see `README.md`.
        """
    )
