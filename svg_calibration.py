"""Extract quantitative centriole geometry from the schematic SVG.

Reads `2x_units_wScale.svg` (two consecutive units + a 500 Angstrom scale
bar) and derives the *unit-internal* geometry: segment lengths, the
pinhead's three attachment points, triplet tubule layout, protofilament
angular positions, and the A-C linker's attachment points.

Everything is returned in the unit's OWN frame (relative to its spoke
direction and its triplet axis), deliberately independent of where the
centriole axis actually is. That matters because the schematic is
hand-drawn and its two units sit 41.4 degrees apart while its two spokes
sit exactly 40 degrees apart -- so a globally-fitted centre is ambiguous,
while the unit-internal shape is consistent between the two units to
within ~1%.

Run directly to print the calibration report:

    ./.venv/bin/python svg_calibration.py
"""

from __future__ import annotations

import base64
import io
import re
from pathlib import Path
from dataclasses import dataclass, field

import numpy as np

SVG_PATH = str(Path(__file__).resolve().parent / "data" / "2x_units_wScale.svg")
SCALE_BAR_ANGSTROM = 500.0


def _flip(p):
    """SVG y-axis points down; flip to a conventional y-up frame."""
    p = np.asarray(p, dtype=float)
    return np.array([p[..., 0], -p[..., 1]]).T if p.ndim > 1 else np.array([p[0], -p[1]])


def _fit_circle(P):
    """Algebraic least-squares circle fit. Returns (centre, radius)."""
    x, y = P[:, 0], P[:, 1]
    A = np.c_[2 * x, 2 * y, np.ones(len(x))]
    b = x**2 + y**2
    cx, cy, c = np.linalg.lstsq(A, b, rcond=None)[0]
    return np.array([cx, cy]), float(np.sqrt(c + cx**2 + cy**2))


def _rot(deg):
    t = np.radians(deg)
    return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])


def angstrom_per_unit(svg: str) -> float:
    """Calibrate from the embedded scale-bar PNG (its solid black bar)."""
    from PIL import Image

    m = re.search(
        r'<image[^>]*transform="translate\([-\d.]+ [-\d.]+\) scale\(([\d.]+)\)"'
        r'[^>]*xlink:href="data:image/png;base64,([^"]+)"',
        svg,
    )
    scale = float(m.group(1))
    img = Image.open(io.BytesIO(base64.b64decode(m.group(2)))).convert("RGBA")
    a = np.array(img)
    dark = (a[..., 3] > 128) & (a[..., :3].mean(axis=2) < 100)
    cols = np.where(dark.any(axis=0))[0]
    bar_px = cols.max() - cols.min() + 1
    return SCALE_BAR_ANGSTROM / (bar_px * scale)


def parse_rings(svg: str):
    """Protofilament labels -> per-unit A/B/C rings (centre, radius, pf angles)."""
    pat = re.compile(
        r'<text class="([^"]+)" transform="translate\(([-\d.]+) ([-\d.]+)\)'
        r'(?: rotate\(([-\d.]+)\))?"><tspan[^>]*>([^<]*)</tspan></text>'
    )
    rows = [
        (m.group(1), float(m.group(2)), float(m.group(3)), float(m.group(4) or 0), m.group(5))
        for m in pat.finditer(svg)
    ]
    units = {}
    for cls in ("cls-7", "cls-8"):
        sub = [r for r in rows if r[0] == cls]
        by_rot = {}
        for r in sub:
            by_rot.setdefault(r[3], []).append(r)
        rings = {}
        for rot, items in by_rot.items():
            P = np.array([[it[1], -it[2]] for it in items])
            centre, radius = _fit_circle(P)
            rings[rot] = dict(
                centre=centre, radius=radius, pts=P, labels=[int(it[4]) for it in items], n=len(items)
            )
        a_rot = [k for k, v in rings.items() if v["n"] == 13][0]
        nines = sorted(
            (float(np.linalg.norm(rings[k]["centre"] - rings[a_rot]["centre"])), k)
            for k in rings
            if rings[k]["n"] == 9
        )
        units[cls] = {"A": rings[a_rot], "B": rings[nines[0][1]], "C": rings[nines[1][1]]}
    return units


def parse_bases(svg: str):
    """Green triplet-base rects -> (inner_end, outer_end, length, thickness)."""
    out = []
    for m in re.finditer(
        r'<rect class="cls-4" x="([-\d.]+)" y="([-\d.]+)" width="([-\d.]+)" height="([-\d.]+)"'
        r' transform="translate\(([-\d.]+) ([-\d.]+)\) rotate\(([-\d.]+)\)"',
        svg,
    ):
        x, y, w, h, tx, ty, a = (float(g) for g in m.groups())
        R = _rot(a)
        p0 = R @ np.array([x, y + h / 2]) + np.array([tx, ty])
        p1 = R @ np.array([x + w, y + h / 2]) + np.array([tx, ty])
        out.append(dict(p0=_flip(p0), p1=_flip(p1), length=w, thickness=h))
    return out


_TOKEN = re.compile(r"([MmLlHhVvCcSsQqTtAaZz])|(-?\d*\.?\d+(?:[eE][-+]?\d+)?)")


def path_vertices(d: str) -> np.ndarray:
    """Minimal SVG path parser -> polyline vertices (curve control points dropped).

    Enough for this schematic, whose shapes are polygons with occasional
    rounded corners; only each command's ENDPOINT is retained, which is
    what the centreline measurements need.
    """
    toks = [(m.group(1), m.group(2)) for m in _TOKEN.finditer(d)]
    pts, nums, cmd = [], [], None
    cur = np.zeros(2)
    start = np.zeros(2)
    # per-command coordinate counts
    need = {"M": 2, "L": 2, "H": 1, "V": 1, "C": 6, "S": 4, "Q": 4, "T": 2, "A": 7, "Z": 0}

    def flush():
        nonlocal cur, start, nums
        if cmd is None:
            return
        u = cmd.upper()
        rel = cmd.islower()
        n = need[u]
        if n == 0:
            cur = start.copy()
            pts.append(cur.copy())
            nums = []
            return
        while len(nums) >= n:
            chunk, nums = nums[:n], nums[n:]
            if u == "M":
                p = np.array(chunk)
                cur = cur + p if rel else p
                start = cur.copy()
            elif u == "L":
                p = np.array(chunk)
                cur = cur + p if rel else p
            elif u == "H":
                cur = np.array([cur[0] + chunk[0] if rel else chunk[0], cur[1]])
            elif u == "V":
                cur = np.array([cur[0], cur[1] + chunk[0] if rel else chunk[0]])
            else:  # curves/arcs: keep only the endpoint
                p = np.array(chunk[-2:])
                cur = cur + p if rel else p
            pts.append(cur.copy())

    for c, num in toks:
        if c:
            flush()
            cmd = c
            nums = []
            if c.upper() == "Z":
                flush()
                cmd = None
        else:
            nums.append(float(num))
    flush()
    return np.array(pts)


def _long_axis_ends(P: np.ndarray):
    """Centreline endpoints of an elongated outline (PCA long axis)."""
    c = P.mean(axis=0)
    X = P - c
    u = np.linalg.svd(X, full_matrices=False)[2][0]
    t = X @ u
    return c + u * t.min(), c + u * t.max()


def parse_spokes(svg: str):
    """Light-blue SAS-6 paths -> spoke centrelines (both ends) + head length.

    Endpoint order is arbitrary here; the caller orients each spoke by
    which end lies nearer its A-tubule. The SAS-6 head is the flared
    arrow at the inner end: it is detected from the half-width profile
    along the long axis (the head flares to ~2.4x the coiled-coil rod's
    half-width), giving the head/rod split.
    """
    spokes = []
    for m in re.finditer(r'<path class="cls-1" d="([^"]+)"', svg):
        P = path_vertices(m.group(1))
        e0, e1 = _long_axis_ends(P)
        centre = P.mean(axis=0)
        X = P - centre
        U = np.linalg.svd(X, full_matrices=False)[2]
        t, w = X @ U[0], np.abs(X @ U[1])
        t = t - t.min()
        rod_half = np.median(w[t > 0.5 * t.max()])  # half-width of the plain rod
        flared = t[w > 1.5 * rod_half]
        head_len = float(flared.max()) if len(flared) else 0.0
        # the flare sits at whichever end the profile puts it
        spokes.append(dict(p0=_flip(e0), p1=_flip(e1), head_length=head_len))
    return spokes


@dataclass
class UnitCalibration:
    """Unit-internal geometry, in Angstrom, measured from the schematic."""

    angstrom_per_unit: float
    spoke_length: float
    spoke_head_length: float
    spoke_rod_length: float
    tubule_radius: float
    pf_pitch_deg: float
    n_pf: dict = field(default_factory=dict)
    ab_spacing: float = 0.0
    bc_spacing: float = 0.0
    base_length: float = 0.0
    base_thickness: float = 0.0
    triplet_axis_vs_spoke_deg: float = 0.0
    pinhead_span: float = 0.0
    pinhead_vs_spoke_deg: float = 0.0
    pf_angle_rel_axis: dict = field(default_factory=dict)
    unit_spacing_deg: float = 0.0
    spoke_spacing_deg: float = 0.0


def calibrate(svg_path: str = SVG_PATH) -> UnitCalibration:
    svg = open(svg_path).read()
    apu = angstrom_per_unit(svg)
    units = parse_rings(svg)
    bases = parse_bases(svg)
    spokes = parse_spokes(svg)

    # Pair each spoke with the unit whose A-tubule it feeds, and orient it:
    # the OUTER end is whichever endpoint lies nearer that A-tubule centre.
    pairing = {}
    for s in spokes:
        best_key, best_d, best_outer = None, np.inf, None
        for k, u in units.items():
            for outer, inner in ((s["p0"], s["p1"]), (s["p1"], s["p0"])):
                d = float(np.linalg.norm(u["A"]["centre"] - outer))
                if d < best_d:
                    best_key, best_d, best_outer = k, d, (outer, inner)
        pairing[best_key] = dict(
            outer=best_outer[0], inner=best_outer[1], head_length=s["head_length"]
        )

    def spoke_dir(s):
        v = s["outer"] - s["inner"]
        return np.degrees(np.arctan2(v[1], v[0]))

    per_unit = []
    for k, u in units.items():
        axis_v = u["C"]["centre"] - u["A"]["centre"]
        axis_deg = np.degrees(np.arctan2(axis_v[1], axis_v[0]))
        s = pairing[k]
        sd = spoke_dir(s)
        # protofilament angles relative to the triplet axis
        pf = {}
        for name in "ABC":
            ring = u[name]
            pf[name] = {}
            for p, lab in zip(ring["pts"], ring["labels"]):
                v = p - ring["centre"]
                pf[name][lab] = ((np.degrees(np.arctan2(v[1], v[0])) - axis_deg + 180) % 360) - 180
        # pinhead: spoke tip -> A pf3/4 contact
        contact_deg = (pf["A"][3] + pf["A"][4]) / 2 + axis_deg
        p_contact = u["A"]["centre"] + u["A"]["radius"] * np.array(
            [np.cos(np.radians(contact_deg)), np.sin(np.radians(contact_deg))]
        )
        pin_v = p_contact - s["outer"]
        per_unit.append(
            dict(
                axis_deg=axis_deg,
                spoke_deg=sd,
                spoke_len=float(np.linalg.norm(s["outer"] - s["inner"])),
                head_len=s["head_length"],
                ab=float(np.linalg.norm(u["B"]["centre"] - u["A"]["centre"])),
                bc=float(np.linalg.norm(u["C"]["centre"] - u["B"]["centre"])),
                r_mt=float(np.mean([u[n]["radius"] for n in "ABC"])),
                pin_span=float(np.linalg.norm(pin_v)),
                pin_deg=np.degrees(np.arctan2(pin_v[1], pin_v[0])),
                pf=pf,
            )
        )

    mean = lambda key: float(np.mean([d[key] for d in per_unit]))
    axis_vs_spoke = float(
        np.mean([((d["axis_deg"] - d["spoke_deg"] + 180) % 360) - 180 for d in per_unit])
    )
    pin_vs_spoke = float(
        np.mean([((d["pin_deg"] - d["spoke_deg"] + 180) % 360) - 180 for d in per_unit])
    )
    pf_avg = {}
    for name in "ABC":
        labs = sorted(per_unit[0]["pf"][name])
        pf_avg[name] = {l: float(np.mean([d["pf"][name][l] for d in per_unit])) for l in labs}

    spoke_spacing = abs(((per_unit[0]["spoke_deg"] - per_unit[1]["spoke_deg"] + 180) % 360) - 180)
    unit_spacing = abs(((per_unit[0]["axis_deg"] - per_unit[1]["axis_deg"] + 180) % 360) - 180)

    # wrap the difference into (-180, 180]; pf5 -> pf6 crosses the +-180 seam
    pitch = float(
        np.mean([abs(((pf_avg["A"][i + 1] - pf_avg["A"][i] + 180) % 360) - 180) for i in range(1, 13)])
    )

    return UnitCalibration(
        angstrom_per_unit=apu,
        spoke_length=mean("spoke_len") * apu,
        spoke_head_length=mean("head_len") * apu,
        spoke_rod_length=(mean("spoke_len") - mean("head_len")) * apu,
        tubule_radius=mean("r_mt") * apu,
        pf_pitch_deg=pitch,
        n_pf={"A": 13, "B": 9, "C": 9},
        ab_spacing=mean("ab") * apu,
        bc_spacing=mean("bc") * apu,
        base_length=float(np.mean([b["length"] for b in bases])) * apu,
        base_thickness=float(np.mean([b["thickness"] for b in bases])) * apu,
        triplet_axis_vs_spoke_deg=axis_vs_spoke,
        pinhead_span=mean("pin_span") * apu,
        pinhead_vs_spoke_deg=pin_vs_spoke,
        pf_angle_rel_axis=pf_avg,
        unit_spacing_deg=unit_spacing,
        spoke_spacing_deg=spoke_spacing,
    )


if __name__ == "__main__":
    c = calibrate()
    print(f"scale                     : {c.angstrom_per_unit:.4f} A per SVG unit")
    print(f"SAS-6 dimer total         : {c.spoke_length:7.1f} A"
          f"   (head {c.spoke_head_length:.1f} + coiled-coil rod {c.spoke_rod_length:.1f})")
    print(f"pinhead span (tip->A pf3/4): {c.pinhead_span:7.1f} A  at {c.pinhead_vs_spoke_deg:+.1f} deg from spoke")
    print(f"triplet base length       : {c.base_length:7.1f} A  (thickness {c.base_thickness:.1f} A)")
    print(f"tubule radius (pf ring)   : {c.tubule_radius:7.1f} A")
    print(f"A->B spacing              : {c.ab_spacing:7.1f} A")
    print(f"B->C spacing              : {c.bc_spacing:7.1f} A")
    print(f"triplet axis vs spoke     : {c.triplet_axis_vs_spoke_deg:+.2f} deg")
    print(f"protofilament pitch       : {c.pf_pitch_deg:.2f} deg  (360/13 = {360/13:.2f})")
    print(f"spoke angular spacing     : {c.spoke_spacing_deg:.2f} deg")
    print(f"triplet angular spacing   : {c.unit_spacing_deg:.2f} deg")
    print()
    print("protofilament angles relative to triplet axis (deg):")
    for name in "ABC":
        s = "  ".join(f"{l}:{v:+7.1f}" for l, v in sorted(c.pf_angle_rel_axis[name].items()))
        print(f"  {name}: {s}")
    print()
    print("connection points:")
    a = c.pf_angle_rel_axis
    print(f"  pinhead -> A pf3&4 : {(a['A'][3]+a['A'][4])/2:+.1f} deg")
    print(f"  linker  -> A pf8   : {a['A'][8]:+.1f} deg")
    print(f"  linker  -> C pf8&9 : {(a['C'][8]+a['C'][9])/2:+.1f} deg")
