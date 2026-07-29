"""Helpers for perturbation studies on top of the `centriole` model.

The core `centriole()` function models a single, fixed parameter set. For
asking "how does perturbing X change the geometry?" (e.g. a SAS-6 mutation
that drops symmetry from 9-fold to 8-fold, or a longer SAS-6 coiled-coil),
what you actually want is:

- `sweep_param`: vary ONE parameter over a range, holding everything else
  fixed, and get back a tidy table of the resulting (success, overlap,
  b_deg, LKm) for each value -- e.g. "as SAS-6 length r increases from 20
  to 60nm, at what point does the model stop finding a non-overlapping
  solution?"
- `plot_grid`: render several named, discrete configurations (e.g.
  "wild-type" vs "8-fold mutant") side by side in one figure, so you can
  visually compare geometries directly instead of opening separate PDFs.

Both are thin wrappers around a centriole model function -- they don't
reimplement any of the geometry, they just call it repeatedly and
organize the results. They default to `centriole()` (the faithful MATLAB
port) but accept `model=centriole_v2` (the cartwheel/pinhead/triplet-base
model in `centriole_v2.py`) or any other function with a compatible
signature: `(..., show_result, ax, return_details) -> CentrioleResult |
(CentrioleResult, details_dict_with_"LKm")`.
"""

from __future__ import annotations

from typing import Any, Callable, Iterable, Mapping, Optional

import matplotlib.pyplot as plt
import pandas as pd

from centriole import centriole


def sweep_param(
    param: str,
    values: Iterable[Any],
    model: Callable = centriole,
    **fixed: Any,
) -> pd.DataFrame:
    """Vary a single model parameter and tabulate the results.

    Parameters
    ----------
    param : str
        Name of the model's keyword argument to vary (e.g. "r", "SYM",
        "GAMMA", or for `centriole_v2`, "SPOKE_LENGTH", "BASE_LENGTH", ...).
    values : iterable
        Values to try for `param`.
    model : callable, default `centriole`
        The model function to sweep -- `centriole` or `centriole_v2`
        (import it and pass e.g. `model=centriole_v2`).
    **fixed
        Any other model keyword arguments, held constant across the sweep
        (e.g. `MTn=3, LK=25, r=77`).

    Returns
    -------
    pandas.DataFrame
        One row per value, columns: the swept param, `success`,
        `overlap`, `b_deg`, `LKm` (achieved linker length).

    Example
    -------
    >>> sweep_param("r", range(20, 90, 5), SYM=9, MTn=3, LK=25)
    >>> from centriole_v2 import centriole_v2
    >>> sweep_param("BASE_LENGTH", range(10, 40, 5), model=centriole_v2, SYM=9, LK=25)
    """

    fixed = dict(fixed)
    fixed.pop("show_result", None)  # sweeps never plot per-call
    fixed.pop("return_details", None)  # sweep always needs LKm -> forced True below

    rows = []
    for value in values:
        kwargs = dict(fixed)
        kwargs[param] = value
        result, details = model(show_result=False, return_details=True, **kwargs)
        rows.append(
            {
                param: value,
                "success": result.success,
                "overlap": result.overlap,
                "b_deg": result.b_deg,
                "LKm": details["LKm"],
            }
        )
    return pd.DataFrame(rows)


def plot_grid(
    configs: Mapping[str, Mapping[str, Any]],
    ncols: int = 3,
    figsize_per_panel: tuple = (4.0, 4.0),
    model: Callable = centriole,
    **fixed: Any,
) -> "tuple[plt.Figure, pd.DataFrame]":
    """Render several named configurations side by side for visual comparison.

    Parameters
    ----------
    configs : mapping of label -> model kwargs
        Each entry is one panel, e.g.
        ``{"WT (9-fold)": {"SYM": 9}, "8-fold mutant": {"SYM": 8}}``.
        Keys given here override any same-named key in `**fixed`.
    ncols : int
        Number of panels per row.
    figsize_per_panel : (width, height)
        Matplotlib inches per panel; the overall figure is sized to fit
        all panels.
    model : callable, default `centriole`
        The model function to render -- `centriole` or `centriole_v2`.
    **fixed
        Model keyword arguments shared by every panel unless a config
        overrides them (e.g. `MTn=3, LK=25, r=77`).

    Returns
    -------
    (matplotlib.figure.Figure, pandas.DataFrame)
        The comparison figure, and a table of (label, success, overlap,
        b_deg, LKm) for all panels -- so you get the numeric comparison
        alongside the visual one.
    """

    fixed = dict(fixed)
    fixed.pop("show_result", None)
    fixed.pop("ax", None)
    fixed.pop("return_details", None)

    labels = list(configs.keys())
    nrows = -(-len(labels) // ncols)  # ceil div
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows), squeeze=False
    )
    axes_flat = axes.flatten()

    rows = []
    for ax, label in zip(axes_flat, labels):
        kwargs = {**fixed, **configs[label]}
        result, details = model(show_result=True, ax=ax, return_details=True, **kwargs)
        ax.set_title(f"{label}\n{ax.get_title()}", fontsize=8)
        rows.append(
            {
                "label": label,
                "success": result.success,
                "overlap": result.overlap,
                "b_deg": result.b_deg,
                "LKm": details["LKm"],
            }
        )

    # Hide any unused trailing panels (e.g. 5 configs in a 2x3 grid).
    for ax in axes_flat[len(labels):]:
        ax.set_visible(False)

    fig.tight_layout()
    return fig, pd.DataFrame(rows)
