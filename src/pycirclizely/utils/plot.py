from __future__ import annotations

import math
from copy import deepcopy
from enum import IntEnum
from typing import Literal

import numpy as np
from plotly.graph_objs import graph_objs as go  # type: ignore[attr-defined]

from pycirclizely import config
from pycirclizely.parser.table import StackedBarTable

from .color import ColorCycler
from .helper import deep_dict_update, precise_position


def get_default_color(kwargs: dict, target: str = "line") -> str:
    """Returns a consistent color based on kwargs or assigns a new one from ColorCycler.

    Parameters
    ----------
    kwargs : dict
        Dictionary of Plotly styling keyword arguments.
    target : str
        The key to check for color (e.g., 'line', 'marker').

    Returns
    -------
    str
        A color string (e.g., "#1f77b4").

    """
    color = kwargs.get(target, {})

    if isinstance(color, dict):
        color = color.get("color")

    if color is None:
        color = ColorCycler.get_color()

    return color


def degrees(rad: float) -> float:
    """Convert radian to positive degree (`0 - 360`)

    Parameters
    ----------
    rad : float
        Target radian

    Returns
    -------
    deg : float
        Positive degree (`0 - 360`)

    """
    # Radian to degree
    deg = math.degrees(rad)
    # Normalize degree in 0 - 360 range
    deg = deg % 360
    # Negative to positive
    if deg < 0:
        deg += 360
    return deg


def is_lower_loc(rad: float) -> bool:
    """Check target radian is lower location or not

    Parameters
    ----------
    rad : float
        Target radian

    Returns
    -------
    result : bool
        Lower location or not

    """
    deg = math.degrees(rad)
    return -270 <= deg < -90 or 90 <= deg < 270


def is_right_loc(rad: float) -> bool:
    """Check target radian is right location or not

    Parameters
    ----------
    rad : float
        Target radian

    Returns
    -------
    result : bool
        Right location or not

    """
    deg = math.degrees(rad)
    return -360 <= deg < -180 or 0 <= deg < 180


def is_ann_rad_shift_target_loc(rad: float) -> bool:
    """Check radian is annotation radian shift target or not

    Parameters
    ----------
    rad : float
        Annotation radian position

    Returns
    -------
    result : bool
        Target or not

    """
    deg = degrees(rad)
    return 30 <= deg <= 150 or 210 <= deg <= 330


def get_loc(
    rad: float,
) -> Literal["upper-right", "lower-right", "lower-left", "upper-left"]:
    """Get location of 4 sections

    Returns
    -------
    loc : str
        Location (`upper-right`|`lower-right`|`lower-left`|`upper-left`)

    """
    deg = degrees(rad)
    if 0 <= deg < 90:
        return "upper-right"
    elif 90 <= deg < 180:
        return "lower-right"
    elif 180 <= deg < 270:
        return "lower-left"
    else:
        return "upper-left"


def get_ann_relpos(rad: float) -> tuple[float, float]:
    """Get relative position for annotate by radian text position

    Parameters
    ----------
    rad : float
        Radian text position

    Returns
    -------
    relpos : tuple[float, float]
        Relative position

    """
    deg = degrees(rad)
    if 0 <= deg <= 180:
        return 0.0, Normalize(0, 180)(deg)
    else:
        return 1.0, 1.0 - Normalize(180, 360)(deg)


def get_plotly_label_params(
    rad: float,
    adjust_rotation: bool,
    orientation: str,
    **kwargs,
) -> dict:
    """Build Plotly label parameters based on radian and orientation."""
    # Start with global defaults
    annotation = deepcopy(config.plotly_annotation_defaults)

    # Override with user-provided kwargs
    annotation = deep_dict_update(annotation, kwargs)

    if adjust_rotation:
        rotation = np.degrees(rad)

        if orientation == "horizontal":
            rotation = rotation % 360
            # Flip if upside-down
            if 90 < rotation <= 270:
                rotation += 180
        elif orientation == "vertical":
            # Point text radially (90° offset from horizontal)
            rotation = (rotation + 90) % 360
            # Flip for vertical text
            if 90 < rotation <= 270:
                rotation += 180

        annotation.update({"textangle": rotation})

    return annotation


def build_plotly_shape(path: str, defaults: dict = {}, **kwargs) -> dict:
    """Build a Plotly shape dictionary with defaults and custom parameters."""
    shape_defaults = deepcopy(defaults)
    shape_defaults = deep_dict_update(shape_defaults, kwargs)
    return {"type": "path", "path": path, **shape_defaults}


def build_scatter_trace(x: list, y: list, mode: str, **kwargs) -> go.Scatter:
    """Build a Plotly Scatter trace with defaults and custom parameters."""
    scatter_config = deepcopy(config.plotly_scatter_defaults)
    scatter_config["mode"] = mode
    scatter_config = deep_dict_update(scatter_config, kwargs)

    return go.Scatter(x=x, y=y, **scatter_config)


def default_hovertext(
    x: list[int] | list[float] | np.ndarray,
    y: list[int] | list[float] | np.ndarray,
    x2: list[int] | list[float] | np.ndarray | None = None,
    sector_name: str | None = None,
    precision_position: int = 2,
) -> list[str]:
    """Generate default hovertext for a Plotly scatter trace."""
    value_format = f".{precision_position}f" if precision_position > 0 else ".0f"

    # Convert numpy arrays to lists if needed
    if isinstance(x, np.ndarray):
        x = x.tolist()
    if isinstance(y, np.ndarray):
        y = y.tolist()
    if x2 is not None and isinstance(x2, np.ndarray):
        x2 = x2.tolist()

    hovertext = []
    for i, (xi, yi) in enumerate(zip(x, y)):
        parts = []
        xi = precise_position(xi, precision_position)
        if sector_name:
            parts.append(f"Sector: {sector_name}")
        if x2 is not None:
            xi2 = precise_position(x2[i], precision_position)
            parts.append(f"Position: {xi}–{xi2}")
        else:
            parts.append(f"Position: {xi}")
        parts.append(f"Value: {format(yi, value_format)}")
        hovertext.append("<br>".join(parts))
    return hovertext


def default_stackedbar_hovertext(
    sb_table: StackedBarTable,
) -> list[str]:
    """Generate default hover text for stacked bars.

    Parameters
    ----------
    sb_table : StackedBarTable
        The stacked bar table object


    Returns
    -------
    list[str]
        Formatted hover text for each segment
    """
    hover_texts = []
    totals = list(sb_table.row_name2sum.values())

    # Get values by column (segment)
    for col_idx, col_name in enumerate(sb_table.col_names):
        col_values = sb_table.stacked_bar_heights[col_idx]
        for row_idx, row_name in enumerate(sb_table.row_names):
            value = col_values[row_idx]
            parts = []

            parts.append(f"<b>{col_name}</b>: {value:.1f}")
            parts.append(f"<b>{row_name}</b> (Total: {totals[row_idx]:.1f})")

            hover_texts.append("<br>".join(parts))

    return hover_texts


class Normalize:
    def __init__(self, vmin, vmax, clip=False):
        if vmin == vmax:
            raise ValueError("vmin and vmax must be different")
        self.vmin = vmin
        self.vmax = vmax
        self.clip = clip

    def __call__(self, value):
        """Normalize a value to the range [0, 1]."""
        normed = (value - self.vmin) / (self.vmax - self.vmin)
        if self.clip:
            return max(0.0, min(1.0, normed))
        return normed


class LinkDirection(IntEnum):
    NONE = 0
    FORWARD = 1
    REVERSE = -1
    BIDIRECTIONAL = 2

    def arrow(self) -> str:
        """Return the arrow representation of the link direction."""
        return {
            LinkDirection.NONE: "-",
            LinkDirection.FORWARD: "->",
            LinkDirection.REVERSE: "<-",
            LinkDirection.BIDIRECTIONAL: "<->",
        }[self]
