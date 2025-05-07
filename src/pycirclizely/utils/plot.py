from __future__ import annotations

import math
from copy import deepcopy
from typing import Literal

import numpy as np
from pycirclizely import config
from plotly.graph_objs import graph_objs as go
from .helper import ColorCycler, deep_dict_update

def get_default_color(kwargs: dict, target: str = "line") -> str:
    """
    Returns a consistent color based on kwargs or assigns a new one from the ColorCycler.

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
    target_dict = kwargs.get(target, {})
    color = target_dict.get("color")

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


def plot_bbox(bbox: Bbox, ax: PolarAxes, **kwargs) -> None:
    """Plot bbox to check bbox area for development

    Parameters
    ----------
    bbox : Bbox
        Bounding box
    ax : PolarAxes
        Polar axes
    **kwargs : dict, optional
        Axes.plot properties (e.g. `color="red", lw=0.5, ls="--", ...`)
        <https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.plot.html>
    """
    trans_bbox = bbox.transformed(ax.transAxes.inverted())
    kwargs.setdefault("clip_on", False)
    x0, y0, x1, y1 = trans_bbox.x0, trans_bbox.y0, trans_bbox.x1, trans_bbox.y1
    x, y = [x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0]
    ax.plot(x, y, transform=ax.transAxes, **kwargs)


def get_plotly_label_params(
    rad: float,
    adjust_rotation: bool,
    orientation: str,
    outer: bool = True,
    only_rotation: bool = False,
    **kwargs,
) -> dict:
    # Start with global defaults
    annotation = deepcopy(config.plotly_annotation_defaults)

    if not only_rotation:
        # Apply orientation-specific defaults
        orientation_defaults = config.plotly_text_orientation_defaults.get(
            orientation, {}
        )
        annotation = deep_dict_update(annotation, orientation_defaults)

        # Handle outer/inner alignment
        if not outer:
            if orientation == "horizontal":
                annotation["yanchor"] = (
                    "bottom" if annotation["yanchor"] == "top" else "top"
                )
            elif orientation == "vertical":
                annotation["xanchor"] = (
                    "right" if annotation["xanchor"] == "left" else "left"
                )
        else:
            annotation["xanchor"] = "center"
            annotation["yanchor"] = "middle"

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
            # ADDED FLIPPING FOR VERTICAL TEXT
            if 90 < rotation <= 270:
                rotation += 180

        annotation.update({"textangle": rotation})

    return annotation


def build_plotly_shape(path: str, defaults: dict = {}, **kwargs) -> dict:
    shape_defaults = deepcopy(defaults)
    shape_defaults = deep_dict_update(shape_defaults, kwargs)
    return {"type": "path", "path": path, **shape_defaults}


def build_scatter_trace(x: list, y: list, mode: str, **kwargs) -> go.Scatter:
    scatter_config = deepcopy(config.plotly_scatter_defaults)
    scatter_config["mode"] = mode
    scatter_config = deep_dict_update(scatter_config, kwargs)
    
    return go.Scatter(x=x, y=y, **scatter_config)


def default_hovertext(
    x: list[float],
    y: list[float],
    sector_name: str | None = None,
    value_format: str = ".2f",
) -> list[str]:
    """
    Generate default hovertext for Plotly traces.

    Parameters
    ----------
    x : list[float]
        List of x-axis (e.g., genomic position or angle) values.
    y : list[float]
        List of y-axis (e.g., data value) values.
    sector_name : str, optional
        Name of the sector to include in the hover text.
    value_format : str, optional
        Format string for the y values (default is ".2f").

    Returns
    -------
    list[str]
        List of formatted hover text strings.
    """
    hovertext = []
    for xi, yi in zip(x, y):
        parts = []
        if sector_name:
            parts.append(f"Sector: {sector_name}")
        parts.append(f"Position: {xi}")
        parts.append(f"Value: {format(yi, value_format)}")
        hovertext.append("<br>".join(parts))
    return hovertext


class Normalize:
    def __init__(self, vmin, vmax, clip=False):
        if vmin == vmax:
            raise ValueError("vmin and vmax must be different")
        self.vmin = vmin
        self.vmax = vmax
        self.clip = clip

    def __call__(self, value):
        normed = (value - self.vmin) / (self.vmax - self.vmin)
        if self.clip:
            return max(0.0, min(1.0, normed))
        return normed