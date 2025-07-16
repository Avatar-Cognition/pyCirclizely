from __future__ import annotations

import math
from copy import deepcopy
from enum import IntEnum
from typing import Any, Literal

import numpy as np
from plotly.basedatatypes import BaseTraceType
from plotly.graph_objs import graph_objs as go  # type: ignore[attr-defined]

from pycirclizely import config

from .color import ColorCycler
from .helper import deep_dict_update

_DEFAULT_CYCLER = ColorCycler("Plotly")


def get_default_color(
    kwargs: dict, target: str = "line", cycler: ColorCycler = _DEFAULT_CYCLER
) -> str:
    """Returns a consistent color based on kwargs or assigns a new one from ColorCycler.

    Args:
        kwargs: Dictionary of Plotly styling keyword arguments.
        target: The key to check for color (e.g., 'line', 'marker').
        color_cycler : ColorCycler, optional
            ColorCycler instance to use. If None, use default one.
    """
    color = (
        kwargs.get(target, {}).get("color")
        if isinstance(kwargs.get(target), dict)
        else kwargs.get(target)
    )
    return color if color is not None else cycler.get_color()


def degrees(rad: float) -> float:
    """Convert radian to positive degree (`0 - 360`)

    Args:
        rad: Target radian
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

    Args:
        rad: Target radian
    """
    deg = math.degrees(rad)
    return -270 <= deg < -90 or 90 <= deg < 270


def is_right_loc(rad: float) -> bool:
    """Check target radian is right location or not

    Args:
        rad: Target radian
    """
    deg = math.degrees(rad)
    return -360 <= deg < -180 or 0 <= deg < 180


def is_ann_rad_shift_target_loc(rad: float) -> bool:
    """Check radian is annotation radian shift target or not

    Args:
        rad: Annotation radian position
    """
    deg = degrees(rad)
    return 30 <= deg <= 150 or 210 <= deg <= 330


def get_loc(
    rad: float,
) -> Literal["upper-right", "lower-right", "lower-left", "upper-left"]:
    """Get location of 4 sections"""
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

    Args:
        rad: Radian text position
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
    **kwargs: Any,
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


def build_plotly_shape(path: str, defaults: dict = {}, **kwargs: Any) -> dict:
    """Build a Plotly shape dictionary with defaults and custom parameters."""
    shape_defaults = deepcopy(defaults)
    shape_defaults = deep_dict_update(shape_defaults, kwargs)
    return {"type": "path", "path": path, **shape_defaults}


def build_scatter_trace(
    x: list | tuple, y: list | tuple, mode: str, **kwargs: Any
) -> BaseTraceType:
    """Build a Plotly Scatter trace with defaults and custom parameters."""
    scatter_config = deepcopy(config.plotly_scatter_defaults)
    scatter_config["mode"] = mode
    scatter_config = deep_dict_update(scatter_config, kwargs)

    return go.Scatter(x=x, y=y, **scatter_config)


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
