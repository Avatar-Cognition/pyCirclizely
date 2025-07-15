from __future__ import annotations

import math
import textwrap
import warnings
from copy import deepcopy
from typing import Any

import numpy as np
import plotly.graph_objects as go
from plotly.basedatatypes import BaseTraceType

from pycirclizely import config, utils
from pycirclizely.patches import PolarSVGPatchBuilder
from pycirclizely.track import Track


class Sector:
    """Circos Sector Class"""

    def __init__(
        self,
        name: str,
        size: float | tuple[float, float],
        rad_lim: tuple[float, float],
        clockwise: bool = True,
    ):
        """
        Args:
            name: Sector name.
            size: Sector size (or range).
            rad_lim: Sector radian limit region.
            clockwise: Sector coordinate direction (clockwise or anti-clockwise).
        """
        self._name = name
        if isinstance(size, (tuple, list)):
            start, end = size[0], size[1]
        else:
            start, end = 0, size
        self._start = start
        self._end = end
        self._size = end - start
        self._rad_lim = rad_lim
        self._clockwise = clockwise
        self._tracks: list[Track] = []

        # Shapes and annotations for Layout
        self._shapes: list[go.layout.Shape] = []
        self._annotations: list[go.layout.Annotation] = []
        self._traces: list[BaseTraceType] = []

    ############################################################
    # Property
    ############################################################

    @property
    def name(self) -> str:
        """Sector name"""
        return self._name

    @property
    def size(self) -> float:
        """Sector size (x coordinate)"""
        return self._size

    @property
    def start(self) -> float:
        """Sector start position (x coordinate)"""
        return self._start

    @property
    def end(self) -> float:
        """Sector end position (x coordinate)"""
        return self._end

    @property
    def center(self) -> float:
        """Sector center position (x coordinate)"""
        return (self.start + self.end) / 2

    @property
    def rad_size(self) -> float:
        """Sector radian size"""
        return max(self.rad_lim) - min(self.rad_lim)

    @property
    def rad_lim(self) -> tuple[float, float]:
        """Sector radian limit"""
        return self._rad_lim

    @property
    def deg_size(self) -> float:
        """Sector degree size"""
        return max(self.deg_lim) - min(self.deg_lim)

    @property
    def deg_lim(self) -> tuple[float, float]:
        """Sector degree limit"""
        return (math.degrees(self.rad_lim[0]), math.degrees(self.rad_lim[1]))

    @property
    def clockwise(self) -> bool:
        """Sector coordinate direction"""
        return self._clockwise

    @property
    def tracks(self) -> list[Track]:
        """Tracks in sector"""
        return self._tracks

    @property
    def shapes(self) -> list[go.layout.Shape]:
        """Layout shapes"""
        return self._shapes

    @property
    def annotations(self) -> list[go.layout.Annotation]:
        """Layout annotations"""
        return self._annotations

    @property
    def traces(self) -> list[BaseTraceType]:
        """Data traces"""
        return self._traces

    ############################################################
    # Public Method
    ############################################################

    def add_track(
        self,
        r_lim: tuple[float, float],
        *,
        r_pad_ratio: float = 0,
        name: str | None = None,
    ) -> Track:
        """
        Args:
            r_lim: Radius limit region.
            r_pad_ratio: Track padding ratio for plot data.
            name: Track name. If None, a default name is generated.
        """
        name = f"Track{len(self.tracks) + 1:02d}" if name is None else name
        if name in [t.name for t in self.tracks]:
            raise ValueError(f"{name=} track is already exists.")
        if not 0 <= min(r_lim) <= max(r_lim) <= 100:
            warn_msg = f"{r_lim=} is unexpected plot range (0 <= r <= 100)."
            warnings.warn(warn_msg, stacklevel=1)
        track = Track(name, r_lim, r_pad_ratio, self)
        self._tracks.append(track)
        return track

    def get_track(self, name: str) -> Track:
        """Get track instance by name."""
        name2track = {t.name: t for t in self.tracks}
        if name not in name2track:
            raise ValueError(f"{name=} track not exists.")
        return name2track[name]

    def get_lowest_r(self) -> float:
        """Get lowest radius position of sector from tracks data."""
        if len(self.tracks) == 0:
            return config.MAX_R
        return min([min(t.r_lim) for t in self.tracks])

    def x_to_rad(self, x: float, ignore_range_error: bool = False) -> float:
        """Convert x coordinate to radian in sector start-end range.

        Args:
            x: X coordinate.
            ignore_range_error: Ignore x coordinate range error.
        """
        # Check target x is in valid sector range
        if not ignore_range_error:
            # Apply relative torelance value to sector range to avoid
            # unexpected invalid range error due to rounding errors (Issue #27, #67)
            min_range = self.start - config.REL_TOL
            max_range = self.end + config.REL_TOL
            if not min_range <= x <= max_range:
                err_msg = f"{x=} is invalid range of '{self.name}' sector.\n{self}"
                raise ValueError(err_msg)

        if not self.clockwise:
            x = (self.start + self.end) - x
        size_ratio = self.rad_size / self.size
        x_from_start = x - self.start
        rad_from_start = x_from_start * size_ratio
        rad = min(self.rad_lim) + rad_from_start
        return rad

    def axis(self, **kwargs: Any) -> None:
        """Plot axis shapes for the sector.

        Args:
            **kwargs: Shape properties.
                <https://plotly.com/python/reference/layout/shapes/>
        """
        kwargs = {} if kwargs is None else kwargs

        # Background shape placed behind other shapes (layer="below")
        fc_behind_kwargs = deepcopy(kwargs)
        fc_behind_kwargs = utils.deep_dict_update(
            fc_behind_kwargs, config.AXIS_FACE_PARAM
        )
        self.rect(self.start, self.end, config.R_LIM, **fc_behind_kwargs)

        # Edge shape placed in front of other shapes (layer="above")
        ec_front_kwargs = deepcopy(kwargs)
        ec_front_kwargs = utils.deep_dict_update(
            ec_front_kwargs, config.AXIS_EDGE_PARAM
        )
        self.rect(self.start, self.end, config.R_LIM, **ec_front_kwargs)

    def text(
        self,
        text: str,
        x: float | None = None,
        r: float = 107,
        *,
        adjust_rotation: bool = True,
        orientation: str = "horizontal",
        ignore_range_error: bool = False,
        **kwargs: Any,
    ) -> None:
        """Plot text within a sector. Uses genomic coordinates (x) mapped to radians.

        Args:
            text: Text content.
            x: Genomic position. If None, sector center is used.
            r: Radius position (default: 105, outer edge).
            adjust_rotation: If True, text rotation auto based on `x` and `orientation`.
            orientation: Text orientation (`horizontal` or `vertical`).
            ignore_range_error: If True, ignores x position outside sector bounds.
            **kwargs: Annotation properties (e.g. `font=dict(size=12, color='red')`).
                <https://plotly.com/python/reference/layout/annotations/>
        """
        x = self.center if x is None else x
        rad = self.x_to_rad(x, ignore_range_error)
        plotly_rad = -(rad - np.pi / 2)
        x_pos = r * np.cos(plotly_rad)
        y_pos = r * np.sin(plotly_rad)

        annotation = utils.plot.get_plotly_label_params(
            rad, adjust_rotation, orientation, **kwargs
        )

        annotation.update(
            {
                "x": x_pos,
                "y": y_pos,
                "text": text,
            }
        )

        annotation_layout = go.layout.Annotation(**annotation)
        self._annotations.append(annotation_layout)

    def line(
        self,
        *,
        r: float | tuple[float, float],
        start: float | None = None,
        end: float | None = None,
        arc: bool = True,
        **kwargs: Any,
    ) -> None:
        """Plot line with sector-relative coordinates.

        Args:
            r: Radius position(s). If float, creates constant-radius line.
            start: Genomic start position. Uses sector start if None.
            end: Genomic end position. Uses sector end if None.
            arc: If True, creates curved arc line (polar projection).
                If False, creates straight chord line.
            **kwargs: Shape properties
                (e.g. `line=dict(color="red", width=2, dash="dash")`).
                <https://plotly.com/python/reference/layout/shapes/>
        """
        # Set default genomic coordinates
        start = self.start if start is None else start
        end = self.end if end is None else end

        # Convert to polar coordinates
        rad_lim = (self.x_to_rad(start), self.x_to_rad(end))
        r_lim = (r, r) if isinstance(r, (float, int)) else r

        # Generate path based on arc preference
        path = (
            PolarSVGPatchBuilder.arc_line(rad_lim, r_lim)
            if arc
            else PolarSVGPatchBuilder.straight_line(rad_lim, r_lim)
        )

        # Create shape with defaults and kwargs
        shape = utils.plot.build_plotly_shape(
            path, config.plotly_shape_defaults, **kwargs
        )
        self._shapes.append(shape)

    def rect(
        self,
        start: int | float | None = None,
        end: int | float | None = None,
        r_lim: tuple[float, float] | None = None,
        **kwargs: Any,
    ) -> None:
        """Plot a rectangle spanning angular and radial ranges.

        Args:
            start: Start position (x coordinate). If None, `sector.start` is set.
            end: End position (x coordinate). If None, `sector.end` is set.
            r_lim: Radius limit region. If None, (0, 100) is set.
            **kwargs: Shape properties
                (e.g. `fillcolor="red", line: {color: "blue", width: 2, ... } ...`)
                <https://plotly.com/python/reference/layout/shapes/>
        """
        start = self.start if start is None else start
        end = self.end if end is None else end
        rad_rect_start = self.x_to_rad(start)
        rad_rect_end = self.x_to_rad(end)

        r_lim = config.R_LIM if r_lim is None else r_lim
        min_rad = min(rad_rect_start, rad_rect_end)
        max_rad = max(rad_rect_start, rad_rect_end)

        radr = (min_rad, min(r_lim))
        width = max_rad - min_rad
        height = max(r_lim) - min(r_lim)

        path = PolarSVGPatchBuilder.arc_rectangle(radr, width, height)
        shape = utils.plot.build_plotly_shape(
            path, config.plotly_shape_defaults, **kwargs
        )
        self._shapes.append(shape)

    ############################################################
    # Private Method
    ############################################################

    def __str__(self):
        min_deg_lim, max_deg_lim = min(self.deg_lim), max(self.deg_lim)
        track_names = [t.name for t in self.tracks]
        return textwrap.dedent(
            f"""
            # Sector = '{self.name}'
            # Size = {self.size} ({self.start} - {self.end})
            # Degree Size = {self.deg_size:.2f} ({min_deg_lim:.2f} - {max_deg_lim:.2f})
            # Track List = {track_names}
            """
        )[1:]
