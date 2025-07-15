from __future__ import annotations

import math
import textwrap
from copy import deepcopy
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from Bio.Phylo.BaseTree import Tree
from Bio.SeqFeature import SeqFeature
from plotly.basedatatypes import BaseTraceType
from plotly.colors import (  # type: ignore[attr-defined]
    get_colorscale,
    sample_colorscale,
)

from pycirclizely import config, utils
from pycirclizely.parser import StackedBarTable
from pycirclizely.patches import PolarSVGPatchBuilder
from pycirclizely.tree import TreeViz
from pycirclizely.types import (
    HoverText,
    LabelFormatter,
    Numeric,
    NumericComponent,
    NumericSequence,
    TextFormatter,
)
from pycirclizely.utils.hover_formatter import HoverFormatter

if TYPE_CHECKING:
    # Avoid Sector <-> Track circular import error at runtime
    from pycirclizely.sector import Sector


class Track:
    """Circos Track Class"""

    def __init__(
        self,
        name: str,
        r_lim: tuple[float, float],
        r_pad_ratio: float,
        parent_sector: Sector,
    ):
        """
        Args:
            name: Track name.
            r_lim: Track radius limit region.
            r_pad_ratio: Track padding ratio for plot data.
            parent_sector: Parent sector of track.
        """
        # Track params
        self._name = name
        self._r_lim = r_lim
        self._r_pad_ratio = r_pad_ratio
        # Inherited from parent sector
        self._parent_sector = parent_sector
        self._rad_lim = parent_sector.rad_lim
        self._start = parent_sector.start
        self._end = parent_sector.end

        # Plotly classes
        self._shapes: list[go.layout.Shape] = []
        self._annotations: list[go.layout.Annotation] = []
        self._traces: list[BaseTraceType] = []
        self._trees: list[TreeViz] = []

    ############################################################
    # Property
    ############################################################

    @property
    def name(self) -> str:
        """Track name"""
        return self._name

    @property
    def size(self) -> float:
        """Track size (x coordinate)"""
        return self.end - self.start

    @property
    def precision_position(self) -> int:
        """Track precision position"""
        return max(1, min(10, int(6 - math.log10(self.size))))

    @property
    def start(self) -> float:
        """Track start position (x coordinate)"""
        return self._start

    @property
    def end(self) -> float:
        """Track end position (x coordinate)"""
        return self._end

    @property
    def center(self) -> float:
        """Track center position (x coordinate)"""
        return (self.start + self.end) / 2

    @property
    def r_size(self) -> float:
        """Track radius size"""
        return max(self.r_lim) - min(self.r_lim)

    @property
    def r_lim(self) -> tuple[float, float]:
        """Track radius limit"""
        return self._r_lim

    @property
    def r_center(self) -> float:
        """Track center radius"""
        return sum(self.r_lim) / 2

    @property
    def r_plot_size(self) -> float:
        """Track radius size for plot data (`r_size` with padding)"""
        return max(self.r_plot_lim) - min(self.r_plot_lim)

    @property
    def r_plot_lim(self) -> tuple[float, float]:
        """Track radius limit for plot data (`r_lim` with padding)"""
        edge_pad_size = (self.r_size * self._r_pad_ratio) / 2
        min_plot_r = min(self.r_lim) + edge_pad_size
        max_plot_r = max(self.r_lim) - edge_pad_size
        return (min_plot_r, max_plot_r)

    @property
    def rad_size(self) -> float:
        """Track radian size"""
        return max(self.rad_lim) - min(self.rad_lim)

    @property
    def rad_lim(self) -> tuple[float, float]:
        """Track radian limit"""
        return self._rad_lim

    @property
    def deg_size(self) -> float:
        """Track degree size"""
        return max(self.deg_lim) - min(self.deg_lim)

    @property
    def deg_lim(self) -> tuple[float, float]:
        """Track degree limit"""
        return (math.degrees(min(self.rad_lim)), math.degrees(max(self.rad_lim)))

    @property
    def parent_sector(self) -> Sector:
        """Parent sector"""
        return self._parent_sector

    @property
    def clockwise(self) -> bool:
        """Track coordinate direction"""
        return self.parent_sector.clockwise

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

    def x_to_rad(self, x: float, ignore_range_error: bool = False) -> float:
        """Convert x coordinate to radian in track start-end range.

        Args:
            x: X coordinate.
            ignore_range_error: Ignore x coordinate range error.
        """
        return self.parent_sector.x_to_rad(x, ignore_range_error)

    def axis(self, **kwargs: Any) -> None:
        """
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
        self.rect(self.start, self.end, ignore_pad=True, **fc_behind_kwargs)

        # Edge shape placed in front of other shapes (layer="above")
        ec_front_kwargs = deepcopy(kwargs)
        ec_front_kwargs = utils.deep_dict_update(
            ec_front_kwargs, config.AXIS_EDGE_PARAM
        )
        self.rect(self.start, self.end, ignore_pad=True, **ec_front_kwargs)

    def text(
        self,
        text: str,
        x: float | None = None,
        r: float | None = None,
        *,
        adjust_rotation: bool = True,
        orientation: str = "horizontal",
        ignore_range_error: bool = False,
        outer: bool = True,
        axis: str | bool = False,
        **kwargs: Any,
    ) -> None:
        """Plot text within a track. Uses genomic coordinates (x) mapped to radians.

        Args:
            text: Text content.
            x: Genomic position. If None, track center is used.
            r: Radius position. If None, track midpoint (`r_center`) is used.
            adjust_rotation: If True, text rotation auto based on `x` and `orientation`.
            orientation: Text orientation (`horizontal` or `vertical`).
            ignore_range_error: If True, ignores x position outside track bounds.
            outer: If True, text aligns outward from center, in horizontal orientation.
            axis: Axis type (`x`, `y`, or False).
            **kwargs: Annotation properties.
                <https://plotly.com/python/reference/layout/annotations/>
        """
        x = self.center if x is None else x
        r = self.r_center if r is None else r
        rad = self.x_to_rad(x, ignore_range_error)
        plotly_rad = -(rad - np.pi / 2)
        x_pos = r * np.cos(plotly_rad)
        y_pos = r * np.sin(plotly_rad)

        annotation = utils.plot.get_plotly_label_params(
            rad, adjust_rotation, orientation, **kwargs
        )

        if axis:
            font_size = annotation["font"]["size"]

            if axis == "x":
                # X-axis labels (circular ticks)
                padding = (font_size * 0.05) + (
                    font_size * len(str(text)) * 0.09
                    if orientation == "vertical"
                    else 0
                )
                padding_angle = (
                    plotly_rad + np.pi / 4 - 0.8
                    if outer
                    else plotly_rad - np.pi / 4 + 0.8
                )
                dx = padding * np.cos(padding_angle)
                dy = padding * np.sin(padding_angle)
                if not outer:
                    dx, dy = -dx, -dy
            else:
                # Y-axis labels (radial ticks)
                padding = (font_size * 0.05) + (font_size * len(str(text)) * 0.09)
                tangent_angle = plotly_rad + (np.pi / 2 if outer else -np.pi / 2 + 0.2)
                dx = padding * np.cos(tangent_angle)
                dy = padding * np.sin(tangent_angle)

            x_pos += dx
            y_pos += dy

        annotation.update({"x": x_pos, "y": y_pos, "text": text})
        annotation_layout = go.layout.Annotation(**annotation)
        self._annotations.append(annotation_layout)

    def rect(
        self,
        start: float,
        end: float,
        *,
        r_lim: tuple[float, float] | None = None,
        ignore_pad: bool = False,
        **kwargs: Any,
    ) -> None:
        """Plot a rectangle within a track, respecting padding settings.

        Args:
            start: Genomic start position (x coordinate).
            end: Genomic end position (x coordinate).
            r_lim: Radial limits (min, max). If None, uses track defaults.
            ignore_pad: If True, ignores track padding when auto-setting `r_lim`.
            **kwargs: Shape properties.
                <https://plotly.com/python/reference/layout/shapes/>
        """
        rad_rect_start = self.x_to_rad(start)
        rad_rect_end = self.x_to_rad(end)
        rad = min(rad_rect_start, rad_rect_end)
        width = abs(rad_rect_end - rad_rect_start)

        if r_lim is not None:
            min_range = min(self.r_lim) - config.REL_TOL
            max_range = max(self.r_lim) + config.REL_TOL
            if not min_range <= min(r_lim) < max(r_lim) <= max_range:
                raise ValueError(f"{r_lim=} is invalid track range.\n{self}")
            radr, height = (rad, min(r_lim)), max(r_lim) - min(r_lim)
        elif ignore_pad:
            radr, height = (rad, min(self.r_lim)), self.r_size
        else:
            radr, height = (rad, min(self.r_plot_lim)), self.r_plot_size

        path = PolarSVGPatchBuilder.arc_rectangle(radr, width, height)
        shape = utils.plot.build_plotly_shape(
            path, config.plotly_shape_defaults, **kwargs
        )
        self._shapes.append(shape)

    def arrow(
        self,
        start: float,
        end: float,
        *,
        r_lim: tuple[float, float] | None = None,
        head_length: float = 2,
        shaft_ratio: float = 0.5,
        **kwargs: Any,
    ) -> None:
        """Plot arrow using SVG path.

        Args:
            start: Start position (x coordinate).
            end: End position (x coordinate).
            r_lim: Radius limit range. If None, `track.r_lim` is set.
            head_length: Arrow head length (Degree unit).
            shaft_ratio: Arrow shaft ratio (0 - 1.0).
            **kwargs: Patch properties.
        """
        rad_arrow_start = self.x_to_rad(start)
        rad_arrow_end = self.x_to_rad(end)

        # Handle radius limits
        if r_lim is None:
            r, dr = min(self.r_plot_lim), self.r_plot_size
        else:
            min_range = min(self.r_lim) - config.REL_TOL
            max_range = max(self.r_lim) + config.REL_TOL
            if not min_range <= min(r_lim) < max(r_lim) <= max_range:
                raise ValueError(f"{r_lim=} is invalid track range.\n{self}")
            r, dr = min(r_lim), max(r_lim) - min(r_lim)

        path = PolarSVGPatchBuilder.arc_arrow(
            rad=rad_arrow_start,
            r=r,
            drad=rad_arrow_end - rad_arrow_start,
            dr=dr,
            head_length=math.radians(head_length),
            shaft_ratio=shaft_ratio,
        )
        shape = utils.plot.build_plotly_shape(
            path, defaults={**config.plotly_arrow_defaults, "opacity": 1}, **kwargs
        )
        self._shapes.append(shape)

    def xticks(
        self,
        x: NumericSequence,
        labels: list[str] | None = None,
        *,
        tick_length: float = 2,
        outer: bool = True,
        show_bottom_line: bool = False,
        label_margin: float = 1,
        label_orientation: str = "horizontal",
        line_kws: dict[str, Any] | None = None,
        text_kws: dict[str, Any] | None = None,
    ) -> None:
        """Plot xticks & labels on user-specified position.

        Args:
            x: X coordinates.
            labels: Labels on xticks. If None, only plot ticks line.
            tick_length: Tick length (Radius unit).
            outer: If True, show ticks on outer. If False, show ticks on inner.
            show_bottom_line: If True, show bottom line.
            label_margin: Label margin size.
            label_orientation: Label orientation (`horizontal` or `vertical`).
            line_kws: Shape properties for ticks/baseline.
                <https://plotly.com/python/reference/layout/shapes/>
            text_kws: Annotation properties for labels.
                <https://plotly.com/python/reference/layout/annotations/>
        """
        line_kws = {} if line_kws is None else deepcopy(line_kws)
        text_kws = {} if text_kws is None else deepcopy(text_kws)

        labels = [""] * len(x) if labels is None else labels
        if len(x) != len(labels):
            err_msg = f"List length is not match ({len(x)=}, {len(labels)=})"
            raise ValueError(err_msg)

        r = max(self.r_lim) if outer else min(self.r_lim)
        tick_r_lim = (r, r + tick_length) if outer else (r - tick_length, r)
        for x_pos, label in zip(x, labels):
            # Plot xticks
            if tick_length > 0:
                self._simpleline((x_pos, x_pos), tick_r_lim, **line_kws)
            # Plot labels
            if label != "":
                if outer:
                    adj_r = max(tick_r_lim) + label_margin
                else:
                    adj_r = min(tick_r_lim) - label_margin

                self.text(
                    label,
                    x_pos,
                    adj_r,
                    orientation=label_orientation,
                    outer=outer,
                    axis="x",
                    **text_kws,
                )

        # Plot bottom line
        if show_bottom_line:
            self._simpleline((self.start, self.end), (r, r), **line_kws)

    def xticks_by_interval(
        self,
        interval: Numeric,
        *,
        tick_length: float = 2,
        outer: bool = True,
        show_bottom_line: bool = False,
        show_label: bool = True,
        show_endlabel: bool = True,
        label_margin: float = 1.7,
        label_orientation: str = "horizontal",
        label_formatter: LabelFormatter = None,
        line_kws: dict[str, Any] | None = None,
        text_kws: dict[str, Any] | None = None,
    ) -> None:
        """Plot xticks & position labels by user-specified interval.

        Args:
            interval: Xticks interval in genomic coordinates.
            tick_length: Tick length in radius units.
            outer: If True, show ticks on outer radius.
            show_bottom_line: If True, show baseline at bottom of ticks.
            show_label: If True, show position labels.
            show_endlabel: If False, hides label at final position to prevent overlap.
            label_margin: Additional radial margin for labels.
            label_orientation: Label orientation ('horizontal' or 'vertical').
            label_formatter: Function to format tick labels.
            line_kws: Shape properties for ticks/baseline.
                <https://plotly.com/python/reference/layout/shapes/>
            text_kws: Annotation properties for labels.
                <https://plotly.com/python/reference/layout/annotations/>
        """
        line_kws = {} if line_kws is None else deepcopy(line_kws)
        text_kws = {} if text_kws is None else deepcopy(text_kws)

        x_list = []
        start_pos, end_pos = self.start - (self.start % interval), self.end + interval
        for x in np.arange(start_pos, end_pos, interval):
            if self.start <= x <= self.end:
                x = int(x) if isinstance(interval, int) else float(x)
                x_list.append(x)

        labels = None
        if show_label:
            map_func = str if label_formatter is None else label_formatter
            labels = list(map(map_func, x_list))
            # No display end xtick label if 'show_endlabel' is False
            if not show_endlabel:
                labels[-1] = ""

        self.xticks(
            x=x_list,
            labels=labels,
            tick_length=tick_length,
            outer=outer,
            show_bottom_line=show_bottom_line,
            label_margin=label_margin,
            label_orientation=label_orientation,
            line_kws=line_kws,
            text_kws=text_kws,
        )

    def yticks(
        self,
        y: NumericSequence,
        labels: list[str] | None = None,
        *,
        vmin: Numeric = 0,
        vmax: Numeric | None = None,
        side: str = "right",
        tick_length: Numeric = 2,
        label_margin: Numeric = 1,
        label_orientation: str = "horizontal",
        line_kws: dict[str, Any] | None = None,
        text_kws: dict[str, Any] | None = None,
    ) -> None:
        """Plot yticks & labels on user-specified position.

        Args:
            y: Y coordinates.
            labels: Labels on yticks. If None, only plot ticks line.
            vmin: Y min value.
            vmax: Y max value. If None, `max(y)` is set.
            side: Ticks side position (`right` or `left`).
            tick_length: Tick length (Degree unit).
            label_margin: Label margin size.
            label_orientation: Label orientation (`horizontal` or `vertical`).
            line_kws: Shape properties for ticks/baseline.
                <https://plotly.com/python/reference/layout/shapes/>
            text_kws: Annotation properties for labels.
                <https://plotly.com/python/reference/layout/annotations/>
        """
        line_kws = {} if line_kws is None else deepcopy(line_kws)
        text_kws = {} if text_kws is None else deepcopy(text_kws)

        labels = [""] * len(y) if labels is None else labels
        if len(y) != len(labels):
            err_msg = f"List length is not match ({len(y)=}, {len(labels)=})"
            raise ValueError(err_msg)

        vmax = max(y) if vmax is None else vmax
        self._check_value_min_max(y, vmin, vmax)

        r = [self._y_to_r(v, vmin, vmax) for v in y]
        for r_pos, label in zip(r, labels):
            # Set plot properties
            x_tick_length = (self.size / self.deg_size) * tick_length
            x_label_margin = (self.size / self.deg_size) * label_margin

            if side == "right":
                x_lim = (self.end, self.end + x_tick_length)
                x_text = self.end + (x_tick_length + x_label_margin)
                outer = False
            elif side == "left":
                x_lim = (self.start, self.start - x_tick_length)
                x_text = self.start - (x_tick_length + x_label_margin)
                outer = True
            else:
                raise ValueError(f"{side=} is invalid ('right' or 'left').")

            # Plot yticks
            if tick_length > 0:
                self._simpleline(
                    x_lim, (r_pos, r_pos), ignore_range_error=True, **line_kws
                )

            if label != "":
                self.text(
                    label,
                    x_text,
                    r_pos,
                    orientation=label_orientation,
                    ignore_range_error=True,
                    outer=outer,
                    axis="y",
                    **text_kws,
                )

    def grid(
        self,
        y_grid_num: int | None = 6,
        x_grid_interval: float | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Args:
            y_grid_num: Y-axis grid line number. If None, y-axis grid line is not shown.
            x_grid_interval: X-axis grid line interval.
                If None, x-axis grid line is not shown.
            **kwargs: Shape properties.
                <https://plotly.com/python/reference/layout/shapes/>
        """
        if y_grid_num is not None and not y_grid_num >= 2:
            raise ValueError(f"{y_grid_num=} is invalid (y_grid_num >= 2).")
        if x_grid_interval is not None and not x_grid_interval > 0:
            raise ValueError(f"{x_grid_interval=} is invalid (x_grid_interval > 0).")

        kwargs = utils.deep_dict_update(config.plotly_grid_defaults, kwargs)

        # Plot y-axis grid line
        if y_grid_num is not None:
            y_vmin, y_vmax = 0.0, float(y_grid_num - 1)
            for y_grid_idx in range(y_grid_num):
                x = [self.start, self.end]
                y: list[float] = [float(y_grid_idx), float(y_grid_idx)]
                self.line(x, y, vmin=y_vmin, vmax=y_vmax, hover_text=None, **kwargs)

        # Plot x-axis grid line
        if x_grid_interval is not None:
            x_vmin, x_vmax = 0.0, 1.0
            x_grid_idx = 0
            while True:
                x_pos = self.start + (x_grid_interval * x_grid_idx)
                if x_pos > self.end:
                    break
                x, y = [x_pos, x_pos], [x_vmin, x_vmax]
                self.line(x, y, vmin=x_vmin, vmax=x_vmax, hover_text=None, **kwargs)
                x_grid_idx += 1

    def line(
        self,
        x: NumericSequence,
        y: NumericSequence,
        *,
        vmin: Numeric = 0,
        vmax: Numeric | None = None,
        arc: bool = True,
        hover_text: HoverText = "default",
        **kwargs: Any,
    ) -> None:
        """
        Args:
            x: Genomic positions along the track.
            y: Data values to plot.
            vmin: Minimum value for radial scaling.
            vmax: Maximum value for radial scaling. If None, uses max(y).
            arc: If True, creates curved arc lines.
                If False, creates straight chord lines.
            hover_text: Hover text for the plot.
            **kwargs: Shape properties.
                <https://plotly.com/python/reference/layout/shapes/>
        """
        if len(x) != len(y):
            raise ValueError(f"x and y lengths must match ({len(x)} vs {len(y)})")

        rad = [self.x_to_rad(pos) for pos in x]
        vmax = max(y) if vmax is None else vmax
        r = [self._y_to_r(val, vmin, vmax) for val in y]

        color = utils.plot.get_default_color(kwargs, target="line")
        kwargs = utils.deep_dict_update(kwargs, {"line": {"color": color}})

        path = PolarSVGPatchBuilder.multi_segment_path(rad, r, arc)
        shape = utils.plot.build_plotly_shape(
            path, defaults=config.plotly_shape_defaults, **kwargs
        )
        self._shapes.append(shape)

        if hover_text is None:
            return

        hover_formatter = HoverFormatter(
            precision_position=self.precision_position,
        )

        hover_formatter.process_hover_text(
            x,
            y,
            colors=color,
            sector_name=self._parent_sector.name,
            hover_text=hover_text,
        )
        hover_formatter.add_hover_positions(rad, r)
        if hover_trace := hover_formatter.create_hover_trace():
            self._traces.append(hover_trace)

    def scatter(
        self,
        x: NumericSequence,
        y: NumericSequence,
        *,
        vmin: Numeric = 0,
        vmax: Numeric | None = None,
        hover_text: HoverText = "default",
        **kwargs: Any,
    ) -> None:
        """
        Args:
            x: X (genomic) positions.
            y: Data values.
            vmin: Minimum value for radial scaling.
            vmax: Maximum value for radial scaling. If None, uses max(y).
            hover_text: Hover text for the plot.
            **kwargs: Scatter trace properties.
        """
        if len(x) != len(y):
            raise ValueError(f"x and y lengths must match ({len(x)} vs {len(y)})")

        color = utils.plot.get_default_color(kwargs, target="marker")
        kwargs = utils.deep_dict_update(kwargs, {"line": {"color": color}})
        kwargs.setdefault("hoverlabel", {"bgcolor": color})

        rad = [self.x_to_rad(pos) for pos in x]
        vmax = max(y) if vmax is None else vmax
        r = [self._y_to_r(val, vmin, vmax) for val in y]
        x_vals, y_vals = zip(
            *[
                PolarSVGPatchBuilder._polar_to_cart(theta, rho)
                for theta, rho in zip(rad, r)
            ]
        )

        trace = utils.plot.build_scatter_trace(x_vals, y_vals, "markers", **kwargs)

        if hover_text is None:
            self._traces.append(trace)
            return

        hover_formatter = HoverFormatter(
            precision_position=self.precision_position,
        )
        hover_formatter.process_hover_text(
            x,
            y,
            colors=color,
            sector_name=self._parent_sector.name,
            hover_text=hover_text,
        )
        trace.update(text=hover_formatter.hover_texts)

        self._traces.append(trace)

    def bar(
        self,
        x: NumericSequence,
        height: NumericSequence,
        width: NumericComponent = 0.8,
        bottom: NumericComponent = 0,
        align: str = "center",
        *,
        vmin: Numeric = 0,
        vmax: Numeric | None = None,
        hover_text: HoverText = "default",
        **kwargs: Any,
    ) -> None:
        """
        Args:
            x: Bar x coordinates (genomic positions).
            height: Bar heights.
            width: Bar widths in genomic coordinates.
            bottom: Bar bottom y-value(s).
            align: Bar alignment ("center" or "edge").
            vmin: Minimum value for radial scaling.
            vmax: Maximum value for radial scaling. If None, uses max(height + bottom).
            hover_text: Hover text for the plot.
            **kwargs: Properties for both shapes and hover text.
        """
        if len(x) != len(height):
            raise ValueError(
                f"x and height lengths must match ({len(x)} vs {len(height)})"
            )
        if align not in ["center", "edge"]:
            raise ValueError(f"{align=} must be either 'center' or 'edge'")

        # Convert inputs and calculate metrics
        x, height, bottom = (
            np.asarray(x),
            np.asarray(height),
            np.asarray(bottom) if not np.isscalar(bottom) else np.full(len(x), bottom),
        )
        top = height + bottom
        vmax = float(np.max(top)) if vmax is None else vmax
        if isinstance(width, (int, float)):
            widths = np.full(len(x), width)
        else:
            widths = np.asarray(width)

        if "colors" in kwargs:
            colors = kwargs.pop("colors")
            if len(colors) != len(x):
                raise ValueError("Length of `colors` must match the number of bars.")
            default_line = 1
        else:
            colors = [utils.plot.get_default_color(kwargs, target="fillcolor")] * len(x)
            default_line = 0

        hover_data: dict[str, list[Any]] = {"rad": [], "r": [], "colors": []}

        for i, (xi, hi, bi, color) in enumerate(zip(x, height, bottom, colors)):
            # Calculate bar geometry
            rad = self.x_to_rad(xi)
            rad_width = self.rad_size * (widths[i] / self.size)

            if align == "center":
                rad_start = rad - rad_width / 2
            else:
                rad_start = rad  # edge alignment

            r_bottom = self._y_to_r(bi, vmin, vmax)
            r_height = self._y_to_r(bi + hi, vmin, vmax) - r_bottom

            path = PolarSVGPatchBuilder.arc_rectangle(
                radr=(rad_start, r_bottom),
                width=rad_width,
                height=r_height,
            )
            shape = utils.plot.build_plotly_shape(
                path,
                defaults=dict(
                    fillcolor=color, line=dict(color=color, width=default_line)
                ),
                **kwargs,
            )
            self._shapes.append(shape)

            if hover_data is not None:
                hover_data["rad"].append(rad)
                hover_data["r"].append(r_bottom + r_height)
                hover_data["colors"].append(color)

        if hover_data is None:
            return

        hover_formatter = HoverFormatter(
            precision_position=self.precision_position,
        )
        hover_formatter.add_hover_positions(hover_data["rad"], hover_data["r"])
        hover_formatter.process_hover_text(
            x,
            height,
            colors=hover_data["colors"],
            sector_name=self._parent_sector.name,
            hover_text=hover_text,
        )

        if hover_trace := hover_formatter.create_hover_trace():
            self._traces.append(hover_trace)

    def stacked_bar(
        self,
        table_data: str | Path | pd.DataFrame | StackedBarTable,
        *,
        delimiter: str = "\t",
        width: float = 0.6,
        cmap: str | dict[str, str] = "T10",
        vmax: float | None = None,
        hover_text: HoverText = "default",
        **kwargs: Any,
    ) -> StackedBarTable:
        """
        Args:
            table_data: Table file or Table DataFrame or StackedBarTable.
            delimiter: Table file delimiter.
            width: Bar width ratio (0.0 - 1.0).
            cmap: Colormap assigned to each stacked bar.
            vmax: Stacked bar max value.
                If None, max value in each row values sum is set.
            hover_text: Hover text for the plot.
            **kwargs: Shape properties.
                <https://plotly.com/python/reference/layout/shapes/>
        """
        if not 0.0 <= width <= 1.0:
            raise ValueError(f"{width=} is invalid (0.0 <= width <= 1.0).")

        sb_table = (
            table_data
            if isinstance(table_data, StackedBarTable)
            else StackedBarTable(table_data, delimiter=delimiter)
        )
        col_name2color = (
            sb_table.get_col_name2color(cmap) if isinstance(cmap, str) else cmap
        )

        x = sb_table.calc_bar_label_x_list(self.size)
        bar_width = (self.size / len(sb_table.row_names)) * width
        vmax = sb_table.row_sum_vmax if vmax is None else vmax

        hover_data: dict[str, list[Any]] = {"rad": [], "r": [], "colors": []}

        for col_idx, (col_name, height, bottom) in enumerate(
            zip(
                sb_table.col_names,
                sb_table.stacked_bar_heights,
                sb_table.stacked_bar_bottoms,
            )
        ):
            color = col_name2color[col_name]
            rad = np.array([self.x_to_rad(pos) for pos in x])
            r_bottom = np.array([self._y_to_r(v, 0, vmax) for v in bottom])
            r_height = (
                np.array([self._y_to_r(v + h, 0, vmax) for v, h in zip(bottom, height)])
                - r_bottom
            )

            for i in range(len(x)):
                # Create bar segment
                path = PolarSVGPatchBuilder.arc_rectangle(
                    radr=(
                        rad[i] - self.rad_size * (bar_width / self.size) / 2,
                        r_bottom[i],
                    ),
                    width=self.rad_size * (bar_width / self.size),
                    height=r_height[i],
                )
                shape = utils.plot.build_plotly_shape(
                    path,
                    defaults=dict(fillcolor=color, line=dict(color=color, width=0)),
                    **kwargs,
                )
                self._shapes.append(shape)

                if hover_data is not None:
                    hover_data["rad"].append(rad[i])
                    hover_data["r"].append(r_bottom[i] + r_height[i] / 2)
                    hover_data["colors"].append(color)

        if hover_text is None:
            return sb_table

        hover_formatter = HoverFormatter(
            precision_position=self.precision_position,
        )
        hover_formatter.add_hover_positions(hover_data["rad"], hover_data["r"])
        hover_formatter.hover_colors = hover_data["colors"]
        hover_formatter.process_stacked_bar_hover_text(
            sb_table, col_name2color, hover_text=hover_text
        )

        if hover_trace := hover_formatter.create_hover_trace():
            self._traces.append(hover_trace)

        return sb_table

    def stacked_barh(
        self,
        table_data: str | Path | pd.DataFrame | StackedBarTable,
        *,
        delimiter: str = "\t",
        width: float = 0.6,
        cmap: str | dict[str, str] = "T10",
        vmax: float | None = None,
        hover_text: HoverText = "default",
        **kwargs: Any,
    ) -> StackedBarTable:
        """
        Args:
            table_data: Table file or Table DataFrame or StackedBarTable.
            delimiter: Table file delimiter.
            width: Bar width ratio (0.0 - 1.0).
            cmap: Colormap assigned to each stacked bar.
            vmax: Stacked bar max value. If None, max row sum is used.
            hover_text: Hover text for the plot.
            **kwargs: Shape properties.
                <https://plotly.com/python/reference/layout/shapes/>
        """
        if not 0.0 <= width <= 1.0:
            raise ValueError(f"{width=} is invalid (0.0 <= width <= 1.0).")

        sb_table = (
            table_data
            if isinstance(table_data, StackedBarTable)
            else StackedBarTable(table_data, delimiter=delimiter)
        )
        col_name2color = (
            sb_table.get_col_name2color(cmap) if isinstance(cmap, str) else cmap
        )

        # Calculate bar positions
        r_lim_list = sb_table.calc_barh_r_lim_list(self.r_plot_lim, width)
        vmax = sb_table.row_sum_vmax if vmax is None else vmax

        hover_data: dict[str, list[Any]] = {
            "x_center": [],
            "r_center": [],
            "colors": [],
        }

        # Plot bars and collect hover data
        for col_idx, (col_name, height, bottom) in enumerate(
            zip(
                sb_table.col_names,
                sb_table.stacked_bar_heights,
                sb_table.stacked_bar_bottoms,
            )
        ):
            color = col_name2color[col_name]
            _kwargs = utils.deep_dict_update(
                dict(fillcolor=color, line=dict(color=color, width=0)), kwargs
            )

            for row_idx, (r_lim, h, b) in enumerate(zip(r_lim_list, height, bottom)):
                self.rect(b, b + h, r_lim=r_lim, **_kwargs)

                if hover_data is not None:
                    hover_data["x_center"].append((b + b + h) / 2)
                    hover_data["r_center"].append((r_lim[0] + r_lim[1]) / 2)
                    hover_data["colors"].append(color)

        if hover_text is None:
            return sb_table

        hover_formatter = HoverFormatter(
            precision_position=self.precision_position,
        )

        rad_positions = [self.x_to_rad(x) for x in hover_data["x_center"]]
        hover_formatter.add_hover_positions(rad_positions, hover_data["r_center"])

        hover_formatter.hover_colors = hover_data["colors"]
        hover_formatter.process_stacked_bar_hover_text(
            sb_table, col_name2color, hover_text=hover_text
        )

        if hover_trace := hover_formatter.create_hover_trace():
            self._traces.append(hover_trace)

        return sb_table

    def fill_between(
        self,
        x: NumericSequence,
        y1: NumericSequence,
        y2: NumericComponent = 0,
        *,
        vmin: Numeric = 0,
        vmax: Numeric | None = None,
        arc: bool = True,
        hover_text: HoverText = "default",
        **kwargs: Any,
    ) -> None:
        """Fill the area between two curves with SVG paths.

        Args:
            x: Genomic positions along the track.
            y1: Upper boundary values to plot.
            y2: Lower boundary values or constant baseline.
            vmin: Minimum value for radial scaling.
            vmax: Maximum value for radial scaling. If None, uses max(y1 + y2).
            arc: If True, creates curved arc fills.
                If False, creates straight chord fills.
            hover_text: Hover text for the plot.
            **kwargs: Shape properties.
                <https://plotly.com/python/reference/layout/shapes/>
        """
        # Input validation
        x = np.asarray(x)
        y1 = np.asarray(y1)
        y2 = np.full_like(y1, y2) if isinstance(y2, (int, float)) else np.asarray(y2)

        if len(x) != len(y1) or len(x) != len(y2):
            raise ValueError(
                f"Input lengths must match ({len(x)}, {len(y1)}, {len(y2)})"
            )

        y_all = np.concatenate([y1, y2])
        vmin = y_all.min() if vmin is None else vmin
        vmax = y_all.max() if vmax is None else vmax
        self._check_value_min_max(y_all, vmin, vmax)

        rad = [self.x_to_rad(pos) for pos in x]
        r1 = [self._y_to_r(v, vmin, vmax) for v in y1]
        r2 = [self._y_to_r(v, vmin, vmax) for v in y2]

        color = utils.plot.get_default_color(kwargs, target="fillcolor")
        kwargs.update({"fillcolor": color})

        path = PolarSVGPatchBuilder.build_filled_path(rad, r1, r2, arc=arc)
        shape = utils.plot.build_plotly_shape(
            path, defaults=dict(fillcolor=color, line=dict(width=0)), **kwargs
        )
        self._shapes.append(shape)

        if hover_text is None:
            return

        hover_formatter = HoverFormatter(
            precision_position=self.precision_position,
        )

        hover_formatter.add_hover_positions(rad, r1)
        hover_formatter.process_hover_text(
            x,
            y1,
            colors=color,
            sector_name=self._parent_sector.name,
            hover_text=hover_text,
        )

        if hover_trace := hover_formatter.create_hover_trace():
            self._traces.append(hover_trace)

    def heatmap(
        self,
        data: list | np.ndarray,
        *,
        vmin: Numeric | None = None,
        vmax: Numeric | None = None,
        start: Numeric | None = None,
        end: Numeric | None = None,
        width: Numeric | None = None,
        cmap: str | list[tuple[float, str]] = "RdBu_r",
        show_value: bool = False,
        hover_text: HoverText = "default",
        coloraxis: str | None = None,
        rect_kws: dict[str, Any] | None = None,
        text_kws: dict[str, Any] | None = None,
    ) -> None:
        """
        Args:
            data: Numerical list, numpy 1d or 2d array.
            vmin: Min value for heatmap plot. If None, `np.min(data)` is set.
            vmax: Max value for heatmap plot. If None, `np.max(data)` is set.
            start: Start position for heatmap plot (x coordinate).
            end: End position for heatmap plot (x coordinate).
            width: Heatmap rectangle x width size.
            cmap: Colormap specification.
            show_value: If True, show data value on heatmap rectangle.
            hover_text: Hover text for the plot.
            coloraxis: Color axis identifier.
            rect_kws: Shape properties.
                <https://plotly.com/python/reference/layout/shapes/>
            text_kws: Text properties.
                <https://plotly.com/python/reference/layout/annotations/>
        """
        rect_kws = {} if rect_kws is None else deepcopy(rect_kws)
        text_kws = {} if text_kws is None else deepcopy(text_kws)

        data = np.array(data)
        if data.ndim == 1:
            data = data.reshape((1, -1))
        elif data.ndim != 2:
            raise ValueError(f"{data=} is not 1d or 2d array!!")

        vmin = np.min(data) if vmin is None else vmin
        vmax = np.max(data) if vmax is None else vmax
        start = self.start if start is None else start
        end = self.end if end is None else end
        self._check_value_min_max(data, vmin, vmax)

        row_num, col_num = data.shape
        unit_r_size = self.r_plot_size / row_num
        unit_x_size = (end - start) / col_num

        if width is not None:
            if (col_num - 1) * width < end - start < col_num * width:
                unit_x_size = width
            else:
                raise ValueError(f"{width=} is invalid ({start=}, {end=})")

        r_range_list = [
            (
                max(self.r_plot_lim) - unit_r_size * (i + 1),
                max(self.r_plot_lim) - unit_r_size * i,
            )
            for i in range(row_num)
        ]

        x_range_list = []
        for i in range(col_num):
            min_x = start + unit_x_size * i
            max_x = min(min_x + unit_x_size, self.end)
            x_range_list.append((min_x, max_x))

        color_scale = (
            get_colorscale(cmap)
            if isinstance(cmap, str)
            else [[pos, utils.parse_color(color)] for pos, color in cmap]
        )
        norm = utils.plot.Normalize(vmin=vmin, vmax=vmax)

        hover_data: dict[str, list[Any]] = {
            "rad": [],
            "r": [],
            "colors": [],
            "start_x": [],
            "end_x": [],
            "values": [],
        }

        # Plot rectangles and collect hover data
        for row_idx, row in enumerate(data):
            for col_idx, v in enumerate(row):
                rect_start, rect_end = x_range_list[col_idx]
                rect_r_lim = r_range_list[row_idx]
                color = sample_colorscale(color_scale, norm(v))[0]

                rect_kws = utils.deep_dict_update(rect_kws, {"fillcolor": color})
                self.rect(rect_start, rect_end, r_lim=rect_r_lim, **rect_kws)

                if hover_data is not None:
                    center_x = (rect_start + rect_end) / 2
                    center_r = sum(rect_r_lim) / 2
                    hover_data["rad"].append(self.x_to_rad(center_x))
                    hover_data["r"].append(center_r)
                    hover_data["colors"].append(color)
                    hover_data["start_x"].append(rect_start)
                    hover_data["end_x"].append(rect_end)
                    hover_data["values"].append(v)

                if show_value:
                    text_value = f"{v:.2f}" if isinstance(v, float) else str(v)
                    self.text(text_value, center_x, center_r, **text_kws)

        if hover_text is None:
            return

        hover_formatter = HoverFormatter(
            precision_position=self.precision_position,
        )

        hover_formatter.add_hover_positions(hover_data["rad"], hover_data["r"])
        hover_formatter.process_hover_text(
            hover_data["start_x"],
            hover_data["values"],
            hover_data["colors"],
            hover_data["end_x"],
            sector_name=self._parent_sector.name,
            hover_text=hover_text,
        )

        if hover_trace := hover_formatter.create_hover_trace(
            marker=dict(
                colorscale=color_scale,
                cmin=vmin,
                cmax=vmax,
                coloraxis=coloraxis if coloraxis else None,
                showscale=False,
            )
        ):
            self._traces.append(hover_trace)

    def tree(
        self,
        tree_data: str | Path | Tree,
        *,
        format: str = "newick",
        outer: bool = True,
        align_leaf_label: bool = True,
        ignore_branch_length: bool = False,
        leaf_label_size: float = 12,
        leaf_label_rmargin: float = 2.0,
        reverse: bool = False,
        ladderize: bool = False,
        line_kws: dict[str, Any] | None = None,
        align_line_kws: dict[str, Any] | None = None,
        label_formatter: TextFormatter = None,
    ) -> TreeViz:
        """
        Args:
            tree_data: Tree data (`File`|`File URL`|`Tree Object`|`Tree String`).
            format: Tree format (`newick`|`phyloxml`|`nexus`|`nexml`|`cdao`).
            outer: If True, plot tree on outer side. If False, plot tree on inner side.
            align_leaf_label: If True, align leaf label.
            ignore_branch_length: If True, ignore branch length for plotting tree.
            leaf_label_size: Leaf label size.
            leaf_label_rmargin: Leaf label radius margin.
            reverse: If True, reverse tree.
            ladderize: If True, ladderize tree.
            line_kws: Shape properties.
                <https://plotly.com/python/reference/layout/shapes/>
            align_line_kws: Shape properties.
                <https://plotly.com/python/reference/layout/shapes/>
            label_formatter:
                User-defined label text format function to change label text content.
        """
        tv = TreeViz(
            tree_data,
            format=format,
            outer=outer,
            align_leaf_label=align_leaf_label,
            ignore_branch_length=ignore_branch_length,
            leaf_label_size=leaf_label_size,
            leaf_label_rmargin=leaf_label_rmargin,
            reverse=reverse,
            ladderize=ladderize,
            line_kws=line_kws,
            align_line_kws=align_line_kws,
            label_formatter=label_formatter,
            track=self,
        )
        self._trees.append(tv)

        return tv

    def genomic_features(
        self,
        features: SeqFeature | list[SeqFeature],
        *,
        plotstyle: str = "box",
        r_lim: tuple[float, float] | None = None,
        hover_text_formatter: Callable[[SeqFeature], str] | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Args:
            features: Biopython's SeqFeature or SeqFeature list.
            plotstyle: Plot style (`box` or `arrow`).
            r_lim: Radius limit range. If None, `track.r_plot_lim` is set.
            hover_text_formatter: User-defined function for hover text format.
            **kwargs: Shape properties.
                <https://plotly.com/python/reference/layout/shapes/>
        """
        if isinstance(features, SeqFeature):
            features = [features]

        if r_lim is None:
            r_lim = self.r_plot_lim
        else:
            if not min(self.r_lim) <= min(r_lim) < max(r_lim) <= max(self.r_lim):
                raise ValueError(f"{r_lim=} is invalid track range.\n{self}")

        # Default feature hover color hierarchy
        bg_color = (
            kwargs.get("fillcolor", None)
            or (kwargs.get("line", {}).get("color") if "line" in kwargs else None)
            or "lightgrey"
        )

        hover_data: dict[str, list[Any]] = {"rad": [], "r": [], "texts": []}

        for feature in features:
            try:
                start = int(str(feature.location.parts[0].start))
                end = int(str(feature.location.parts[-1].end))
            except ValueError:
                print(f"Failed to parse feature's start-end position.\n{feature}")
                continue

            # Reverse coordinates if negative strand
            if feature.location.strand == -1:
                start, end = end, start

            # Plot the feature with style
            if plotstyle == "box":
                self.rect(start, end, r_lim=r_lim, **kwargs)
            elif plotstyle == "arrow":
                self.arrow(start, end, r_lim=r_lim, **kwargs)
            else:
                raise ValueError(f"{plotstyle=} is invalid ('box' or 'arrow').")

            if hover_text_formatter:
                midpoint = (start + end) / 2
                hover_data["rad"].append(self.x_to_rad(midpoint))
                hover_data["r"].append(sum(r_lim) / 2)
                hover_data["texts"].append(hover_text_formatter(feature))

        if hover_text_formatter is None:
            return

        hover_formatter = HoverFormatter(
            precision_position=self.precision_position,
        )

        hover_formatter.add_hover_positions(hover_data["rad"], hover_data["r"])
        hover_formatter.hover_texts = hover_data["texts"]
        hover_formatter.hover_colors = [bg_color] * len(hover_data["rad"])

        if hover_trace := hover_formatter.create_hover_trace():
            self._traces.append(hover_trace)

    ############################################################
    # Private Method
    ############################################################

    def _y_to_r(self, y: float, vmin: float, vmax: float) -> float:
        """Convert y coordinate to radius in track."""
        norm = utils.plot.Normalize(vmin, vmax)
        r = min(self.r_plot_lim) + (self.r_plot_size * norm(y))
        return r

    def _simpleline(
        self,
        x_lim: tuple[float, float],
        r_lim: tuple[float, float],
        ignore_range_error=False,
        **kwargs: Any,
    ) -> None:
        """Plot a line between two positions at specified radial distances."""
        # Convert genomic positions to radians
        rad_start = self.x_to_rad(x_lim[0], ignore_range_error=ignore_range_error)
        rad_end = self.x_to_rad(x_lim[1], ignore_range_error=ignore_range_error)

        # Convert to Plotly's coordinate system (0=up, clockwise)
        plotly_rad_start = -(rad_start - np.pi / 2)
        plotly_rad_end = -(rad_end - np.pi / 2)

        if np.isclose(plotly_rad_start, plotly_rad_end):
            # Special case: straight radial line (for ticks in axis)
            x0 = r_lim[0] * np.cos(plotly_rad_start)
            y0 = r_lim[0] * np.sin(plotly_rad_start)
            x1 = r_lim[1] * np.cos(plotly_rad_end)
            y1 = r_lim[1] * np.sin(plotly_rad_end)

            path = f"M {x0} {y0} L {x1} {y1}"
        else:
            # General arc line
            path = PolarSVGPatchBuilder.arc_line(
                rad_lim=(rad_start, rad_end), r_lim=r_lim
            )

        # Build and add shape
        shape = utils.plot.build_plotly_shape(
            path, config.plotly_shape_defaults, **kwargs
        )
        self._shapes.append(shape)

    def _check_value_min_max(
        self,
        value: float | NumericSequence,
        vmin: float,
        vmax: float,
    ) -> None:
        """Check if value(s) is in valid min-max range."""
        if isinstance(value, (list, tuple, np.ndarray)):
            if isinstance(value, np.ndarray):
                value = list(value.flatten())
            for v in value:
                if not vmin <= v <= vmax:
                    err_msg = f"value={v} is not in valid range ({vmin=}, {vmax=})"
                    raise ValueError(err_msg)
        else:
            float_value: float = value  # type: ignore[assignment]
            if not vmin <= float_value <= vmax:
                err_msg = f"{value=} is not in valid range ({vmin=}, {vmax=})"
                raise ValueError(err_msg)

    def __str__(self):
        min_deg_lim, max_deg_lim = min(self.deg_lim), max(self.deg_lim)
        min_r_lim, max_r_lim = min(self.r_lim), max(self.r_lim)
        return textwrap.dedent(
            f"""
            # Track = '{self.name}' (Parent Sector = '{self.parent_sector.name}')
            # Size = {self.size} ({self.start} - {self.end})
            # Degree Size = {self.deg_size:.2f} ({min_deg_lim:.2f} - {max_deg_lim:.2f})
            # Radius Size = {self.r_size:.2f} ({min_r_lim:.2f} - {max_r_lim:.2f})
            """
        )[1:]
