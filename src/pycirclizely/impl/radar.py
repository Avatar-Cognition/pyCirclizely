from __future__ import annotations

import math
from copy import deepcopy
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Optional, Union

import numpy as np
import pandas as pd

from pycirclizely import config, utils
from pycirclizely.parser.table import RadarTable
from pycirclizely.types import LabelFormatter, Numeric

if TYPE_CHECKING:
    from ..circos import Circos


def radar_chart_impl(
    table: Union[str, Path, pd.DataFrame, RadarTable],
    *,
    r_lim: tuple[float, float] = (0, 100),
    vmin: Numeric = 0,
    vmax: Numeric = 100,
    fill: bool = True,
    marker_size: int = 0,
    bg_color: Optional[str] = "#eeeeee80",
    circular: bool = False,
    cmap: Union[str, dict[str, str]] = "Set2",
    show_grid_label: bool = True,
    grid_interval_ratio: Optional[float] = 0.2,
    grid_line_kws: Optional[dict[str, Any]] = None,
    grid_label_kws: Optional[dict[str, Any]] = None,
    grid_label_formatter: LabelFormatter = None,
    label_kws_handler: Optional[Callable[[str], dict[str, Any]]] = None,
    line_kws_handler: Optional[Callable[[str], dict[str, Any]]] = None,
    marker_kws_handler: Optional[Callable[[str], dict[str, Any]]] = None,
) -> Circos:
    """Radar chart implementation for visualizing multi-dimensional data."""
    from ..circos import Circos  # Local import to avoid circular import

    if not vmin < vmax:
        raise ValueError(f"vmax must be larger than vmin ({vmin=}, {vmax=})")
    size = vmax - vmin

    grid_line_kws = {} if grid_line_kws is None else deepcopy(grid_line_kws)
    grid_line_kws = utils.deep_dict_update(config.radar_grid_defaults, grid_line_kws)

    grid_label_kws = {} if grid_label_kws is None else deepcopy(grid_label_kws)
    grid_label_kws = utils.deep_dict_update(
        config.radar_annotation_defaults, grid_label_kws
    )

    radar_table = table if isinstance(table, RadarTable) else RadarTable(table)
    circos = Circos(dict(radar=radar_table.col_num))
    sector = circos.sectors[0]
    track = sector.add_track(r_lim)
    x = np.arange(radar_table.col_num + 1)

    # Plot background color
    if bg_color:
        track.fill_between(
            x,
            [vmax] * len(x),
            arc=circular,
            hover_text=None,
            fillcolor=utils.parse_color(bg_color),
        )

    # Plot grid line
    if grid_interval_ratio:
        if not 0 < grid_interval_ratio <= 1.0:
            raise ValueError(f"{grid_interval_ratio=} is invalid.")
        # Plot horizontal grid line & label
        stop, step = vmax + (size / 1000), size * grid_interval_ratio
        for v in np.arange(vmin, stop, step):
            y = [v] * len(x)
            track.line(
                x,
                y,
                vmin=vmin,
                vmax=vmax,
                arc=circular,
                hover_text=None,
                **grid_line_kws,
            )
            if show_grid_label:
                r = track._y_to_r(v, vmin, vmax)
                # Format grid label
                if grid_label_formatter:
                    text = grid_label_formatter(v)
                else:
                    v = float(f"{v:.9f}")  # Correct rounding error
                    text = f"{v:.0f}" if math.isclose(int(v), float(v)) else str(v)
                track.text(text, 0, r, **grid_label_kws)
        # Plot vertical grid line
        for p in x[:-1]:
            track.line(
                [p, p],
                [vmin, vmax],
                vmin=vmin,
                vmax=vmax,
                hover_text=None,
                **grid_line_kws,
            )

    # Plot radar charts
    row_name2color = (
        radar_table.get_row_name2color(cmap) if isinstance(cmap, str) else cmap
    )

    for row_name, values in radar_table.row_name2values.items():
        y = values + [values[0]]
        color = row_name2color[row_name]

        # Create hover_text
        hover_texts = []
        for idx, (col_name, value) in enumerate(zip(radar_table.col_names, values)):
            hover_texts.append(
                f"Class: {row_name}<br>Feature: {col_name}<br>Value: {value:.2f}"
            )
        # Add the first point again to close the polygon
        hover_texts.append(hover_texts[0])

        # Plot line
        line_kws = line_kws_handler(row_name) if line_kws_handler else {}
        defaults = utils.deep_dict_update(
            config.plotly_shape_defaults, dict(line=dict(color=color))
        )
        line_kws = utils.deep_dict_update(defaults, line_kws)

        track.line(
            x,
            y,
            vmin=vmin,
            vmax=vmax,
            arc=False,
            hover_text=hover_texts,
            **line_kws,
        )

        # Plot markers
        if marker_size > 0:
            marker_kws = marker_kws_handler(row_name) if marker_kws_handler else {}
            defaults = utils.deep_dict_update(
                config.plotly_scatter_defaults,
                dict(marker=dict(size=marker_size, color=color)),
            )
            marker_kws = utils.deep_dict_update(defaults, marker_kws)
            track.scatter(
                x,
                y,
                vmin=vmin,
                vmax=vmax,
                hovertext=None,
                fillcolor=color,
                **marker_kws,
            )

        # Fill area under the radar chart
        if fill:
            track.fill_between(
                x,
                y,
                y2=vmin,
                vmin=vmin,
                vmax=vmax,
                arc=False,
                fillcolor=color,
                hover_text=None,
                opacity=0.3,
            )

    # Plot column names
    for idx, col_name in enumerate(radar_table.col_names):
        deg = 360 * (idx / sector.size)
        label_kws = label_kws_handler(col_name) if label_kws_handler else {}
        label_kws = utils.deep_dict_update(config.plotly_annotation_defaults, label_kws)
        if math.isclose(deg, 0):
            label_kws.update(yanchor="bottom")
        elif math.isclose(deg, 180):
            label_kws.update(yanchor="top")
        elif 0 < deg < 180:
            label_kws.update(xanchor="left")
        elif 180 < deg < 360:
            label_kws.update(xanchor="right")
        track.text(col_name, idx, r=105, adjust_rotation=False, **label_kws)

    return circos
