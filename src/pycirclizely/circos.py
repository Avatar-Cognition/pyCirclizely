from __future__ import annotations

import itertools
import math
import textwrap
from collections import defaultdict
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from Bio.Phylo.BaseTree import Tree
from plotly.basedatatypes import BaseTraceType

from pycirclizely import config, utils
from pycirclizely.parser import Bed
from pycirclizely.parser.matrix import Matrix
from pycirclizely.parser.table import RadarTable
from pycirclizely.patches import PolarSVGPatchBuilder
from pycirclizely.sector import Sector
from pycirclizely.track import Track
from pycirclizely.tree import TreeViz
from pycirclizely.utils.plot import LinkDirection


class Circos:
    """Circos Visualization Class"""

    def __init__(
        self,
        sectors: Mapping[str, int | float | tuple[float, float]],
        start: float = 0,
        end: float = 360,
        *,
        space: float | list[float] = 0,
        endspace: bool = True,
        sector2clockwise: dict[str, bool] | None = None,
        show_axis_for_debug: bool = False,
    ):
        """Parameters
        ----------
        sectors : Mapping[str, int | float | tuple[float, float]]
            Sector name & size (or range) dict
        start : float, optional
            Plot start degree (`-360 <= start < end <= 360`)
        end : float, optional
            Plot end degree (`-360 <= start < end <= 360`)
        space : float | list[float], optional
            Space degree(s) between sector
        endspace : bool, optional
            If True, insert space after the end sector
        sector2clockwise : dict[str, bool] | None, optional
            Sector name & clockwise bool dict. By default, `clockwise=True`.
        show_axis_for_debug : bool, optional
            Show axis for position check debugging (Developer option)

        """
        sector2clockwise = {} if sector2clockwise is None else sector2clockwise

        # Check start-end degree range
        self._check_degree_range(start, end)

        # Calculate sector region & add sector
        whole_deg_size = end - start
        space_num = len(sectors) if endspace else len(sectors) - 1
        if isinstance(space, (list, tuple)):
            if len(space) != space_num:
                err_msg = f"{space=} is invalid.\n"
                err_msg += f"Length of space list must be {space_num}."
                raise ValueError(err_msg)
            space_list = list(space) + [0]
            space_deg_size = sum(space)
        else:
            space_list = [space] * space_num + [0]
            space_deg_size = space * space_num
        whole_deg_size_without_space = whole_deg_size - space_deg_size
        if whole_deg_size_without_space < 0:
            err_msg = textwrap.dedent(
                f"""
                Too large sector space size is set!!
                Circos Degree Size = {whole_deg_size} ({start} - {end})
                """
                # Total Sector Space Size = {space_deg_size}
                # List of Sector Space Size = {space_list}
            )[1:-1]
            raise ValueError(err_msg)

        sector2range = self._to_sector2range(sectors)
        sector_total_size = sum([max(r) - min(r) for r in sector2range.values()])

        rad_pos = math.radians(start)

        self._sectors: list[Sector] = []
        for idx, (sector_name, sector_range) in enumerate(sector2range.items()):
            sector_size = max(sector_range) - min(sector_range)
            sector_size_ratio = sector_size / sector_total_size
            deg_size = whole_deg_size_without_space * sector_size_ratio
            rad_size = math.radians(deg_size)
            rad_lim = (rad_pos, rad_pos + rad_size)
            rad_pos += rad_size + math.radians(space_list[idx])
            clockwise = sector2clockwise.get(sector_name, True)
            sector = Sector(sector_name, sector_range, rad_lim, clockwise)
            self._sectors.append(sector)

        self._deg_lim = (start, end)
        self._rad_lim = (math.radians(start), math.radians(end))
        self._show_axis_for_debug = show_axis_for_debug

        # Plotly classes
        self._shapes: list[go.layout.Shape] = []
        self._annotations: list[go.layout.Annotation] = []
        self._traces: list[BaseTraceType] = []
        self._coloraxes: list[go.layout.Coloraxis] = []

    ############################################################
    # Property
    ############################################################

    @property
    def rad_size(self) -> float:
        """Circos radian size"""
        return max(self.rad_lim) - min(self.rad_lim)

    @property
    def rad_lim(self) -> tuple[float, float]:
        """Circos radian limit"""
        return self._rad_lim

    @property
    def deg_size(self) -> float:
        """Circos degree size"""
        return max(self.deg_lim) - min(self.deg_lim)

    @property
    def deg_lim(self) -> tuple[float, float]:
        """Circos degree limit"""
        return self._deg_lim

    @property
    def sectors(self) -> list[Sector]:
        """Sectors"""
        return self._sectors

    @property
    def tracks(self) -> list[Track]:
        """Tracks (from sectors)"""
        tracks = []
        for sector in self.sectors:
            for track in sector.tracks:
                tracks.append(track)
        return tracks

    ############################################################
    # Public Method
    ############################################################

    @staticmethod
    def radar_chart(
        table: str | Path | pd.DataFrame | RadarTable,
        *,
        r_lim: tuple[float, float] = (0, 100),
        vmin: float = 0,
        vmax: float = 100,
        fill: bool = True,
        marker_size: int = 0,
        bg_color: str | None = "#eeeeee80",
        circular: bool = False,
        cmap: str | dict[str, str] = "Set2",
        show_grid_label: bool = True,
        grid_interval_ratio: float | None = 0.2,
        grid_line_kws: dict[str, Any] | None = None,
        grid_label_kws: dict[str, Any] | None = None,
        grid_label_formatter: Callable[[float], str] | None = None,
        label_kws_handler: Callable[[str], dict[str, Any]] | None = None,
        line_kws_handler: Callable[[str], dict[str, Any]] | None = None,
        marker_kws_handler: Callable[[str], dict[str, Any]] | None = None,
    ) -> Circos:
        """Plot radar chart

        Parameters
        ----------
        table : str | Path | pd.DataFrame | RadarTable
            Table file or Table dataframe or RadarTable instance
        r_lim : tuple[float, float], optional
            Radar chart radius limit region (0 - 100)
        vmin : float, optional
            Min value
        vmax : float, optional
            Max value
        fill : bool, optional
            If True, fill color of radar chart.
        marker_size : int, optional
            Marker size
        bg_color : str | None, optional
            Background color
        circular : bool, optional
            If True, plot with circular style.
        cmap : str | dict[str, str], optional
            Colormap assigned to each target row(index) in table.
            User can set plotly's built-in palettes (e.g. `T10`, `Set2`) or
            target_name -> color dict (e.g. `dict(A="red", B="blue", C="green", ...)`)
        show_grid_label : bool, optional
            If True, show grid label.
        grid_interval_ratio : float | None, optional
            Grid interval ratio (0.0 - 1.0)
        grid_line_kws : dict[str, Any] | None, optional
            Keyword arguments passed to `track.line()` method
            (e.g. `dict(line=dict(color="black", dash="dot", width=1.0), ...)`)
        grid_label_kws : dict[str, Any] | None, optional
            Keyword arguments passed to `track.text()` method
            (e.g. `dict(font=dict(size=12, color="red"), ...)`)
        grid_label_formatter : Callable[[float], str] | None, optional
            User-defined function to format grid label (e.g. `lambda v: f"{v:.1f}%"`).
        label_kws_handler : Callable[[str], dict[str, Any]] | None, optional
            Handler function for keyword arguments passed to `track.text()` method.
            Handler function takes each column name of table as an argument.
        line_kws_handler : Callable[[str], dict[str, Any]] | None, optional
            Handler function for keyword arguments passed to `track.line()` method.
            Handler function takes each row(index) name of table as an argument.
        marker_kws_handler : Callable[[str], dict[str, Any]] | None, optional
            Handler function for keyword arguments passed to `track.scatter()` method.
            Handler function takes each row(index) name of table as an argument.

        Returns
        -------
        circos : Circos
            Circos instance initialized for radar chart
        """
        if not vmin < vmax:
            raise ValueError(f"vmax must be larger than vmin ({vmin=}, {vmax=})")
        size = vmax - vmin

        # Setup default properties
        grid_line_kws = {} if grid_line_kws is None else deepcopy(grid_line_kws)
        grid_line_kws = utils.deep_dict_update(
            config.radar_grid_defaults, grid_line_kws
        )

        grid_label_kws = {} if grid_label_kws is None else deepcopy(grid_label_kws)
        grid_label_kws = utils.deep_dict_update(
            config.radar_annotation_defaults, grid_label_kws
        )

        # Initialize circos for radar chart
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
        if isinstance(cmap, str):
            row_name2color = radar_table.get_row_name2color(cmap)
        else:
            row_name2color = cmap

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
            label_kws = utils.deep_dict_update(
                config.plotly_annotation_defaults, label_kws
            )
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

    @staticmethod
    def chord_diagram(
        matrix: str | Path | pd.DataFrame | Matrix,
        *,
        start: float = 0,
        end: float = 360,
        space: float | list[float] = 0,
        endspace: bool = True,
        r_lim: tuple[float, float] = (97, 100),
        cmap: str | dict[str, str] = "Viridis",
        link_cmap: list[tuple[str, str, str]] | None = None,
        ticks_interval: int | None = None,
        order: str | list[str] | None = None,
        label_kws: dict[str, Any] | None = None,
        ticks_kws: dict[str, Any] | None = None,
        link_kws: dict[str, Any] | None = None,
        link_kws_handler: Callable[[str, str], dict[str, Any] | None] | None = None,
    ) -> Circos:
        """Plot chord diagram

        Circos tracks and links are auto-defined from Matrix

        Parameters
        ----------
        matrix : str | Path | pd.DataFrame | Matrix
            Matrix file or Matrix dataframe or Matrix instance
        start : float, optional
            Plot start degree (-360 <= start < end <= 360)
        end : float, optional
            Plot end degree (-360 <= start < end <= 360)
        space : float | list[float], optional
            Space degree(s) between sector
        endspace : bool, optional
            If True, insert space after the end sector
        r_lim : tuple[float, float], optional
            Outer track radius limit region (0 - 100)
        cmap : str | dict[str, str], optional
            Colormap assigned to each outer track and link.
            User can set plotly's colormap (e.g. `Viridis`, `T10`) or
            label_name -> color dict (e.g. `dict(A="red", B="blue", C="green", ...)`)
        link_cmap : list[tuple[str, str, str]] | None, optional
            Link colormap to overwrite link colors automatically set by cmap.
            User can set list of `from_label`, `to_label`, `color` tuple
            (e.g. `[("A", "B", "red"), ("A", "C", "#ffff00"), ...]`)
        ticks_interval : int | None, optional
            Ticks interval. If None, ticks are not plotted.
        order : str | list[str] | None, optional
            Sort order of matrix for plotting Chord Diagram. If `None`, no sorting.
            If `asc`|`desc`, sort in ascending(or descending) order by node size.
            If node name list is set, sort in user specified node order.
        label_kws : dict[str, Any] | None, optional
            Keyword arguments passed to `sector.text()` method
            (e.g. `dict(r=110, orientation="vertical", font=dict(size=15), ...)`)
        ticks_kws : dict[str, Any] | None, optional
            Keyword arguments passed to `track.xticks_by_interval()` method
            (e.g. `dict(text_kws=dict(font=dict(size=10)),
                label_orientation="vertical", ...)`)
        link_kws : dict[str, Any] | None, optional
            Keyword arguments passed to `circos.link()` method
            (e.g. `dict(direction=1, line=dict(color="black", width=0.5),
                opacity=0.8, ...)`)
        link_kws_handler : Callable[[str, str], dict[str, Any] | None] | None, optional
            User-defined function to handle keyword arguments for each link.
            This option allows user to set or override properties such as
            `fillcolor`, `opacity`, `layer`, etc... on each link.
            Handler function arguments `[str, str]` means `[from_label, to_label]`.

        Returns
        -------
        circos : Circos
            Circos instance initialized from Matrix
        """
        link_cmap = [] if link_cmap is None else deepcopy(link_cmap)
        label_kws = {} if label_kws is None else deepcopy(label_kws)
        ticks_kws = {} if ticks_kws is None else deepcopy(ticks_kws)
        link_kws = {} if link_kws is None else deepcopy(link_kws)

        # If input matrix is file path, convert to Matrix instance
        if isinstance(matrix, (str, Path, pd.DataFrame)):
            matrix = Matrix(matrix)

        # Sort matrix if order is set
        if order is not None:
            matrix = matrix.sort(order)

        # Get name2color dict from user-specified colormap
        names = matrix.all_names
        name2color: dict[str, str]
        if isinstance(cmap, str):
            utils.ColorCycler.set_palette(cmap)
            colors = utils.ColorCycler.get_color_list(len(names))
            name2color = dict(zip(names, colors))
        else:
            if isinstance(cmap, defaultdict):
                name2color = cmap
            else:
                name2color = defaultdict(lambda: "grey")
                name2color.update(cmap)

        # Initialize circos sectors
        circos = Circos(matrix.to_sectors(), start, end, space=space, endspace=endspace)
        for sector in circos.sectors:
            # Plot label, outer track axis & xticks
            sector.text(sector.name, **label_kws)
            outer_track = sector.add_track(r_lim)
            color = name2color[sector.name]
            outer_track.axis(fillcolor=color)
            if ticks_interval is not None:
                outer_track.xticks_by_interval(ticks_interval, **ticks_kws)

        # Plot links
        fromto_label2color = {f"{t[0]}-->{t[1]}": t[2] for t in link_cmap}
        for link in matrix.to_links():
            from_label, to_label = link[0][0], link[1][0]
            fromto_label = f"{from_label}-->{to_label}"
            # Set link color
            if fromto_label in fromto_label2color:
                color = fromto_label2color[fromto_label]
            else:
                color = name2color[from_label]
            # Update link properties by user-defined handler function
            _link_kws = deepcopy(link_kws)
            _link_kws.update(fillcolor=color)
            if link_kws_handler is not None:
                handle_link_kws = link_kws_handler(from_label, to_label)
                if handle_link_kws is not None:
                    _link_kws.update(handle_link_kws)
            circos.link(*link, **_link_kws)

        return circos

    @staticmethod
    def initialize_from_tree(
        tree_data: str | Path | Tree,
        *,
        start: float = 0,
        end: float = 360,
        r_lim: tuple[float, float] = (50, 100),
        format: str = "newick",
        outer: bool = True,
        align_leaf_label: bool = True,
        ignore_branch_length: bool = False,
        leaf_label_size: float = 14,
        leaf_label_rmargin: float = 2.0,
        reverse: bool = False,
        ladderize: bool = False,
        line_kws: dict[str, Any] | None = None,
        label_formatter: Callable[[str], str] | None = None,
        align_line_kws: dict[str, Any] | None = None,
    ) -> tuple[Circos, TreeViz]:
        """Initialize Circos instance from phylogenetic tree

        Circos sector and track are auto-defined by phylogenetic tree

        Parameters
        ----------
        tree_data : str | Path | Tree
            Tree data (`File`|`File URL`|`Tree Object`|`Tree String`)
        start : float, optional
            Plot start degree (-360 <= start < end <= 360)
        end : float, optional
            Plot end degree (-360 <= start < end <= 360)
        r_lim : tuple[float, float], optional
            Tree track radius limit region (0 - 100)
        format : str, optional
            Tree format (`newick`|`phyloxml`|`nexus`|`nexml`|`cdao`)
        outer : bool, optional
            If True, plot tree on outer side. If False, plot tree on inner side.
        align_leaf_label: bool, optional
            If True, align leaf label.
        ignore_branch_length : bool, optional
            If True, ignore branch length for plotting tree.
        leaf_label_size : float, optional
            Leaf label size
        leaf_label_rmargin : float, optional
            Leaf label radius margin
        reverse : bool, optional
            If True, reverse tree
        ladderize : bool, optional
            If True, ladderize tree
        line_kws : dict[str, Any] | None, optional
            Shape properties (default: None)
            e.g. `dict(line=dict(color="red", width=1, dash="dash"))`
            See: <https://plotly.com/python/reference/layout/shapes/>
        align_line_kws : dict[str, Any] | None, optional
            Shape properties (default: None)
            e.g. `dict(line=dict(color="black", dash="dot"), opacity=0.5)`
            See: <https://plotly.com/python/reference/layout/shapes/>
        label_formatter : Callable[[str], str] | None, optional
            User-defined label text format function to change plot label text content.
            For example, if you want to change underscore of the label to space,
            set `lambda t: t.replace("_", " ")`.

        Returns
        -------
        circos : Circos
            Circos instance
        tv : TreeViz
            TreeViz instance
        """
        # Initialize circos sector with tree size
        tree = TreeViz.load_tree(tree_data, format=format)
        leaf_num = tree.count_terminals()
        circos = Circos(dict(tree=leaf_num), start=start, end=end)
        sector = circos.sectors[0]

        # Plot tree on track
        track = sector.add_track(r_lim)
        tv = track.tree(
            tree,
            format=format,
            outer=outer,
            align_leaf_label=align_leaf_label,
            ignore_branch_length=ignore_branch_length,
            leaf_label_size=leaf_label_size,
            leaf_label_rmargin=leaf_label_rmargin,
            reverse=reverse,
            ladderize=ladderize,
            line_kws=line_kws,
            label_formatter=label_formatter,
            align_line_kws=align_line_kws,
        )
        return circos, tv

    @staticmethod
    def initialize_from_bed(
        bed_file: str | Path,
        start: float = 0,
        end: float = 360,
        *,
        space: float | list[float] = 0,
        endspace: bool = True,
        sector2clockwise: dict[str, bool] | None = None,
    ) -> Circos:
        """Initialize Circos instance from BED file

        Circos sectors are auto-defined by BED chromosomes

        Parameters
        ----------
        bed_file : str | Path
            Chromosome BED format file (zero-based coordinate)
        start : float, optional
            Plot start degree (-360 <= start < end <= 360)
        end : float, optional
            Plot end degree (-360 <= start < end <= 360)
        space : float | list[float], optional
            Space degree(s) between sector
        endspace : bool, optional
            If True, insert space after the end sector
        sector2clockwise : dict[str, bool] | None, optional
            Sector name & clockwise bool dict. By default, `clockwise=True`.

        Returns
        -------
        circos : Circos
            Circos instance initialized from BED file

        """
        records = Bed(bed_file).records
        sectors = {rec.chr: (rec.start, rec.end) for rec in records}
        return Circos(
            sectors,
            start,
            end,
            space=space,
            endspace=endspace,
            sector2clockwise=sector2clockwise,
        )

    def add_cytoband_tracks(
        self,
        r_lim: tuple[float, float],
        cytoband_file: str | Path,
        *,
        track_name: str = "cytoband",
        cytoband_cmap: dict[str, str] | None = None,
        show_hovertext: bool = False,
    ) -> None:
        """Add track & plot chromosome cytoband on each sector

        Parameters
        ----------
        r_lim : tuple[float, float]
            Radius limit region (0 - 100)
        cytoband_file : str | Path
            Cytoband tsv file (UCSC format)
        track_name : str, optional
            Cytoband track name. By default, `cytoband`.
        cytoband_cmap : dict[str, str] | None, optional
            User-defined cytoband colormap. If None, use Circos style colormap.
            (e.g. `{"gpos100": "#000000", "gneg": "#FFFFFF", ...}`)
        show_hovertext : bool, optional
            If True, shows hovertext with band information. Default is False.
        """
        if cytoband_cmap is None:
            cytoband_cmap = config.CYTOBAND_COLORMAP
        cytoband_records = Bed(cytoband_file).records

        for sector in self.sectors:
            track = sector.add_track(r_lim, name=track_name)
            track.axis()

            # Prepare hover data if needed
            hover_x, hover_y, hover_texts, colors = [], [], [], []

            for rec in cytoband_records:
                if sector.name == rec.chr:
                    color = cytoband_cmap.get(str(rec.score), "white")
                    kwargs = utils.deep_dict_update(
                        config.cytoband_defaults, {"fillcolor": color}
                    )
                    track.rect(rec.start, rec.end, **kwargs)

                    if show_hovertext:
                        # Calculate midpoint for hover point
                        midpoint = (rec.start + rec.end) / 2
                        rad = track.x_to_rad(midpoint)
                        r = sum(r_lim) / 2
                        cx, cy = PolarSVGPatchBuilder._polar_to_cart(rad, r)
                        hover_x.append(cx)
                        hover_y.append(cy)
                        colors.append(color)
                        hover_texts.append(
                            f"Chromosome: {rec.chr}<br>"
                            f"Start: {rec.start:,}<br>"
                            f"End: {rec.end:,}<br>"
                            f"Band: {rec.name}<br>"
                            f"Type: {rec.score}"
                        )

            if show_hovertext and hover_x:
                hover_trace = utils.plot.build_scatter_trace(
                    hover_x,
                    hover_y,
                    mode="markers",
                    text=hover_texts,
                    marker=dict(
                        size=20,
                        opacity=0,
                    ),
                    hoverlabel={"bgcolor": colors},
                )
                track._traces.append(hover_trace)

    def get_sector(self, name: str) -> Sector:
        """Get sector by name

        Parameters
        ----------
        name : str
            Sector name

        Returns
        -------
        sector : Sector
            Sector

        """
        name2sector = {s.name: s for s in self.sectors}
        if name not in name2sector:
            raise ValueError(f"{name=} sector not found.")
        return name2sector[name]

    def get_group_sectors_deg_lim(
        self,
        group_sector_names: list[str],
    ) -> tuple[float, float]:
        """Get degree min-max limit in target group sectors

        Parameters
        ----------
        group_sector_names : list[str]
            Group sector names

        Returns
        -------
        group_sectors_deg_lim : tuple[float, float]
            Degree limit in group sectors

        """
        group_sectors = [self.get_sector(name) for name in group_sector_names]
        min_deg = min([min(s.deg_lim) for s in group_sectors])
        max_deg = max([max(s.deg_lim) for s in group_sectors])
        return min_deg, max_deg

    def axis(self, **kwargs) -> None:
        """Plot axis

        Parameters
        ----------
        **kwargs : dict, optional
            Shape properties
            (e.g. `fillcolor="red", line=dict(color="green", width=2, dash="dash")`)
            <https://plotly.com/python/reference/layout/shapes/>

        """
        kwargs = {} if kwargs is None else kwargs

        # Background shape placed behind other shapes (layer="below")
        fc_behind_kwargs = deepcopy(kwargs)
        fc_behind_kwargs = utils.deep_dict_update(
            fc_behind_kwargs, config.AXIS_FACE_PARAM
        )
        self.rect(**fc_behind_kwargs)

        # Edge shape placed in front of other shapes (layer="above")
        ec_front_kwargs = deepcopy(kwargs)
        ec_front_kwargs = utils.deep_dict_update(
            ec_front_kwargs, config.AXIS_EDGE_PARAM
        )
        self.rect(**ec_front_kwargs)

    def text(
        self,
        text: str,
        *,
        r: float = 0,
        deg: float = 0,
        adjust_rotation: bool = False,
        orientation: str = "horizontal",
        **kwargs,
    ) -> None:
        """Plot text on the entire circos plot. Uses angular positioning (0-360°).
        Angle is adjusted to Plotly's coordinate system:
            - 0° points upward (Plotly's default)
            - Angles increase clockwise

        Parameters
        ----------
        text : str
            Text content
        r : float, optional
            Radius position (default: 0, centered).
        deg : float, optional
            Degree position (0-360). 0° points upward.
        adjust_rotation : bool, optional
            If True, text rotation is auto set based on `deg` and `orientation`.
        orientation : str, optional
            Text orientation (`horizontal` or `vertical`).
        **kwargs : dict, optional
            Annotation properties (e.g. `font=dict(size=12, color='red')`).
            See: <https://plotly.com/python/reference/layout/annotations/>

        """
        rad = np.radians(deg)
        plotly_rad = -(rad - np.pi / 2)  # Convert to Plotly's polar coordinates
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
        deg_lim: tuple[float, float] | None = None,
        arc: bool = True,
        **kwargs,
    ) -> None:
        """Plot line

        Parameters
        ----------
        r : float | tuple[float, float]
            Line radius position (0 - 100). If r is float, (r, r) is set.
        deg_lim : tuple[float, float] | None, optional
            Degree limit region (-360 - 360). If None, `circos.deg_lim` is set.
        arc : bool, optional
            If True, plot arc style line for polar projection.
            If False, simply plot linear style line.
        **kwargs : dict, optional
            Line properties (e.g. `line=dict(color="red", width=2, dash="dash")`)
            See: <https://plotly.com/python/reference/layout/shapes/>

        """
        deg_lim = self.deg_lim if deg_lim is None else deg_lim
        start_deg, end_deg = min(deg_lim), max(deg_lim)
        rad_lim = (math.radians(start_deg), math.radians(end_deg))

        # Convert radius to tuple if needed
        r_lim = (r, r) if isinstance(r, (float, int)) else r

        path = (
            PolarSVGPatchBuilder.arc_line(rad_lim, r_lim)
            if arc
            else PolarSVGPatchBuilder.straight_line(rad_lim, r_lim)
        )

        shape = utils.plot.build_plotly_shape(
            path, config.plotly_shape_defaults, **kwargs
        )
        self._shapes.append(shape)

    def rect(
        self,
        r_lim: tuple[float, float] = (0, 100),
        deg_lim: tuple[float, float] | None = None,
        **kwargs,
    ) -> None:
        """Plot a rectangle spanning angular and radial ranges.
        Angle is adjusted to Plotly's coordinate system:
            - 0° points upward (Plotly's default)
            - Angles increase clockwise

        Parameters
        ----------
        r_lim : tuple[float, float]
            Radial limits (min, max) between 0 (center) and 100 (outer edge).
        deg_lim : tuple[float, float] | None, optional
            Angular limits in degrees (-360 to 360). If None, uses `circos.deg_lim`.
        **kwargs : dict, optional
            Shape properties (e.g. `fillcolor="red", line=dict(color="blue", width=2)`)
            See: <https://plotly.com/python/reference/layout/shapes/>

        """
        deg_lim = self.deg_lim if deg_lim is None else deg_lim
        rad_start = math.radians(deg_lim[0])
        rad_end = math.radians(deg_lim[1])

        min_rad, max_rad = min(rad_start, rad_end), max(rad_start, rad_end)

        # Build rectangle path
        radr = (min_rad, min(r_lim))
        width = max_rad - min_rad
        height = max(r_lim) - min(r_lim)

        path = PolarSVGPatchBuilder.arc_rectangle(radr, width, height)
        shape = utils.plot.build_plotly_shape(
            path, config.plotly_shape_defaults, **kwargs
        )
        self._shapes.append(shape)

    def link(
        self,
        sector_region1: tuple[str, float, float],
        sector_region2: tuple[str, float, float],
        r1: float | None = None,
        r2: float | None = None,
        *,
        height_ratio: float = 0.5,
        direction: int = 0,
        arrow_length_ratio: float = 0.05,
        allow_twist: bool = True,
        **kwargs,
    ) -> None:
        """Plot curved links between genomic regions using SVG paths.

        Parameters
        ----------
        sector_region1 : tuple[str, float, float]
            First region (sector_name, start, end)
        sector_region2 : tuple[str, float, float]
            Second region (sector_name, start, end)
        r1 : float | None, optional
            Radius for first region (None uses track bottom)
        r2 : float | None, optional
            Radius for second region (None uses track bottom)
        height_ratio : float, optional
            Controls curve height (default: 0.5)
        direction : int, optional
            0=no arrow, 1=forward, -1=reverse, 2=bidirectional
        arrow_length_ratio : float, optional
            Arrow size relative to link length
        allow_twist : bool, optional
            Whether to allow twisted ribbons
        **kwargs : dict, optional
            Shape properties (e.g. `fillcolor="red", line=dict(color="blue", width=2)`)
            Hover text for link (e.g. `hovertext="Link: ..."`).
            See: <https://plotly.com/python/reference/layout/shapes/>

        """
        # Extract regions
        name1, start1, end1 = sector_region1
        name2, start2, end2 = sector_region2

        # Get default hovertext or pop from kwargs
        arrow_symbol = LinkDirection(direction).arrow()
        hovertext = kwargs.pop(
            "hovertext",
            f"Link: {name1}:{start1}-{end1} {arrow_symbol} {name2}:{start2}-{end2}",
        )

        # Get sectors and calculate positions
        sector1, sector2 = self.get_sector(name1), self.get_sector(name2)
        r1 = sector1.get_lowest_r() if r1 is None else r1
        r2 = sector2.get_lowest_r() if r2 is None else r2
        rad_start1, rad_end1 = sector1.x_to_rad(start1), sector1.x_to_rad(end1)
        rad_start2, rad_end2 = sector2.x_to_rad(start2), sector2.x_to_rad(end2)

        # Handle twist resolution
        if not allow_twist:
            if (rad_end1 - rad_start1) * (rad_end2 - rad_start2) > 0:
                rad_start2, rad_end2 = rad_end2, rad_start2

        # Create Bezier curve path
        path = PolarSVGPatchBuilder.bezier_ribbon_path(
            rad_start1,
            rad_end1,
            r1,
            rad_start2,
            rad_end2,
            r2,
            height_ratio,
            direction,
            arrow_length_ratio,
        )

        shape = utils.plot.build_plotly_shape(
            path, defaults=config.plotly_arrow_defaults, **kwargs
        )
        self._shapes.append(shape)

        # Add invisible scatter points for hovertext at link positions
        hover_x, hover_y = zip(
            *[
                PolarSVGPatchBuilder._polar_to_cart((rad_start1 + rad_end1) / 2, r1),
                PolarSVGPatchBuilder._polar_to_cart((rad_start2 + rad_end2) / 2, r2),
            ]
        )
        hover_trace = utils.plot.build_scatter_trace(
            list(hover_x),
            list(hover_y),
            mode="markers",
            text=hovertext,
            marker=dict(size=20, opacity=0),
            hoverlabel={"bgcolor": shape["fillcolor"]},
        )
        self._traces.append(hover_trace)

    def link_line(
        self,
        sector_pos1: tuple[str, float],
        sector_pos2: tuple[str, float],
        r1: float | None = None,
        r2: float | None = None,
        *,
        height_ratio: float = 0.5,
        direction: int = 0,
        arrow_height: float = 3.0,
        arrow_width: float = 0.05,
        **kwargs,
    ) -> None:
        """Plot link line to specified position within or between sectors

        Parameters
        ----------
        sector_pos1 : tuple[str, float]
            Link line sector position1 (name, position)
        sector_pos2 : tuple[str, float]
            Link line sector position2 (name, position)
        r1 : float | None, optional
            Link line radius end position for sector_pos1.
            If None, lowest radius position of track in target sector is set.
        r2 : float | None, optional
            Link line radius end position for sector_pos2.
            If None, lowest radius position of track in target sector is set.
        height_ratio : float, optional
            Bezier curve height ratio
        direction : int, optional
            `0`: No direction edge shape (Default)
            `1`: Forward direction arrow edge shape (pos1 -> pos2)
            `-1`: Reverse direction arrow edge shape (pos1 <- pos2)
            `2`: Bidirectional arrow edge shape (pos1 <-> pos2)
        arrow_height : float, optional
            Arrow height size (Radius unit)
        arrow_width : float, optional
            Arrow width size (Degree unit)
        **kwargs : dict, optional
            Shape properties (e.g. `fillcolor="red", line=dict(color="blue", width=2)`)
            Hover text for link (e.g. `hovertext="Link: ..."`).
            See: <https://plotly.com/python/reference/layout/shapes/>

        """
        # Set data for plot link
        name1, pos1 = sector_pos1
        name2, pos2 = sector_pos2

        # Get default hovertext or pop from kwargs
        arrow_symbol = LinkDirection(direction).arrow()
        hovertext = kwargs.pop(
            "hovertext", f"Link: {name1}:{pos1} {arrow_symbol} {name2}:{pos2}"
        )

        # Get coordinates
        sector1, sector2 = self.get_sector(name1), self.get_sector(name2)
        r1 = sector1.get_lowest_r() if r1 is None else r1
        r2 = sector2.get_lowest_r() if r2 is None else r2
        rad_pos1, rad_pos2 = sector1.x_to_rad(pos1), sector2.x_to_rad(pos2)

        # Create Bezier curve path
        path = PolarSVGPatchBuilder.bezier_line_path(
            rad_pos1,
            r1,
            rad_pos2,
            r2,
            height_ratio,
            direction,
            arrow_height,
            arrow_width,
        )

        shape = utils.plot.build_plotly_shape(
            path, defaults=config.plotly_linelink_defaults, **kwargs
        )
        self._shapes.append(shape)

        # Add invisible scatter points for hovertext at link positions
        hover_x, hover_y = zip(
            *[
                PolarSVGPatchBuilder._polar_to_cart(rad_pos1, r1),
                PolarSVGPatchBuilder._polar_to_cart(rad_pos2, r2),
            ]
        )
        hover_trace = utils.plot.build_scatter_trace(
            list(hover_x),
            list(hover_y),
            mode="markers",
            text=hovertext,
            marker=dict(size=20, opacity=0),
            hoverlabel={"bgcolor": shape["line"]["color"]},
        )
        self._traces.append(hover_trace)

    def colorbar(
        self,
        *,
        vmin: int | float = 0,
        vmax: int | float = 1,
        cmap: str = "RdBu_r",
        **kwargs,
    ) -> str:
        """Plot colorbar using Plotly's coloraxis system.

        Parameters
        ----------
        vmin : int | float, optional
            Colorbar min value
        vmax : int | float, optional
            Colorbar max value
        cmap : str, optional
            Colormap name
        **kwargs : dict, optional
            Colorbar properties (e.g. `orientation="v", tickfont=dict(size=12)`)
            See: <https://plotly.com/python/reference/layout/coloraxis/>
        """
        # Create and store coloraxis config
        coloraxis_config = {
            "cmin": vmin,
            "cmax": vmax,
            "colorscale": cmap,
            "colorbar": kwargs,
        }

        coloraxis_name = (
            "coloraxis"
            if len(self._coloraxes) == 0
            else f"coloraxis{len(self._coloraxes) + 1}"
        )
        self._coloraxes.append(coloraxis_config)

        return coloraxis_name

    def plotfig(
        self,
        dpi: int = 100,
        figsize: tuple[float, float] = (8, 8),
        **kwargs,
    ) -> go.Figure:
        """Create the Plotly Circos-style figure.

        Parameters
        ----------
        dpi : int, optional
            Dots per inch (used to scale figsize)
        figsize : tuple[float, float], optional
            Size of figure in inches (width, height)
        **kwargs : dict
            Additional layout settings to override defaults

        Returns
        -------
        fig : go.Figure
            Plotly figure object

        """
        layout_dict = self._initialize_plotly_layout(figsize=figsize, dpi=dpi)
        layout_dict = utils.deep_dict_update(layout_dict, kwargs)

        # Plot trees (call to generate shapes, annotations and traces)
        for tv in self._get_all_treeviz_list():
            tv._plot_tree_line()
            tv._plot_tree_label()

        layout_dict["shapes"] = self._get_all_shapes()
        layout_dict["annotations"] = self._get_all_annotations()

        for i, coloraxis in enumerate(self._coloraxes):
            axis_key = "coloraxis" if i == 0 else f"coloraxis{i+1}"
            layout_dict[axis_key] = coloraxis

        data_dict = self._get_all_traces()

        return go.Figure(data=data_dict, layout=go.Layout(layout_dict))

    ############################################################
    # Private Method
    ############################################################

    def _check_degree_range(self, start: float, end: float) -> None:
        """Check start-end degree range (`-360 <= start < end <= 360`)

        Parameters
        ----------
        start : float
            Start degree range
        end : float
            End degree range

        """
        min_deg, max_deg = -360, 360
        if not min_deg <= start < end <= max_deg:
            err_msg = "start-end must be "
            err_msg += f"'{min_deg} <= start < end <= {max_deg}' ({start=}, {end=})"
            raise ValueError(err_msg)
        if end - start > max_deg:
            err_msg = f"'end - start' must be less than {max_deg} ({start=}, {end=})"
            raise ValueError(err_msg)

    def _to_sector2range(
        self,
        sectors: Mapping[str, int | float | tuple[float, float]],
    ) -> dict[str, tuple[float, float]]:
        """Convert sectors to sector2range"""
        sector2range: dict[str, tuple[float, float]] = {}
        for name, value in sectors.items():
            if isinstance(value, (tuple, list)):
                sector_start, sector_end = value
                if not sector_start < sector_end:
                    err_msg = f"{sector_end=} must be larger than {sector_start=}."
                    raise ValueError(err_msg)
                sector2range[name] = (sector_start, sector_end)
            else:
                sector2range[name] = (0, value)
        return sector2range

    @staticmethod
    def _initialize_plotly_layout(
        figsize: tuple[float, float] = (8, 8),
        dpi: int = 100,
    ) -> dict:
        """Initialize default Plotly layout based on config and figure size."""
        width = int(figsize[0] * dpi)
        height = int(figsize[1] * dpi)

        layout: dict = deepcopy(config.plotly_layout_defaults)

        layout["width"] = width
        layout["height"] = height
        layout["xaxis"]["range"] = config.AXIS_RANGE
        layout["yaxis"]["range"] = config.AXIS_RANGE

        return layout

    def _get_all_shapes(self) -> list[go.layout.Shape]:
        """Gather all shape dictionaries from self, sectors, and tracks."""
        circos_shapes = self._shapes
        sector_shapes = list(itertools.chain(*[s._shapes for s in self.sectors]))
        track_shapes = list(itertools.chain(*[t._shapes for t in self.tracks]))
        return circos_shapes + sector_shapes + track_shapes

    def _get_all_annotations(self) -> list[go.layout.Annotation]:
        """Gather all annotation dictionaries from self, sectors, and tracks."""
        circos_ann = self._annotations
        sector_ann = list(itertools.chain(*[s._annotations for s in self.sectors]))
        track_ann = list(itertools.chain(*[t._annotations for t in self.tracks]))
        return circos_ann + sector_ann + track_ann

    def _get_all_traces(self) -> list[BaseTraceType]:
        """Gather all traces from self, sectors, and tracks.

        Returns
        -------
        List[BaseTraceType]
            Combined list of all trace objects (scatter, bar, etc.)

        """
        # Get traces from main Circos object
        circos_traces = self._traces

        # Get traces from all sectors (flatten nested lists)
        sector_traces = list(
            itertools.chain(*[s._traces for s in self.sectors if hasattr(s, "_traces")])
        )

        # Get traces from all tracks (flatten nested lists)
        track_traces = list(
            itertools.chain(*[t._traces for t in self.tracks if hasattr(t, "_traces")])
        )

        return circos_traces + sector_traces + track_traces

    def _get_all_treeviz_list(self) -> list[TreeViz]:
        """Get all tree visualization instance list from tracks

        Returns
        -------
        all_treeviz_list : list[TreeViz]
            All tree visualization instance list
        """
        return list(itertools.chain(*[t._trees for t in self.tracks]))
