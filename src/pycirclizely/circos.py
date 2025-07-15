from __future__ import annotations

import itertools
import math
import textwrap
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any, Callable, Optional, Union

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from Bio.Phylo.BaseTree import Tree
from plotly.basedatatypes import BaseTraceType

from pycirclizely import config, utils
from pycirclizely.impl.chord import chord_diagram_impl
from pycirclizely.impl.radar import radar_chart_impl
from pycirclizely.parser import Bed
from pycirclizely.parser.matrix import Matrix
from pycirclizely.parser.table import RadarTable
from pycirclizely.patches import PolarSVGPatchBuilder
from pycirclizely.sector import Sector
from pycirclizely.track import Track
from pycirclizely.tree import TreeViz
from pycirclizely.types import HoverText, LabelFormatter, Numeric, TextFormatter
from pycirclizely.utils.hover_formatter import HoverFormatter


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
        """
        Args:
            sectors: Sector name & size (or range) dict.
            start: Plot start degree (`-360 <= start < end <= 360`).
            end: Plot end degree (`-360 <= start < end <= 360`).
            space: Space degree(s) between sectors.
            endspace: If True, insert space after the end sector.
            sector2clockwise: Sector name & clockwise bool dict.
                By default, `clockwise=True`.
            show_axis_for_debug: Show axis for position check debugging (Dev option).
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
        """
        Args:
            table: Table file or Table dataframe or RadarTable instance.
            r_lim: Radar chart radius limit region (0 - 100).
            vmin: Min value.
            vmax: Max value.
            fill: If True, fill color of radar chart.
            marker_size: Marker size.
            bg_color: Background color.
            circular: If True, plot with circular style.
            cmap: Colormap assigned to each target row(index) in table.
            show_grid_label: If True, show grid label.
            grid_interval_ratio: Grid interval ratio (0.0 - 1.0).
            grid_line_kws: Keyword arguments passed to `track.line()` method.
            grid_label_kws: Keyword arguments passed to `track.text()` method.
            grid_label_formatter: User-defined function to format grid label.
            label_kws_handler: Handler function for keyword arguments passed
                to `track.text()` method.
            line_kws_handler: Handler function for keyword arguments passed
                to `track.line()` method.
            marker_kws_handler: Handler function for keyword arguments passed
                to `track.scatter()` method.
        """
        return radar_chart_impl(
            table=table,
            r_lim=r_lim,
            vmin=vmin,
            vmax=vmax,
            fill=fill,
            marker_size=marker_size,
            bg_color=bg_color,
            circular=circular,
            cmap=cmap,
            show_grid_label=show_grid_label,
            grid_interval_ratio=grid_interval_ratio,
            grid_line_kws=grid_line_kws,
            grid_label_kws=grid_label_kws,
            grid_label_formatter=grid_label_formatter,
            label_kws_handler=label_kws_handler,
            line_kws_handler=line_kws_handler,
            marker_kws_handler=marker_kws_handler,
        )

    @staticmethod
    def chord_diagram(
        matrix: Union[str, Path, pd.DataFrame, Matrix],
        *,
        start: float = 0,
        end: float = 360,
        space: Union[float, list[float]] = 0,
        endspace: bool = True,
        r_lim: tuple[float, float] = (97, 100),
        cmap: Union[str, dict[str, str]] = "Viridis",
        link_cmap: Optional[list[tuple[str, str, str]]] = None,
        ticks_interval: Optional[int] = None,
        order: Union[str, list[str], None] = None,
        label_kws: Optional[dict[str, Any]] = None,
        ticks_kws: Optional[dict[str, Any]] = None,
        link_kws: Optional[dict[str, Any]] = None,
        link_kws_handler: Optional[
            Callable[[str, str], Optional[dict[str, Any]]]
        ] = None,
    ) -> Circos:
        """Plot chord diagram from matrix data.

        Args:
            matrix: Matrix file or Matrix dataframe or Matrix instance.
            start: Plot start degree (-360 <= start < end <= 360).
            end: Plot end degree (-360 <= start < end <= 360).
            space: Space degree(s) between sectors.
            endspace: If True, insert space after the end sector.
            r_lim: Outer track radius limit region (0 - 100).
            cmap: Colormap assigned to each outer track and link.
            link_cmap: Link colormap to overwrite link colors automatically set by cmap.
            ticks_interval: Ticks interval. If None, ticks are not plotted.
            order: Sort order of matrix for plotting Chord Diagram.
            label_kws: Keyword arguments passed to `sector.text()` method.
            ticks_kws: Keyword arguments passed to `track.xticks_by_interval()` method.
            link_kws: Keyword arguments passed to `circos.link()` method.
            link_kws_handler: User-defined function to handle keyword arguments
                for each link.
        """
        return chord_diagram_impl(
            matrix=matrix,
            start=start,
            end=end,
            space=space,
            endspace=endspace,
            r_lim=r_lim,
            cmap=cmap,
            link_cmap=link_cmap,
            ticks_interval=ticks_interval,
            order=order,
            label_kws=label_kws,
            ticks_kws=ticks_kws,
            link_kws=link_kws,
            link_kws_handler=link_kws_handler,
        )

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
        label_formatter: TextFormatter = None,
        align_line_kws: dict[str, Any] | None = None,
    ) -> tuple[Circos, TreeViz]:
        """Circos sector and track are auto-defined by phylogenetic tree.

        Args:
            tree_data: Tree data (file, file URL, tree object, or tree string).
            start: Plot start degree (-360 <= start < end <= 360).
            end: Plot end degree (-360 <= start < end <= 360).
            r_lim: Tree track radius limit region (0 - 100).
            format: Tree format (newick, phyloxml, nexus, nexml, or cdao).
            outer: If True, plot tree on outer side. If False, plot tree on inner side.
            align_leaf_label: If True, align leaf label.
            ignore_branch_length: If True, ignore branch length for plotting tree.
            leaf_label_size: Leaf label size.
            leaf_label_rmargin: Leaf label radius margin.
            reverse: If True, reverse tree.
            ladderize: If True, ladderize tree.
            line_kws: Shape properties
                (e.g., dict(line=dict(color="red", width=1, dash="dash"))).
                <https://plotly.com/python/reference/layout/shapes/>
            align_line_kws: Shape properties
                (e.g., dict(line=dict(color="black", dash="dot"), opacity=0.5)).
                <https://plotly.com/python/reference/layout/shapes/>
            label_formatter: User-defined label text format function
                to change plot label text content. For example, if you want to change
                underscores in the label to spaces, set `lambda t: t.replace("_", " ")`.
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
        """Circos sectors are auto-defined by BED chromosomes.

        Args:
            bed_file: Chromosome BED format file (zero-based coordinate).
            start: Plot start degree (-360 <= start < end <= 360).
            end: Plot end degree (-360 <= start < end <= 360).
            space: Space degree(s) between sectors.
            endspace: If True, insert space after the end sector.
            sector2clockwise: Sector name & clockwise bool dict.

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
        """
        Args:
            r_lim: Radius limit region (0 - 100).
            cytoband_file: Cytoband tsv file (UCSC format).
            track_name: Cytoband track name. By default, `cytoband`.
            cytoband_cmap: User-defined cytoband colormap.
                If None, use Circos style colormap.
            show_hovertext: If True, shows hovertext with band information.
        """
        if cytoband_cmap is None:
            cytoband_cmap = config.CYTOBAND_COLORMAP
        cytoband_records = Bed(cytoband_file).records

        for sector in self.sectors:
            track = sector.add_track(r_lim, name=track_name)
            track.axis()

            if show_hovertext:
                hover_formatter = HoverFormatter(
                    precision_position=0,
                )

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

                        # Add hover position and text
                        hover_formatter.add_hover_positions([rad], [r])
                        hover_formatter.hover_colors.append(color)
                        hover_formatter.hover_texts.append(
                            f"Chromosome: {rec.chr}<br>"
                            f"Start: {rec.start:,}<br>"
                            f"End: {rec.end:,}<br>"
                            f"Band: {rec.name}<br>"
                            f"Type: {rec.score}"
                        )
            if show_hovertext:
                if hover_trace := hover_formatter.create_hover_trace():
                    track._traces.append(hover_trace)

    def get_sector(self, name: str) -> Sector:
        """Get sector by name"""
        name2sector = {s.name: s for s in self.sectors}
        if name not in name2sector:
            raise ValueError(f"{name=} sector not found.")
        return name2sector[name]

    def get_group_sectors_deg_lim(
        self,
        group_sector_names: list[str],
    ) -> tuple[float, float]:
        """Get degree min-max limit in target group sectors."""
        group_sectors = [self.get_sector(name) for name in group_sector_names]
        min_deg = min([min(s.deg_lim) for s in group_sectors])
        max_deg = max([max(s.deg_lim) for s in group_sectors])
        return min_deg, max_deg

    def axis(self, **kwargs: Any) -> None:
        """
        Args:
            **kwargs: Shape properties
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
        **kwargs: Any,
    ) -> None:
        """Plot text using angular positioning (0-360°).
        Angle is adjusted to Plotly's coordinate system:
            - 0° points upward (Plotly's default)
            - Angles increase clockwise

        Args:
            text: Text content.
            r: Radius position (default: 0, centered).
            deg: Degree position (0-360). 0° points upward.
            adjust_rotation: If True, text rotation is auto set based on
                `deg` and `orientation`.
            orientation: Text orientation (`horizontal` or `vertical`).
            **kwargs: Annotation properties (e.g. `font=dict(size=12, color='red')`).
                <https://plotly.com/python/reference/layout/annotations/>
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
        **kwargs: Any,
    ) -> None:
        """
        Args:
            r: Line radius position (0 - 100). If r is float, (r, r) is set.
            deg_lim: Degree limit region (-360 - 360).
                If None, `circos.deg_lim` is set.
            arc: If True, plot arc style line for polar projection.
                If False, simply plot linear style line.
            **kwargs: Line properties
                (e.g. `line=dict(color="red", width=2, dash="dash")`)
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
        **kwargs: Any,
    ) -> None:
        """Plot a rectangle spanning angular and radial ranges.
        Angle is adjusted to Plotly's coordinate system:
            - 0° points upward (Plotly's default)
            - Angles increase clockwise

        Args:
            r_lim: Radial limits (min, max) between 0 (center) and 100 (outer edge).
            deg_lim: Angular limits in degrees (-360 to 360).
                If None, uses `circos.deg_lim`.
            **kwargs: Shape properties
                (e.g. `fillcolor="red", line=dict(color="blue", width=2)`)
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
        hover_text: HoverText = "default",
        **kwargs: Any,
    ) -> None:
        """Plot curved link between two regions within or between sectors.

        Args:
            sector_region1: First region (sector_name, start, end).
            sector_region2: Second region (sector_name, start, end).
            r1: Radius for first region (None uses track bottom).
            r2: Radius for second region (None uses track bottom).
            height_ratio: Controls curve height (default: 0.5).
            direction: 0=no arrow, 1=forward, -1=reverse, 2=bidirectional.
            arrow_length_ratio: Arrow size relative to link length.
            allow_twist: Whether to allow twisted ribbons.
            hover_text: Custom hover text or default.
            **kwargs: Shape properties
                (e.g. `fillcolor="red", line=dict(color="blue", width=2)`).
                <https://plotly.com/python/reference/layout/shapes/>
        """
        name1, start1, end1 = sector_region1
        name2, start2, end2 = sector_region2

        sector1, sector2 = self.get_sector(name1), self.get_sector(name2)
        r1 = sector1.get_lowest_r() if r1 is None else r1
        r2 = sector2.get_lowest_r() if r2 is None else r2
        rad_start1, rad_end1 = sector1.x_to_rad(start1), sector1.x_to_rad(end1)
        rad_start2, rad_end2 = sector2.x_to_rad(start2), sector2.x_to_rad(end2)

        if not allow_twist:
            if (rad_end1 - rad_start1) * (rad_end2 - rad_start2) > 0:
                rad_start2, rad_end2 = rad_end2, rad_start2

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

        if hover_text is None:
            return

        hover_formatter = HoverFormatter()

        # Add hover positions at both ends of the link
        hover_rad = [(rad_start1 + rad_end1) / 2, (rad_start2 + rad_end2) / 2]
        hover_r = [r1, r2]
        hover_formatter.add_hover_positions(hover_rad, hover_r)
        color = shape.get("fillcolor") or "grey"
        hover_formatter.hover_colors = [color, color]
        arrow_symbol = utils.plot.LinkDirection(direction).arrow()
        hover_formatter.hover_texts = [
            f"Link: {name1}:{start1}-{end1} {arrow_symbol} {name2}:{start2}-{end2}"
        ] * 2

        # Create and add hover trace
        if hover_trace := hover_formatter.create_hover_trace():
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
        hover_text: HoverText = "default",
        **kwargs: Any,
    ) -> None:
        """Plot link line to specified position within or between sectors

        Args:
            sector_pos1: Link line sector position1 (name, position).
            sector_pos2: Link line sector position2 (name, position).
            r1: Link line radius end position for sector_pos1.
                If None, lowest radius position of track in target sector is set.
            r2: Link line radius end position for sector_pos2.
                If None, lowest radius position of track in target sector is set.
            height_ratio: Bezier curve height ratio.
            direction: Direction of arrow edge shape.
                0: No direction edge shape (Default).
                1: Forward direction arrow edge shape (pos1 -> pos2).
                -1: Reverse direction arrow edge shape (pos1 <- pos2).
                2: Bidirectional arrow edge shape (pos1 <-> pos2).
            arrow_height: Arrow height size (Radius unit).
            arrow_width: Arrow width size (Degree unit).
            hover_text: Custom hover text or default.
            **kwargs: Shape properties
                (e.g. `fillcolor="red", line=dict(color="blue", width=2)`).
                See: <https://plotly.com/python/reference/layout/shapes/>
        """
        name1, pos1 = sector_pos1
        name2, pos2 = sector_pos2

        sector1, sector2 = self.get_sector(name1), self.get_sector(name2)
        r1 = sector1.get_lowest_r() if r1 is None else r1
        r2 = sector2.get_lowest_r() if r2 is None else r2
        rad_pos1, rad_pos2 = sector1.x_to_rad(pos1), sector2.x_to_rad(pos2)

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

        if hover_text is None:
            return

        hover_formatter = HoverFormatter()

        # Add hover positions at both ends of the link
        hover_rad = [rad_pos1, rad_pos2]
        hover_r = [r1, r2]
        hover_formatter.add_hover_positions(hover_rad, hover_r)
        color = shape.get("line", {}).get("color") or "grey"
        hover_formatter.hover_colors = [color, color]
        arrow_symbol = utils.plot.LinkDirection(direction).arrow()
        hover_formatter.hover_texts = [
            f"Link: {name1}:{pos1} {arrow_symbol} {name2}:{pos2}"
        ] * 2

        # Create and add hover trace
        if hover_trace := hover_formatter.create_hover_trace():
            self._traces.append(hover_trace)

    def colorbar(
        self,
        *,
        vmin: Numeric = 0,
        vmax: Numeric = 1,
        cmap: str = "RdBu_r",
        **kwargs: Any,
    ) -> str:
        """Plot colorbar using Plotly's coloraxis system.

        Args:
            vmin: Colorbar min value.
            vmax: Colorbar max value.
            cmap: Colormap name.
            **kwargs: Colorbar properties
                (e.g., `orientation="v", tickfont=dict(size=12)`).
                See: <https://plotly.com/python/reference/layout/coloraxis/>.
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
        figsize: tuple[float, float] = (7, 7),
        **kwargs: Any,
    ) -> go.Figure:
        """Create Plotly figure with all shapes, annotations, and traces.

        Args:
            dpi: Dots per inch (used to scale figsize).
            figsize: Size of figure in inches (width, height).
            **kwargs: Additional layout settings to override defaults.
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
        """Check start-end degree range (`-360 <= start < end <= 360`)."""
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
        figsize: tuple[float, float] = (7, 7),
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
        sector_shapes = list(
            itertools.chain.from_iterable(s._shapes for s in self.sectors)
        )
        track_shapes = list(
            itertools.chain.from_iterable(t._shapes for t in self.tracks)
        )
        return circos_shapes + sector_shapes + track_shapes

    def _get_all_annotations(self) -> list[go.layout.Annotation]:
        """Gather all annotation dictionaries from self, sectors, and tracks."""
        circos_ann = self._annotations
        sector_ann = list(
            itertools.chain.from_iterable(s._annotations for s in self.sectors)
        )
        track_ann = list(
            itertools.chain.from_iterable(t._annotations for t in self.tracks)
        )
        return circos_ann + sector_ann + track_ann

    def _get_all_traces(self) -> list[BaseTraceType]:
        """Gather all traces from self, sectors, and tracks."""
        # Get traces from main Circos object
        circos_traces = self._traces

        # Get traces from all sectors (flatten nested lists)
        sector_traces = list(
            itertools.chain.from_iterable(
                s._traces for s in self.sectors if hasattr(s, "_traces")
            )
        )

        # Get traces from all tracks (flatten nested lists)
        track_traces = list(
            itertools.chain.from_iterable(
                t._traces for t in self.tracks if hasattr(t, "_traces")
            )
        )

        return circos_traces + sector_traces + track_traces

    def _get_all_treeviz_list(self) -> list[TreeViz]:
        """Get all tree visualization instance list from tracks"""
        return list(itertools.chain.from_iterable(t._trees for t in self.tracks))
