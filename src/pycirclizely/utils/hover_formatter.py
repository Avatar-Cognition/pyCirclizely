import math
from dataclasses import dataclass, field
from typing import Any, Optional, Sequence, Union

import numpy as np
from plotly.basedatatypes import BaseTraceType

from pycirclizely import config
from pycirclizely.parser.table import StackedBarTable
from pycirclizely.patches import PolarSVGPatchBuilder
from pycirclizely.types import HoverText, NumericSequence

from .helper import deep_dict_update, precise_position
from .plot import build_scatter_trace


@dataclass
class HoverFormatter:
    hover_texts: list[str] = field(default_factory=list)
    hover_x: list[float] = field(default_factory=list)
    hover_y: list[float] = field(default_factory=list)
    hover_colors: list[str] = field(default_factory=list)
    precision_position: int = 2
    # Number of hover samples generated for each item (set by add_hover_regions).
    # When set, hover_texts / hover_colors (one per item) are repeated accordingly.
    region_counts: list[int] = field(default_factory=list)

    def process_hover_text(
        self,
        x: NumericSequence,
        y: NumericSequence,
        colors: Union[str, list[str]],
        x2: Optional[NumericSequence] = None,
        sector_name: Optional[str] = None,
        hover_text: HoverText = "default",
    ) -> None:
        """Process hover text for standard plots."""
        if hover_text == "default":
            self.hover_texts = self._default_hovertext(
                x, y, x2=x2, sector_name=sector_name
            )
        elif isinstance(hover_text, list):
            if len(hover_text) != len(x):
                raise ValueError(
                    f"hover_text length ({len(hover_text)}) must match data length"
                    f"{len(x)})"
                )
            self.hover_texts = hover_text
        else:
            raise TypeError("hover_text must be 'default', list[str], or None")

        # Handle colors - accept either single color or list
        if isinstance(colors, str):
            self.hover_colors = [colors] * len(x)
        else:
            if len(colors) != len(x):
                raise ValueError(
                    f"colors length ({len(colors)}) must match data length ({len(x)})"
                )
            self.hover_colors = colors

    def process_stacked_bar_hover_text(
        self,
        sb_table: StackedBarTable,
        col_name2color: dict[str, str],
        hover_text: HoverText = "default",
    ) -> None:
        """
        Process hover text for stacked bar plots.

        Args:
            hover_text: Either "default", a list of strings, or None
            sb_table: StackedBarTable object
            col_name2color: Dictionary mapping column names to colors
        """
        if hover_text is None:
            return

        if hover_text == "default":
            self.hover_texts = self._default_stackedbar_hovertext(sb_table)
        elif isinstance(hover_text, list):
            expected_length = len(sb_table.row_names) * len(sb_table.col_names)
            if len(hover_text) != expected_length:
                raise ValueError(
                    f"""hover_text length ({len(hover_text)}) must match number of
                     segments ({expected_length})"""
                )
            self.hover_texts = hover_text
        else:
            raise TypeError("hover_text must be 'default', list[str], or None")

        self.hover_colors = []
        for col_name in sb_table.col_names:
            color = col_name2color[col_name]
            self.hover_colors.extend([color] * len(sb_table.row_names))

    def add_hover_positions(
        self,
        rad: NumericSequence,
        r: NumericSequence,
    ) -> None:
        """
        Add hover positions from polar coordinates.

        Args:
            rad: Angular coordinates in radians
            r: Radial coordinates
        """
        for theta, rho in zip(rad, r):
            cx, cy = PolarSVGPatchBuilder._polar_to_cart(theta, rho)
            self.hover_x.append(cx)
            self.hover_y.append(cy)

    def add_hover_regions(
        self,
        rad_lims: Sequence[tuple[float, float]],
        r_lims: Sequence[tuple[float, float]],
        spacing: float | None = None,
    ) -> None:
        """
        Add hover positions covering polar rectangles (heatmap cells, bars...).

        Plotly only knows how to hover *points*, so each region is filled with a
        grid of invisible samples (one per sub-cell of about `spacing` data units).
        Plotly then shows the region whose sample is closest to the cursor, which
        is the region under the cursor (up to half a sub-cell at the borders),
        regardless of the region shape or how small/elongated it is.

        Args:
            rad_lims: (start, end) angle of each region in radians.
            r_lims: (bottom, top) radius of each region.
            spacing: Target distance between samples, in data (radius) units.
        """
        spacing = config.HOVER_SAMPLE_SPACING if spacing is None else spacing
        for (rad0, rad1), (r0, r1) in zip(rad_lims, r_lims):
            r_lo, r_hi = min(r0, r1), max(r0, r1)
            rad_lo, rad_hi = min(rad0, rad1), max(rad0, rad1)
            arc = (rad_hi - rad_lo) * (r_lo + r_hi) / 2
            n_r = max(1, math.ceil((r_hi - r_lo) / spacing))
            n_t = max(1, math.ceil(arc / spacing))
            n_t = min(n_t, config.HOVER_MAX_SAMPLES_PER_REGION // n_r or 1)
            rs = r_lo + (r_hi - r_lo) * (np.arange(n_r) + 0.5) / n_r
            rads = rad_lo + (rad_hi - rad_lo) * (np.arange(n_t) + 0.5) / n_t
            grid_rad, grid_r = np.meshgrid(rads, rs)
            self.add_hover_positions(grid_rad.ravel(), grid_r.ravel())
            self.region_counts.append(n_r * n_t)

    def create_hover_trace(self, **kwargs: Any) -> BaseTraceType:
        """Create a Plotly dummy hover trace."""
        texts, colors = self.hover_texts, self.hover_colors
        if self.region_counts:
            if len(texts) != len(self.region_counts):
                raise ValueError(
                    f"hover_text length ({len(texts)}) must match number of "
                    f"regions ({len(self.region_counts)})"
                )
            texts = list(np.repeat(np.asarray(texts, dtype=object), self.region_counts))
            if isinstance(colors, (list, tuple)) and len(colors) == len(
                self.region_counts
            ):
                colors = list(
                    np.repeat(np.asarray(colors, dtype=object), self.region_counts)
                )
        args = dict(
            text=texts,
            marker=dict(color=colors),
            hoverlabel=dict(bgcolor=colors),
        )

        # Start from defaults, then update with args, then with kwargs
        all_args = deep_dict_update(
            deep_dict_update(config.hovertext_dummy_defaults, args), kwargs
        )

        mode = all_args.pop("mode", "lines")
        return build_scatter_trace(
            x=self.hover_x, y=self.hover_y, mode=mode, **all_args
        )

    def clear(self) -> None:
        """Clear all stored hover data."""
        self.hover_texts.clear()
        self.hover_x.clear()
        self.hover_y.clear()
        self.hover_colors.clear()
        self.region_counts.clear()

    ############################################################
    # Private Method
    ############################################################

    def _default_hovertext(
        self,
        x: NumericSequence,
        y: NumericSequence,
        x2: Optional[NumericSequence] = None,
        sector_name: Optional[str] = None,
    ) -> list[str]:
        value_format = (
            f".{self.precision_position}f" if self.precision_position > 0 else ".0f"
        )

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
            xi = precise_position(xi, self.precision_position)
            if sector_name:
                parts.append(f"Sector: {sector_name}")
            if x2 is not None:
                xi2 = precise_position(x2[i], self.precision_position)
                parts.append(f"Position: {xi}–{xi2}")
            else:
                parts.append(f"Position: {xi}")
            parts.append(f"Value: {format(yi, value_format)}")
            hovertext.append("<br>".join(parts))
        return hovertext

    def _default_stackedbar_hovertext(
        self,
        sb_table: StackedBarTable,
    ) -> list[str]:
        """
        Generate default hover text for stacked bar plots.

        Args:
            sb_table: The stacked bar table object

        Returns
        -------
            List of formatted hover text strings
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
