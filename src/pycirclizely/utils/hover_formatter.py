from dataclasses import dataclass, field
from typing import Any, Optional, Union

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

    def create_hover_trace(self, **kwargs: Any) -> BaseTraceType:
        """Create a Plotly dummy hover trace."""
        args = dict(
            text=self.hover_texts,
            marker=dict(color=self.hover_colors),
            hoverlabel=dict(bgcolor=self.hover_colors),
        )

        # Start from defaults, then update with args, then with kwargs
        all_args = deep_dict_update(
            deep_dict_update(config.hovertext_dummy_defaults, args), kwargs
        )

        return build_scatter_trace(
            x=self.hover_x, y=self.hover_y, mode="markers", **all_args
        )

    def clear(self) -> None:
        """Clear all stored hover data."""
        self.hover_texts.clear()
        self.hover_x.clear()
        self.hover_y.clear()
        self.hover_colors.clear()

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
