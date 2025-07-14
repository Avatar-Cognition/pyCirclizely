from dataclasses import dataclass
from enum import Enum
from typing import Optional, Union

from plotly.colors import qualitative, sequential  # type: ignore[attr-defined]
from webcolors import hex_to_rgb, name_to_rgb


@dataclass
class ColorPalette:
    """Container for color palette configuration."""

    class Type(Enum):
        """Type of ploty color palettes."""

        QUALITATIVE = "qualitative"
        SEQUENTIAL = "sequential"

    name: str
    colors: list[str]
    palette_type: Type
    counter: int = 0


class ColorCycler:
    """Color cycler that manages color selection from palettes.

    The cycler can work in two modes:
    1. Automatic cycling: When no index is provided, it remembers the last position
       in the palette and returns the next color each time (reset with reset_cycle())
    2. Indexed access: When an index is provided, it returns the color at that position
       (modulo the palette length) without affecting the cycle position.

    For sequential palettes, colors are distributed evenly across the available range.
    """

    def __init__(
        self,
        palette: Union[str, list[str]] = "Plotly",
        palette_type: Optional[ColorPalette.Type] = None,
    ):
        """
        Args:
            palette: Name of the palette to use or list of color strings.
            palette_type: Either QUALITATIVE or SEQUENTIAL. Required for custom
                         color lists, optional for built-in palettes.
        """
        self._palette = self._initialize_palette(palette, palette_type)

    def _initialize_palette(
        self,
        palette: Union[str, list[str]],
        palette_type: Optional[ColorPalette.Type] = None,
    ) -> ColorPalette:
        """Initialize a color palette from name or list of colors."""
        if isinstance(palette, str):
            # Handle built-in plotly palettes
            palette_type = palette_type or self._detect_palette_type(palette)

            if palette_type == ColorPalette.Type.QUALITATIVE:
                if not hasattr(qualitative, palette):
                    raise ValueError(f"Qualitative palette '{palette}' not found")
                colors = getattr(qualitative, palette)
            elif palette_type == ColorPalette.Type.SEQUENTIAL:
                if not hasattr(sequential, palette):
                    raise ValueError(f"Sequential palette '{palette}' not found")
                colors = getattr(sequential, palette)
        else:
            # Handle custom color lists
            if not isinstance(palette_type, ColorPalette.Type):
                raise ValueError(
                    "palette_type must be specified for custom color lists"
                )
            colors = palette

        return ColorPalette(
            name=palette if isinstance(palette, str) else "custom",
            colors=colors,
            palette_type=palette_type,
        )

    def _detect_palette_type(self, palette_name: str) -> ColorPalette.Type:
        """Detect if a palette is qualitative or sequential."""
        if hasattr(qualitative, palette_name):
            return ColorPalette.Type.QUALITATIVE
        elif hasattr(sequential, palette_name):
            return ColorPalette.Type.SEQUENTIAL
        raise ValueError(f"Could not detect palette type for '{palette_name}'")

    def get_color(self, index: Optional[int] = None) -> str:
        """Get a color from the palette.

        Args:
            index: Optional position in the palette. If None, returns the next color
                   in the automatic cycle sequence.

        """
        if index is None:
            index = self._palette.counter
            self._palette.counter += 1

        if self._palette.palette_type == ColorPalette.Type.SEQUENTIAL:
            # For sequential palettes, distribute colors evenly
            idx = int(
                (index % len(self._palette.colors))
                * (len(self._palette.colors) - 1)
                / max(1, len(self._palette.colors) - 1)
            )
            return self._palette.colors[idx]
        else:
            # For qualitative palettes, cycle through colors
            return self._palette.colors[index % len(self._palette.colors)]

    def get_colors(self, count: int) -> list[str]:
        """Get multiple colors from the palette.

        Args:
            count: Number of colors to return. Must be positive.

        """
        if count <= 0:
            raise ValueError(f"Count must be positive, got {count}")

        if self._palette.palette_type == ColorPalette.Type.SEQUENTIAL:
            if count == 1:
                return [self._palette.colors[len(self._palette.colors) // 2]]
            return [
                self._palette.colors[
                    int(i * (len(self._palette.colors) - 1) / (count - 1))
                ]
                for i in range(count)
            ]
        else:
            return [
                self._palette.colors[i % len(self._palette.colors)]
                for i in range(count)
            ]

    def reset_cycle(self) -> None:
        """Reset the automatic color cycle to the beginning of the palette."""
        self._palette.counter = 0

    def set_palette(
        self,
        palette: Union[str, list[str]],
        palette_type: Optional[ColorPalette.Type] = None,
    ) -> None:
        """Change the current color palette.

        Args:
            palette: Either a palette name or list of color strings.
            palette_type: Either QUALITATIVE or SEQUENTIAL. Required for custom
                         color lists, optional for built-in palettes.
        """
        self._palette = self._initialize_palette(palette, palette_type)
        self.reset_cycle()

    @property
    def palette_name(self) -> str:
        """Get the name of the current palette."""
        return self._palette.name

    @property
    def palette_type(self) -> ColorPalette.Type:
        """Get the type of the current palette."""
        return self._palette.palette_type


def parse_color(color: str) -> str:
    """
    Convert css, hex (including 8-digit with alpha),
    or named colors into an rgb/rgba string.

    Args:
        color: A string representing the color in CSS, hex, or named format.
    """
    if isinstance(color, str):
        color = color.strip().lower()

        # Handle rgb() or rgba() strings directly
        if color.startswith(("rgb(", "rgba(")):
            return color

        # Handle hex colors (including 8-digit with alpha)
        elif color.startswith("#"):
            hex_code = color.lstrip("#")

            # 8-digit hex (#RRGGBBAA)
            if len(hex_code) == 8:
                r = int(hex_code[0:2], 16)
                g = int(hex_code[2:4], 16)
                b = int(hex_code[4:6], 16)
                a = round(int(hex_code[6:8], 16) / 255, 2)  # Convert to 0-1 float
                return f"rgba({r}, {g}, {b}, {a})"

            # Standard 6-digit hex (#RRGGBB)
            elif len(hex_code) == 6:
                rgb = hex_to_rgb(color)
                return f"rgb({rgb.red}, {rgb.green}, {rgb.blue})"

            # 3-digit hex (#RGB)
            elif len(hex_code) == 3:
                # Expand to 6-digit and process
                expanded = f"#{hex_code[0] * 2}{hex_code[1] * 2}{hex_code[2] * 2}"
                rgb = hex_to_rgb(expanded)
                return f"rgb({rgb.red}, {rgb.green}, {rgb.blue})"

            else:
                raise ValueError(f"Invalid hex color code: {color}")

    # Handle named colors
    try:
        rgb = name_to_rgb(color)
        return f"rgb({rgb.red}, {rgb.green}, {rgb.blue})"
    except ValueError:
        raise ValueError(f"Could not parse color: {color}")
