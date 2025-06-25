from __future__ import annotations

import collections.abc
from typing import Any, Dict, List, Mapping, Optional, cast

from Bio.SeqFeature import SeqFeature
from plotly.colors import (  # type: ignore[attr-defined]
    qualitative,
    sequential,
)
from webcolors import hex_to_rgb, name_to_rgb


class ColorCycler:
    """Color cycler class supporting both qualitative and sequential color scales."""

    _palette_counters: Dict[str, int] = {}
    _current_palette: str = "Plotly"
    _palette_colors: Dict[str, List[str]] = {
        "Plotly": qualitative.Plotly,
        "default_qualitative": qualitative.Plotly,
        "default_sequential": sequential.Viridis,
    }
    _palette_type: Dict[str, str] = {
        "Plotly": "qualitative",
        "default_qualitative": "qualitative",
        "default_sequential": "sequential",
    }

    def __new__(cls, n: Optional[int] = None) -> str:  # type: ignore[misc]
        """Return a color from the cycle when the class is called."""
        return cls.get_color(n)

    @classmethod
    def get_color(cls, n: Optional[int] = None) -> str:
        """Get a color from the current palette."""
        colors = cls._palette_colors[cls._current_palette]

        if cls._palette_type[cls._current_palette] == "sequential":
            if n is None:
                if cls._current_palette not in cls._palette_counters:
                    cls._palette_counters[cls._current_palette] = 0
                n = cls._palette_counters[cls._current_palette]
                cls._palette_counters[cls._current_palette] += 1

            # Evenly distribute colors across sequential scale
            idx = int((n % len(colors)) * (len(colors) - 1) / max(1, len(colors) - 1))
            return colors[idx]
        else:
            # Original qualitative palette behavior
            if n is None:
                if cls._current_palette not in cls._palette_counters:
                    cls._palette_counters[cls._current_palette] = 0
                n = cls._palette_counters[cls._current_palette]
                cls._palette_counters[cls._current_palette] += 1
            return colors[n % len(colors)]

    @classmethod
    def get_color_list(cls, n: Optional[int] = None) -> List[str]:
        """Get a list of `n` colors from the current palette."""
        colors = cls._palette_colors[cls._current_palette]

        if n is None:
            return colors.copy()
        if n <= 0:
            raise ValueError(f"{n=} is invalid number (Must be 'n > 0').")

        if cls._palette_type[cls._current_palette] == "sequential":
            # For sequential palettes, evenly sample colors across the scale
            if n == 1:  # Return middle color for single value
                return [colors[len(colors) // 2]]
            return [colors[int(i * (len(colors) - 1) / (n - 1))] for i in range(n)]
        else:
            # For qualitative palettes, cycle through colors
            return [colors[i % len(colors)] for i in range(n)]

    @classmethod
    def get_sequential_color(
        cls, value: float, min_val: float = 0.0, max_val: float = 1.0
    ) -> str:
        """Get a color from a sequential palette based on a value range."""
        if cls._palette_type[cls._current_palette] != "sequential":
            raise ValueError("Current palette is not sequential.")
        colors = cls._palette_colors[cls._current_palette]
        normalized = (value - min_val) / (max_val - min_val)
        idx = int(normalized * (len(colors) - 1))
        return colors[max(0, min(idx, len(colors) - 1))]

    @classmethod
    def get_sequential_colors(
        cls,
        values: List[float],
        min_val: Optional[float] = None,
        max_val: Optional[float] = None,
    ) -> List[str]:
        """Get a list of colors from a sequential palette based on value ranges."""
        if cls._palette_type[cls._current_palette] != "sequential":
            raise ValueError("Current palette is not sequential.")
        if min_val is None:
            min_val = min(values)
        if max_val is None:
            max_val = max(values)
        return [cls.get_sequential_color(v, min_val, max_val) for v in values]

    @classmethod
    def reset_cycle(cls, palette_name: Optional[str] = None) -> None:
        """Reset the color cycle counter."""
        if palette_name is None:
            cls._palette_counters = {}
        else:
            cls._palette_counters[palette_name] = 0

    @classmethod
    def set_palette(cls, name: str, palette_type: Optional[str] = None) -> None:
        """Set the current color palette by name.

        Args:
            name: Name of the palette to set as current.
            palette_type: Type of palette ('qualitative' or 'sequential').
                If None, will try to detect automatically.
        """
        # Try to detect palette type if not specified
        if palette_type is None:
            if hasattr(qualitative, name):
                palette_type = "qualitative"
            elif hasattr(sequential, name):
                palette_type = "sequential"
            else:
                raise ValueError(f"Invalid built-in palette: '{name}'")

        # Set the palette based on detected/confirmed type
        if palette_type == "qualitative":
            if not hasattr(qualitative, name):
                raise ValueError(
                    f"Palette '{name}' not found in plotly.colors.qualitative"
                )
            if name not in cls._palette_colors:
                cls._palette_colors[name] = getattr(qualitative, name)
        elif palette_type == "sequential":
            if not hasattr(sequential, name):
                raise ValueError(
                    f"Palette '{name}' not found in plotly.colors.sequential"
                )
            if name not in cls._palette_colors:
                cls._palette_colors[name] = getattr(sequential, name)
        else:
            raise ValueError(
                f"Invalid palette_type: {palette_type}. "
                "Must be 'qualitative' or 'sequential'"
            )

        cls._current_palette = name
        cls._palette_type[name] = palette_type
        cls.reset_cycle(name)  # Reset counter when changing palette

    @classmethod
    def current_palette_name(cls) -> str:
        """Get the name of the current color palette."""
        return cls._current_palette

    @classmethod
    def current_palette_type(cls) -> str:
        """Get the type of the current color palette (qualitative or sequential)."""
        return cls._palette_type[cls._current_palette]


def deep_dict_update(
    orig_dict: Dict[str, Any], new_dict: Mapping[str, Any]
) -> Dict[str, Any]:
    """From deep-dict-update package https://pypi.org/project/deep-dict-update/
    Recursively updates a nested dictionary with the content of another dictionary.

    Parameters
    ----------
    - orig_dict (Dict[str, Any]): The original dictionary to be updated.
    - new_dict (Mapping[str, Any]): The dictionary containing updates.

    Returns
    -------
    - Dict[str, Any]: The updated dictionary.

    Notes
    -----
    - If a key in `new_dict` is not present in `orig_dict`,
      it will be added to `orig_dict`.
    - If a value in `new_dict` is a dictionary,
      the corresponding value in `orig_dict` will be updated recursively.
    - If a value in `new_dict` is a list of dictionaries,
      each dictionary in the list will be updated recursively.
    - For non-dictionary and non-list values,
      the value in `orig_dict` will be updated directly.

    """
    orig_dict = dict(orig_dict)
    for key, val in dict(new_dict).items():
        if key not in orig_dict:
            # If key is not present in orig_dict, initialize with an empty dictionary
            orig_dict[key] = {}

        if isinstance(val, collections.abc.Mapping):
            # If both orig_dict[key] and val are dictionaries, recursively update
            tmp = deep_dict_update(orig_dict[key], val)
            orig_dict[key] = tmp
        elif isinstance(val, list):
            # If the value is a list, iterate through the items
            # and apply dict_update for each dictionary in the list
            orig_dict[key] = [
                (
                    deep_dict_update(
                        orig_dict[key][i] if i < len(orig_dict[key]) else {},
                        cast(Dict[str, Any], item),  # Cast to Dict for type checker
                    )
                    if isinstance(item, collections.abc.Mapping)
                    else item
                )
                for i, item in enumerate(val)
            ]
        else:
            # For non-dictionary and non-list values, update directly
            orig_dict[key] = val

    return orig_dict


def calc_group_spaces(
    groups: list[int],
    *,
    space_bw_group: float = 15,
    space_in_group: float = 2,
    endspace: bool = True,
) -> list[float]:
    """Calculate spaces between/within groups

    This function can be used to easily calculate the space size
    when data is separated into multiple groups for display.
    For example, to set up a space to divide `[A, B, C, D, E, F, G, H, I, J]`
    into three groups such as `[(A, B, C, D), (E, F, G), (H, I, J)]`,
    set `groups=[4, 3, 3]`.

    Parameters
    ----------
    groups : list[int]
        List of each group number (e.g. `[4, 3, 3]`)
    space_bw_group : float, optional
        Space size between group
    space_in_group : float, optional
        Space size within group
    endspace : bool, optional
        If True, insert space after the end group

    Returns
    -------
    spaces : list[float]
        Spaces between/within groups

    """
    if len(groups) == 0:
        raise ValueError(f"{len(groups)=} is invalid.")
    elif len(groups) == 1:
        group_num = groups[0]
        return (
            [space_in_group] * (group_num - 1)
            if not endspace
            else [space_in_group] * group_num
        )
    else:
        spaces: list[float] = []
        for group_num in groups:
            group_spaces = [space_in_group] * (group_num - 1)
            group_spaces.append(space_bw_group)
            spaces.extend(group_spaces)
        return spaces[:-1] if not endspace else spaces


def is_pseudo_feature(feature: SeqFeature) -> bool:
    """Check target feature is pseudo or not from qualifiers tag

    Parameters
    ----------
    feature : SeqFeature
        Target feature

    Returns
    -------
    check_result : bool
        pseudo check result

    """
    quals = feature.qualifiers
    return True if "pseudo" in quals or "pseudogene" in quals else False


def precise_position(val: float, position_precision: int) -> float:
    """Round positions while preserving important decimals."""
    # First round to handle floating-point artifacts
    rounded = round(val, position_precision + 2)
    # Then round to target precision
    return round(rounded, position_precision)


def parse_color(color):
    """
    Convert css, hex (including 8-digit with alpha),
    or named colors into an rgb/rgba string.
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
