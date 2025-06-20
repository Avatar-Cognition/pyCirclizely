from __future__ import annotations

import collections.abc
from typing import Any, Dict, Mapping, Optional, cast

from Bio.SeqFeature import SeqFeature
from plotly.colors import label_rgb, qualitative  # type: ignore[attr-defined]
from webcolors import hex_to_rgb, name_to_rgb


class ColorCycler:
    """Color cycler class using Plotly qualitative palettes."""

    _palette_counters: Dict[str, int] = {}
    _current_palette: str = "Plotly"
    _palette_colors: Dict[str, list[str]] = {"Plotly": qualitative.Plotly}

    def __new__(cls, n: Optional[int] = None) -> str:  # type: ignore[misc]
        """Return a color from the cycle when the class is called."""
        return cls.get_color(n)

    @classmethod
    def get_color(cls, n: Optional[int] = None) -> str:
        """Get a color from the current palette."""
        colors = cls._palette_colors[cls._current_palette]

        if n is None:
            if cls._current_palette not in cls._palette_counters:
                cls._palette_counters[cls._current_palette] = 0

            n = cls._palette_counters[cls._current_palette]
            cls._palette_counters[cls._current_palette] += 1

        return colors[n % len(colors)]

    @classmethod
    def get_color_list(cls, n: Optional[int] = None) -> list[str]:
        """Get a list of `n` colors from the current palette."""
        colors = cls._palette_colors[cls._current_palette]
        if n is None:
            return colors.copy()
        if n <= 0:
            raise ValueError(f"{n=} is invalid number (Must be 'n > 0').")
        return [colors[i % len(colors)] for i in range(n)]

    @classmethod
    def reset_cycle(cls, palette_name: Optional[str] = None) -> None:
        """Reset the color cycle counter."""
        if palette_name is None:
            cls._palette_counters = {}
        else:
            cls._palette_counters[palette_name] = 0

    @classmethod
    def set_palette(cls, name: str) -> None:
        """Set the current color palette by name."""
        if not hasattr(qualitative, name):
            raise ValueError(f"Palette '{name}' not found in plotly.colors.qualitative")
        if name not in cls._palette_colors:
            cls._palette_colors[name] = getattr(qualitative, name)
        cls._current_palette = name
        cls.reset_cycle(name)  # Reset counter when changing palette

    @classmethod
    def current_palette_name(cls) -> str:
        """Get the name of the current color palette."""
        return cls._current_palette


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
    """Convert css, hex or other color codes into a rgb coded string."""
    if isinstance(color, str) and color.strip().lower().startswith("rgb("):
        return color
    elif isinstance(color, str) and color.startswith("#"):
        rgb = hex_to_rgb(color)
    else:
        rgb = name_to_rgb(color)
    return label_rgb(rgb)
