from __future__ import annotations

from pathlib import Path
from urllib.parse import urlparse
from urllib.request import urlopen

from Bio.SeqFeature import SeqFeature
from PIL import Image

from plotly.colors import qualitative
from typing import Optional, Dict, Any
import collections.abc

class ColorCycler:
    """Color cycler class using Plotly qualitative palettes."""

    _counter: int = 0
    _palette_name: str = "Plotly"
    _colors: list[str] = qualitative.Plotly  # Default palette

    def __new__(cls, n: Optional[int] = None) -> str:
        """Return a color from the cycle, same as get_color."""
        return cls.get_color(n)

    @classmethod
    def reset_cycle(cls) -> None:
        """Reset the color cycle counter."""
        cls._counter = 0

    @classmethod
    def set_palette(cls, name: str) -> None:
        """Set the color palette by name (e.g. 'Plotly', 'D3', 'Dark24')."""
        if not hasattr(qualitative, name):
            raise ValueError(f"Palette '{name}' not found in plotly.colors.qualitative")
        cls._colors = getattr(qualitative, name)
        cls._palette_name = name
        cls._counter = 0

    @classmethod
    def get_color(cls, n: Optional[int] = None) -> str:
        """Get a color from the palette, either by index or the next in the cycle."""
        if n is None:
            n = cls._counter
            cls._counter += 1
        return cls._colors[n % len(cls._colors)]

    @classmethod
    def get_color_list(cls, n: Optional[int] = None) -> list[str]:
        """Get a list of `n` colors from the palette (cycled if n > palette size)."""
        if n is None:
            return cls._colors.copy()
        if n <= 0:
            raise ValueError(f"{n=} is invalid number (Must be 'n > 0').")
        return [cls._colors[i % len(cls._colors)] for i in range(n)]

    @classmethod
    def current_palette_name(cls) -> str:
        """Get the name of the current color palette."""
        return cls._palette_name

def deep_dict_update(orig_dict: Dict[str, Any], new_dict: Dict[str, Any]) -> Dict[str, Any]:
    """ From deep-dict-update package https://pypi.org/project/deep-dict-update/
    Recursively updates a nested dictionary with the content of another dictionary.

    Parameters:
    - orig_dict (Dict[str, Any]): The original dictionary to be updated.
    - new_dict (Dict[str, Any]): The dictionary containing updates.

    Returns:
    - Dict[str, Any]: The updated dictionary.

    Notes:
    - If a key in `new_dict` is not present in `orig_dict`, it will be added to `orig_dict`.
    - If a value in `new_dict` is a dictionary, the corresponding value in `orig_dict` will be updated recursively.
    - If a value in `new_dict` is a list of dictionaries, each dictionary in the list will be updated recursively.
    - For non-dictionary and non-list values, the value in `orig_dict` will be updated directly.
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
                deep_dict_update(orig_dict[key][i] if i < len(orig_dict[key]) else {}, item) if isinstance(item, collections.abc.Mapping) else item
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
        spaces = [space_in_group] * group_num
    else:
        spaces: list[float] = []
        for group_num in groups:
            group_spaces = [space_in_group] * (group_num - 1)
            group_spaces.extend([space_bw_group])
            spaces.extend(group_spaces)
    if endspace:
        return spaces
    else:
        return spaces[:-1]


def load_image(img: str | Path | Image.Image) -> Image.Image:
    """Load target image as PIL Image

    Parameters
    ----------
    img : str | Path | Image.Image
        Load target image (`File Path`|`URL`|`PIL Image`)

    Returns
    -------
    im : Image.Image
        PIL Image (mode=`RGBA`)
    """
    if isinstance(img, str) and urlparse(img).scheme in ("http", "https"):
        im = Image.open(urlopen(img))
    elif isinstance(img, (str, Path)):
        im = Image.open(str(img))
    else:
        im = img
    return im.convert("RGBA")


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
