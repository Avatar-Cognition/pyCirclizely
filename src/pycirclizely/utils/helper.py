import collections.abc
from typing import Any, Mapping, cast

from Bio.SeqFeature import SeqFeature


def deep_dict_update(
    orig_dict: dict[str, Any], new_dict: Mapping[str, Any]
) -> dict[str, Any]:
    """Recursively updates a nested dictionary with the content of another dictionary.

    Args:
        orig_dict: The original dictionary to be updated.
        new_dict: The dictionary containing updates.
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
                        cast(dict[str, Any], item),  # Cast to dict for type checker
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
    """Calculate spaces between/within groups.

    Args:
        groups: List of each group number (e.g., `[4, 3, 3]`).
        space_bw_group: Space size between groups.
        space_in_group: Space size within groups.
        endspace: If True, insert space after the end group.
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
    """Check target feature is pseudo or not from qualifiers tag.

    Args:
        feature: Target feature.
    """
    quals = feature.qualifiers
    return True if "pseudo" in quals or "pseudogene" in quals else False


def precise_position(val: float, position_precision: int) -> float:
    """Round positions while preserving important decimals.

    Args:
        val: Value to round.
        position_precision: Precision for rounding.
    """
    # First round to handle floating-point artifacts
    rounded = round(val, position_precision + 2)
    # Then round to target precision
    return round(rounded, position_precision)
