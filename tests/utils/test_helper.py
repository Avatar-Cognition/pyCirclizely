import pytest
from plotly.colors import qualitative  # type: ignore[attr-defined]

from pycirclizely.utils import ColorCycler, calc_group_spaces


def test_color_cycler():
    """Test color cycler"""
    # Check get color list length
    ColorCycler.set_palette("Plotly")
    assert len(ColorCycler.get_color_list()) == len(qualitative.Plotly)
    assert len(ColorCycler.get_color_list(5)) == 5
    assert len(ColorCycler.get_color_list(20)) == 20

    # Check cycle index, color
    assert ColorCycler.get_color(0) != ColorCycler.get_color(1)
    assert ColorCycler.get_color(0) == ColorCycler.get_color(len(qualitative.Plotly))
    assert ColorCycler.get_color(15) == ColorCycler.get_color(
        15 + len(qualitative.Plotly)
    )

    # Check cycle counter
    assert ColorCycler.get_color() != ColorCycler.get_color()
    assert ColorCycler._palette_counters.get(ColorCycler._current_palette, 0) == 2

    # Check reset cycle
    ColorCycler.reset_cycle()
    assert ColorCycler._palette_counters.get(ColorCycler._current_palette, 0) == 0

    # Check palette change
    ColorCycler.set_palette("Alphabet")
    with pytest.raises(ValueError):
        ColorCycler.set_palette("invalid name")
    assert len(ColorCycler.get_color_list()) == len(qualitative.Alphabet)


def test_calc_group_spaces():
    """Test `calc_group_spaces`"""
    # Case1. Blank list (error)
    with pytest.raises(ValueError):
        calc_group_spaces([])

    # Case2. List length = 1 (endspace=True)
    spaces = calc_group_spaces([5])
    expected_spaces = [2, 2, 2, 2, 2]
    assert spaces == expected_spaces

    # Case3. List length = 1 (endspace=False)
    spaces = calc_group_spaces([5], space_in_group=3, endspace=False)
    expected_spaces = [3, 3, 3, 3]
    assert spaces == expected_spaces

    # Case4. List length > 1 (endspace=True)
    spaces = calc_group_spaces([4, 3, 3])
    expected_spaces = [2, 2, 2, 15, 2, 2, 15, 2, 2, 15]
    assert spaces == expected_spaces

    # Case5. List length > 1 (endspace=False)
    spaces = calc_group_spaces(
        [4, 3, 3], space_bw_group=8, space_in_group=1, endspace=False
    )
    expected_spaces = [1, 1, 1, 8, 1, 1, 8, 1, 1]
    assert spaces == expected_spaces
