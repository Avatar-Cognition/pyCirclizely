import pytest
from plotly.colors import qualitative  # type: ignore[attr-defined]

from pycirclizely.utils import ColorCycler, calc_group_spaces


def test_color_cycler():
    """Test color cycler"""
    # Check get colors length
    plotly_cycler = ColorCycler("Plotly")
    assert len(plotly_cycler.get_colors(len(qualitative.Plotly))) == len(
        qualitative.Plotly
    )
    assert len(plotly_cycler.get_colors(5)) == 5
    assert len(plotly_cycler.get_colors(20)) == 20

    # Check cycle index, color
    assert plotly_cycler.get_color(0) != plotly_cycler.get_color(1)
    assert plotly_cycler.get_color(0) == plotly_cycler.get_color(
        len(qualitative.Plotly)
    )
    assert plotly_cycler.get_color(15) == plotly_cycler.get_color(
        15 + len(qualitative.Plotly)
    )

    # Check cycle counter
    cycler1 = ColorCycler("Plotly")
    first_color = cycler1.get_color()
    second_color = cycler1.get_color()
    assert first_color != second_color
    assert cycler1._palette.counter == 2  # Check internal counter

    # Check reset cycle
    cycler1.reset_cycle()
    assert cycler1._palette.counter == 0

    # Check palette change
    alphabet_cycler = ColorCycler("Alphabet")
    assert len(alphabet_cycler.get_colors(len(qualitative.Alphabet))) == len(
        qualitative.Alphabet
    )

    # Test invalid palette
    with pytest.raises(ValueError):
        invalid_cycler = ColorCycler("invalid name")
        invalid_cycler.get_colors(5)

    # Test sequential palette
    viridis_cycler = ColorCycler("Viridis")
    colors = viridis_cycler.get_colors(5)
    assert len(colors) == 5
    assert colors[0] != colors[1]  # Should be different colors from sequential scale


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
