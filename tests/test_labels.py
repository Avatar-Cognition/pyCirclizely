import numpy as np
import pytest

from pycirclizely import Circos
from pycirclizely.utils.plot import GHOST_STYLE, pad_label_away_from_anchor


def _ann(circos):
    return list(circos.plotfig().layout.annotations)


@pytest.mark.parametrize("orientation", ["horizontal", "vertical"])
@pytest.mark.parametrize("outer", [True, False])
def test_xtick_label_anchored_at_tick_end(orientation, outer):
    """Label anchor sits exactly at tick end + margin, whatever the text length
    (the pixel extent is handled by a transparent copy of the text).
    """
    circos = Circos({"A": 10})
    track = circos.sectors[0].add_track((60, 80))
    labels = ["1", "a much longer label"]
    track.xticks(
        [3, 7],
        labels,
        outer=outer,
        tick_length=2,
        label_margin=1,
        label_orientation=orientation,
    )
    r_expected = 83 if outer else 57
    for ann, x in zip(_ann(circos), [3, 7]):
        assert np.hypot(ann.x, ann.y) == pytest.approx(r_expected)
        assert GHOST_STYLE in ann.text


def test_ghost_side_follows_reading_direction():
    """Ghost side follows reading direction."""
    # baseline pointing along the outward direction -> visible text second
    assert pad_label_away_from_anchor("ab", (1, 0), 0).endswith(">ab</span>ab")
    # baseline pointing against it -> visible text first
    assert pad_label_away_from_anchor("ab", (-1, 0), 0).startswith("ab<span")
    # baseline perpendicular, outward = text top -> visible line on top
    assert pad_label_away_from_anchor("ab", (0, 1), 0).startswith("ab<br>")
    assert pad_label_away_from_anchor("ab", (0, -1), 0).endswith("<br>ab")
