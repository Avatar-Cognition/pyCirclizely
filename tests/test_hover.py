import math

import numpy as np
import pytest

from pycirclizely import Circos, config


def _polar(x: float, y: float) -> tuple[float, float]:
    """Data cartesian -> (pycirclizely rad in [0, 2pi), r)."""
    return math.atan2(x, y) % (2 * math.pi), math.hypot(x, y)


def _hover_traces(fig):
    """Hover-only traces of ``fig`` (text hoverinfo, zero-width line)."""
    return [t for t in fig.data if t.hoverinfo == "text" and t.line.width == 0]


def _inside(rad, r, rad0, rad1, r0, r1, tol=1e-9):
    """Whether polar point ``(rad, r)`` lies in the annular sector given."""
    lo, hi = sorted((rad0 % (2 * math.pi), rad1 % (2 * math.pi)))
    return r0 - tol <= r <= r1 + tol and lo - tol <= rad <= hi + tol


def test_hover_trace_is_not_marker_based():
    """plotly.js treats every marker within max(3px, radius) of the cursor as
    equally close, so the *first* point wins instead of the closest one.
    Line-only traces use a 0.5px radius.
    """
    assert "markers" not in config.hovertext_dummy_defaults["mode"]


def test_heatmap_coloraxis_kept_with_line_hover_trace():
    """The colorbar still renders: plotly only draws a coloraxis referenced by
    a marker trace, and hover traces are now line-only.
    """
    circos = Circos({"A": 10})
    ca = circos.colorbar(vmin=0, vmax=1)
    track = circos.sectors[0].add_track((50, 90))
    track.heatmap(np.random.rand(2, 10), vmin=0, vmax=1, coloraxis=ca, hover_text=None)
    fig = circos.plotfig()
    assert any(
        "markers" in (t.mode or "") and t.marker.coloraxis == ca for t in fig.data
    )


def test_heatmap_hover_samples_cover_each_cell():
    """Every heatmap hover sample lies inside its own cell, and every cell gets
    several samples (one per cell let neighbours win on large cells).
    """
    circos = Circos({"A": 10}, start=0, end=90)
    track = circos.sectors[0].add_track((50, 90))
    data = np.arange(40).reshape(4, 10)
    ids = [f"{i}|{j}" for i in range(4) for j in range(10)]
    track.heatmap(data, hover_text=ids)
    (trace,) = _hover_traces(circos.plotfig())

    unit_r = track.r_plot_size / 4
    seen: dict[str, int] = {}
    for x, y, text in zip(trace.x, trace.y, trace.text):
        i, j = map(int, text.split("|"))
        rad, r = _polar(x, y)
        r1 = max(track.r_plot_lim) - unit_r * i
        assert _inside(
            rad, r, track.x_to_rad(j), track.x_to_rad(j + 1), r1 - unit_r, r1
        ), f"hover sample for cell {text} lies outside the cell"
        seen[text] = seen.get(text, 0) + 1
    assert set(seen) == set(ids)
    # cells are bigger than the sampling step -> several samples per cell
    assert min(seen.values()) > 1


@pytest.mark.parametrize("align", ["center", "edge"])
def test_bar_hover_samples_inside_bar(align):
    """Every bar hover sample lies inside its own bar, for both alignments
    (``align="edge"`` used to put the single sample on the start edge).
    """
    circos = Circos({"A": 10})
    track = circos.sectors[0].add_track((50, 90))
    x = np.arange(1, 9)
    h = np.arange(1, 9, dtype=float)
    track.bar(x, h, width=1.0, align=align, hover_text=[str(v) for v in x])
    (trace,) = _hover_traces(circos.plotfig())
    for px, py, text in zip(trace.x, trace.y, trace.text):
        k = int(text)
        x0 = k - 0.5 if align == "center" else k
        rad, r = _polar(px, py)
        assert _inside(
            rad,
            r,
            track.x_to_rad(x0),
            track.x_to_rad(x0 + 1),
            min(track.r_plot_lim),
            track._y_to_r(k, 0, 8),
        )


def test_bar_hover_text_none_disables_hover():
    """``bar(hover_text=None)`` emits no hover trace (it used to raise
    ``TypeError``).
    """
    circos = Circos({"A": 10})
    track = circos.sectors[0].add_track((50, 90))
    track.bar([1, 2], [1, 2], hover_text=None)
    assert _hover_traces(circos.plotfig()) == []


def test_layout_has_equal_aspect_ratio():
    """1 data unit spans the same pixels in x and y, so the circle is not
    stretched on non-square figures.
    """
    fig = Circos({"A": 10}).plotfig(figsize=(12, 6))
    assert fig.layout.yaxis.scaleanchor == "x"
    assert fig.layout.yaxis.scaleratio == 1
