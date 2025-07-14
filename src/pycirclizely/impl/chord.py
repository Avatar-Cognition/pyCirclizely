from __future__ import annotations

from collections import defaultdict
from copy import deepcopy
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Optional, Union

import pandas as pd

from pycirclizely import utils
from pycirclizely.parser.matrix import Matrix

if TYPE_CHECKING:
    from ..circos import Circos


def chord_diagram_impl(
    matrix: Union[str, Path, pd.DataFrame, Matrix],
    *,
    start: float = 0,
    end: float = 360,
    space: Union[float, list[float]] = 0,
    endspace: bool = True,
    r_lim: tuple[float, float] = (97, 100),
    cmap: Union[str, dict[str, str]] = "Viridis",
    link_cmap: Optional[list[tuple[str, str, str]]] = None,
    ticks_interval: Optional[int] = None,
    order: Union[str, list[str], None] = None,
    label_kws: Optional[dict[str, Any]] = None,
    ticks_kws: Optional[dict[str, Any]] = None,
    link_kws: Optional[dict[str, Any]] = None,
    link_kws_handler: Optional[Callable[[str, str], Optional[dict[str, Any]]]] = None,
) -> Circos:
    """Chord diagram implementation for visualizing relationships in a matrix."""
    from ..circos import Circos  # Local import to avoid circular import

    link_cmap = [] if link_cmap is None else deepcopy(link_cmap)
    label_kws = {} if label_kws is None else deepcopy(label_kws)
    ticks_kws = {} if ticks_kws is None else deepcopy(ticks_kws)
    link_kws = {} if link_kws is None else deepcopy(link_kws)

    if isinstance(matrix, (str, Path, pd.DataFrame)):
        matrix = Matrix(matrix)

    if order is not None:
        matrix = matrix.sort(order)

    # Get name2color dict from user-specified colormap
    names = matrix.all_names
    name2color: dict[str, str]
    if isinstance(cmap, str):
        color_cycler = utils.ColorCycler(cmap)
        colors = [color_cycler.get_color(i) for i in range(len(names))]
        name2color = dict(zip(names, colors))
    else:
        if isinstance(cmap, defaultdict):
            name2color = cmap
        else:
            name2color = defaultdict(lambda: "grey")
            name2color.update(cmap)

    # Initialize circos sectors
    circos = Circos(matrix.to_sectors(), start, end, space=space, endspace=endspace)
    for sector in circos.sectors:
        # Plot label, outer track axis & xticks
        sector.text(sector.name, **label_kws)
        outer_track = sector.add_track(r_lim)
        color = name2color[sector.name]
        outer_track.axis(fillcolor=color)
        if ticks_interval is not None:
            outer_track.xticks_by_interval(ticks_interval, **ticks_kws)

    # Plot links
    fromto_label2color = {f"{t[0]}-->{t[1]}": t[2] for t in link_cmap}
    for link in matrix.to_links():
        from_label, to_label = link[0][0], link[1][0]
        fromto_label = f"{from_label}-->{to_label}"
        # Set link color
        if fromto_label in fromto_label2color:
            color = fromto_label2color[fromto_label]
        else:
            color = name2color[from_label]
        # Update link properties by user-defined handler function
        _link_kws = deepcopy(link_kws)
        _link_kws.update(fillcolor=color)
        if link_kws_handler is not None:
            handle_link_kws = link_kws_handler(from_label, to_label)
            if handle_link_kws is not None:
                _link_kws.update(handle_link_kws)
        circos.link(*link, **_link_kws)

    return circos
