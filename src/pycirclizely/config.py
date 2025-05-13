# from __future__ import annotations

import math
from enum import IntEnum
from typing import ClassVar

###########################################################
# Constant Value Config
###########################################################

# Fundamental Plot Parameters
MIN_R = 0
MAX_R = 100
R_PLOT_MARGIN = 15
ARC_POINTS = 100
R_LIM = (MIN_R, MAX_R)
AXIS_FACE_PARAM = dict(layer="below", line=dict(color="rgba(0,0,0,0)"))
AXIS_EDGE_PARAM = dict(layer="above", fillcolor=None)
REL_TOL = 1e-10  # Relative Tolerance
AXIS_RANGE = [-MAX_R - R_PLOT_MARGIN, MAX_R + R_PLOT_MARGIN]

# Circos Color Scheme
# http://circos.ca/tutorials/lessons/configuration/colors/
CYTOBAND_COLORMAP = {
    "gpos100": "#000000",  # 0,0,0
    "gpos": "#000000",  # 0,0,0
    "gpos75": "#828282",  # 130,130,130
    "gpos66": "#A0A0A0",  # 160,160,160
    "gpos50": "#C8C8C8",  # 200,200,200
    "gpos33": "#D2D2D2",  # 210,210,210
    "gpos25": "#C8C8C8",  # 200,200,200
    "gvar": "#DCDCDC",  # 220,220,220
    "gneg": "#FFFFFF",  # 255,255,255
    "acen": "#D92F27",  # 217,47,39
    "stalk": "#647FA4",  # 100,127,164
}


class Direction(IntEnum):
    """Link BezierCurve Direction Enum"""

    REVERSE = -1
    NONE = 0
    FORWARD = 1
    BIDIRECTIONAL = 2


###########################################################
# Plotly Default Configuration
###########################################################

plotly_layout_defaults = {
    "title": {
        "font": {"color": "black", "family": "Times New Roman", "size": 18},
        "text": None,
    },
    "hovermode": "closest",
    "showlegend": False,
    "xaxis": {
        "autorange": True,
        "showgrid": False,
        "zeroline": False,
        "showticklabels": False,
    },
    "yaxis": {
        "autorange": True,
        "showgrid": False,
        "zeroline": False,
        "showticklabels": False,
    },
    "plot_bgcolor": "rgba(0,0,0,0)",  # Transparent background inside the axes
}

# Plotly annotation defaults
plotly_annotation_defaults = {
    "font": {
        "size": 10,
        "color": "black",
    },
    "showarrow": False,
}

# Plotly shape defaults
plotly_shape_defaults = {
    "fillcolor": None,
    "line": {"color": "black", "width": 1},
    "layer": "above",
}

# Plotly link defaults
plotly_link_defaults = {
    "fillcolor": "grey",
    "opacity": 0.5,
    "line": {"width": 0.2, "color": "white"},
    "layer": "above",
}

# Text orientation-specific overrides
plotly_text_orientation_defaults = {
    "horizontal": {
        "yanchor": "top",  # Default for outer horizontal text
    },
    "vertical": {
        "xanchor": "left",  # Default for outer vertical text
    },
}

plotly_scatter_defaults = {
    'marker': {
        'size': 6,
        'opacity': 1.0,
    },
    'line': {
        'width': 2,
    },
    'hoverinfo': 'text',
    'showlegend': False
}

###########################################################
# GitHub Eukaryote & Prokaryote Dataset Config
###########################################################

# GITHUB_DATA_URL = "https://raw.githubusercontent.com/moshi4/pycirclizely-data/master/"

EUKARYOTE_DATASET = {
    "hg38": [
        "hg38_chr.bed",
        "hg38_cytoband.tsv",
        "hg38_genomic_link.tsv",
    ],
    "hs1": [
        "hs1_chr.bed",
        "hs1_cytoband.tsv",
        "hs1_genomic_link.tsv",
    ],
    "mm10": [
        "mm10_chr.bed",
        "mm10_cytoband.tsv",
        "mm10_genomic_link.tsv",
    ],
    "mm39": [
        "mm39_chr.bed",
        "mm39_cytoband.tsv",
        "mm39_genomic_link.tsv",
    ],
}

PROKARYOTE_FILES = [
    "enterobacteria_phage.gbk",
    "enterobacteria_phage.gff",
    "mycoplasma_alvi.gbk",
    "mycoplasma_alvi.gff",
    "escherichia_coli.gbk.gz",
    "escherichia_coli.gff.gz",
]
