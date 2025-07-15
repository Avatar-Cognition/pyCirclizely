from __future__ import annotations

from enum import IntEnum

###########################################################
# Constant Value Config
###########################################################

# Fundamental Plot Parameters
MIN_R = 0
MAX_R = 100
R_PLOT_MARGIN = 13
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
    "showlegend": True,
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

plotly_annotation_defaults = {
    "font": {
        "size": 11,
        "color": "black",
    },
    "showarrow": False,
    "xanchor": "center",
    "yanchor": "middle",
}

plotly_shape_defaults = {
    "fillcolor": None,
    "line": {"color": "black", "width": 1},
    "layer": "between",
}

plotly_grid_defaults = {
    "fillcolor": None,
    "line": {"color": "grey", "width": 0.5},
    "opacity": 0.5,
    "layer": "below",
}

plotly_arrow_defaults = {
    "fillcolor": "grey",
    "opacity": 0.5,
    "line": {"width": 0.2, "color": "white"},
    "layer": "between",
}

plotly_linelink_defaults = {
    "line": {"width": 1.2, "color": "grey"},
    "layer": "above",
}

plotly_scatter_defaults = {
    "marker": {
        "size": 4,
        "opacity": 1.0,
        "symbol": "circle",
    },
    "line": {
        "width": 2,
    },
    "hoverinfo": "text",
    "showlegend": False,
}

tree_alignline_defaults = {
    "line": {
        "dash": "dot",
    },
    "opacity": 0.5,
    "layer": "between",
}

tree_highlight_defaults = {
    "line": {
        "width": 0,
    },
    "opacity": 0.6,
    "fillcolor": "lightyellow",
    "layer": "below",
}

tree_hovertext_defaults = {
    "hoverlavel": {
        "bgcolor": "lightgrey",
    },
}

cytoband_defaults = {
    "line": {
        "width": 0,
    },
    "layer": "below",
}

radar_line_defaults = {
    "line": {
        "width": 2,
    },
    "layer": "between",
}

radar_grid_defaults = {
    "line": {
        "color": "grey",
        "width": 0.5,
        "dash": "dash",
    },
    "layer": "below",
}

radar_annotation_defaults = {
    "font": {
        "color": "dimgrey",
        "size": 10,
    },
    "xanchor": "left",
    "yanchor": "top",
}

hovertext_dummy_defaults = {
    "hoverinfo": "text",
    "marker": {
        "size": 20,
        "opacity": 0,
    },
}


###########################################################
# GitHub Eukaryote & Prokaryote Dataset Config
###########################################################

GITHUB_DATA_URL = "https://raw.githubusercontent.com/moshi4/pycirclize-data/master/"

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
