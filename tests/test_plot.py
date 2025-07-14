import random
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objs as go
from Bio import Phylo

from pycirclizely import Circos
from pycirclizely.parser import Genbank, Gff, StackedBarTable
from pycirclizely.utils import (
    ColorCycler,
    load_eukaryote_example_dataset,
    load_example_tree_file,
    load_prokaryote_example_file,
)

np.random.seed(0)
random.seed(0)


###########################################################
# Circos Class Plot
###########################################################


class TestCircosPlots:
    """Test class for Circos-level plotting methods"""

    def test_circos_axis_plot(self, fig: go.Figure):
        """Test `circos.axis()`"""
        assert isinstance(fig, go.Figure)

    def test_circos_text_plot(self, circos: Circos):
        """Test `circos.text()`"""
        circos.text("center")
        circos.text("top", r=100)
        circos.text("right", r=100, deg=90)
        circos.text("right-middle", r=50, deg=90)
        circos.text("bottom", r=100, deg=180)
        circos.text("left", r=100, deg=270)
        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_circos_line_plot(self, circos: Circos):
        """Test `circos.line()`"""
        circos.line(r=100)
        circos.line(r=80, deg_lim=(0, 270), line=dict(color="red"))
        circos.line(
            r=60, deg_lim=(90, 360), line=dict(color="blue", width=3, dash="dot")
        )

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_circos_rect_plot(self, circos: Circos):
        """Test `circos.rect()`"""
        circos.rect((80, 100))
        circos.rect((60, 80), deg_lim=(0, 270), fillcolor="tomato")
        circos.rect(
            (30, 50),
            deg_lim=(90, 360),
            fillcolor="lime",
            line=dict(color="grey", width=2),
        )
        circos.rect((30, 100), deg_lim=(0, 90), fillcolor="orange", opacity=0.2)

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_circos_link_plot(self, circos: Circos):
        """Test `circos.link()`"""
        # Plot links in various styles
        circos.link(("A", 0, 1), ("A", 7, 8))
        circos.link(("A", 1, 2), ("A", 7, 6), fillcolor="skyblue")
        circos.link(("A", 9, 10), ("B", 4, 3), direction=1, fillcolor="tomato")
        circos.link(
            ("B", 5, 7), ("C", 6, 8), direction=1, line=dict(color="black", width=1)
        )
        circos.link(
            ("B", 18, 16),
            ("B", 11, 13),
            r1=90,
            r2=90,
            fillcolor="violet",
            line=dict(color="red"),
        )
        circos.link(("C", 1, 3), ("B", 2, 0), direction=1, fillcolor="limegreen")
        circos.link(
            ("C", 11.5, 14),
            ("A", 4, 3),
            direction=2,
            fillcolor="olive",
            line=dict(color="black"),
        )

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_circos_link_line_plot(self, circos: Circos):
        """Test `circos.link_line()`"""
        # Plot link lines in various style
        circos.link_line(("A", 0), ("A", 5))
        circos.link_line(
            ("A", 5), ("B", 5), direction=1, line=dict(width=2, dash="dot")
        )
        circos.link_line(
            ("B", 20), ("C", 0), r1=80, r2=90, direction=-1, line=dict(color="blue")
        )
        circos.link_line(
            ("C", 5), ("C", 15), direction=2, arrow_height=6, arrow_width=4
        )

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_circos_colorbar_plot(self, circos: Circos):
        """Test `circos.colorbar()`"""
        # Plot colorbar in various style
        vmin1, vmax1 = 0, 100
        circos.colorbar(vmin=vmin1, vmax=vmax1)
        circos.colorbar(vmin=vmin1, vmax=vmax1, cmap="RdBu")

        vmin2, vmax2 = -200, 200
        circos.colorbar(vmin=vmin2, vmax=vmax2, cmap="viridis")
        circos.colorbar(
            vmin=vmin2,
            vmax=vmax2,
            cmap="viridis",
            orientation="v",
            title=dict(
                text="Vertical colorbar", side="bottom", font=dict(size=15, color="red")
            ),
            tickcolor="red",
            tickfont=dict(color="red", size=10),
        )

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_radar_chart_plot(self, tsv_radar_table_file: Path):
        """Test radar chart plot"""
        circos = Circos.radar_chart(tsv_radar_table_file, vmax=100, marker_size=6)

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_chord_diagram_plot(self, tsv_matrix_file: pd.DataFrame):
        """Test chord diagram plot"""
        circos = Circos.chord_diagram(tsv_matrix_file)

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_cytoband_plot(self, hg38_testdata_dir: Path):
        """Test hg38 cytoband plot"""
        # Add tracks for cytoband plot
        chr_bed_file, cytoband_file, _ = load_eukaryote_example_dataset(
            "hg38", cache_dir=hg38_testdata_dir
        )
        circos = Circos.initialize_from_bed(chr_bed_file, space=2)
        circos.add_cytoband_tracks((95, 100), cytoband_file)

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_phylogenetic_tree_plot(self):
        """Test phylogenetic tree plot"""
        tree_file = load_example_tree_file("alphabet.nwk")
        circos, tv = Circos.initialize_from_tree(tree_file)

        tv.highlight("A", fillcolor="red")
        tv.highlight(["D", "E", "F"], fillcolor="blue")

        tv.marker(["G", "H"], marker=dict(color="green"))
        tv.marker(["J", "K"], marker=dict(color="magenta"), descendent=False)

        tv.set_node_label_props("A", font=dict(color="red"))

        tv.set_node_line_props(["P", "O", "N"], line=dict(color="orange"))
        tv.set_node_line_props(["S", "R"], line=dict(color="lime"), descendent=False)
        tv.set_node_line_props(
            ["X", "Y", "Z"], line=dict(color="purple"), apply_label_color=True
        )

        tv.show_node_info()

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)


###########################################################
# Sector Class Plot
###########################################################


class TestSectorPlots:
    """Test class for Sector-level plotting methods"""

    def test_sector_axis_plot(self, circos: Circos):
        """Test `sector.axis()`"""
        sector_a = circos.get_sector("A")
        sector_a.axis()
        sector_b = circos.get_sector("B")
        sector_b.axis(fillcolor="lightgrey", line=dict(width=2, dash="dot"))
        sector_c = circos.get_sector("C")
        sector_c.axis(fillcolor="tomato", line=dict(color="blue"))

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_sector_text_plot(self, circos: Circos):
        """Test `sector.text()`"""
        name2color = {"A": "red", "B": "blue", "C": "green"}
        for sector in circos.sectors:
            sector.axis()
            sector.text(sector.name)
            color = name2color[sector.name]
            sector.text(f"Center\n{sector.name}", r=50, font=dict(color=color, size=12))

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_sector_line_plot(self, circos: Circos):
        """Test `sector_line()`"""
        for sector in circos.sectors:
            sector.axis()
            sector.line(r=90)
            sector_center = (sector.start + sector.end) / 2
            sector.line(r=80, end=sector_center, line=dict(color="red"))
            sector.line(r=80, start=sector_center, line=dict(color="blue"))
            sector.line(r=60, line=dict(color="green", width=2, dash="dot"))

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_sector_rect_plot(self, circos: Circos):
        """Test `sector.rect()`"""
        color_cycler = ColorCycler("T10")
        for sector in circos.sectors:
            sector.axis()
            sector.rect(r_lim=(90, 100), fillcolor="tomato")
            sector_center = (sector.start + sector.end) / 2
            sector.rect(end=sector_center, r_lim=(70, 80), line=dict(color="skyblue"))
            sector.rect(
                start=sector_center, r_lim=(70, 80), line=dict(color="limegreen")
            )
            for i in range(int(sector.size)):
                sector.rect(
                    i,
                    i + 1,
                    (50, 60),
                    fillcolor=color_cycler.get_color(),
                    line=dict(color="black", width=1),
                )
            start, end = sector.start + 3, sector.end - 3
            sector.rect(start, end, (30, 100), line=dict(color="orange"), opacity=0.2)

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)


###########################################################
# Track Class Plot
###########################################################


class TestTrackPlots:
    """Test class for Track-level plotting methods"""

    def test_track_axis_plot(self, circos: Circos):
        """Test `track.axis()`"""
        for sector in circos.sectors:
            track1 = sector.add_track((90, 100))
            track1.axis()
            track2 = sector.add_track((70, 80))
            track2.axis(fillcolor="lightgrey")

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_text_plot(self, circos: Circos):
        """Test `track.text()`"""
        for sector in circos.sectors:
            track1 = sector.add_track((90, 100))
            track1.axis()
            track1.text(sector.name)
            track2 = sector.add_track((70, 80))
            track2.axis(fillcolor="lightgrey")
            track2.text(
                sector.name, orientation="vertical", font=dict(color="red", size=15)
            )

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_rect_plot(self, circos: Circos):
        """Test `track.rect()`"""
        for sector in circos.sectors:
            track1 = sector.add_track((90, 100))
            track1.axis()
            # Plot rect & text (style1)
            color_cycler = ColorCycler("T10")
            for i in range(int(track1.size)):
                start, end = i, i + 1
                track1.rect(start, end, fillcolor=color_cycler.get_color())
                track1.text(str(end), (end + start) / 2)
            # Plot rect & text (style2)
            track2 = sector.add_track((70, 80))
            for i in range(int(track2.size)):
                start, end = i, i + 1
                track2.rect(
                    start,
                    end,
                    fillcolor=color_cycler.get_color(),
                    line=dict(color="white", width=1),
                )
                track2.text(
                    str(end),
                    (end + start) / 2,
                    font=dict(color="white"),
                    orientation="vertical",
                )

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_arrow_plot(self, circos: Circos):
        """Test `track.arrow()`"""
        color_cycler = ColorCycler("T10")
        sectors = {"A": 10, "B": 20, "C": 15}
        circos = Circos(sectors, space=5)
        for sector in circos.sectors:
            sector.axis()
            # Plot forward arrow with default style
            track1 = sector.add_track((90, 100))
            for i in range(int(track1.size)):
                start, end = i, i + 1
                track1.arrow(start, end, fillcolor=color_cycler.get_color())
            # Plot reverse arrow with user-specified style
            track2 = sector.add_track((70, 80))
            for i in range(int(track2.size)):
                start, end = i, i + 1
                track2.arrow(
                    end,
                    start,
                    head_length=4,
                    shaft_ratio=1.0,
                    fillcolor=color_cycler.get_color(),
                )

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_xticks_plot(self, circos: Circos):
        """Test `track.xticks()`"""
        for sector in circos.sectors:
            track1 = sector.add_track((90, 100))
            track1.axis()
            # Plot outer xticks
            pos_list = list(range(0, int(track1.size) + 1))
            labels = [f"{i:02d}" for i in pos_list]
            track1.xticks(pos_list, labels)
            # Plot inner xticks label
            labels = [f"Label{i:02d}" for i in pos_list]
            track1.xticks(
                pos_list,
                labels,
                outer=False,
                tick_length=0,
                label_margin=2,
                label_orientation="vertical",
            )

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_xticks_by_interval_plot(self, circos: Circos):
        """Test `track.xticks_by_interval()`"""
        mb_size = 10
        for sector in circos.sectors:
            # Major & Minor xticks
            track1 = sector.add_track((90, 100))
            track1.axis()
            track1.xticks_by_interval(mb_size, label_orientation="vertical")
            track1.xticks_by_interval(mb_size / 5, tick_length=1, show_label=False)
            # Mb formatted xticks
            track2 = sector.add_track((80, 90))
            track2.xticks_by_interval(
                mb_size,
                outer=False,
                show_bottom_line=True,
                label_orientation="vertical",
                label_formatter=lambda v: f"{v / mb_size:.1f} Mb",
            )

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_yticks_plot(self, circos: Circos):
        """Test `track.yticks()`"""
        for sector in circos.sectors:
            # Plot yticks
            track1 = sector.add_track((80, 100))
            track1.axis()
            y = [0, 5, 10, 15, 20]
            y_labels = list(map(str, y))
            track1.yticks(y, y_labels)
            # Plto yticks label
            track2 = sector.add_track((50, 70), r_pad_ratio=0.1)
            track2.axis()
            y = [10, 15, 20]
            y_labels = ["Label1", "Label2", "Label3"]
            track2.yticks(
                y,
                y_labels,
                vmin=10,
                vmax=25,
                side="left",
                line_kws=dict(line=dict(color="red", width=1)),
                text_kws=dict(font=dict(color="blue")),
            )

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_grid_plot(self, circos: Circos):
        """Test `track.grid()`"""
        for sector in circos.sectors:
            # Plot Y-axis grid line (Default: 6 grid line)
            track1 = sector.add_track((80, 100), r_pad_ratio=0.1)
            track1.axis()
            track1.grid()
            # Plot X-axis grid line (interval=1)
            track2 = sector.add_track((55, 75))
            track2.axis()
            track2.grid(y_grid_num=None, x_grid_interval=1, line=dict(color="red"))
            # Plot both XY-axis grid line
            track3 = sector.add_track((30, 50))
            track3.grid(
                y_grid_num=11, x_grid_interval=0.5, line=dict(color="blue", dash="dash")
            )

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_line_plot(self, circos: Circos):
        """Test `track.line()`"""
        for sector in circos.sectors:
            track = sector.add_track((80, 100), r_pad_ratio=0.1)
            track.axis()
            track.xticks_by_interval(1)
            vmin, vmax = 0, 10
            # Line between start-end two points
            track.line(
                [track.start, track.end], [vmin, vmax], line=dict(width=1.5, dash="dot")
            )
            # Line of random value points
            x = np.linspace(track.start, track.end, int(track.size) * 5 + 1)
            y = np.random.randint(vmin, vmax, len(x))
            track.line(x, y)

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_scatter_plot(self, circos: Circos):
        """Test `track.scatter()`"""
        for sector in circos.sectors:
            track = sector.add_track((80, 100), r_pad_ratio=0.1)
            track.axis()
            track.xticks_by_interval(1)
            vmin, vmax = 0, 10
            x = np.linspace(track.start, track.end, int(track.size) * 5 + 1)
            y = np.random.randint(vmin, vmax, len(x))
            track.scatter(x, y)

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_bar_plot(self, circos: Circos):
        """Test `track.bar()`"""
        for sector in circos.sectors:
            vmin, vmax = 1, 10
            x = np.linspace(sector.start + 0.5, sector.end - 0.5, int(sector.size))
            y = np.random.randint(vmin, vmax, len(x))
            # Plot bar (default)
            track1 = sector.add_track((80, 100), r_pad_ratio=0.1)
            track1.axis()
            track1.xticks_by_interval(1)
            track1.xticks_by_interval(0.1, tick_length=1, show_label=False)
            track1.bar(x, y)
            # Plot stacked bar with user-specified params
            track2 = sector.add_track((50, 70))
            track2.axis()
            track2.xticks_by_interval(1, outer=False)

            color_cycler = ColorCycler("T10")
            tab10_colors = [color_cycler.get_color() for _ in range(len(x))]
            track2.bar(
                x,
                y,
                width=1.0,
                colors=tab10_colors,
                line=dict(color="grey", width=0.5),
                vmax=vmax * 2,
            )

            color_cycler = ColorCycler("Pastel1")
            pastel_colors = [color_cycler.get_color() for _ in range(len(x))]
            y2 = np.random.randint(vmin, vmax, len(x))
            track2.bar(
                x,
                y2,
                width=1.0,
                bottom=y,
                colors=pastel_colors,
                line=dict(
                    color="grey",
                    width=0.5,
                ),
                vmax=vmax * 2,
            )

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_stacked_bar_plot(self):
        """Test `track.stacked_bar()`"""
        # Generate matrix data for stacked bar plot
        row_num, col_num = 12, 6
        matrix = np.random.randint(5, 20, (row_num, col_num))
        row_names = [f"R{i}" for i in range(row_num)]
        col_names = [f"group{i}" for i in range(col_num)]
        table_df = pd.DataFrame(matrix, index=row_names, columns=col_names)

        # Initialize Circos sector & track
        circos = Circos(sectors=dict(bar=len(table_df.index)))
        sector = circos.sectors[0]
        track = sector.add_track((50, 100))

        # Plot stacked bar
        sb_table = track.stacked_bar(
            table_df,
            width=0.6,
            cmap="Set3",
        )
        x_list = sb_table.calc_bar_label_x_list(track.size)
        track.xticks(
            x=x_list,
            labels=sb_table.row_names,
            outer=False,
            tick_length=0,
            label_margin=2,
            label_orientation="horizontal",
            text_kws=dict(font=dict(size=12)),
        )

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_stacked_barh_plot(self):
        """Test `track.stacked_barh()`"""
        # Generate & load matrix data for horizontal stacked bar plot
        row_names = list("ABCDEF")
        col_names = ["group1", "group2", "group3", "group4", "group5", "group6"]
        matrix = np.random.randint(5, 20, (len(row_names), len(col_names)))
        table_df = pd.DataFrame(matrix, index=row_names, columns=col_names)
        sb_table = StackedBarTable(table_df)

        # Initialize Circos sector & track (0 <= range <= 270)
        circos = Circos(sectors=dict(bar=sb_table.row_sum_vmax), start=0, end=270)
        sector = circos.sectors[0]
        track = sector.add_track((30, 100))
        track.axis(fillcolor="lightgrey", line=dict(color="black"), opacity=0.5)

        # Plot horizontal stacked bar & label & xticks
        track.stacked_barh(sb_table.dataframe, cmap="T10", width=0.6)
        label_r_list = sb_table.calc_barh_label_r_list(track.r_plot_lim)
        for label_r, row_name in zip(label_r_list, sb_table.row_names):
            track.text(f"{row_name} ", x=0, r=label_r, xanchor="right")
        track.xticks_by_interval(interval=5)
        track.xticks_by_interval(interval=1, tick_length=1, show_label=False)

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_fill_between_plot(self, circos: Circos):
        """Test `track.fill_between()`"""
        for sector in circos.sectors:
            vmin, vmax = 0, 10
            # Plot fill_between with simple lines
            track1 = sector.add_track((80, 100), r_pad_ratio=0.1)
            track1.axis()
            track1.xticks_by_interval(1)
            track1.fill_between(
                x=[track1.start, track1.end], y1=[vmin, vmax], y2=[vmin, vmax / 2]
            )
            # Plot fill_between with random points line
            track2 = sector.add_track((50, 70), r_pad_ratio=0.1)
            track2.axis()
            x = np.linspace(track2.start, track2.end, int(track2.size) * 5 + 1)
            y = np.random.randint(vmin, vmax, len(x))
            track2.fill_between(x, y, line=dict(color="black", width=1, dash="dash"))

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_heatmap_plot(self, circos: Circos):
        """Test `track.heatmap()`"""
        for sector in circos.sectors:
            # Plot heatmap
            track1 = sector.add_track((80, 100))
            track1.axis()
            track1.xticks_by_interval(1)
            data = np.random.randint(0, 10, (4, int(sector.size)))
            track1.heatmap(data, show_value=True)
            # Plot heatmap with labels
            track2 = sector.add_track((50, 70))
            track2.axis()
            x = np.linspace(1, int(track2.size), int(track2.size)) - 0.5
            xlabels = [str(int(v + 1)) for v in x]
            track2.xticks(x, xlabels, outer=False)
            track2.yticks([0.5, 1.5, 2.5, 3.5, 4.5], list("ABCDE"), vmin=0, vmax=5)
            data = np.random.randint(0, 100, (5, int(sector.size)))
            track2.heatmap(
                data,
                cmap="viridis",
                rect_kws=dict(fillcolor="white", line=dict(width=1)),
            )

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_tree_plot(self):
        """Test `track.heatmap()`"""
        # Load newick tree
        tree_text = "((((A:1,B:1)100:1,(C:1,D:1)100:1)100:1,(E:2,F:2)90:1):1,G:6)100;"
        tree = Phylo.read(StringIO(tree_text), "newick")
        # Initialize circos sector by tree size
        circos = Circos(sectors={"Tree": tree.count_terminals()})
        sector = circos.sectors[0]
        # Plot tree
        track = sector.add_track((50, 100))
        track.axis(line=dict(color="lightgrey"))
        track.tree(tree, leaf_label_size=12)

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_genomic_features_genbank_plot(
        self,
        prokaryote_testdata_dir: Path,
    ):
        """Test `track.genomic_features()` with genbank file"""
        # Load Genbank file
        gbk_file = load_prokaryote_example_file(
            "enterobacteria_phage.gbk",
            cache_dir=prokaryote_testdata_dir,
        )
        gbk = Genbank(gbk_file)
        # Initialize circos sector by genome size
        circos = Circos(sectors={gbk.name: gbk.range_size})
        circos.text("Enterobacteria phage\n(NC_000902)", font=dict(size=15))
        sector = circos.sectors[0]
        # Outer track
        outer_track = sector.add_track((98, 100))
        outer_track.axis(fillcolor="lightgrey")
        outer_track.xticks_by_interval(
            5000, label_formatter=lambda v: f"{v / 1000:.0f} Kb"
        )
        outer_track.xticks_by_interval(1000, tick_length=1, show_label=False)
        # Plot forward & reverse CDS genomic features
        cds_track = sector.add_track((90, 95))
        cds_track.genomic_features(
            gbk.extract_features("CDS", target_strand=1),
            plotstyle="arrow",
            fillcolor="salmon",
        )
        cds_track.genomic_features(
            gbk.extract_features("CDS", target_strand=-1),
            plotstyle="arrow",
            fillcolor="skyblue",
        )

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)

    def test_track_genomic_features_gff_plot(
        self,
        prokaryote_testdata_dir: Path,
    ):
        """Test `track.genomic_features()` with gff file"""
        # Load Genbank file
        gff_file = load_prokaryote_example_file(
            "enterobacteria_phage.gff",
            cache_dir=prokaryote_testdata_dir,
        )
        gff = Gff(gff_file)
        # Initialize circos sector by genome size
        circos = Circos(sectors={gff.name: gff.range_size})
        circos.text("Enterobacteria phage\n(NC_000902)", font=dict(size=15))
        sector = circos.sectors[0]
        # Outer track
        outer_track = sector.add_track((98, 100))
        outer_track.axis(fillcolor="lightgrey")
        outer_track.xticks_by_interval(
            5000, label_formatter=lambda v: f"{v / 1000:.0f} Kb"
        )
        outer_track.xticks_by_interval(1000, tick_length=1, show_label=False)
        # Plot forward & reverse CDS genomic features
        cds_track = sector.add_track((90, 95))
        cds_track.genomic_features(
            gff.extract_features("CDS", target_strand=1),
            plotstyle="arrow",
            fillcolor="salmon",
        )
        cds_track.genomic_features(
            gff.extract_features("CDS", target_strand=-1),
            plotstyle="arrow",
            fillcolor="skyblue",
        )

        fig = circos.plotfig()
        assert isinstance(fig, go.Figure)
