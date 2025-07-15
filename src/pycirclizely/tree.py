from __future__ import annotations

import io
import os
from collections import Counter, defaultdict
from copy import deepcopy
from functools import cached_property
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Literal
from urllib.parse import urlparse
from urllib.request import urlopen

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from pycirclizely import config, utils
from pycirclizely.types import TextFormatter

if TYPE_CHECKING:
    from pycirclizely.track import Track


class TreeViz:
    """Phylogenetic Tree Visualization Class

    Interface for changing tree properties and adding tree annotations in a track
    """

    def __init__(
        self,
        tree_data: str | Path | Tree,  # type: ignore
        *,
        format: str = "newick",
        outer: bool = True,
        align_leaf_label: bool = True,
        ignore_branch_length: bool = False,
        leaf_label_size: float = 12,
        leaf_label_rmargin: float = 2.0,
        reverse: bool = False,
        ladderize: bool = False,
        line_kws: dict[str, Any] | None = None,
        align_line_kws: dict[str, Any] | None = None,
        label_formatter: TextFormatter = None,
        track: Track,
    ):
        """
        Args:
            tree_data: Tree data (`File`|`File URL`|`Tree Object`|`Tree String`).
            format: Tree format (`newick`|`phyloxml`|`nexus`|`nexml`|`cdao`).
            outer: If True, plot tree on outer side. If False, plot tree on inner side.
            align_leaf_label: If True, align leaf label.
            ignore_branch_length: If True, ignore branch length for plotting tree.
            leaf_label_size: Leaf label size.
            leaf_label_rmargin: Leaf label radius margin.
            reverse: If True, reverse tree.
            ladderize: If True, ladderize tree.
            line_kws: Shape properties.
                <https://plotly.com/python/reference/layout/shapes/>
            align_line_kws: Shape properties.
                <https://plotly.com/python/reference/layout/shapes/>
            label_formatter: User-defined label text format function
                to change label text content.
            track: Track for tree visualization.
        """
        tree = self.load_tree(tree_data, format=format)

        # Set unique node name and branch length if not exists
        tree, _ = self._set_uniq_innode_name(tree)
        self._check_node_name_dup(tree)
        max_tree_depth = max(tree.depths().values())
        if ignore_branch_length or max_tree_depth == 0:
            tree = self._to_ultrametric_tree(tree)
        if ladderize:
            tree.ladderize()
        if reverse:
            for clade in tree.find_clades():
                clade.clades = clade.clades[::-1]
        self._tree = tree

        # Set plot parameters
        self._outer = outer
        self._align_leaf_label = align_leaf_label
        self._leaf_label_size = leaf_label_size
        self._leaf_label_rmargin = leaf_label_rmargin
        self._track = track

        line_kws = {} if line_kws is None else deepcopy(line_kws)
        line_kws = utils.deep_dict_update(config.plotly_shape_defaults, line_kws)
        self._line_kws = line_kws

        align_line_kws = {} if align_line_kws is None else deepcopy(align_line_kws)
        align_line_kws = utils.deep_dict_update(
            config.tree_alignline_defaults, align_line_kws
        )
        self._align_line_kws = align_line_kws

        self._node2label_props: dict[str, dict[str, Any]] = defaultdict(lambda: {})
        self._node2line_props: dict[str, dict[str, Any]] = defaultdict(lambda: {})

        self._label_formatter = label_formatter

    ############################################################
    # Properties
    ############################################################

    @property
    def track(self) -> Track:
        """Track for tree visualization"""
        return self._track

    @property
    def tree(self) -> Tree:
        """BioPython's Tree Object"""
        return self._tree

    @cached_property
    def leaf_num(self) -> int:
        """Leaf number"""
        return len(self.leaf_labels)

    @cached_property
    def leaf_labels(self) -> list[str]:
        """Leaf labels"""
        return [str(n.name) for n in self.tree.get_terminals()]

    @cached_property
    def innode_labels(self) -> list[str]:
        """Internal node labels"""
        return [str(n.name) for n in self.tree.get_nonterminals()]

    @cached_property
    def all_node_labels(self) -> list[str]:
        """All node labels"""
        return self.leaf_labels + self.innode_labels

    @cached_property
    def max_tree_depth(self) -> float:
        """Max tree depth (root -> leaf max branch length)"""
        return max(self.tree.depths().values())

    @cached_property
    def name2xr(self) -> dict[str, tuple[float, float]]:
        """Tree node name & node xr coordinate dict"""
        return self._calc_name2xr()

    @cached_property
    def name2rect(self) -> dict[str, dict[str, float]]:
        """Tree node name & rectangle dict"""
        return self._calc_name2rect()

    ############################################################
    # Public Method
    ############################################################

    @staticmethod
    def load_tree(data: str | Path | Tree, format: str) -> Tree:
        """
        Args:
            data: Tree data.
            format: Tree format.
        """
        if isinstance(data, str) and urlparse(data).scheme in ("http", "https"):
            # Load tree file from URL
            treeio = io.StringIO(urlopen(data).read().decode(encoding="utf-8"))
            return Phylo.read(treeio, format=format)
        elif isinstance(data, (str, Path)) and os.path.isfile(data):
            # Load tree file
            with open(data, encoding="utf-8") as f:
                return Phylo.read(f, format=format)
        elif isinstance(data, str):
            # Load tree string
            return Phylo.read(io.StringIO(data), format=format)
        elif isinstance(data, Tree):
            return deepcopy(data)
        else:
            raise ValueError(f"{data=} is invalid input tree data!!")

    def search_target_node_name(
        self,
        query: str | list[str] | tuple[str],
    ) -> str:
        """
        Args:
            query: Search query node name(s). If multiple node names are set,
                MRCA(Most Recent Common Ancestor) node is set.
        """
        self._check_node_name_exist(query)
        if isinstance(query, (list, tuple)):
            target_node_name = self.tree.common_ancestor(*query).name
        else:
            target_node_name = query
        return target_node_name

    def get_target_xlim(
        self,
        query: str | list[str] | tuple[str],
    ) -> tuple[float, float]:
        """
        Args:
            query: Search query node name(s) for getting x limit.
                If multiple node names are set,
                MRCA(Most Recent Common Ancestor) node is set.
        """
        target_node_name = self.search_target_node_name(query)
        target_rect = self.name2rect[target_node_name]
        target_xlim = (
            target_rect["x0"],
            target_rect["x0"] + target_rect["width"],
        )
        return target_xlim

    def show_node_info(
        self,
        *,
        node_type: Literal["all", "internal", "leaf"] = "all",
        hover_text_formatter: Callable[[list[dict]], list[str]] | None = None,
        **kwargs: Any,
    ) -> None:
        """Show node information as hovertext.

        Args:
            node_type: Which nodes to show information for:
                - "all": All nodes
                - "internal": Only internal nodes (default)
                - "leaf": Only leaf nodes
            hover_text_formatter: User-defined function for hover text format.
            **kwargs: Scatter trace properties that override defaults.
                Note: The scatter points will be invisible by default (opacity=0).
        """
        kwargs = utils.deep_dict_update(dict(opacity=0), kwargs)
        default_bgcolor = config.tree_hovertext_defaults.get("hoverlabel", {}).get(
            "bgcolor", "lightgrey"
        )

        node: Clade
        x: list[float] = []
        r: list[float] = []
        node_info: list[dict] = []
        hover_colors: list[str] = []
        rmin, rmax = self.track.r_plot_lim

        # Determine which nodes to process based on node_type
        if node_type == "all":
            nodes_to_process = list(self.tree.find_clades())
        elif node_type == "internal":
            nodes_to_process = list(self.tree.get_nonterminals())
        elif node_type == "leaf":
            nodes_to_process = list(self.tree.get_terminals())
        else:
            raise ValueError(
                f"Invalid node_type: {node_type}. Must be 'all', 'internal', or 'leaf'"
            )

        for node in nodes_to_process:
            node_name = str(node.name)
            is_terminal = node.is_terminal()

            if node_type == "internal" and is_terminal:
                continue
            if node_type == "leaf" and not is_terminal:
                continue

            # Get target node x, r coordinate
            node_x, node_r = self.name2xr[node_name]
            if is_terminal and self._align_leaf_label:
                node_r = rmax if self._outer else rmin
            x.append(node_x)
            r.append(node_r)

            # Create node info dictionary
            node_info.append(
                {
                    "name": node_name,
                    "branch_length": getattr(node, "branch_length", None),
                    "is_terminal": is_terminal,
                    "confidence": getattr(node, "confidence", None),
                    "n_descendents": (
                        0 if is_terminal else len(list(node.find_clades())) - 1
                    ),
                }
            )

            if is_terminal:
                label_props = self._node2label_props[node_name]
                line_props = self._node2line_props[node_name]
                if "font" in label_props and "color" in label_props["font"]:
                    hover_color = label_props["font"]["color"]
                else:
                    hover_color = default_bgcolor
            else:
                if node_name in self._node2line_props:
                    line_props = self._node2line_props[node_name]
                    if "line" in line_props and "color" in line_props["line"]:
                        hover_color = line_props["line"]["color"]
                else:
                    hover_color = default_bgcolor

            hover_colors.append(hover_color)

        # Generate hover text
        if hover_text_formatter:
            hover_text = hover_text_formatter(node_info)
        else:
            hover_text = self._generate_phylogeny_hovertext(node_info)

        user_bg_color = kwargs.get("hoverlabel", {}).get("bgcolor", None)
        user_marker_color = kwargs.get("marker", {}).get("color", None)

        bg_colors = user_bg_color or user_marker_color or hover_colors
        kwargs = utils.deep_dict_update(kwargs, {"hoverlabel": {"bgcolor": bg_colors}})

        # Plot invisible scatter point with hover text
        self.track.scatter(x, r, vmin=rmin, vmax=rmax, hovertext=hover_text, **kwargs)

    def highlight(
        self,
        query: str | list[str] | tuple[str],
        **kwargs: Any,
    ) -> None:
        """
        Args:
            query: Search query node name(s) for highlight. If multiple node names
                are set, MRCA(Most Recent Common Ancestor) node is set.
            **kwargs: Shape properties
                (e.g. `fillcolor="red", line=dict(color="blue", width=2), opacity=0.5`)
                See: <https://plotly.com/python/reference/layout/shapes/>
        """
        # Get target rectangle for highlight
        target_node_name = self.search_target_node_name(query)
        rect = self.name2rect[target_node_name]

        # Set defaults for highligth
        kwargs = utils.deep_dict_update(config.tree_highlight_defaults, kwargs)

        # Setup track.rect() parameters
        start, end = rect["x0"], rect["x0"] + rect["width"]
        r_lim = (rect["y0"], rect["y0"] + rect["height"])

        self.track.rect(start, end, r_lim=r_lim, **kwargs)

    def marker(
        self,
        query: str | list[str] | tuple[str],
        *,
        descendent: bool = True,
        hover_text: list[str] | None = None,
        **kwargs: Any,
    ) -> None:
        """
        Args:
            query: Search query node name(s) for plotting marker.
                If multiple node names are set,
                MRCA(Most Recent Common Ancestor) node is set.
            descendent: If True, plot markers on target node's descendent as well.
            hover_text: Custom hover text for each marker.
                If None, generates default phylogenetic info.
            **kwargs: Scatter trace properties that override defaults.
                Common options include:
                - marker: dict with properties like size, color, symbol,..
                - mode: 'markers', 'lines', 'markers+lines'
                - name: legend name for the trace
        """
        target_node_name = self.search_target_node_name(query)

        # Set markers (x, r) coordinates (include descendent nodes)
        x: list[float] = []
        r: list[float] = []
        node_info: list[dict] = []
        rmin, rmax = self.track.r_plot_lim
        clade: Clade = next(self.tree.find_clades(target_node_name))
        descendent_nodes: list[Clade] = list(clade.find_clades())
        for descendent_node in descendent_nodes:
            node_x, node_r = self.name2xr[str(descendent_node.name)]
            if descendent_node.is_terminal() and self._align_leaf_label:
                node_r = rmax if self._outer else rmin
            x.append(node_x)
            r.append(node_r)

            node_info.append(
                {
                    "name": descendent_node.name,
                    "branch_length": getattr(descendent_node, "branch_length", None),
                    "is_terminal": descendent_node.is_terminal(),
                    "confidence": getattr(descendent_node, "confidence", None),
                    "n_descendents": (
                        len(list(descendent_node.find_clades())) - 1
                        if not descendent_node.is_terminal()
                        else 0
                    ),
                }
            )

        # If `descendent=False`, remove descendent nodes (x, r) coordinate
        if not descendent:
            x, r = [x[0]], [r[0]]
            node_info = [node_info[0]]

        # Generate default hover text if not provided
        if hover_text is None:
            hover_text = self._generate_phylogeny_hovertext(node_info)

        self.track.scatter(x, r, vmin=rmin, vmax=rmax, hover_text=hover_text, **kwargs)

    def set_node_label_props(self, target_node_label: str, **kwargs: Any) -> None:
        """
        Args:
            target_node_label: Target node label name.
            **kwargs: Annotation properties
                (e.g. `dict(font=dict(size=12, color="black"))`).
                <https://plotly.com/python/reference/layout/annotations/>.
        """
        self.search_target_node_name(target_node_label)
        self._node2label_props[target_node_label] = utils.deep_dict_update(
            self._node2label_props[target_node_label], kwargs
        )

    def set_node_line_props(
        self,
        query: str | list[str] | tuple[str],
        *,
        descendent: bool = True,
        apply_label_color: bool = False,
        **kwargs: Any,
    ) -> None:
        """
        Args:
            query: Search query node name(s) for coloring tree node line.
                If multiple node names are set,
                MRCA(Most Recent Common Ancestor) node is set.
            descendent: If True, set properties on target node's descendent as well.
            apply_label_color: If True & `descendent=True` & kwargs contain color
                keyword, apply node line color to node label color as well.
            **kwargs: Shape properties
                (e.g. `dict(line=dict(color="red", width=1, dash="dash"))`)
                See: <https://plotly.com/python/reference/layout/shapes/>
        """
        target_node_name = self.search_target_node_name(query)

        clade: Clade = next(self.tree.find_clades(target_node_name))
        if descendent:
            descendent_nodes: list[Clade] = list(clade.find_clades())
            for descendent_node in descendent_nodes:
                node_name = str(descendent_node.name)
                self._node2line_props[node_name] = kwargs
                if apply_label_color and "color" in kwargs["line"]:
                    self._node2label_props[node_name] = utils.deep_dict_update(
                        self._node2label_props[node_name],
                        dict(font=dict(color=kwargs["line"]["color"])),
                    )
        else:
            self._node2line_props[str(clade.name)] = kwargs

    ############################################################
    # Private Method
    ############################################################

    def _set_uniq_innode_name(self, tree: Tree) -> tuple[Tree, list[str]]:
        """Set unique internal node name (N_1, N_2, ..., N_XXX)."""
        tree = deepcopy(tree)
        uniq_innode_names: list[str] = []
        for idx, node in enumerate(tree.get_nonterminals(), 1):
            uniq_innode_name = f"N_{idx}"
            if node.name is None:
                node.name = uniq_innode_name
                uniq_innode_names.append(uniq_innode_name)
        return tree, uniq_innode_names

    def _to_ultrametric_tree(self, tree: Tree) -> Tree:
        """Convert to ultrametric tree."""
        tree = deepcopy(tree)
        # Get unit branch depth info
        name2depth = {str(n.name): float(d) for n, d in tree.depths(True).items()}
        name2depth = dict(sorted(name2depth.items(), key=lambda t: t[1], reverse=True))
        max_tree_depth = max(name2depth.values())
        # Reset node branch length
        for node in tree.find_clades():
            node.branch_length = None
        tree.root.branch_length = 0
        # Calculate appropriate ultrametric tree branch length
        for name, depth in name2depth.items():
            node = next(tree.find_clades(name))
            if not node.is_terminal():
                continue
            path: list[Clade] | None = tree.get_path(node)
            if path is None:
                raise ValueError(f"{name=} node not exists?")
            if depth == max_tree_depth:
                for path_node in path:
                    path_node.branch_length = 1
            else:
                # Collect nodes info which has branch length
                bl_sum, bl_exist_node_count = 0, 0
                for path_node in path:
                    if path_node.branch_length is not None:
                        bl_sum += path_node.branch_length
                        bl_exist_node_count += 1
                # Set branch length to no branch length nodes
                other_bl = (max_tree_depth - bl_sum) / (len(path) - bl_exist_node_count)
                for path_node in path:
                    if path_node.branch_length is None:
                        path_node.branch_length = other_bl
        return tree

    def _check_node_name_dup(self, tree: Tree) -> None:
        """Check node name duplication in tree."""
        all_node_names = [str(n.name) for n in tree.find_clades()]
        err_msg = ""
        for node_name, count in Counter(all_node_names).items():
            if count > 1:
                err_msg += f"{node_name=} is duplicated in tree ({count=}).\n"
        if err_msg != "":
            err_msg += "\nTreeViz cannot handle tree with duplicate node names!!"
            raise ValueError("\n" + err_msg)

    def _check_node_name_exist(
        self,
        query: str | list[str] | tuple[str],
    ) -> None:
        """Check node name exist in tree."""
        if isinstance(query, str):
            query = [query]
        err_msg = ""
        for node_name in query:
            if node_name not in self.all_node_labels:
                err_msg += f"{node_name=} is not found in tree.\n"
        if err_msg != "":
            err_msg = f"\n{err_msg}\nAvailable node names:\n{self.all_node_labels}"
            raise ValueError(err_msg)

    def _calc_name2xr(self) -> dict[str, tuple[float, float]]:
        """Calculate node name & xr coordinate."""
        track = self.track
        name2depth = {str(n.name): float(d) for n, d in self.tree.depths().items()}
        # Calculate x, r unit size of depth
        x_unit_size = track.size / self.tree.count_terminals()
        r_unit_size = track.r_plot_size / self.max_tree_depth
        # Calculate leaf node (x, r) coordinate
        name2xr: dict[str, tuple[float, float]] = {}
        node: Clade
        for idx, node in enumerate(self.tree.get_terminals()):
            x = track.start + (x_unit_size * idx) + (x_unit_size / 2)
            if self._outer:
                r = min(track.r_plot_lim) + r_unit_size * name2depth[str(node.name)]
            else:
                r = max(track.r_plot_lim) - r_unit_size * name2depth[str(node.name)]
            name2xr[str(node.name)] = (x, r)
        # Calculate internal node (x, r) coordinate
        for node in self.tree.get_nonterminals(order="postorder"):
            x = sum([name2xr[n.name][0] for n in node.clades]) / len(node.clades)
            if self._outer:
                r = min(track.r_plot_lim) + r_unit_size * name2depth[str(node.name)]
            else:
                r = max(track.r_plot_lim) - r_unit_size * name2depth[str(node.name)]
            name2xr[str(node.name)] = (x, r)
        return name2xr

    def _calc_name2rect(self) -> dict[str, dict[str, float]]:
        """Calculate tree node name & rectangle."""
        name2rect: dict[str, dict[str, float]] = {}
        for name, xr in self.name2xr.items():
            # Get parent node
            node: Clade = next(self.tree.find_clades(name))
            if node == self.tree.root:
                parent_node = node
            else:
                tree_path = self.tree.get_path(node.name)
                tree_path = [self.tree.root] + tree_path  # type: ignore
                parent_node = tree_path[-2]

            # Get child node xr coordinates
            child_node_names = [str(n.name) for n in node.find_clades()]
            x_list: list[float] = []
            r_list: list[float] = []
            for child_node_name in child_node_names:
                x, r = self.name2xr[child_node_name]
                x_list.append(x)
                r_list.append(r)

            # Calculate rectangle min-max xr coordinate
            x_unit_size = self.track.size / len(self.leaf_labels)
            xmin = min(x_list) - (x_unit_size / 2)
            xmax = max(x_list) + (x_unit_size / 2)

            parent_xr = self.name2xr[str(parent_node.name)]
            upper_r = (xr[1] + parent_xr[1]) / 2
            r_plot_lim = self.track.r_plot_lim
            lower_r = max(r_plot_lim) if self._outer else min(r_plot_lim)
            rmin, rmax = min(upper_r, lower_r), max(upper_r, lower_r)

            # Set rectangle
            rect = {
                "x0": xmin,
                "y0": rmin,
                "width": xmax - xmin,
                "height": rmax - rmin,
            }

            name2rect[name] = rect

        return name2rect

    def _plot_tree_line(self) -> None:
        """Plot tree line."""
        # Plot tree line by node (x, r) coordinate
        for node in self.tree.get_nonterminals():
            parent_x, parent_r = self.name2xr[node.name]
            child_node: Clade
            for child_node in node.clades:
                child_x, child_r = self.name2xr[str(child_node.name)]
                # Set node color if exists
                _line_kws = deepcopy(self._line_kws)
                _line_kws = utils.deep_dict_update(
                    _line_kws,
                    self._node2line_props[str(child_node.name)],
                )
                # Plot horizontal line
                h_line_points = (parent_x, child_x), (parent_r, parent_r)
                self.track._simpleline(*h_line_points, **_line_kws)
                # Plot vertical line
                v_line_points = (child_x, child_x), (parent_r, child_r)
                self.track._simpleline(*v_line_points, **_line_kws)
                # Plot vertical line for leaf node alignment
                if self._align_leaf_label and child_node.is_terminal():
                    r_plot_lim = self.track.r_plot_lim
                    end_r = max(r_plot_lim) if self._outer else min(r_plot_lim)
                    v_align_line_points = (child_x, child_x), (child_r, end_r)
                    _align_line_kws = deepcopy(self._align_line_kws)
                    _align_line_kws = utils.deep_dict_update(
                        dict(line=dict(color=_line_kws["line"]["color"])),
                        _align_line_kws,
                    )
                    self.track._simpleline(*v_align_line_points, **_align_line_kws)

    def _plot_tree_label(self) -> None:
        """Plot tree label."""
        text_kws: dict[str, Any] = dict(
            font=dict(size=self._leaf_label_size), orientation="vertical"
        )
        for node in self.tree.get_terminals():
            # Set label text (x, r) position
            label = str(node.name)
            x, r = self.name2xr[label]
            if self._align_leaf_label:
                r_plot_lim = self.track.r_plot_lim
                r = max(r_plot_lim) if self._outer else min(r_plot_lim)
            rmargin = self._leaf_label_rmargin
            r = r + rmargin if self._outer else r - rmargin

            # Set label text properties
            text_kws = utils.deep_dict_update(text_kws, self._node2label_props[label])

            # Apply label text format function if defined
            if self._label_formatter is not None:
                label = self._label_formatter(label)

            # Plot label if text size > 0
            if float(text_kws["font"]["size"]) > 0:
                self.track.text(label, x, r, outer=self._outer, axis="x", **text_kws)

    def _generate_phylogeny_hovertext(self, node_info: list[dict]) -> list[str]:
        """Generate phylogenetic hover text for nodes."""
        hovertext = []
        for info in node_info:
            parts = []
            parts.append(f"Node: {info['name']}")

            if info["is_terminal"]:
                parts.append("Type: Terminal (Leaf)")
            else:
                parts.append("Type: Internal")
                parts.append(f"Descendents: {info['n_descendents']}")

            if info["branch_length"] is not None:
                parts.append(f"Branch length: {info['branch_length']:.4f}")

            if info["confidence"] is not None:
                parts.append(f"Confidence: {info['confidence']:.2f}")

            hovertext.append("<br>".join(parts))
        return hovertext
