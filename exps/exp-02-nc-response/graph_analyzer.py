"""Interactive graph analysis tool for fitness landscape graphs.

Provides a GraphAnalyzer class that loads GraphML files and supports:
- Summary of macroscopic structure (peaks, top connection nodes)
- Interactive plotly visualization with hover info
- Node inspection and genotype search
"""

import json
import math

import networkx as nx
import numpy as np
import plotly.graph_objects as go


class GraphAnalyzer:
    """Analyze and visualize a fitness landscape graph from a GraphML file."""

    def __init__(self, graphml_path: str) -> None:
        self.path = graphml_path
        self.graph = self._load_graph(graphml_path)
        self._peaks: list[tuple[str, dict]] = []
        self._connections: list[tuple[str, dict]] = []
        self._classify_nodes()

    @staticmethod
    def _load_graph(graphml_path: str) -> nx.DiGraph:
        """Load GraphML and deserialize JSON-encoded attributes."""
        graph = nx.read_graphml(graphml_path)
        for node in graph.nodes():
            data = graph.nodes[node]
            if isinstance(data.get("group_mutants"), str):
                data["group_mutants"] = json.loads(data["group_mutants"])
            for key in ("fitness", "group_size", "is_peak", "contain_wildtype"):
                if key in data and isinstance(data[key], str):
                    try:
                        data[key] = float(data[key])
                    except ValueError:
                        pass
            if "group_size" in data:
                data["group_size"] = int(data["group_size"])
            if "is_peak" in data:
                data["is_peak"] = int(data["is_peak"])
        return graph

    def _classify_nodes(self) -> None:
        """Sort nodes into peaks and connections."""
        self._peaks = [
            (n, d) for n, d in self.graph.nodes(data=True) if d.get("is_peak") == 1
        ]
        self._connections = [
            (n, d) for n, d in self.graph.nodes(data=True) if d.get("is_peak") == 0
        ]
        self._peaks.sort(key=lambda x: x[1]["fitness"], reverse=True)
        self._connections.sort(key=lambda x: x[1]["group_size"], reverse=True)

    def summary(self, top_k: int = 10) -> None:
        """Print overview of the graph structure.

        Args:
            top_k: Number of top connection nodes to display.
        """
        g = self.graph
        print(f"Graph: {self.path}")
        print(f"Nodes: {g.number_of_nodes()}  Edges: {g.number_of_edges()}")
        print(f"Peaks: {len(self._peaks)}  Connection nodes: {len(self._connections)}")
        print()

        print("=== Peak nodes (sorted by fitness) ===")
        print(f"  {'Node':<22s} {'Fitness':>8s} {'Size':>6s}  Logo")
        print(f"  {'-'*22} {'-'*8} {'-'*6}  {'-'*30}")
        for n, d in self._peaks:
            print(
                f"  {n:<22s} {d['fitness']:>8.3f} {d['group_size']:>6d}  "
                f"{d.get('logo_string', '')}"
            )

        print()
        print(f"=== Top {top_k} connection nodes (sorted by group_size) ===")
        print(f"  {'Node':<22s} {'Fitness':>8s} {'Size':>6s}  Logo")
        print(f"  {'-'*22} {'-'*8} {'-'*6}  {'-'*30}")
        for n, d in self._connections[:top_k]:
            print(
                f"  {n:<22s} {d['fitness']:>8.3f} {d['group_size']:>6d}  "
                f"{d.get('logo_string', '')}"
            )

    def plot(
        self,
        highlight_nodes: list[str] | None = None,
        width: int = 900,
        height: int = 700,
    ) -> go.Figure:
        """Create interactive plotly visualization of the graph.

        X position from spring layout, Y position from fitness (height).
        Node size proportional to log(group_size). Color: red=peak, blue=connection.
        Hover shows node details.

        Args:
            highlight_nodes: Nodes to highlight with a green border.
            width: Figure width in pixels.
            height: Figure height in pixels.

        Returns:
            Plotly Figure object.
        """
        g = self.graph
        highlight_set = set(highlight_nodes or [])

        # Layout: spring for x, fitness for y
        pos_layout = nx.spring_layout(g, seed=42, k=2.0 / math.sqrt(g.number_of_nodes()))
        pos = {}
        for node in g.nodes():
            pos[node] = (pos_layout[node][0], g.nodes[node]["fitness"])

        # Edge traces
        edge_x, edge_y = [], []
        for u, v in g.edges():
            x0, y0 = pos[u]
            x1, y1 = pos[v]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])

        edge_trace = go.Scatter(
            x=edge_x,
            y=edge_y,
            mode="lines",
            line={"width": 0.5, "color": "rgba(150,150,150,0.4)"},
            hoverinfo="none",
        )

        # Node traces — separate for peaks, connections, and highlighted
        def _make_node_trace(
            nodes: list[str], color: str, name: str, symbol: str = "circle"
        ) -> go.Scatter:
            x_vals, y_vals, sizes, hovers, borders = [], [], [], [], []
            for n in nodes:
                d = g.nodes[n]
                x_vals.append(pos[n][0])
                y_vals.append(pos[n][1])
                sz = max(8, 5 * math.log1p(d["group_size"]))
                sizes.append(sz)
                hover = (
                    f"<b>{n}</b><br>"
                    f"Fitness: {d['fitness']:.3f}<br>"
                    f"Group size: {d['group_size']}<br>"
                    f"Peak: {'Yes' if d.get('is_peak') == 1 else 'No'}<br>"
                    f"Logo: {d.get('logo_string', 'N/A')}"
                )
                hovers.append(hover)
                borders.append(
                    "green" if n in highlight_set else "rgba(0,0,0,0.3)"
                )

            border_width = [3 if n in highlight_set else 1 for n in nodes]

            return go.Scatter(
                x=x_vals,
                y=y_vals,
                mode="markers",
                name=name,
                marker={
                    "size": sizes,
                    "color": color,
                    "symbol": symbol,
                    "line": {"width": border_width, "color": borders},
                    "opacity": 0.85,
                },
                text=hovers,
                hoverinfo="text",
            )

        peak_nodes = [n for n, _ in self._peaks]
        conn_nodes = [n for n, _ in self._connections]

        fig = go.Figure()
        fig.add_trace(edge_trace)
        fig.add_trace(
            _make_node_trace(conn_nodes, "rgba(70,130,210,0.8)", "Connection")
        )
        fig.add_trace(
            _make_node_trace(peak_nodes, "rgba(220,60,60,0.8)", "Peak")
        )

        fig.update_layout(
            width=width,
            height=height,
            showlegend=True,
            xaxis={"visible": False},
            yaxis={"title": "Fitness", "side": "left"},
            plot_bgcolor="white",
            hovermode="closest",
            margin={"l": 60, "r": 20, "t": 40, "b": 40},
        )

        return fig

    def inspect(self, node: str) -> None:
        """Print detailed information about a node.

        Args:
            node: Node label to inspect.
        """
        if node not in self.graph:
            print(f"Node '{node}' not found in graph.")
            return

        d = self.graph.nodes[node]
        print(f"Node: {node}")
        print(f"  Fitness:          {d['fitness']:.4f}")
        print(f"  Group size:       {d['group_size']}")
        print(f"  Is peak:          {'Yes' if d.get('is_peak') == 1 else 'No'}")
        print(f"  Contains WT:      {'Yes' if d.get('contain_wildtype') else 'No'}")
        print(f"  Logo:             {d.get('logo_string', 'N/A')}")

        # Predecessors (nodes pointing TO this node)
        preds = list(self.graph.predecessors(node))
        if preds:
            print(f"  Predecessors ({len(preds)}):")
            for p in sorted(preds, key=lambda n: self.graph.nodes[n]["fitness"], reverse=True):
                w = self.graph[p][node].get("weight", "?")
                pd = self.graph.nodes[p]
                print(
                    f"    {p:<22s}  fit={pd['fitness']:.3f}  "
                    f"size={pd['group_size']:>5d}  w={w}"
                )

        # Successors (nodes this node points TO)
        succs = list(self.graph.successors(node))
        if succs:
            print(f"  Successors ({len(succs)}):")
            for s in sorted(succs, key=lambda n: self.graph.nodes[n]["fitness"], reverse=True):
                w = self.graph[node][s].get("weight", "?")
                sd = self.graph.nodes[s]
                print(
                    f"    {s:<22s}  fit={sd['fitness']:.3f}  "
                    f"size={sd['group_size']:>5d}  w={w}"
                )

        if not preds and not succs:
            print("  No neighbors (isolated node).")

    def get_group_mutants(self, node: str) -> dict[str, float]:
        """Return the group_mutants dict for a node.

        Args:
            node: Node label.

        Returns:
            Dict mapping genotype strings to fitness values.
        """
        if node not in self.graph:
            raise KeyError(f"Node '{node}' not found in graph.")
        return self.graph.nodes[node].get("group_mutants", {})

    def find_genotype(self, genotype: str) -> str | None:
        """Find which node contains a specific genotype string.

        Args:
            genotype: A 13-character mutant profile string.

        Returns:
            The node label containing the genotype, or None.
        """
        for node, data in self.graph.nodes(data=True):
            if genotype in data.get("group_mutants", {}):
                return node
        return None
