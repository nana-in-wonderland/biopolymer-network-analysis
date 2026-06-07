"""
Biopolymer Inner Structure Analysis using Graph Theory
======================================================
Models the molecular network of a gelatin-starch biopolymer composite
using graph-theoretic methods.

Components modelled:
    - Starch          : glucose backbone with amylopectin-like branching
    - Gelatin         : peptide-bonded protein fragments
    - Calcium carbonate (CaCO₃) : inorganic filler (stone powder)
    - Acetic acid     : protonation modifier (vinegar)

Nodes  = individual monomer / molecule units (with 3D positions)
Edges  = chemical bonds (glycosidic, peptide, hydrogen, ionic)

Author : Nafis Yousefi Rad
Thesis : Master's Thesis — Biopolymer Inner Structure Analysis Based on Graph Theory
"""

import random
from collections import defaultdict

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


# ---------------------------------------------------------------------------
# Colour maps
# ---------------------------------------------------------------------------

COMPONENT_COLORS = {
    'starch':  '#FFD700',
    'gelatin': '#FF4500',
    'filler':  '#808080',
}

BOND_COLORS = {
    'glycosidic':        '#00C853',
    'glycosidic_branch': '#69F0AE',
    'peptide':           '#2979FF',
    'hydrogen':          '#FF1744',
    'ionic':             '#AA00FF',
}

NODE_SIZES = {'starch': 50, 'gelatin': 70, 'filler': 90}


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class ChemicalPolymerNetwork:
    """
    Represents a biopolymer composite as a molecular graph.

    Spatial (3D) positions are assigned stochastically to every node
    and used to determine which non-covalent interactions form,
    mimicking the distance-dependence of hydrogen and ionic bonds.
    """

    def __init__(self, seed: int = 42):
        """
        Args:
            seed: Random seed for full reproducibility.
        """
        self.seed = seed
        random.seed(seed)
        np.random.seed(seed)

    # ------------------------------------------------------------------
    # Network construction
    # ------------------------------------------------------------------

    def build_network(self, composition: dict) -> tuple:
        """
        Build a polymer network from composition.

        Args:
            composition: dict with keys 'starch', 'gelatin', 'stone', 'vinegar'
                         and numeric amounts (mass in g or volume in mL).

        Returns:
            (G, positions) where G is a NetworkX Graph and positions is a
            dict mapping node id → np.ndarray of shape (3,).
        """
        random.seed(self.seed)
        np.random.seed(self.seed)

        G         = nx.Graph()
        positions = {}
        node_id   = 0

        # --- Starch backbone with amylopectin branching -------------------
        n_starch      = int(composition['starch'] * 10)
        branch_prob   = 0.2
        prev_node     = None
        branch_offset = n_starch

        for i in range(n_starch):
            G.add_node(node_id, type='glucose', component='starch')
            positions[node_id] = np.random.rand(3)

            if prev_node is not None:
                G.add_edge(prev_node, node_id, bond_type='glycosidic')

                if random.random() < branch_prob:
                    b = branch_offset
                    branch_offset += 1
                    G.add_node(b, type='glucose', component='starch')
                    positions[b] = positions[node_id] + np.random.rand(3) * 0.1
                    G.add_edge(node_id, b, bond_type='glycosidic_branch')

            prev_node = node_id
            node_id  += 1

        node_id = branch_offset

        # --- Gelatin chain ------------------------------------------------
        n_gelatin = int(composition['gelatin'] * 10)
        for i in range(n_gelatin):
            G.add_node(node_id, type='gelatin_fragment', component='gelatin')
            positions[node_id] = np.random.rand(3)
            if i > 0:
                G.add_edge(node_id - 1, node_id, bond_type='peptide')
            node_id += 1

        # --- CaCO₃ filler -------------------------------------------------
        n_stone = int(composition['stone'] * 10)
        for i in range(n_stone):
            G.add_node(node_id, type='calcium_carbonate', component='filler')
            positions[node_id] = np.random.rand(3)
            node_id += 1

        # --- Non-covalent interactions ------------------------------------
        self._add_hydrogen_bonds(G, positions)
        self._add_ionic_interactions(G, positions)
        self._add_acetic_acid_effects(G, composition['vinegar'])

        return G, positions

    def _add_hydrogen_bonds(self, G: nx.Graph, positions: dict,
                            cutoff: float = 0.3):
        """
        Hydrogen bonds between starch (OH) and gelatin (NH / C=O) units.
        Bond forms when spatial distance < cutoff.
        """
        starch_nodes  = [n for n, d in G.nodes(data=True) if d['component'] == 'starch']
        gelatin_nodes = [n for n, d in G.nodes(data=True) if d['component'] == 'gelatin']

        for n1 in starch_nodes:
            for n2 in gelatin_nodes:
                if np.linalg.norm(positions[n1] - positions[n2]) < cutoff:
                    G.add_edge(n1, n2, bond_type='hydrogen')

    def _add_ionic_interactions(self, G: nx.Graph, positions: dict,
                                cutoff: float = 0.4):
        """
        Ionic interactions between Ca²⁺ (filler) and polymer OH/carboxylate groups.
        Bond forms when spatial distance < cutoff.
        """
        filler_nodes  = [n for n, d in G.nodes(data=True) if d['component'] == 'filler']
        polymer_nodes = [n for n, d in G.nodes(data=True)
                         if d['component'] in ('starch', 'gelatin')]

        for f in filler_nodes:
            for p in polymer_nodes:
                if np.linalg.norm(positions[f] - positions[p]) < cutoff:
                    G.add_edge(f, p, bond_type='ionic')

    def _add_acetic_acid_effects(self, G: nx.Graph, vinegar_amount: float):
        """
        Acetic acid protonates OH groups on starch and gelatin.
        Protonation probability scales with vinegar concentration.
        """
        prob = vinegar_amount / 1000.0
        for node, data in G.nodes(data=True):
            if data['component'] in ('starch', 'gelatin'):
                if random.random() < prob:
                    G.nodes[node]['protonated'] = True

    # ------------------------------------------------------------------
    # Analysis
    # ------------------------------------------------------------------

    def analyze_path_lengths(self, G: nx.Graph,
                             sample_size: int = 1000) -> tuple:
        """
        Estimate path-length distribution via random pair sampling.

        Args:
            G           : polymer network graph
            sample_size : number of node pairs to sample

        Returns:
            (path_lengths, component_paths)
        """
        nodes           = list(G.nodes())
        path_lengths    = []
        component_paths = defaultdict(list)

        for _ in range(min(sample_size, len(nodes) * (len(nodes) - 1) // 2)):
            n1, n2 = random.sample(nodes, 2)
            try:
                length = nx.shortest_path_length(G, n1, n2)
                path_lengths.append(length)
                key = f"{G.nodes[n1]['component']}-{G.nodes[n2]['component']}"
                component_paths[key].append(length)
            except nx.NetworkXNoPath:
                pass

        return path_lengths, component_paths

    def analyze_network_properties(self, G: nx.Graph,
                                   centrality_sample: int = 200) -> dict:
        """
        Compute key graph-theoretic descriptors.

        Betweenness centrality is estimated on a sampled subgraph
        to keep computation tractable on large networks.

        Args:
            G                 : polymer network graph
            centrality_sample : number of nodes for centrality estimation

        Returns:
            dict of network metrics
        """
        sampled = random.sample(list(G.nodes()), min(centrality_sample, G.number_of_nodes()))
        subG    = G.subgraph(sampled)

        return {
            'clustering':          nx.average_clustering(G),
            'density':             nx.density(G),
            'n_components':        nx.number_connected_components(G),
            'avg_degree':          np.mean([d for _, d in G.degree()]),
            'degree_distribution': [d for _, d in G.degree()],
            'betweenness':         nx.betweenness_centrality(subG, normalized=True),
        }

    # ------------------------------------------------------------------
    # Visualisation
    # ------------------------------------------------------------------

    def plot_all(self, G: nx.Graph, positions: dict,
                 path_lengths: list, properties: dict,
                 save_path: str = 'analysis.png'):
        """
        Generate a 5-panel analysis figure and save to disk.

        Panels:
            1. 3D network visualisation
            2. Path length distribution with KDE
            3. Degree distribution
            4. Betweenness centrality distribution
            5. Network statistics table
        """
        fig = plt.figure(figsize=(20, 14))
        fig.suptitle('Biopolymer Network Analysis', fontsize=18, fontweight='bold', y=0.98)

        ax1 = fig.add_subplot(231, projection='3d')
        self._draw_3d_network(ax1, G, positions)
        ax1.set_title('3D Chemical Network Structure', fontweight='bold')

        ax2 = fig.add_subplot(232)
        self._draw_path_distribution(ax2, path_lengths)
        ax2.set_title('Path Length Distribution', fontweight='bold')

        ax3 = fig.add_subplot(233)
        sns.histplot(properties['degree_distribution'], kde=True, ax=ax3,
                     color='#4C9BE8', edgecolor='white')
        ax3.set_xlabel('Node Degree')
        ax3.set_ylabel('Count')
        ax3.set_title('Degree Distribution', fontweight='bold')

        ax4 = fig.add_subplot(234)
        sns.histplot(list(properties['betweenness'].values()), kde=True, ax=ax4,
                     color='#E8A04C', edgecolor='white')
        ax4.set_xlabel('Betweenness Centrality')
        ax4.set_ylabel('Count')
        ax4.set_title('Betweenness Centrality Distribution', fontweight='bold')

        ax5 = fig.add_subplot(235)
        self._draw_stats_table(ax5, G, path_lengths, properties)

        fig.add_subplot(236).axis('off')

        plt.tight_layout()
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Saved: {save_path}")

    def plot_interaction_heatmap(self, component_paths: dict,
                                 save_path: str = 'interaction_heatmap.png'):
        """
        Heatmap of average shortest-path lengths between component pairs.
        """
        components = sorted({c.split('-')[0] for c in component_paths})
        n          = len(components)
        matrix     = np.full((n, n), np.nan)

        for i, c1 in enumerate(components):
            for j, c2 in enumerate(components):
                paths = component_paths.get(f"{c1}-{c2}", []) + \
                        component_paths.get(f"{c2}-{c1}", [])
                if paths:
                    matrix[i, j] = np.mean(paths)

        fig, ax = plt.subplots(figsize=(8, 6))
        sns.heatmap(matrix, xticklabels=components, yticklabels=components,
                    annot=True, fmt='.2f', cmap='viridis', ax=ax,
                    linewidths=0.5, linecolor='white')
        ax.set_title('Average Path Lengths Between Component Types',
                     fontsize=14, fontweight='bold')
        fig.tight_layout()
        fig.savefig(save_path, dpi=150)
        plt.close(fig)
        print(f"Saved: {save_path}")

    # ------------------------------------------------------------------
    # Private drawing helpers
    # ------------------------------------------------------------------

    def _draw_3d_network(self, ax, G: nx.Graph, positions: dict):
        for component, color in COMPONENT_COLORS.items():
            nodes = [n for n, d in G.nodes(data=True) if d['component'] == component]
            if nodes:
                pos = np.array([positions[n] for n in nodes])
                ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2],
                           c=color, s=NODE_SIZES[component],
                           label=component.title(), alpha=0.7,
                           edgecolors='white', linewidth=0.4)

        for n1, n2, data in G.edges(data=True):
            p1, p2 = positions[n1], positions[n2]
            color  = BOND_COLORS.get(data.get('bond_type', ''), '#999999')
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], [p1[2], p2[2]],
                    c=color, alpha=0.25, linewidth=0.5)

        ax.legend(fontsize=8, loc='upper left')
        ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')

    def _draw_path_distribution(self, ax, path_lengths: list):
        sns.histplot(path_lengths, kde=True, stat='density', ax=ax,
                     color='#4CE87A', edgecolor='white')
        mean_pl   = np.mean(path_lengths)
        median_pl = np.median(path_lengths)
        ax.axvline(mean_pl,   color='#E84C4C', linestyle='--',
                   label=f'Mean: {mean_pl:.2f}')
        ax.axvline(median_pl, color='#4C9BE8', linestyle='--',
                   label=f'Median: {median_pl:.2f}')
        ax.set_xlabel('Path Length')
        ax.set_ylabel('Density')
        ax.legend(fontsize=9)

    def _draw_stats_table(self, ax, G: nx.Graph,
                          path_lengths: list, properties: dict):
        ax.axis('off')
        stats_text = (
            f"Network Statistics\n"
            f"{'─' * 28}\n"
            f"Nodes               {G.number_of_nodes():>8}\n"
            f"Edges               {G.number_of_edges():>8}\n"
            f"Connected components{properties['n_components']:>8}\n"
            f"Average degree      {properties['avg_degree']:>8.2f}\n"
            f"Network density     {properties['density']:>8.4f}\n"
            f"Clustering coeff.   {properties['clustering']:>8.3f}\n"
            f"Mean path length    {np.mean(path_lengths):>8.2f}\n"
            f"Median path length  {np.median(path_lengths):>8.2f}\n"
        )
        ax.text(0.05, 0.95, stats_text,
                fontsize=10, verticalalignment='top',
                fontfamily='monospace',
                transform=ax.transAxes,
                bbox=dict(boxstyle='round', facecolor='#F5F5F5', alpha=0.9))


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    composition = {
        'gelatin': 46,    # g
        'starch':  36,    # g
        'stone':   18,    # g  (stone powder ≈ CaCO₃)
        'vinegar': 34,    # mL (~5 % acetic acid)
    }

    print("Building biopolymer network...")
    polymer = ChemicalPolymerNetwork(seed=42)
    G, positions = polymer.build_network(composition)

    print("Analysing path lengths...")
    path_lengths, component_paths = polymer.analyze_path_lengths(G, sample_size=1000)

    print("Computing network properties...")
    properties = polymer.analyze_network_properties(G, centrality_sample=200)

    print("\n=== Network Statistics ===")
    for key, val in properties.items():
        if key not in ('degree_distribution', 'betweenness'):
            print(f"  {key:<25}: {val:.4f}" if isinstance(val, float) else f"  {key:<25}: {val}")

    print("\nGenerating plots...")
    polymer.plot_all(G, positions, path_lengths, properties)
    polymer.plot_interaction_heatmap(component_paths)

    print("\nDone.")


if __name__ == "__main__":
    main()
