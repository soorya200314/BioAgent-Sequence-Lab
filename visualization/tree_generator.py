"""Phylogenetic tree visualization module using SciPy and Matplotlib."""

import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform


def generate_tree(
    distance_matrix: list[list[float]],
    labels: list[str],
    output_path: str = None,
) -> plt.Figure:
    """Generate a phylogenetic tree (dendrogram) from a distance matrix and save/return it.

    Args:
        distance_matrix: Square, symmetric 2D list of pairwise distances.
        labels: List of species/sequence names corresponding to matrix indices.
        output_path: Optional filepath where the resulting PNG image should be saved.

    Returns:
        The Matplotlib figure object.
    """
    # Convert distance matrix to a numpy array
    dist_arr = np.array(distance_matrix)

    # Enforce symmetry and diagonal zeros to avoid numeric precision issues
    dist_arr = (dist_arr + dist_arr.T) / 2.0
    np.fill_diagonal(dist_arr, 0.0)

    # Convert square matrix to condensed 1D distance vector required by SciPy
    condensed_dist = squareform(dist_arr)

    # Perform UPGMA (Average linkage) clustering
    linkage_matrix = linkage(condensed_dist, method="average")

    # Set up matplotlib style (clean, modern aesthetics)
    fig = plt.figure(figsize=(10, 6), facecolor="#FAFAFA")
    ax = plt.gca()
    ax.set_facecolor("#FAFAFA")

    # Draw the dendrogram
    # 'left' orientation places root on the left, leaves on the right
    dendrogram_data = dendrogram(
        linkage_matrix,
        labels=labels,
        orientation="left",
        leaf_font_size=12,
        color_threshold=None,  # Do not auto-color based on threshold
        link_color_func=lambda x: "#2C3E50",  # Sleek dark blue-grey links
    )

    # Aesthetic enhancements
    plt.title(
        "Phylogenetic Tree (UPGMA Hierarchical Clustering)",
        fontsize=14,
        fontweight="bold",
        color="#2C3E50",
        pad=20,
    )
    plt.xlabel("Evolutionary Distance (1 - Sequence Identity)", fontsize=11, labelpad=10, color="#34495E")

    # Clean axes
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("#BDC3C7")
    ax.spines["bottom"].set_linewidth(1.2)
    ax.tick_params(axis="x", colors="#34495E", labelsize=10)
    ax.tick_params(axis="y", left=False)  # Remove y tick marks since labels are text

    # Add a horizontal grid line for clarity
    ax.grid(axis="x", linestyle="--", alpha=0.5, color="#BDC3C7")

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=300, facecolor="#FAFAFA")
        plt.close(fig)
    return fig

