import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from typing import Iterable, Tuple, List



def plot_inertia_silhouette_vs_k(
    Z: np.ndarray,
    *,
    ks: Iterable[int] = range(2, 15),
    n_init: int = 20,
    random_state: int = 42,
    figsize: Tuple[int, int] = (7, 4),
):
    """
    Plota Inertia (elbow) e Silhouette score vs k.
    """
    inertias: List[float] = []
    silhouettes: List[float] = []

    for k in ks:
        km = KMeans(
            n_clusters=k,
            n_init=n_init,
            random_state=random_state,
        )
        labels = km.fit_predict(Z)
        inertias.append(km.inertia_)
        silhouettes.append(silhouette_score(Z, labels))

    fig, ax1 = plt.subplots(figsize=figsize)

    # Inertia
    ax1.plot(ks, inertias, marker="o", color="tab:blue")
    ax1.set_xlabel("k")
    ax1.set_ylabel("Inertia", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")

    # Silhouette
    ax2 = ax1.twinx()
    ax2.plot(ks, silhouettes, linestyle="--", color="tab:orange")
    ax2.set_ylabel("Silhouette score", color="tab:orange")
    ax2.tick_params(axis="y", labelcolor="tab:orange")

    ax1.set_title("Elbow (Inertia) + Silhouette vs k")
    fig.tight_layout()

    return fig, (ax1, ax2)


def plot_clusters_pca(
    Z: np.ndarray,
    labels: np.ndarray,
    pca_idxs: Tuple[int, int] = (0, 1),
    *,
    alpha: float = 0.8,
    cmap: str = "tab10",
    figsize: Tuple[int, int] = (6, 5),
    show_centroids: bool = True,
):
    """
    Scatter PC_{i} vs PC{i+1} colorido por cluster, com centróides.
    """
    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(
        Z[:, pca_idxs[0]],
        Z[:, pca_idxs[1]],
        c=labels,
        cmap=cmap,
        alpha=alpha,
        s=100,
        edgecolors='k',
        linewidths=0.3,
    )

    unique_labels = np.unique(labels)
    cmap_obj = plt.get_cmap(cmap)
    norm = Normalize(vmin=labels.min(), vmax=labels.max())

    legend_handles = [
        mpatches.Patch(color=cmap_obj(norm(k)), label=f'Cluster {k}')
        for k in unique_labels
    ]

    if show_centroids:
        centroids = np.vstack([
            Z[labels == k].mean(axis=0) for k in unique_labels
        ])
        ax.scatter(
            centroids[:, pca_idxs[0]],
            centroids[:, pca_idxs[1]],
            c="black",
            marker="X",
            s=100,
            linewidths=1.5,
            zorder=5,
            alpha=0.5
        )
        legend_handles.append(
            Line2D([0], [0], marker='X', color='black', linestyle='None',
                   markersize=8, label='Centroids')
        )

    ax.legend(handles=legend_handles)

    ax.set_xlabel(f"PC{pca_idxs[0] + 1}")
    ax.set_ylabel(f"PC{pca_idxs[1] + 1}")
    ax.set_title("Cluster centroids in PCA space")

    return fig, ax