import numpy as np
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from typing import Tuple, List


def plot_cumulative_variance(
    pca: PCA,
    *,
    thresholds: List[float] = [0.90, 0.95],
    colors: List[str] = ["green", "red"],
    xlim_max: int = 10,
    figsize: Tuple[int, int] = (6, 4),
    dpi: int = 120,
):
    """
    Plots cumulative explained variance vs number of components.
    Draws horizontal + vertical reference lines for each threshold.
    """
    cumvar = np.cumsum(pca.explained_variance_ratio_)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.plot(range(1, len(cumvar) + 1), cumvar, marker="o", markersize=5)

    for threshold, color in zip(thresholds, colors):
        n = np.argmax(cumvar >= threshold) + 1
        ax.axhline(threshold, ls="--", color=color, label=f">={int(threshold*100)}% variance")
        ax.axvline(n, ls="--", color=color, label=f"{n} components")

    ax.set(
        xlabel="Number of components",
        ylabel="Cumulative Explained Variance",
        xlim=[1, xlim_max],
    )
    ax.legend()
    fig.tight_layout()

    return fig, ax


def plot_explained_variance_ratio(
    pca: PCA,
    *,
    n_components: int = 5,
    figsize: Tuple[int, int] = (6, 4),
    dpi: int = 120,
):
    """
    Bar chart of explained variance ratio per principal component.
    """
    evr = pca.explained_variance_ratio_

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.bar(range(1, len(evr) + 1), evr)
    ax.set(
        xlabel="Principal Component",
        ylabel="Explained Variance Ratio",
        xlim=[0.5, n_components + 0.5],
    )
    fig.tight_layout()

    return fig, ax


def plot_pca_scatter(
    X_pca,
    y,
    *,
    n_plots: int = 3,
    figsize_per_plot: Tuple[int, int] = (5, 6),
    dpi: int = 120,
    s: int = 20,
):
    """
    Scatter plots of PC1 vs PC2, PC1 vs PC3, ..., PC1 vs PC{n_plots+1},
    coloured by label (y).
    """
    fig, axes = plt.subplots(
        1, n_plots,
        figsize=(figsize_per_plot[0] * n_plots, figsize_per_plot[1]),
        dpi=dpi,
    )

    for i, ax in enumerate(axes):
        for label in np.unique(y):
            mask = y == label
            ax.scatter(X_pca.loc[mask, "PC1"], X_pca.loc[mask, f"PC{i + 2}"],
                       label=label, s=s)
        ax.set(xlabel="PC1", ylabel=f"PC{i + 2}", title=f"PC1 vs PC{i + 2}")

    axes[0].legend()
    fig.tight_layout()

    return fig, axes
