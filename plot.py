# plot.py

import numpy as np
import matplotlib.pyplot as plt


def plot_volcano(logFC, neg_log_P_Value, colors, output, prefix):
    legend_labels = {
        "gray": "Non-Significant",
        "red": "Positive FC",
        "blue": "Negative FC",
    }

    plt.figure(figsize=(8, 6))
    for color in np.unique(colors):
        plt.scatter(
            logFC[colors == color],
            neg_log_P_Value[colors == color],
            color=color,
            alpha=0.7,
            label=legend_labels[color],
        )

    plt.title("Volcano Plot for all samples")
    plt.xlabel("logFC (Fold Change)")
    #plt.ylabel("-log10(P.Value)")
    plt.ylabel("-log10 FDR")
    plt.legend()
    plt.grid(True)
    plt.savefig(f"{output}/{prefix}volcano_grouped.png")
