"""Module 02 step D: reliability diagram from LOOCV predictions (run
train_loocv.py first)."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results"


def main():
    per_genome = pd.read_csv(RESULTS_DIR / "loocv_per_genome_predictions.tsv", sep="\t")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    antibiotics = per_genome["antibiotic"].unique()

    for ax, ab in zip(axes, antibiotics):
        sub = per_genome[per_genome["antibiotic"] == ab].copy()
        sub["bin"] = pd.cut(sub["y_prob_resistant"], bins=[0, 0.2, 0.4, 0.6, 0.8, 1.0], include_lowest=True)
        grouped = sub.groupby("bin", observed=True).agg(
            n=("y_true", "size"),
            mean_pred=("y_prob_resistant", "mean"),
            observed=("y_true", "mean"),
        ).dropna()

        ax.plot([0, 1], [0, 1], linestyle="--", color="gray", label="Perfect calibration")
        ax.scatter(grouped["mean_pred"], grouped["observed"], s=grouped["n"] * 30, alpha=0.7)
        ax.set_title(f"{ab}\n(n={len(sub)}, LOOCV)")
        ax.set_xlabel("Mean predicted P(resistant)")
        ax.set_ylabel("Observed resistant rate")
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
        ax.legend(fontsize=8)

    fig.suptitle("Reliability plots: cross-antibiotic co-resistance LOOCV baseline "
                 "(NOT a genome-derived model -- see README)")
    fig.tight_layout()
    out_path = RESULTS_DIR / "reliability_plot.png"
    fig.savefig(out_path, dpi=150)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
