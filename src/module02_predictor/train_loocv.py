"""
Module 02 step B: LOOCV logistic regression baseline, evaluated on real
cross-antibiotic co-resistance features (see data_prep.py docstring for
why this baseline exists instead of a genome-derived-feature model).

Why LOOCV: with n<=25 per antibiotic, a held-out train/test split would
leave single-digit test sets with unstable metrics; Leave-One-Out
Cross-Validation uses every sample as a test point exactly once and is
the standard choice for cohorts this small.

For each target antibiotic:
  - Drop genomes with an "Intermediate" call for that antibiotic from the
    binary R/S modeling task (too few I calls to model as a third class;
    they are reported separately, not silently dropped from the report).
  - Features = every OTHER antibiotic's real phenotype for that genome,
    encoded Resistant=1, Susceptible=0, Intermediate=0.5, missing=NaN.
  - Only features with >=15/25 non-missing real values across the cohort
    are used, to avoid one or two-observation columns dominating a
    5-parameter-ish logistic fit.
  - Missing values are imputed using the TRAINING FOLD's mean only
    (recomputed inside each LOOCV split) to avoid test-set leakage.
  - LogisticRegression(penalty="l2", class_weight="balanced").
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    balanced_accuracy_score, f1_score, roc_auc_score,
    average_precision_score, brier_score_loss, recall_score,
    confusion_matrix,
)

from data_prep import load_labeled_cohort, get_all_antibiotics, TARGET_ANTIBIOTICS

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results"

ENCODE = {"Resistant": 1.0, "Susceptible": 0.0, "Intermediate": 0.5}
MIN_NON_MISSING = 15
NO_CALL_CONFIDENCE_THRESHOLD = 0.70  # matches src/module02_predictor/predict.py


def build_feature_frame(table: dict) -> pd.DataFrame:
    genome_ids = sorted(table.keys())
    antibiotics = get_all_antibiotics(table)
    df = pd.DataFrame(index=genome_ids, columns=antibiotics, dtype=float)
    for gid, labels in table.items():
        for ab, phenotype in labels.items():
            df.loc[gid, ab] = ENCODE.get(phenotype, np.nan)
    return df


def loocv_for_antibiotic(df: pd.DataFrame, target_ab: str) -> dict:
    eligible = df.index[df[target_ab].isin([0.0, 1.0])]
    excluded_intermediate = df.index[df[target_ab] == 0.5].tolist()

    feature_cols = [c for c in df.columns if c != target_ab]
    non_missing_counts = df.loc[eligible, feature_cols].notna().sum()
    feature_cols = non_missing_counts[non_missing_counts >= MIN_NON_MISSING].index.tolist()

    y_true, y_prob, y_pred, genome_order = [], [], [], []

    for held_out in eligible:
        train_idx = [g for g in eligible if g != held_out]
        X_train_raw = df.loc[train_idx, feature_cols]
        col_means = X_train_raw.mean()
        X_train = X_train_raw.fillna(col_means).values
        y_train = df.loc[train_idx, target_ab].values.astype(int)

        X_test = df.loc[[held_out], feature_cols].fillna(col_means).values

        clf = LogisticRegression(penalty="l2", class_weight="balanced", max_iter=1000)
        clf.fit(X_train, y_train)
        prob_resistant = clf.predict_proba(X_test)[0, 1]

        y_true.append(int(df.loc[held_out, target_ab]))
        y_prob.append(float(prob_resistant))
        y_pred.append(int(prob_resistant >= 0.5))
        genome_order.append(held_out)

    y_true_arr = np.array(y_true)
    y_pred_arr = np.array(y_pred)
    y_prob_arr = np.array(y_prob)

    tn, fp, fn, tp = confusion_matrix(y_true_arr, y_pred_arr, labels=[0, 1]).ravel()

    confidence = np.maximum(y_prob_arr, 1 - y_prob_arr)
    would_be_no_call = confidence < NO_CALL_CONFIDENCE_THRESHOLD
    no_call_rate = float(would_be_no_call.mean())
    if (~would_be_no_call).sum() > 0:
        called_balanced_acc = float(balanced_accuracy_score(
            y_true_arr[~would_be_no_call], y_pred_arr[~would_be_no_call]))
    else:
        called_balanced_acc = None

    metrics = {
        "antibiotic": target_ab,
        "n_eligible_RS": int(len(eligible)),
        "n_excluded_intermediate": len(excluded_intermediate),
        "excluded_intermediate_genomes": excluded_intermediate,
        "n_features_used": len(feature_cols),
        "features_used": feature_cols,
        "balanced_accuracy": float(balanced_accuracy_score(y_true_arr, y_pred_arr)),
        "f1": float(f1_score(y_true_arr, y_pred_arr, zero_division=0)),
        "resistant_recall": float(recall_score(y_true_arr, y_pred_arr, pos_label=1, zero_division=0)),
        "susceptible_recall": float(recall_score(y_true_arr, y_pred_arr, pos_label=0, zero_division=0)),
        "brier_score": float(brier_score_loss(y_true_arr, y_prob_arr)),
        "confusion_matrix": {"tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp)},
        "no_call_rate_at_threshold_0.70": no_call_rate,
        "balanced_accuracy_on_called_only": called_balanced_acc,
    }

    if len(set(y_true_arr.tolist())) == 2:
        metrics["auroc"] = float(roc_auc_score(y_true_arr, y_prob_arr))
        metrics["pr_auc"] = float(average_precision_score(y_true_arr, y_prob_arr))
    else:
        metrics["auroc"] = None
        metrics["pr_auc"] = None
        metrics["note"] = "AUROC/PR-AUC undefined: only one class present among eligible genomes."

    per_genome = pd.DataFrame({
        "genome_id": genome_order,
        "antibiotic": target_ab,
        "y_true": y_true_arr,
        "y_prob_resistant": y_prob_arr,
        "y_pred": y_pred_arr,
    })
    return metrics, per_genome


def reliability_bins(per_genome_all: pd.DataFrame, n_bins: int = 5) -> pd.DataFrame:
    df = per_genome_all.copy()
    df["bin"] = pd.cut(df["y_prob_resistant"], bins=np.linspace(0, 1, n_bins + 1), include_lowest=True)
    grouped = df.groupby("bin", observed=True).agg(
        n=("y_true", "size"),
        mean_predicted=("y_prob_resistant", "mean"),
        observed_resistant_rate=("y_true", "mean"),
    ).reset_index()
    grouped["bin"] = grouped["bin"].astype(str)
    return grouped


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    table = load_labeled_cohort()
    df = build_feature_frame(table)

    all_metrics = {}
    per_genome_frames = []
    for ab in TARGET_ANTIBIOTICS:
        metrics, per_genome = loocv_for_antibiotic(df, ab)
        all_metrics[ab] = metrics
        per_genome_frames.append(per_genome)
        print(f"\n=== {ab} ===")
        print(f"  n (R/S, Intermediate excluded)={metrics['n_eligible_RS']}  "
              f"excluded_intermediate={metrics['n_excluded_intermediate']}")
        print(f"  features used ({metrics['n_features_used']}): {metrics['features_used']}")
        print(f"  balanced_accuracy={metrics['balanced_accuracy']:.3f}  f1={metrics['f1']:.3f}  "
              f"R_recall={metrics['resistant_recall']:.3f}  S_recall={metrics['susceptible_recall']:.3f}")
        print(f"  AUROC={metrics['auroc']}  PR-AUC={metrics['pr_auc']}  Brier={metrics['brier_score']:.3f}")
        print(f"  confusion_matrix={metrics['confusion_matrix']}")
        print(f"  no_call_rate(conf<0.70)={metrics['no_call_rate_at_threshold_0.70']:.3f}  "
              f"balanced_acc_on_called_only={metrics['balanced_accuracy_on_called_only']}")

    per_genome_all = pd.concat(per_genome_frames, ignore_index=True)
    per_genome_all.to_csv(RESULTS_DIR / "loocv_per_genome_predictions.tsv", sep="\t", index=False)

    with open(RESULTS_DIR / "loocv_metrics.json", "w") as f:
        json.dump(all_metrics, f, indent=2)

    reliability = reliability_bins(per_genome_all)
    reliability.to_csv(RESULTS_DIR / "reliability_bins.tsv", sep="\t", index=False)

    print(f"\nWrote per-genome predictions, metrics, and reliability bins to {RESULTS_DIR}")


if __name__ == "__main__":
    main()
