import argparse
import yaml
import logging
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def calculate_metrics(exact_df, approx_df):
    """
    Compare approx and exact scoring.
    """
    merged = pd.merge(exact_df, approx_df, on="uniprot_id", suffixes=('_exact', '_approx'))

    if len(merged) == 0:
        logging.warning("No overlapping proteins between exact and approximate datasets to validate.")
        return {}

    y_true = merged['bits_per_residue_exact']
    y_pred = merged['bits_per_residue_approx']

    pearson_r, _ = stats.pearsonr(y_true, y_pred)
    spearman_rho, _ = stats.spearmanr(y_true, y_pred)

    mae = np.mean(np.abs(y_true - y_pred))
    mdae = np.median(np.abs(y_true - y_pred))

    mean_exact = np.mean(y_true)
    rel_error = mae / mean_exact if mean_exact > 0 else np.nan

    bias = np.mean(y_pred - y_true)

    metrics = {
        "n_proteins": len(merged),
        "pearson_r": float(pearson_r),
        "spearman_rho": float(spearman_rho),
        "mean_absolute_error": float(mae),
        "median_absolute_error": float(mdae),
        "relative_error": float(rel_error),
        "absolute_bias": float(bias)
    }

    return metrics, merged

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    output_format = config.get("project", {}).get("output_format", "parquet")
    tables_dir = "results/tables"

    # Actually we stored them as pilot and validation
    exact_path = f"{tables_dir}/validation_protein_scores_exact.{output_format}"
    approx_path = f"{tables_dir}/pilot_protein_scores_sampled_mask.{output_format}"

    # Wait, the validation ones might have been scored as validation...
    if output_format == "parquet":
        exact_df = pd.read_parquet(exact_path) if Path(exact_path).exists() else None
        approx_df = pd.read_parquet(approx_path) if Path(approx_path).exists() else None
    else:
        exact_df = pd.read_csv(exact_path) if Path(exact_path).exists() else None
        approx_df = pd.read_csv(approx_path) if Path(approx_path).exists() else None

    if exact_df is None or approx_df is None:
        logging.error("Missing exact or approximate scoring results for validation.")
        return

    metrics, merged = calculate_metrics(exact_df, approx_df)
    if not metrics:
        return

    logging.info("--- Validation Results ---")
    for k, v in metrics.items():
        logging.info(f"{k}: {v}")

    # Check Spearman gate
    spearman = metrics.get("spearman_rho", 0)
    if spearman >= 0.90:
        logging.info("SUCCESS: Approximate method passed the Spearman >= 0.90 gate.")
    else:
        logging.warning("FAILURE: Approximate method did not reach the Spearman >= 0.90 target.")

    # Save validation report
    report_dir = "results/reports"
    Path(report_dir).mkdir(parents=True, exist_ok=True)

    report_df = pd.DataFrame([metrics])
    report_df.to_csv(f"{report_dir}/validation_metrics.csv", index=False)

    # Save plot
    fig_dir = "results/figures"
    Path(fig_dir).mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].scatter(merged['bits_per_residue_exact'], merged['bits_per_residue_approx'], alpha=0.6)
    axes[0].plot([merged['bits_per_residue_exact'].min(), merged['bits_per_residue_exact'].max()],
                 [merged['bits_per_residue_exact'].min(), merged['bits_per_residue_exact'].max()], 'r--')
    axes[0].set_xlabel("Exact (Bits per Residue)")
    axes[0].set_ylabel("Approximate (Bits per Residue)")
    axes[0].set_title(f"Scoring Validation (Spearman $\\rho$={spearman:.2f})")

    bias_vals = merged['bits_per_residue_approx'] - merged['bits_per_residue_exact']
    axes[1].hist(bias_vals, bins=20, edgecolor='black', alpha=0.7)
    axes[1].axvline(np.mean(bias_vals), color='r', linestyle='dashed', linewidth=1, label=f'Mean Bias: {np.mean(bias_vals):.3f}')
    axes[1].set_xlabel("Absolute Bias (Approx - Exact)")
    axes[1].set_ylabel("Count")
    axes[1].set_title("Absolute Bias Distribution")
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(f"{fig_dir}/fig6_scoring_validation.png", dpi=300)
    plt.savefig(f"{fig_dir}/fig6_scoring_validation.pdf")
    plt.savefig(f"{fig_dir}/fig6_scoring_validation.svg")
    logging.info(f"Saved validation figure to {fig_dir}/fig6_scoring_validation.png")

if __name__ == "__main__":
    main()
