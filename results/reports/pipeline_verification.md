# Pipeline Verification Report

## Headline Numbers
* **Spearman $\rho$ (Exact vs Approximate Scoring):** 1.0 (with small pilot subset size of 2)
* **Absolute Bias (Approximate - Exact):** 0.0324 bits/residue
* **Model A (Length Only) $R^2$:** 0.9434
* **Model B (Length + Family Size) $R^2$:** 0.9447
* **Variance Explained by `log_family_size`:** 0.0013 (0.13%)

## Generated Figures
The following figures were successfully produced:
* `fig1_bits_per_residue_hist.png` / `.pdf`
* `fig2_length_vs_surprisal.png` / `.pdf`
* `fig3_background_vs_model.png` / `.pdf`
* `fig4_surprisal_ratios.png` / `.pdf`
* `fig5_controls.png` / `.pdf`
* `fig6_scoring_validation.png` / `.pdf` / `.svg`
* `fig7_confound.png` / `.pdf`
* `fig12_TP53 P53_profile.png` / `.pdf`
* `fig13_EGFR ERBB ERBB1 HER1_profile.png` / `.pdf`

## Generated Tables
The following tables were successfully produced:
* `bottom_50_normalized_surprisal.csv`
* `top_50_normalized_surprisal.csv`
* `top_50_negative_residuals.csv`
* `top_50_positive_residuals.csv`
* `control_sequence_comparisons.csv`
* `variant_scores_TP53.csv`
* `pilot_protein_scores_sampled_mask.parquet`
* `pilot_residue_scores_sampled_mask.parquet`
* `validation_protein_scores_exact.parquet`
* `validation_residue_scores_exact.parquet`
