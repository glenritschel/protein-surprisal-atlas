import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from typing import Tuple
import scipy.stats
import logging

def fit_distribution_models(df: pd.DataFrame, cluster_var: str = "uniref50_cluster_id", add_residuals: bool = True) -> Tuple[pd.DataFrame, dict]:
    """
    Fits Model A (Length) and Model B (Length + Family Size).
    Calculates family-clustered sensitivity.
    """

    d = df.dropna(subset=['bits_per_residue', 'sequence_length', 'log_family_size']).copy()

    # Primary Models: bits_per_residue
    model_a = smf.ols("bits_per_residue ~ sequence_length", data=d).fit()
    model_b = smf.ols("bits_per_residue ~ sequence_length + log_family_size", data=d).fit()

    # Secondary Models: total_surprisal_bits (diagnostic)
    total_model_a = smf.ols("total_surprisal_bits ~ sequence_length", data=d).fit()
    total_model_b = smf.ols("total_surprisal_bits ~ sequence_length + log_family_size", data=d).fit()

    # Clustered Model B
    try:
        model_b_clustered = smf.ols("bits_per_residue ~ sequence_length + log_family_size", data=d).fit(cov_type='cluster', cov_kwds={'groups': d[cluster_var]})
        clustered_pvalues = model_b_clustered.pvalues.to_dict()
    except Exception as e:
        logging.warning(f"Clustered model fit failed: {e}")
        clustered_pvalues = None

    # Multiple testing correction using Benjamini-Hochberg FDR
    pvalues_list = list(model_b.pvalues)
    _, pvals_corrected, _, _ = multipletests(pvalues_list, alpha=0.05, method='fdr_bh')

    adjusted_pvalues = dict(zip(model_b.pvalues.index, pvals_corrected))

    if add_residuals:
        df.loc[d.index, 'length_adjusted_surprisal'] = model_a.resid
        df.loc[d.index, 'length_and_family_adjusted_surprisal'] = model_b.resid

    # Spearman correlation
    spearman_corr = None
    if 'bits_per_residue' in df.columns and 'log_family_size' in df.columns:
        valid_df = df[['bits_per_residue', 'log_family_size']].dropna()
        if len(valid_df) > 0:
            spearman_corr, _ = scipy.stats.spearmanr(valid_df['bits_per_residue'], valid_df['log_family_size'])

    results = {
        "model_a_r2": model_a.rsquared,
        "model_b_r2": model_b.rsquared,
        "variance_explained_by_family": model_b.rsquared - model_a.rsquared,
        "total_model_a_r2": total_model_a.rsquared,
        "total_model_b_r2": total_model_b.rsquared,
        "spearman_corr_bits_vs_family": spearman_corr,
        "model_b_params": model_b.params.to_dict(),
        "model_b_pvalues": model_b.pvalues.to_dict(),
        "model_b_adjusted_pvalues_fdr_bh": adjusted_pvalues,
        "model_b_clustered_pvalues": clustered_pvalues
    }

    return df, results
