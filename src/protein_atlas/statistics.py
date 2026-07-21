import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
from typing import Tuple

def fit_distribution_models(df: pd.DataFrame, cluster_var: str = "uniref50_cluster_id", add_residuals: bool = True) -> Tuple[pd.DataFrame, dict]:
    """
    Fits Model A (Length) and Model B (Length + Family Size).
    Calculates family-clustered sensitivity.
    """

    # Model A
    model_a = smf.ols("total_surprisal_bits ~ sequence_length", data=df).fit()

    # Model B
    model_b = smf.ols("total_surprisal_bits ~ sequence_length + log_family_size", data=df).fit()

    # Clustered Model B
    try:
        model_b_clustered = smf.ols("total_surprisal_bits ~ sequence_length + log_family_size", data=df).fit(cov_type='cluster', cov_kwds={'groups': df[cluster_var]})
        clustered_pvalues = model_b_clustered.pvalues.to_dict()
    except Exception as e:
        clustered_pvalues = None

    # Multiple testing correction using Benjamini-Hochberg FDR
    pvalues_list = list(model_b.pvalues)
    _, pvals_corrected, _, _ = multipletests(pvalues_list, alpha=0.05, method='fdr_bh')

    adjusted_pvalues = dict(zip(model_b.pvalues.index, pvals_corrected))

    if add_residuals:
        df['length_adjusted_surprisal'] = model_a.resid
        df['length_and_family_adjusted_surprisal'] = model_b.resid

    results = {
        "model_a_r2": model_a.rsquared,
        "model_b_r2": model_b.rsquared,
        "variance_explained_by_family": model_b.rsquared - model_a.rsquared,
        "model_b_params": model_b.params.to_dict(),
        "model_b_pvalues": model_b.pvalues.to_dict(),
        "model_b_adjusted_pvalues_fdr_bh": adjusted_pvalues,
        "model_b_clustered_pvalues": clustered_pvalues
    }

    return df, results
