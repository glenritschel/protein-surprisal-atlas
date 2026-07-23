import pandas as pd
import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

def fit_association_models(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    baseline_vars = ['sequence_length', 'log_family_size', 'low_complexity_fraction']
    target = 'bits_per_residue'

    annotations = ['gnomad_loeuf', 'gnomad_pli', 'disorder_fraction']

    results = []

    for A in annotations:
        if A not in df.columns:
            continue

        # Complete case subset for this specific annotation + baseline vars
        subset = df[[target, A, 'uniref50_cluster_id'] + baseline_vars].dropna().copy()

        n = len(subset)
        if n < 10:
            continue

        # Fit baseline model
        baseline_formula = f"{target} ~ sequence_length + log_family_size + low_complexity_fraction"
        baseline_model = smf.ols(baseline_formula, data=subset).fit()
        r2_baseline = baseline_model.rsquared

        # Fit full model with all covariates together: target ~ covariates + A
        full_formula = f"{baseline_formula} + {A}"
        full_model = smf.ols(full_formula, data=subset).fit()
        r2_full = full_model.rsquared

        delta_r2 = r2_full - r2_baseline

        coef = full_model.params[A]
        ols_p = full_model.pvalues[A]

        # Cluster robust p-value
        try:
            robust_model = full_model.get_robustcov_results(cov_type='cluster', groups=subset['uniref50_cluster_id'])
            robust_p = robust_model.pvalues[robust_model.model.exog_names.index(A)]
        except Exception:
            robust_p = np.nan

        # Robust check: Spearman(length_and_family_adjusted_surprisal, A)
        # Note: length_and_family_adjusted_surprisal must be pre-computed before passing into this function.
        if 'length_and_family_adjusted_surprisal' in df.columns:
            subset_spearman = df[['length_and_family_adjusted_surprisal', A]].dropna()
            if len(subset_spearman) > 0:
                spearman_corr, spearman_p = spearmanr(subset_spearman['length_and_family_adjusted_surprisal'], subset_spearman[A])
            else:
                spearman_corr, spearman_p = np.nan, np.nan
        else:
            spearman_corr, spearman_p = np.nan, np.nan

        results.append({
            'annotation': A,
            'n': n,
            'coefficient': coef,
            'ols_p': ols_p,
            'cluster_robust_p': robust_p,
            'delta_r2': delta_r2,
            'spearman_corr': spearman_corr,
            'spearman_p': spearman_p
        })

    res_df = pd.DataFrame(results)

    if not res_df.empty:
        valid_p = res_df['cluster_robust_p'].notna()
        res_df['adjusted_p'] = np.nan
        if valid_p.any():
            _, pvals_corrected, _, _ = multipletests(res_df.loc[valid_p, 'cluster_robust_p'], method='fdr_bh')
            res_df.loc[valid_p, 'adjusted_p'] = pvals_corrected

    return res_df
