import pandas as pd
import numpy as np
from src.protein_atlas.statistics import fit_distribution_models

def test_family_size_per_residue_model():
    # Build synthetic DataFrame where bits_per_residue depends negatively on log_family_size
    np.random.seed(42)
    n = 1000

    sequence_length = np.random.randint(100, 1000, n)
    log_family_size = np.random.uniform(0, 10, n)

    # Base bits_per_residue is 4.0, goes down by 0.05 per unit of log_family_size
    # sequence length also has some effect just to test model A vs model B
    true_bpr = 4.0 - 0.05 * log_family_size + 0.0001 * sequence_length

    # Add some noise
    bits_per_residue = true_bpr + np.random.normal(0, 0.1, n)

    # total surprisal bits (for the secondary models to not crash)
    total_surprisal_bits = bits_per_residue * sequence_length

    df = pd.DataFrame({
        'sequence_length': sequence_length,
        'log_family_size': log_family_size,
        'bits_per_residue': bits_per_residue,
        'total_surprisal_bits': total_surprisal_bits,
        'uniref50_cluster_id': np.random.randint(1, 50, n) # dummy cluster ids
    })

    df_out, results = fit_distribution_models(df)

    # Check that variance explained by family is > 0
    assert results['variance_explained_by_family'] > 0

    # Check that family coefficient is negative
    assert 'log_family_size' in results['model_b_params']
    assert results['model_b_params']['log_family_size'] < 0

    # Check that it's significant
    assert results['model_b_pvalues']['log_family_size'] < 0.05

    # Check that spearman correlation is calculated and negative
    assert results['spearman_corr_bits_vs_family'] is not None
    assert results['spearman_corr_bits_vs_family'] < 0

    # Check that secondary total models are returned
    assert 'total_model_a_r2' in results
    assert 'total_model_b_r2' in results
