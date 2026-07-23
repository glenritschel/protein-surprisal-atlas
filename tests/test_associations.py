import pandas as pd
import numpy as np
from src.protein_atlas.associations import fit_association_models

def test_fit_association_models():
    # Create synthetic data
    np.random.seed(42)
    n = 100
    df = pd.DataFrame({
        'bits_per_residue': np.random.randn(n) * 0.1 + 2.0,
        'sequence_length': np.random.randint(100, 500, n),
        'log_family_size': np.random.uniform(0, 10, n),
        'low_complexity_fraction': np.random.uniform(0, 0.5, n),
        'uniref50_cluster_id': [f"cluster_{i%10}" for i in range(n)],
        'length_and_family_adjusted_surprisal': np.random.randn(n) * 0.1,
    })

    # Introduce a strong negative association with disorder_fraction
    df['disorder_fraction'] = np.random.uniform(0, 1, n)
    df['bits_per_residue'] -= 2.0 * df['disorder_fraction']

    # Introduce a strong positive association with gnomad_loeuf
    df['gnomad_loeuf'] = np.random.uniform(0, 2, n)
    df['bits_per_residue'] += 1.5 * df['gnomad_loeuf']

    # Introduce a negative association with depmap_gene_effect
    df['depmap_gene_effect'] = np.random.uniform(-2, 0, n)
    df['bits_per_residue'] += 1.0 * df['depmap_gene_effect'] # depmap effect is negative, this reduces bits_per_residue

    # We leave gnomad_pli missing for testing

    res = fit_association_models(df)

    # Assert keys / columns
    expected_cols = [
        'annotation', 'n', 'coefficient', 'ols_p', 'cluster_robust_p',
        'delta_r2', 'spearman_corr', 'spearman_p', 'adjusted_p'
    ]
    for c in expected_cols:
        assert c in res.columns

    # Check expected coefficient signs
    disorder_row = res[res['annotation'] == 'disorder_fraction'].iloc[0]
    assert disorder_row['coefficient'] < 0
    assert disorder_row['n'] == 100

    loeuf_row = res[res['annotation'] == 'gnomad_loeuf'].iloc[0]
    assert loeuf_row['coefficient'] > 0
    assert loeuf_row['n'] == 100

    depmap_row = res[res['annotation'] == 'depmap_gene_effect'].iloc[0]
    assert depmap_row['coefficient'] > 0
    assert depmap_row['n'] == 100

    # Delta R2 should be positive
    assert disorder_row['delta_r2'] > 0
    assert loeuf_row['delta_r2'] > 0
    assert depmap_row['delta_r2'] > 0
