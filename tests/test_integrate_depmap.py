import pandas as pd
import numpy as np
from src.protein_atlas.integrate_depmap import parse_depmap_label, join_depmap_to_pilot

def test_parse_depmap_label():
    assert parse_depmap_label("TP53 (7157)") == "TP53"
    assert parse_depmap_label("EGFR (1956)") == "EGFR"
    assert parse_depmap_label("MYC") == "MYC"

def test_join_depmap_to_pilot():
    pilot_df = pd.DataFrame({
        'gene_symbol': ["TP53", "EGFR", "SMIM26 LINC00493", "UNKNOWN_GENE"],
        'sequence': ["M...", "M...", "M...", "M..."]
    })

    depmap_df = pd.DataFrame({
        'gene_symbol': ["TP53", "EGFR", "SMIM26"],
        'depmap_gene_effect': [-1.5, -0.8, -0.2],
        'depmap_common_essential': [True, False, False]
    })

    merged = join_depmap_to_pilot(pilot_df, depmap_df)

    assert 'depmap_gene_effect' in merged.columns
    assert 'depmap_common_essential' in merged.columns

    assert len(merged) == 4

    tp53_row = merged[merged['gene_symbol'] == 'TP53'].iloc[0]
    assert tp53_row['depmap_gene_effect'] == -1.5
    assert tp53_row['depmap_common_essential'] == True

    smim_row = merged[merged['gene_symbol'] == 'SMIM26 LINC00493'].iloc[0]
    assert smim_row['depmap_gene_effect'] == -0.2
    assert smim_row['depmap_common_essential'] == False

    unknown_row = merged[merged['gene_symbol'] == 'UNKNOWN_GENE'].iloc[0]
    assert pd.isna(unknown_row['depmap_gene_effect'])
    assert pd.isna(unknown_row['depmap_common_essential'])
