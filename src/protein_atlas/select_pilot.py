import pandas as pd
import numpy as np

VALIDATION_GENES = ["TP53", "EGFR", "BRCA1", "CFTR", "COL1A1", "ACTB", "GAPDH", "INS", "HBB", "APOE"]

def select_pilot_set(df: pd.DataFrame, target_size: int = 500, seed: int = 42) -> pd.DataFrame:
    np.random.seed(seed)

    def is_validation(gene_str):
        if not isinstance(gene_str, str):
            return False
        parts = gene_str.split()
        for v in VALIDATION_GENES:
            if v in parts:
                return True
        return False

    df['is_validation'] = df['gene_symbol'].apply(is_validation)

    val_df = df[df['is_validation']].copy()
    val_df = val_df.drop_duplicates(subset=['gene_symbol']).head(len(VALIDATION_GENES) * 2)

    remaining_size = target_size - len(val_df)
    if remaining_size <= 0:
        return val_df.head(target_size)

    other_df = df[~df['is_validation']].copy()
    other_df['length_bin'] = pd.qcut(other_df['sequence_length'], q=5, labels=False, duplicates='drop')

    sampled_dfs = []
    bins = other_df['length_bin'].unique()
    per_bin = remaining_size // len(bins)

    for b in bins:
        bin_df = other_df[other_df['length_bin'] == b]
        sampled_dfs.append(bin_df.sample(n=min(len(bin_df), per_bin), random_state=seed))

    sampled_df = pd.concat(sampled_dfs)

    shortfall = remaining_size - len(sampled_df)
    if shortfall > 0:
        leftover = other_df[~other_df['uniprot_id'].isin(sampled_df['uniprot_id'])]
        sampled_df = pd.concat([sampled_df, leftover.sample(n=shortfall, random_state=seed)])

    final_df = pd.concat([val_df, sampled_df]).reset_index(drop=True)
    return final_df
