import hashlib
import pandas as pd
from typing import Tuple

VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

def normalize_sequence(seq: str) -> str:
    if not isinstance(seq, str):
        return ""
    return "".join(seq.upper().split())

def is_valid_sequence(seq: str) -> bool:
    if not seq:
        return False
    return set(seq).issubset(VALID_AMINO_ACIDS)

def get_sequence_checksum(seq: str) -> str:
    return hashlib.sha256(seq.encode("utf-8")).hexdigest()

def count_ambiguous(seq: str) -> int:
    return len([c for c in seq if c not in VALID_AMINO_ACIDS])

def run_qc(df: pd.DataFrame, config: dict) -> Tuple[pd.DataFrame, dict]:
    total_downloaded = len(df)

    if 'Sequence' not in df.columns:
        raise ValueError("Sequence column missing from dataframe.")

    df['Sequence'] = df['Sequence'].apply(normalize_sequence)
    df['Sequence length'] = df['Sequence'].apply(len)

    mask_empty = df['Sequence length'] == 0
    empty_count = mask_empty.sum()
    df = df[~mask_empty]

    df['ambiguous_residue_count'] = df['Sequence'].apply(count_ambiguous)

    ambig_policy = config.get("sequence", {}).get("ambiguous_residue_policy", "flag")
    rejected_ambig_count = 0
    if ambig_policy == "reject":
        mask_valid = df['ambiguous_residue_count'] == 0
        rejected_ambig_count = (~mask_valid).sum()
        df = df[mask_valid]

    df['sequence_sha256'] = df['Sequence'].apply(get_sequence_checksum)

    initial_len = len(df)
    df = df.drop_duplicates(subset=["Entry"])
    dupe_accessions = initial_len - len(df)

    total_retained = len(df)

    qc_report = {
        "total_downloaded": int(total_downloaded),
        "retained": int(total_retained),
        "excluded": int(total_downloaded - total_retained),
        "exclusion_reasons": {
            "empty_sequences": int(empty_count),
            "rejected_ambiguous": int(rejected_ambig_count),
            "duplicate_accessions": int(dupe_accessions)
        },
        "length_distribution": {
            "min": int(df['Sequence length'].min()) if total_retained > 0 else 0,
            "max": int(df['Sequence length'].max()) if total_retained > 0 else 0,
            "mean": float(df['Sequence length'].mean()) if total_retained > 0 else 0,
            "median": float(df['Sequence length'].median()) if total_retained > 0 else 0,
        },
        "ambiguous_residue_total_count": int(df['ambiguous_residue_count'].sum())
    }

    return df, qc_report
