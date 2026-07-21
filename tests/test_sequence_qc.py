import pytest
import pandas as pd
from src.protein_atlas.sequence_qc import normalize_sequence, is_valid_sequence, run_qc

def test_normalize_sequence():
    assert normalize_sequence(" aBc\nDef\t") == "ABCDEF"
    assert normalize_sequence(None) == ""

def test_is_valid_sequence():
    assert is_valid_sequence("ACDEFGHIKLMNPQRSTVWY") is True
    assert is_valid_sequence("ACDX") is False
    assert is_valid_sequence("") is False

def test_run_qc():
    df = pd.DataFrame({
        "Entry": ["P1", "P2", "P3", "P4", "P1"],
        "Sequence": ["ACD", "  a c ", "ACDX", "", "ACD"],
    })
    config = {"sequence": {"ambiguous_residue_policy": "reject"}}
    qc_df, report = run_qc(df, config)

    assert report["total_downloaded"] == 5
    assert report["exclusion_reasons"]["empty_sequences"] == 1
    assert report["exclusion_reasons"]["rejected_ambiguous"] == 1
    assert report["exclusion_reasons"]["duplicate_accessions"] == 1
    assert report["retained"] == 2
    assert set(qc_df["Entry"]) == {"P1", "P2"}
