import pytest
import pandas as pd
from src.protein_atlas.select_pilot import select_pilot_set

def test_deterministic_pilot():
    df = pd.DataFrame({
        "uniprot_id": [f"P{i}" for i in range(100)],
        "gene_symbol": [f"G{i}" for i in range(100)],
        "sequence_length": list(range(100, 200))
    })

    pilot1 = select_pilot_set(df, target_size=10, seed=42)
    pilot2 = select_pilot_set(df, target_size=10, seed=42)

    assert pilot1.equals(pilot2)

    pilot3 = select_pilot_set(df, target_size=10, seed=43)
    assert not pilot1.equals(pilot3)
