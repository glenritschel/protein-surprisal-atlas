import pytest
import os
import shutil
import pandas as pd
import json
import yaml
from pathlib import Path
import subprocess

def test_score_proteome_scaling(tmp_path):
    # Setup mock config
    config = {
        "project": {"output_format": "parquet"},
        "model": {
            "name": "facebook/esm2_t12_35M_UR50D",
            "device": "cpu",
            "precision": "float32",
            "scoring_method": "naive_onepass"
        },
        "scoring": {},
        "sequence": {"maximum_length": 100},
        "baselines": {"use_order0": True}
    }
    config_file = tmp_path / "config.yaml"
    with open(config_file, "w") as f:
        yaml.dump(config, f)

    # Setup mock input
    interim_dir = tmp_path / "data" / "interim"
    interim_dir.mkdir(parents=True)
    mock_df = pd.DataFrame({
        "uniprot_id": [f"P{i:04d}" for i in range(25)],
        "gene_symbol": [f"GENE{i}" for i in range(25)],
        "sequence": ["M" + "A" * 9 for _ in range(25)],
        "sequence_length": [10 for _ in range(25)]
    })
    mock_df.to_parquet(interim_dir / "human_proteome_qc.parquet")

    # We will patch the paths in score_proteome.py by running it with a custom python script
    # that overrides the paths, or we can just run it in the tmp_path as CWD.

    script = """
import sys
import os
from unittest.mock import patch

# Mock out interim dir
def mock_interim():
    return "data/interim"

if __name__ == '__main__':
    from scripts.score_proteome import main
    main()
"""
    runner_script = tmp_path / "runner.py"
    with open(runner_script, "w") as f:
        f.write(script)

    # Create results dirs
    (tmp_path / "results" / "tables").mkdir(parents=True)
    (tmp_path / "results" / "logs").mkdir(parents=True)

    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).parent.parent.absolute())

    # Run 1: process first 15 records
    # We can use --limit 15
    cmd1 = ["python", str(Path(__file__).parent.parent / "scripts" / "score_proteome.py"), "--config", str(config_file), "--limit", "15"]
    subprocess.run(cmd1, cwd=tmp_path, env=env, check=True)

    protein_out_dir = tmp_path / "results" / "tables" / "proteome_protein_scores_naive_onepass"
    residue_out_dir = tmp_path / "results" / "tables" / "proteome_residue_scores_naive_onepass"

    assert protein_out_dir.exists()
    assert residue_out_dir.exists()

    parts1 = list(protein_out_dir.glob("part-*.parquet"))
    # Should write parts of size 10, so 1 part of 10, 1 part of 5
    assert len(parts1) == 2

    # Combine outputs
    df1 = pd.concat([pd.read_parquet(p) for p in parts1])
    assert len(df1) == 15
    assert set(df1["uniprot_id"]) == set([f"P{i:04d}" for i in range(15)])

    # Run 2: process all 25 records. It should skip the first 15 and process the remaining 10.
    cmd2 = ["python", str(Path(__file__).parent.parent / "scripts" / "score_proteome.py"), "--config", str(config_file), "--limit", "25"]
    subprocess.run(cmd2, cwd=tmp_path, env=env, check=True)

    parts2 = list(protein_out_dir.glob("part-*.parquet"))
    # Should add 1 part of 10
    assert len(parts2) == 3

    # Combine outputs
    df2 = pd.concat([pd.read_parquet(p) for p in parts2])
    assert len(df2) == 25
    assert len(df2["uniprot_id"].unique()) == 25

    # Verify no duplicates
    assert len(df2) == len(df2["uniprot_id"].unique())
