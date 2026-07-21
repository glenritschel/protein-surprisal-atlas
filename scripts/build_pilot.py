import argparse
import yaml
import logging
import pandas as pd
from pathlib import Path

from src.protein_atlas.select_pilot import select_pilot_set
from src.protein_atlas.family_size import add_family_sizes

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    interim_dir = "data/interim"
    output_format = config.get("project", {}).get("output_format", "parquet")
    input_path = f"{interim_dir}/human_proteome_qc.{output_format}"

    logging.info(f"Loading QC'd data from {input_path}")
    if output_format == "parquet":
        df = pd.read_parquet(input_path)
    else:
        df = pd.read_csv(input_path)

    target_size = config.get("pilot", {}).get("target_size", 500)
    seed = config.get("project", {}).get("seed", 42)

    logging.info(f"Selecting {target_size} pilot proteins (seed={seed})...")
    pilot_df = select_pilot_set(df, target_size=target_size, seed=seed)

    logging.info("Retrieving UniRef cluster sizes...")
    pilot_df = add_family_sizes(pilot_df)

    output_path = f"{interim_dir}/pilot_set.{output_format}"
    if output_format == "parquet":
        pilot_df.to_parquet(output_path, index=False)
    else:
        pilot_df.to_csv(output_path, index=False)

    logging.info(f"Pilot set saved to {output_path}")

if __name__ == "__main__":
    main()
