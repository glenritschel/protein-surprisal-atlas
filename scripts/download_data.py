import argparse
import yaml
import json
import logging
from pathlib import Path
from src.protein_atlas.download_uniprot import download_human_proteome
from src.protein_atlas.sequence_qc import run_qc

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    raw_dir = "data/raw"
    interim_dir = "data/interim"
    Path(interim_dir).mkdir(parents=True, exist_ok=True)

    logging.info("Downloading data from UniProt...")
    df, metadata = download_human_proteome(raw_dir)

    logging.info("Running QC...")
    qc_df, qc_report = run_qc(df, config)

    with open(f"{raw_dir}/download_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    with open(f"{interim_dir}/qc_report.json", "w") as f:
        json.dump(qc_report, f, indent=2)

    output_format = config.get("project", {}).get("output_format", "parquet")
    output_path = f"{interim_dir}/human_proteome_qc.{output_format}"

    rename_cols = {
        "Entry": "uniprot_id",
        "Gene Names": "gene_symbol",
        "Protein names": "protein_name",
        "Sequence length": "sequence_length",
        "Sequence": "sequence"
    }
    qc_df = qc_df.rename(columns=rename_cols)

    if output_format == "parquet":
        qc_df.to_parquet(output_path, index=False)
    else:
        qc_df.to_csv(output_path, index=False)

    logging.info(f"Data saved to {output_path}")

if __name__ == "__main__":
    main()
