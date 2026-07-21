import argparse
import yaml
import logging
import pandas as pd
import os
import json
from pathlib import Path
from tqdm import tqdm

from src.protein_atlas.esm_model import ESM2Scorer
from src.protein_atlas.score_exact import score_sequence_exact
from src.protein_atlas.score_approximate import score_sequence_approximate, score_sequence_naive
from src.protein_atlas.score_windows import score_sequence_windowed
from src.protein_atlas.baselines import order0_baseline, uniform_baseline

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def score_protein(row, scorer, config, method, validation_only):
    seq = row['sequence']
    seq_len = row['sequence_length']

    max_len = config.get("sequence", {}).get("maximum_length", 900)
    overlap = config.get("sequence", {}).get("window_overlap", 100)
    merge_rule = config.get("sequence", {}).get("window_merge_rule", "central")

    kwargs = {}
    if method == "exact":
        score_fn = score_sequence_exact
        kwargs["batch_size"] = config.get("scoring", {}).get("exact_batch_size", 32)
    elif method == "sampled_mask":
        score_fn = score_sequence_approximate
        kwargs["mask_fraction"] = config.get("scoring", {}).get("sampled_mask_fraction", 0.15)
        kwargs["passes"] = config.get("scoring", {}).get("sampled_mask_passes", 7)
    elif method == "naive_onepass":
        score_fn = score_sequence_naive
    else:
        raise ValueError(f"Unknown scoring method: {method}")

    res = score_sequence_windowed(scorer, seq, score_fn, max_len=max_len, overlap=overlap, merge_rule=merge_rule, **kwargs)

    o0_bits = order0_baseline(seq)
    uni_bits = uniform_baseline(seq_len)

    res["uniform_baseline_bits"] = uni_bits
    res["order0_baseline_bits"] = o0_bits
    res["order1_baseline_bits"] = o0_bits

    res["surprisal_ratio_order0"] = res["total_surprisal_bits"] / o0_bits if o0_bits > 0 else 0
    res["surprisal_ratio_uniform"] = res["total_surprisal_bits"] / uni_bits if uni_bits > 0 else 0

    surps = pd.Series(res["residue_surprisals"])
    res["mean_residue_surprise"] = surps.mean()
    res["median_residue_surprise"] = surps.median()
    res["maximum_residue_surprise"] = surps.max()
    res["minimum_residue_surprise"] = surps.min()
    res["residue_surprise_std"] = surps.std()

    res["uniprot_id"] = row["uniprot_id"]
    res["gene_symbol"] = row["gene_symbol"]

    return res

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument("--method", default=None, help="Override config scoring_method")
    parser.add_argument("--validation-only", action="store_true", help="Only score the validation set")
    args = parser.parse_args()

    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    method = args.method or config.get("model", {}).get("scoring_method", "sampled_mask")

    interim_dir = "data/interim"
    output_format = config.get("project", {}).get("output_format", "parquet")
    input_path = f"{interim_dir}/pilot_set.{output_format}"

    if output_format == "parquet":
        df = pd.read_parquet(input_path)
    else:
        df = pd.read_csv(input_path)

    if args.validation_only:
        val_size = config.get("pilot", {}).get("exact_validation_size", 50)
        df = df.head(val_size)
        suffix = "validation"
    else:
        suffix = "pilot"

    model_name = config.get("model", {}).get("name")
    device = config.get("model", {}).get("device")
    precision = config.get("model", {}).get("precision")

    scorer = ESM2Scorer(model_name=model_name, device=device, precision=precision)

    output_dir = "results/tables"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    protein_out = f"{output_dir}/{suffix}_protein_scores_{method}.{output_format}"
    residue_out = f"{output_dir}/{suffix}_residue_scores_{method}.{output_format}"
    error_out = f"results/logs/scoring_errors_{suffix}_{method}.jsonl"
    Path("results/logs").mkdir(exist_ok=True, parents=True)

    protein_records = []
    residue_records = []

    completed_ids = set()
    if os.path.exists(protein_out):
        try:
            if output_format == "parquet":
                existing = pd.read_parquet(protein_out)
            else:
                existing = pd.read_csv(protein_out)
            completed_ids = set(existing["uniprot_id"].unique())
            protein_records = existing.to_dict('records')
            logging.info(f"Resuming. Found {len(completed_ids)} completed proteins.")

            if os.path.exists(residue_out):
                if output_format == "parquet":
                    ex_res = pd.read_parquet(residue_out)
                else:
                    ex_res = pd.read_csv(residue_out)
                residue_records = ex_res.to_dict('records')
        except Exception as e:
            logging.warning(f"Failed to read checkpoints: {e}. Starting fresh.")
            completed_ids = set()
            protein_records = []
            residue_records = []

    for _, row in tqdm(df.iterrows(), total=len(df)):
        uid = row["uniprot_id"]
        if uid in completed_ids:
            continue

        try:
            res = score_protein(row, scorer, config, method, args.validation_only)

            p_record = row.to_dict()
            for k, v in res.items():
                if k not in ["residue_surprisals", "residue_probabilities"]:
                    p_record[k] = v
            protein_records.append(p_record)

            for i in range(row["sequence_length"]):
                residue_records.append({
                    "uniprot_id": uid,
                    "gene_symbol": row["gene_symbol"],
                    "position": i + 1,
                    "reference_residue": row["sequence"][i],
                    "surprisal_bits": res["residue_surprisals"][i],
                    "probability": res["residue_probabilities"][i],
                    "window_id": i // config.get("sequence", {}).get("maximum_length", 900) if res["windowed"] else 0,
                    "distance_from_window_edge": -1
                })

            if len(protein_records) % 10 == 0:
                pd.DataFrame(protein_records).to_parquet(protein_out, index=False) if output_format=="parquet" else pd.DataFrame(protein_records).to_csv(protein_out, index=False)
                pd.DataFrame(residue_records).to_parquet(residue_out, index=False) if output_format=="parquet" else pd.DataFrame(residue_records).to_csv(residue_out, index=False)

        except Exception as e:
            logging.exception(f"Failed on {uid}: {e}")
            with open(error_out, "a") as ef:
                ef.write(json.dumps({"uniprot_id": uid, "error": str(e)}) + "\n")

    if protein_records:
        pd.DataFrame(protein_records).to_parquet(protein_out, index=False) if output_format=="parquet" else pd.DataFrame(protein_records).to_csv(protein_out, index=False)
        pd.DataFrame(residue_records).to_parquet(residue_out, index=False) if output_format=="parquet" else pd.DataFrame(residue_records).to_csv(residue_out, index=False)
        logging.info(f"Finished. Saved to {protein_out}")

if __name__ == "__main__":
    main()
