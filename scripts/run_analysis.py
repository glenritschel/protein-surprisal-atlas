import argparse
import yaml
import logging
import pandas as pd
import json
from pathlib import Path

from src.protein_atlas.statistics import fit_distribution_models
from src.protein_atlas.plots import (
    plot_surprisal_distribution,
    plot_length_vs_surprisal,
    plot_background_vs_model,
    plot_surprisal_ratio,
    plot_family_size_confound,
    plot_residue_profile,
    plot_controls
)
from src.protein_atlas.score_controls import shuffle_sequence, generate_random_uniform, generate_random_background, reverse_sequence
from src.protein_atlas.variant_scoring import score_variants_batch
from src.protein_atlas.esm_model import ESM2Scorer

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def get_controls_scores(df: pd.DataFrame, scorer: ESM2Scorer) -> pd.DataFrame:
    """Generate and score controls for a subset of proteins (e.g. top 10) to save time."""
    sub_df = df.head(10).copy()
    results = []

    import torch
    from src.protein_atlas.score_approximate import score_sequence_approximate

    for _, row in sub_df.iterrows():
        seq = row['sequence']
        uid = row['uniprot_id']

        # Natural
        nat_res = score_sequence_approximate(scorer, seq)
        results.append({'uniprot_id': uid, 'type': 'natural', 'total_surprisal_bits': nat_res['total_surprisal_bits']})

        # Shuffled
        shuf = shuffle_sequence(seq)
        shuf_res = score_sequence_approximate(scorer, shuf)
        results.append({'uniprot_id': uid, 'type': 'shuffled', 'total_surprisal_bits': shuf_res['total_surprisal_bits']})

        # Reversed
        rev = reverse_sequence(seq)
        rev_res = score_sequence_approximate(scorer, rev)
        results.append({'uniprot_id': uid, 'type': 'reversed', 'total_surprisal_bits': rev_res['total_surprisal_bits']})

        # Random Uniform
        rand_uni = generate_random_uniform(len(seq))
        rand_uni_res = score_sequence_approximate(scorer, rand_uni)
        results.append({'uniprot_id': uid, 'type': 'random_uniform', 'total_surprisal_bits': rand_uni_res['total_surprisal_bits']})

        # Random Background
        rand_bg = generate_random_background(len(seq))
        rand_bg_res = score_sequence_approximate(scorer, rand_bg)
        results.append({'uniprot_id': uid, 'type': 'random_background', 'total_surprisal_bits': rand_bg_res['total_surprisal_bits']})

    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    output_format = config.get("project", {}).get("output_format", "parquet")
    method = config.get("model", {}).get("scoring_method", "sampled_mask")

    tables_dir = "results/tables"
    fig_dir = "results/figures"
    rep_dir = "results/reports"

    Path(fig_dir).mkdir(parents=True, exist_ok=True)
    Path(rep_dir).mkdir(parents=True, exist_ok=True)

    protein_path = f"{tables_dir}/pilot_protein_scores_{method}.{output_format}"
    residue_path = f"{tables_dir}/pilot_residue_scores_{method}.{output_format}"

    if output_format == "parquet":
        df = pd.read_parquet(protein_path)
        res_df = pd.read_parquet(residue_path)
    else:
        df = pd.read_csv(protein_path)
        res_df = pd.read_csv(residue_path)

    # Statistical Models (FDR included in function now)
    df, stats_results = fit_distribution_models(df)

    with open(f"{rep_dir}/regression_results.json", "w") as f:
        json.dump(stats_results, f, indent=2)

    top_50 = df.nlargest(50, 'bits_per_residue')
    bottom_50 = df.nsmallest(50, 'bits_per_residue')
    top_pos_resid = df.nlargest(50, 'length_and_family_adjusted_surprisal')
    top_neg_resid = df.nsmallest(50, 'length_and_family_adjusted_surprisal')

    top_50.to_csv(f"{tables_dir}/top_50_normalized_surprisal.csv", index=False)
    bottom_50.to_csv(f"{tables_dir}/bottom_50_normalized_surprisal.csv", index=False)
    top_pos_resid.to_csv(f"{tables_dir}/top_50_positive_residuals.csv", index=False)
    top_neg_resid.to_csv(f"{tables_dir}/top_50_negative_residuals.csv", index=False)

    if output_format == "parquet":
        df.to_parquet(protein_path, index=False)
    else:
        df.to_csv(protein_path, index=False)

    plot_surprisal_distribution(df, fig_dir)
    plot_length_vs_surprisal(df, fig_dir)
    plot_background_vs_model(df, fig_dir)
    plot_surprisal_ratio(df, fig_dir)
    plot_family_size_confound(df, fig_dir)

    def find_exact(gene):
        for g in df['gene_symbol'].dropna().unique():
            if gene in g.split():
                return g
        return None

    tp53_gene = find_exact("TP53")
    if tp53_gene:
        plot_residue_profile(res_df, tp53_gene, fig_dir, 12)

    egfr_gene = find_exact("EGFR")
    if egfr_gene:
        plot_residue_profile(res_df, egfr_gene, fig_dir, 13)

    logging.info("Scoring controls & variants for a subset...")
    model_name = config.get("model", {}).get("name")
    device = config.get("model", {}).get("device")
    precision = config.get("model", {}).get("precision")

    scorer = ESM2Scorer(model_name=model_name, device=device, precision=precision)

    controls_df = get_controls_scores(df, scorer)
    controls_df.to_csv(f"{tables_dir}/control_sequence_comparisons.csv", index=False)
    plot_controls(controls_df, fig_dir)

    # Variant scoring placeholder (TP53 example)
    if tp53_gene:
        tp53_seq = df[df['gene_symbol'] == tp53_gene].iloc[0]['sequence']
        # Mute first residue to Ala if not already, or something
        ref_res = tp53_seq[0]
        alt_res = 'A' if ref_res != 'A' else 'G'
        mock_variants = [{'position': 1, 'reference': ref_res, 'alternate': alt_res}]
        var_results = score_variants_batch(scorer, tp53_seq, mock_variants)
        pd.DataFrame(var_results).to_csv(f"{tables_dir}/variant_scores_TP53.csv", index=False)

if __name__ == "__main__":
    main()
