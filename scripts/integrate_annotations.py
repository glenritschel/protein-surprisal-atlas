import os
import yaml
import pandas as pd
from src.protein_atlas.integrate_gnomad import fetch_gnomad_metrics, process_gnomad_metrics, join_gnomad_to_pilot
from src.protein_atlas.integrate_depmap import fetch_depmap_data, process_depmap_data, join_depmap_to_pilot
from src.protein_atlas.annotate_disorder import get_disorder_fractions
from src.protein_atlas.annotate_lowcomplexity import compute_low_complexity_fraction

def main():
    with open('config.yaml', 'r') as f:
        config = yaml.safe_load(f)

    input_file = "results/tables/pilot_protein_scores_sampled_mask.parquet"
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"{input_file} not found.")

    print(f"Loading {input_file}...")
    df = pd.read_parquet(input_file)

    # 1. Low complexity
    print("Computing low complexity fraction...")
    window_size = config.get('annotations', {}).get('low_complexity_window_size', 20)
    entropy_threshold = config.get('annotations', {}).get('low_complexity_entropy_threshold', 3.0)

    df['low_complexity_fraction'] = df['sequence'].apply(
        lambda seq: compute_low_complexity_fraction(seq, window_size, entropy_threshold)
    )

    # 2. MobiDB Disorder
    print("Fetching disorder fraction...")
    df = get_disorder_fractions(df, max_workers=4)

    # 3. gnomAD LOEUF & pLI
    print("Integrating gnomAD constraint metrics...")
    gnomad_fp, sha256, version, url = fetch_gnomad_metrics()
    gnomad_ens_gene, gnomad_sym_gene, gnomad_transcript = process_gnomad_metrics(gnomad_fp)
    df = join_gnomad_to_pilot(df, gnomad_ens_gene, gnomad_sym_gene, gnomad_transcript)

    # 4. DepMap Gene Effect
    print("Integrating DepMap gene effect...")
    depmap_effect_path, depmap_effect_sha256, depmap_essential_path, depmap_essential_sha256, depmap_label, depmap_doi = fetch_depmap_data()
    depmap_df = process_depmap_data(depmap_effect_path, depmap_essential_path)
    df = join_depmap_to_pilot(df, depmap_df)

    # Write output
    output_file = "results/tables/pilot_annotated.parquet"
    print(f"Writing annotated dataset to {output_file}...")
    df.to_parquet(output_file)

    # Generate coverage report
    print("Generating coverage report...")
    os.makedirs("results/reports", exist_ok=True)
    report_file = "results/reports/annotation_coverage.md"
    n_total = len(df)

    lowcomp_present = df['low_complexity_fraction'].notna().sum()
    disorder_present = df['disorder_fraction'].notna().sum()
    loeuf_present = df['gnomad_loeuf'].notna().sum()
    pli_present = df['gnomad_pli'].notna().sum()
    depmap_present = df['depmap_gene_effect'].notna().sum()

    join_counts = df['gnomad_join_method'].value_counts().to_dict() if 'gnomad_join_method' in df.columns else {}

    with open(report_file, 'w') as f:
        f.write("# Biological Annotation Coverage Report\n\n")
        f.write(f"Total sequences: {n_total}\n\n")

        f.write("## Coverage\n")
        f.write(f"- **Low Complexity Fraction**: {lowcomp_present} ({lowcomp_present/n_total*100:.1f}%)\n")
        f.write(f"- **Disorder Fraction**: {disorder_present} ({disorder_present/n_total*100:.1f}%)\n")
        f.write(f"- **gnomAD LOEUF**: {loeuf_present} ({loeuf_present/n_total*100:.1f}%)\n")
        f.write(f"- **gnomAD pLI**: {pli_present} ({pli_present/n_total*100:.1f}%)\n")
        f.write(f"- **DepMap Gene Effect**: {depmap_present} ({depmap_present/n_total*100:.1f}%)\n\n")

        f.write("## gnomAD Join Methods\n")
        for method, count in join_counts.items():
            f.write(f"- **{method}**: {count}\n")

        f.write("\n## Provenance\n")
        f.write(f"- **gnomAD URL**: {url}\n")
        f.write(f"- **gnomAD Version**: {version}\n")
        f.write(f"- **gnomAD SHA256**: {sha256}\n")
        f.write(f"- **DepMap URL (Effect)**: {DEPMAP_RELEASE_URL_EFFECT if 'DEPMAP_RELEASE_URL_EFFECT' in globals() else 'https://ndownloader.figshare.com/files/51064667'}\n")
        f.write(f"- **DepMap URL (Essentials)**: {DEPMAP_RELEASE_URL_ESSENTIAL if 'DEPMAP_RELEASE_URL_ESSENTIAL' in globals() else 'https://ndownloader.figshare.com/files/51064916'}\n")
        f.write(f"- **DepMap Label**: {depmap_label}\n")
        f.write(f"- **DepMap DOI**: {depmap_doi}\n")
        f.write(f"- **DepMap SHA256 (Effect)**: {depmap_effect_sha256}\n")
        f.write(f"- **DepMap SHA256 (Essentials)**: {depmap_essential_sha256}\n")

    print("Integration complete.")

if __name__ == "__main__":
    main()
