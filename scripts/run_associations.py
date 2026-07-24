import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from src.protein_atlas.associations import fit_association_models

import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scope", choices=["pilot", "proteome"], default="pilot", help="Scope of the associations")
    args = parser.parse_args()

    if args.scope == "pilot":
        input_file = "results/tables/pilot_annotated.parquet"
        out_table = "results/tables/association_results.csv"
        out_md = "results/reports/associations.md"
        out_fig_prefix = "results/figures/fig"
    else:
        input_file = "results/tables/proteome_annotated.parquet"
        out_table = "results/tables/proteome_association_results.csv"
        out_md = "results/reports/proteome_associations.md"
        out_fig_prefix = "results/figures/proteome_fig"

    if not os.path.exists(input_file):
        raise FileNotFoundError(f"{input_file} not found.")

    print(f"Loading {input_file}...")
    df = pd.read_parquet(input_file)

    # 1. Pre-calculate length_and_family_adjusted_surprisal BEFORE passing to fit_association_models
    if 'length_and_family_adjusted_surprisal' not in df.columns:
        print("Calculating 'length_and_family_adjusted_surprisal' as residual of bits_per_residue ~ sequence_length + log_family_size")
        import statsmodels.formula.api as smf
        # We need a complete case subset for the residual
        subset_res = df[['bits_per_residue', 'sequence_length', 'log_family_size']].dropna()
        if not subset_res.empty:
            model = smf.ols("bits_per_residue ~ sequence_length + log_family_size", data=subset_res).fit()
            df['length_and_family_adjusted_surprisal'] = pd.Series(model.resid, index=subset_res.index)
        else:
            df['length_and_family_adjusted_surprisal'] = float('nan')

    # 2. Run associations
    print("Fitting association models...")
    res_df = fit_association_models(df)

    # Save results table
    os.makedirs(os.path.dirname(out_table), exist_ok=True)
    res_df.to_csv(out_table, index=False)
    print(f"Saved association table to {out_table}")

    # Write reports
    os.makedirs(os.path.dirname(out_md), exist_ok=True)
    with open(out_md, 'w') as f:
        f.write("# Biological Association Analyses\n\n")
        f.write("## Overview\n")
        f.write("This report presents the association between biological constraints (gnomAD metrics and intrinsic disorder) and protein-language-model surprisal (bits per residue), after adjusting for sequence length, homolog abundance (log_family_size), and low-complexity fraction.\n\n")

        f.write("## Interpretation Caveats\n")
        f.write("- **Coverage Limitations**: Biological annotations have partial coverage across the proteome. The sample size ($n$) varies across tests, and associations are evaluated only on complete-case subsets.\n")
        f.write("- **Collinearity**: Disorder and low-complexity partly overlap (disorder regions are often low-complexity), so interpret their coefficients jointly, not in isolation.\n")
        f.write("- **Causality**: Associations are cross-sectional and do not establish biological mechanism or causality.\n\n")

        f.write("## Results\n\n")

        f.write("| Annotation | n | Effect Size (Coef) | Delta R² | Robust p-value | Adjusted p-value (FDR) | Spearman Rank Corr |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for _, row in res_df.iterrows():
            f.write(f"| {row['annotation']} | {row['n']} | {row['coefficient']:.4f} | {row['delta_r2']:.4f} | {row['cluster_robust_p']:.2e} | {row['adjusted_p']:.2e} | {row['spearman_corr']:.4f} |\n")
    print(f"Saved report to {out_md}")

    # Generate Figures
    print("Generating figures...")
    os.makedirs("results/figures", exist_ok=True)

    # fig8: residual vs gnomAD LOEUF
    if 'gnomad_loeuf' in df.columns and df['gnomad_loeuf'].notna().any():
        plt.figure(figsize=(6, 5))
        sns.regplot(data=df, x='gnomad_loeuf', y='length_and_family_adjusted_surprisal',
                    scatter_kws={'alpha':0.5, 'color':'#2c3e50'},
                    line_kws={'color':'#e74c3c'})
        plt.title('Adjusted Surprisal vs. gnomAD LOEUF')
        plt.xlabel('gnomAD LOEUF')
        plt.ylabel('Residual (length+family adjusted bits/res)')
        plt.tight_layout()
        plt.savefig(f'{out_fig_prefix}8_residual_vs_gnomad_loeuf.png', dpi=300)
        plt.savefig(f'{out_fig_prefix}8_residual_vs_gnomad_loeuf.pdf')
        plt.close()

    # fig9: residual vs disorder_fraction
    if 'disorder_fraction' in df.columns and df['disorder_fraction'].notna().any():
        plt.figure(figsize=(6, 5))
        sns.regplot(data=df, x='disorder_fraction', y='length_and_family_adjusted_surprisal',
                    scatter_kws={'alpha':0.5, 'color':'#2980b9'},
                    line_kws={'color':'#e67e22'})
        plt.title('Adjusted Surprisal vs. Intrinsic Disorder Fraction')
        plt.xlabel('Disorder Fraction')
        plt.ylabel('Residual (length+family adjusted bits/res)')
        plt.tight_layout()
        plt.savefig(f'{out_fig_prefix}9_residual_vs_disorder_fraction.png', dpi=300)
        plt.savefig(f'{out_fig_prefix}9_residual_vs_disorder_fraction.pdf')
        plt.close()

    # fig10: correlation matrix
    cols_for_corr = ['bits_per_residue', 'length_and_family_adjusted_surprisal',
                     'log_family_size', 'low_complexity_fraction',
                     'gnomad_loeuf', 'gnomad_pli', 'disorder_fraction']

    cols_for_corr = [c for c in cols_for_corr if c in df.columns]

    if len(cols_for_corr) > 1:
        corr = df[cols_for_corr].corr(method='spearman')
        plt.figure(figsize=(8, 6))
        sns.heatmap(corr, annot=True, cmap='coolwarm', fmt=".2f", vmin=-1, vmax=1, center=0)
        plt.title('Spearman Correlation Matrix')
        plt.tight_layout()
        plt.savefig(f'{out_fig_prefix}10_correlation_matrix.png', dpi=300)
        plt.savefig(f'{out_fig_prefix}10_correlation_matrix.pdf')
        plt.close()

    print("Done.")

if __name__ == "__main__":
    main()
