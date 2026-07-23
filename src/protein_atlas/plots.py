import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import seaborn as sns
import re

def plot_surprisal_distribution(df: pd.DataFrame, output_dir: str):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(df['bits_per_residue'], bins=30, edgecolor='black', alpha=0.7)
    ax.set_xlabel("Bits per Residue")
    ax.set_ylabel("Number of Proteins")
    ax.set_title("Distribution of Bits per Residue across Pilot")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/fig1_bits_per_residue_hist.png", dpi=300)
    plt.savefig(f"{output_dir}/fig1_bits_per_residue_hist.pdf")
    plt.close()

def plot_length_vs_surprisal(df: pd.DataFrame, output_dir: str):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df['sequence_length'], df['total_surprisal_bits'], alpha=0.6, s=20)
    ax.set_xlabel("Sequence Length")
    ax.set_ylabel("Total Surprisal (Bits)")
    ax.set_title("Protein Length vs Total Surprisal")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/fig2_length_vs_surprisal.png", dpi=300)
    plt.savefig(f"{output_dir}/fig2_length_vs_surprisal.pdf")
    plt.close()

def plot_background_vs_model(df: pd.DataFrame, output_dir: str):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df['order0_baseline_bits'], df['total_surprisal_bits'], alpha=0.6, s=20)
    max_val = max(df['order0_baseline_bits'].max(), df['total_surprisal_bits'].max())
    ax.plot([0, max_val], [0, max_val], 'r--', label='y=x (No compression)')
    ax.set_xlabel("Order-0 Background Bits")
    ax.set_ylabel("Model Surprisal Bits")
    ax.set_title("Background Bits vs Model Surprisal")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"{output_dir}/fig3_background_vs_model.png", dpi=300)
    plt.savefig(f"{output_dir}/fig3_background_vs_model.pdf")
    plt.close()

def plot_surprisal_ratio(df: pd.DataFrame, output_dir: str):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(df['surprisal_ratio_order0'], bins=30, alpha=0.7, label='Order-0 Ratio', density=True)
    ax.hist(df['surprisal_ratio_uniform'], bins=30, alpha=0.5, label='Uniform Ratio', density=True)
    ax.set_xlabel("Surprisal Ratio")
    ax.set_ylabel("Density")
    ax.set_title("Distribution of Surprisal Ratios")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"{output_dir}/fig4_surprisal_ratios.png", dpi=300)
    plt.savefig(f"{output_dir}/fig4_surprisal_ratios.pdf")
    plt.close()

def plot_family_size_confound(df: pd.DataFrame, output_dir: str):
    fig, ax = plt.subplots(figsize=(8, 6))
    if 'log_family_size' in df.columns:
        ax.scatter(df['log_family_size'], df['bits_per_residue'], alpha=0.6, s=20)
        ax.set_xlabel("Log Family Size (UniRef50)")
        ax.set_ylabel("Bits per Residue")
        ax.set_title("Homolog Abundance Confound")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/fig7_confound.png", dpi=300)
        plt.savefig(f"{output_dir}/fig7_confound.pdf")
    plt.close()

def plot_residual_vs_depmap_effect(df: pd.DataFrame, output_dir: str, fignum: int = 11):
    fig, ax = plt.subplots(figsize=(8, 6))
    if 'depmap_gene_effect' in df.columns and 'length_and_family_adjusted_surprisal' in df.columns:
        subset = df[['depmap_gene_effect', 'length_and_family_adjusted_surprisal']].dropna()
        ax.scatter(subset['depmap_gene_effect'], subset['length_and_family_adjusted_surprisal'], alpha=0.6, s=20)
        ax.set_xlabel("DepMap Gene Effect (Mean)")
        ax.set_ylabel("Residual Surprisal (Bits)")
        ax.set_title("Residual Surprisal vs DepMap Gene Effect")
        plt.tight_layout()
        plt.savefig(f"{output_dir}/fig{fignum}_residual_vs_depmap_effect.png", dpi=300)
        plt.savefig(f"{output_dir}/fig{fignum}_residual_vs_depmap_effect.pdf")
    plt.close()

def plot_residue_profile(res_df: pd.DataFrame, gene: str, output_dir: str, fignum: int):
    gene_df = res_df[res_df['gene_symbol'] == gene]
    if len(gene_df) == 0:
        return

    gene_df = gene_df.sort_values('position')

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(gene_df['position'], gene_df['surprisal_bits'], label='Residue Surprisal')
    mean_val = gene_df['surprisal_bits'].mean()
    ax.axhline(mean_val, color='r', linestyle='--', label=f'Mean = {mean_val:.2f}')

    ax.set_xlabel("Position")
    ax.set_ylabel("Surprisal (bits)")
    ax.set_title(f"Residue-Level Surprisal Profile: {gene}")
    ax.legend()
    plt.tight_layout()

    primary = str(gene).split()[0]
    safe = re.sub(r'[^A-Za-z0-9_.-]', '_', primary)

    plt.savefig(f"{output_dir}/fig{fignum}_{safe}_profile.png", dpi=300)
    plt.savefig(f"{output_dir}/fig{fignum}_{safe}_profile.pdf")
    plt.close()

def plot_controls(controls_df: pd.DataFrame, output_dir: str):
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.boxplot(data=controls_df, x='type', y='total_surprisal_bits', ax=ax)
    ax.set_title("Control Sequence Comparison")
    ax.set_xlabel("Sequence Type")
    ax.set_ylabel("Total Surprisal (Bits)")
    plt.tight_layout()
    plt.savefig(f"{output_dir}/fig5_controls.png", dpi=300)
    plt.savefig(f"{output_dir}/fig5_controls.pdf")
    plt.close()
