# protein-surprisal-atlas

Estimating protein language-model surprisal (pseudo-log-likelihood in bits) across the reviewed human proteome, and testing its biological signal beyond sequence length and homolog abundance.

## Terminology and Math Definitions
The core quantity evaluated in this project is **model surprisal** or **pseudo-log-likelihood (PLL) in bits**.

For a sequence $x$ of length $L$ and a masked language model $p_\theta$, the model surprisal is calculated as the sum of the negative log-probabilities of each true residue when masked:
$$ S(x) = \sum_{i=1}^{L} -\log_2 p_\theta(x_i \mid x_{\setminus i}) $$

This value is relative to baselines such as uniform probability (where each of the 20 amino acids is equally likely) and order-0 background probability (based on human proteome amino acid frequencies).

## Scientific Warning

**This project estimates model-relative sequence surprisal (pseudo-log-likelihood in bits) under a trained protein language model. It does not calculate a true Shannon code length or achievable compression size, does not calculate Kolmogorov complexity, does not prove protein functionality, and does not provide clinical interpretation of variants. Surprisal is strongly influenced by homolog abundance in the model's training data.**

The quantity `S(x)` is a pseudo-log-likelihood: the per-position conditionals do not multiply into a proper joint distribution and `S(x)` is not an achievable code length. Treat it as a model surprisal score, useful primarily in relative and rank terms.

## Installation
Dependencies are managed via conda and pip. The project supports both CPU and GPU execution.

```bash
# Create the conda environment and install dependencies
make env
conda activate protein-surprisal
```

## Dataset Acquisition
The pipeline starts by downloading the reviewed human proteome directly from UniProt via their REST API.

```bash
make download
```

## Execution (Pilot vs Full)
The pipeline is designed to be configurable via `config.yaml` for either a small pilot or a full run. A pilot run (`target_size` < total proteins) samples proteins across sequence length distributions to give representative results quickly.

To build the pilot set and extract Uniref50 cluster sizes for homolog abundance:
```bash
make build_pilot
```

To run the scoring approximation (Sampled Masking) on the dataset:
```bash
make score_pilot
```

To run exact scoring on a small validation subset and validate the approximation metrics:
```bash
make validate_scoring
```

Finally, to fit statistical models and emit plots and tables:
```bash
make run_analysis
```

## Expected Outputs
All outputs are directed to the `results/` folder:
- **`results/tables/`**: Contains raw Parquet files with protein and residue-level surprisals, top/bottom ranked proteins, and control sequence comparisons.
- **`results/figures/`**: Contains distribution histograms, regression plots of length and family size versus surprisal, and specific residue profile plots for validation genes (e.g., TP53, EGFR).
- **`results/reports/`**: Contains statistical outputs including multiple-testing corrected (FDR) P-values and $R^2$ values for the models.

## Reproducibility
Random seeds used for sequence sampling and control sequence generation (shuffling, random uniform, and random background) are drawn deterministically based on `config.project.seed`, the sequence's UniProt ID, and the replicate number, guaranteeing deterministic outputs for controls.

## Limitations
- **Confounding Variables:** `log_family_size` (homolog abundance) is a major confounder for model surprisal.
- **Not True Joint Probability:** Because we use a Masked Language Model, $S(x)$ is a pseudo-log-likelihood and not a true factorized joint probability.

## Citation
If you use this pipeline or data, please cite the project repository.
