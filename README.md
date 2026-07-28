# protein-surprisal-atlas

Estimating **protein language-model surprisal** — how predictable each residue of a protein is to a
pretrained protein language model (ESM-2) — across the full reviewed human proteome, and testing what
biological properties that predictability reflects.

[![Open pilot in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/glenritschel/protein-surprisal-atlas/blob/main/notebooks/01_pilot_gpu.ipynb)

## What this is

For every residue of a protein we compute ESM-2's **masked-marginal surprisal** — the number of bits
needed to encode the observed amino acid given its context, `−log₂ P(xᵢ | x₋ᵢ)`. Summed and normalized
by length, this "bits per residue" measures how *surprising*, or how *unpredictable*, a sequence is to
the model. The project scores this across the reviewed human proteome (~20,000 proteins) and asks a
simple question: after accounting for the obvious nuisances (protein length, how many homologs a
protein has, and low-complexity content), does per-residue surprisal carry biological signal?

> **Important — what this quantity is *not*.** The masked-marginal sum is a *pseudo-log-likelihood*,
> not a normalized code length. It is **not** an achievable compression size, and it is **not**
> Kolmogorov complexity (which is uncomputable). Treat it as a *model-relative surprisal score*,
> meaningful in relative and rank terms. See [`docs/METHODS.md`](docs/METHODS.md).

## Key result

Across the full proteome (n ≈ 20,000), per-residue surprisal is **dominated by intrinsic disorder**,
with **weak but robust and mutually consistent** associations to gene-level constraint that only become
detectable at scale:

| Annotation | n | coef | ΔR² | Spearman ρ |
|---|---|---|---|---|
| Intrinsic disorder (MobiDB) | 20,361 | **+0.59** | **4.6%** | +0.22 |
| gnomAD LOEUF | 17,798 | −0.12 | 1.1% | −0.09 |
| DepMap CRISPR essentiality | 17,834 | −0.08 | 0.4% | −0.06 |
| gnomAD pLI | 17,798 | +0.05 | 0.2% | +0.05 |

Model: `bits_per_residue ~ length + log_family_size + low_complexity_fraction + annotation`, with
UniRef50 cluster-robust standard errors and Benjamini–Hochberg FDR. Disorder's effect is ~4–25× larger
than any constraint effect; all three constraint measures agree in direction (more-constrained genes
are marginally *less* predictable). Notably, these constraint associations are **null in a 500-protein
pilot** and appear only at proteome scale — a reminder that at large n, significance must be read
through effect size, not the p-value. Full interpretation and the pilot comparison are in
[`docs/RESULTS.md`](docs/RESULTS.md).

## Companion analysis: delta-surprisal for variant effect

The result above concerns *wild-type* per-residue surprisal. A separate question is whether the
*mutation-effect* form — **delta-surprisal**, `surprisal(mut) − surprisal(wt)` at a masked position —
carries signal the wild-type quantity does not. It does, with a qualified, honest answer:

- **It predicts pathogenicity, but only at scale.** On 2,145 balanced ClinVar missense variants across
  eight cancer-predisposition genes, delta-surprisal separates pathogenic from benign at
  **AUROC 0.889** with ESM-2 650M, far above a −BLOSUM62 baseline (0.694). With the 35M model used for
  the main atlas it is no better than that baseline (pooled 0.70; TP53 0.48). The small-model result is
  a false negative.
- **It is dominated by AlphaMissense.** On the identical variants AlphaMissense scores **0.955** and
  beats delta-surprisal on every gene; a rank-average of the two (0.940) falls *below* AlphaMissense
  alone, so delta-surprisal adds no orthogonal signal.

Delta-surprisal reproduces known ESM zero-shot variant-effect performance and confirms the signal is
real, but it is not a new or competitive predictor. Full method, per-gene table, and the reproducible
script are in [`variant_effect/`](variant_effect/).

## Quickstart

The pipeline runs as three Colab notebooks (see [`notebooks/README.md`](notebooks/README.md)):

1. **`01_pilot_gpu.ipynb`** — 500-protein pilot end-to-end: scoring, exact-vs-approximate validation
   (the fast estimator matches exact masked-marginal scoring at Spearman 0.94), controls, analysis.
2. **`02_proteome_scoring_gpu.ipynb`** — score the full reviewed proteome (partitioned, resumable).
3. **`03_proteome_annotations_cpu.ipynb`** — annotate the scored proteome and run the associations.

Or run the modules directly (see `docs/METHODS.md` for the CLI). Python 3.11; CPU or CUDA.

## Repository layout

```
src/protein_atlas/   core modules (scoring, windowing, baselines, annotations, statistics)
scripts/             CLI entry points (download, score, integrate annotations, associations)
notebooks/           Colab notebooks for the full pipeline
docs/                METHODS.md (methodology) and RESULTS.md (findings + limitations)
tests/               unit + integration tests
results/             figures / reports / summary tables (large data tables are gitignored)
variant_effect/      delta-surprisal variant-effect companion analysis (script + results)
```

## Data sources

UniProtKB/Swiss-Prot (sequences), UniRef (homolog/family size), MobiDB (intrinsic disorder),
gnomAD v4.1 (loss-of-function constraint), DepMap Public 24Q4 (CRISPR gene essentiality). Download
URLs, versions, and checksums are recorded in the pipeline's provenance. See `docs/METHODS.md`.

## Limitations (brief)

The proteome atlas uses a single model (ESM-2 35M); the `variant_effect/` companion additionally uses
650M. The disorder association is partly convergent validity (disorder predictors and protein language
models both learn from evolutionary conservation); associations are cross-sectional; annotation
coverage is partial (~87–100%). Full list in `docs/RESULTS.md`.

## Citation

See [`CITATION.cff`](CITATION.cff).

## License

Code is licensed under the MIT License (see [`LICENSE`](LICENSE)). Data and figures derived from
UniProt and the annotation sources above are made available under CC BY 4.0.
