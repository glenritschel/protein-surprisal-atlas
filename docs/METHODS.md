# Methods

## 1. The quantity: model surprisal (pseudo-log-likelihood in bits)

For a protein sequence `x = (x₁, …, x_L)` we compute the **masked-marginal surprisal** under a
pretrained masked protein language model (ESM-2):

```
S(x) = − Σᵢ log₂ P(xᵢ | x₋ᵢ)
```

where `P(xᵢ | x₋ᵢ)` is the model's probability for the observed residue at position `i` when that
position is replaced by the mask token. We report the length-normalized value, **bits per residue**,
`S(x) / L`.

**What this is not.** Because ESM-2 is a *masked* (not autoregressive) model, the per-position
conditionals do not multiply into a normalized joint distribution — `S(x)` is a *pseudo*-log-likelihood
(pseudo-perplexity), not a Shannon code length. It is therefore **not** an achievable compression size,
and it is **not** Kolmogorov complexity (which is uncomputable). It is a *model-relative surprisal
score*, meaningful in relative and rank terms. All interpretation in this project respects that
distinction.

## 2. Scoring

- **Exact masked-marginal** — one masked forward pass per residue (the `L` masked variants of a protein
  are batched). Used to validate the fast estimator and for the pilot validation subset.
- **Sampled masking (default, fast)** — positions are partitioned deterministically into `k` folds
  (default `k = 7`); each fold is masked in turn and the observed residues in that fold are scored.
  This costs ~`k` forward passes per protein instead of `L`. Deterministic (no RNG), so reproducible.
  Validated against exact scoring at **Spearman ρ = 0.94** (n = 50; see `RESULTS.md`).
- **Windowing** — proteins longer than the model context are scored in sliding windows (default length
  900, overlap 100); each residue's score is taken from the window in which it sits farthest from an
  edge. Long proteins are never silently truncated.

## 3. Baselines and controls

- **Baselines.** Surprisal is reported against an order-0 human amino-acid-frequency background and a
  uniform 20-letter baseline (`L·log₂20`). The masked-marginal sum is *not* a code length, so ratios
  against these baselines are descriptive, not achievable compression ratios.
- **Controls.** For pilot proteins we score composition-preserving shuffles, the reversed sequence, and
  length-matched random sequences (uniform and background frequencies). Natural sequences are markedly
  less surprising than any control (see `RESULTS.md`), confirming the metric captures real structure.

## 4. Annotations

Each protein is annotated with the following, each downloaded with recorded provenance
(URL, version, sha256) and joined non-destructively:

| Annotation | Source | Field / definition | Join key |
|---|---|---|---|
| `log_family_size` | UniRef | log of UniRef50 cluster member count | UniProt accession |
| `low_complexity_fraction` | (computed) | fraction of residues in ≥1 low-entropy window (20-aa window, 3.0-bit threshold) | — |
| `disorder_fraction` | MobiDB | `prediction-disorder-th_50` content fraction | UniProt accession |
| `gnomad_loeuf`, `gnomad_pli` | gnomAD v4.1 | LOEUF / pLI, MANE-select transcript | gene symbol → gene |
| `depmap_gene_effect` | DepMap 24Q4 | mean CRISPR gene effect across cell lines (more negative = more essential) | HUGO symbol |

Homolog abundance (family size) is treated as the **primary nuisance covariate**: protein language-model
surprisal is strongly influenced by how many homologs a protein has in the training data, so it must be
adjusted for before interpreting any biological signal.

## 5. Statistical model

For each biological annotation `A` we fit a single full-covariate ordinary least-squares model on the
complete-case subset:

```
bits_per_residue ~ sequence_length + log_family_size + low_complexity_fraction + A
```

and report the coefficient of `A`, its ordinary and **UniRef50 cluster-robust** p-values (proteins in
the same family are non-independent), the **incremental ΔR²** over the nuisance-only model, and the
Spearman correlation of the length+family-adjusted residual with `A`. Benjamini–Hochberg FDR is applied
across the annotation set. A separate diagnostic model on *total* surprisal is reported but not
interpreted: length trivially explains ~95% of total bits, which is why the primary analysis is on
bits-per-residue.

**At large n, report effect sizes.** At proteome scale (n ≈ 20,000) even a ~0.2% effect is highly
significant; the ΔR² and Spearman columns, not the p-values, carry the interpretation.

## 6. Reproducibility

Deterministic scoring (fixed seed; no RNG in sampled masking), recorded provenance and checksums for
every external dataset, resumable caches for the per-accession annotation fetches, and partitioned,
restart-safe outputs for the proteome run. Exact bitwise reproducibility across CPU/GPU is not expected
(floating-point non-associativity); reproducibility is asserted within numerical tolerance.
