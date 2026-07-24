# Results

## 1. The fast estimator is faithful

Sampled-masking vs exact masked-marginal scoring, on a 50-protein validation set spanning length,
membrane, disorder, enzyme, and disease categories:

- Spearman ρ = **0.94**, Pearson r = 0.90
- median absolute error 0.024 bits/residue, relative error ~2.4%
- small positive bias (+0.077 bits/residue)

Rank agreement is more than sufficient to score the full proteome with the fast estimator rather than
exact scoring.

## 2. Homolog abundance is the primary nuisance, and it is modest

Per-residue surprisal is essentially uncorrelated with protein length (R² ≈ 0.001) — confirming that
per-residue normalization is the right frame. Adding homolog abundance (`log_family_size`) gives a
significant negative coefficient (more homologs → more predictable), as expected from training-data
redundancy, but it explains only ~3% of per-residue variance. It is a genuine confounder that must be
adjusted for, but it does not dominate — leaving substantial variance for biology.

## 3. Controls: natural sequences are more predictable than any scramble

Mean total surprisal rises monotonically across control types:

```
natural  <  reversed  <  shuffled  <  random (background freq.)  <  random (uniform)
```

Real sequences are markedly less surprising than composition-matched shuffles or random sequences, and
uniform-random (ignoring amino-acid frequencies) is the most surprising — confirming the metric
captures learned sequence structure.

## 4. Proteome-scale associations (n ≈ 20,000) — the primary result

Model: `bits_per_residue ~ length + log_family_size + low_complexity_fraction + A`, cluster-robust SEs,
FDR across the four annotations.

| Annotation | n | coef | ΔR² | cluster-robust p | Spearman ρ |
|---|---|---|---|---|---|
| Intrinsic disorder (MobiDB) | 20,361 | **+0.593** | **0.046** | 3.3×10⁻¹¹⁸ | +0.217 |
| gnomAD LOEUF | 17,798 | −0.121 | 0.011 | 1.9×10⁻²⁵ | −0.086 |
| DepMap CRISPR essentiality | 17,834 | −0.079 | 0.0042 | 3.7×10⁻¹⁷ | −0.062 |
| gnomAD pLI | 17,798 | +0.051 | 0.0017 | 3.8×10⁻⁸ | +0.052 |

Coverage: disorder 99.7%, gnomAD 87%, DepMap 87%, low-complexity 100%, family size ~100%.

**Interpretation.**

- **Disorder dominates.** More disordered proteins are less predictable to the model. The effect
  survives the low-complexity adjustment (ΔR² ≈ 4.6%, ρ = 0.22) and is ~4–25× larger than any
  constraint effect.
- **Constraint is weak but real and consistent.** All three gene-constraint measures are highly
  significant yet have small effect sizes (ΔR² 1.1% / 0.4% / 0.2%). Crucially they **agree in
  direction**: lower LOEUF, higher pLI, and more-negative DepMap gene-effect all denote greater
  constraint/essentiality, and the coefficient signs (−, +, −) all indicate that more-constrained
  genes are marginally *less* predictable. Three independent constraint axes — loss-of-function
  depletion, dosage sensitivity, and CRISPR fitness dependency — converging on the same weak effect
  argues it is real rather than dataset-specific.
- **Read effect sizes, not p-values.** At n ≈ 18,000 a ~0.2% effect reaches p ≈ 10⁻⁸. Significance
  here means "non-zero," not "important." The ΔR² and Spearman columns carry the meaning.

## 5. The pilot, and a power caveat

A 500-protein stratified pilot found the disorder association already present and robust
(coef +0.52, ΔR² ≈ 3.4%, FDR p = 1.7×10⁻⁴), but **all three constraint measures non-significant**
(LOEUF p = 0.66, pLI p = 0.22, DepMap p = 0.16). The proteome analysis shows those constraint "nulls"
were an artifact of limited power: a ~1% effect is undetectable at n = 441 but trivially significant at
n ≈ 18,000. Underpowered nulls should not be read as absence of effect — the pilot's clean "disorder
yes, constraint no" was, in part, a sample-size story.

## 6. Extremes (illustrative)

The most predictable proteins fall into two textbook classes: repeat-rich / low-complexity (e.g. LPA,
kringle repeats, ~1.0 bits/residue; zinc-finger proteins) and members of very large families (olfactory
receptors; actin). The least predictable are non-canonical or unusual (selenoproteins; short
mitochondrial micro-peptides; pseudogene/antisense products). This motivated including a low-complexity
covariate so that "repetitive → predictable" is separated from the disorder signal.

## 7. Limitations

1. Single model (ESM-2 35M); larger models not yet tested.
2. The disorder association is partly **convergent validity** — MobiDB disorder predictions and ESM
   surprisal both derive from sequence/evolutionary signal, so their agreement is partly two
   conservation measures agreeing, not an independent mechanism.
3. Associations are cross-sectional and do not establish mechanism or causality.
4. Annotation coverage is partial (~87–100%); tests use complete-case subsets.
5. Low-complexity and disorder partially overlap and should be interpreted jointly.
6. Canonical isoforms only.
7. Surprisal is model-relative and influenced by training-data composition (see `METHODS.md`).
