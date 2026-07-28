# Delta-surprisal variant-effect analysis (companion to the wild-type surprisal atlas)

Reproducible script: `delta_surprisal_full_comparison.py` (self-contained, GPU, one cell).

## Motivation

The main atlas measures *wild-type* per-residue surprisal across the human proteome and finds it
is dominated by intrinsic disorder, with only weak associations to gene-level constraint. This
raises a natural question: does the *mutation-effect* form of surprisal behave differently?

Delta-surprisal is `surprisal(mut) - surprisal(wt)` at a masked position, i.e.
`log p(wt | context) - log p(mut | context)`. It is mathematically the negative of the standard
ESM zero-shot variant-effect score, so a positive prior exists that it should separate pathogenic
from benign missense. This analysis tests that directly, with a trivial baseline and the current
state of the art as reference points.

## Data and method

ESM-2 masked-marginal scoring, windowed to 1022 residues for long proteins. Variants: ClinVar
pathogenic/likely-pathogenic vs benign/likely-benign missense, balanced per gene, across eight
cancer-predisposition genes (TP53, BRCA1, BRCA2, PTEN, MLH1, MSH2, VHL, RET). n = 2145
(1269 pathogenic / 876 benign). Baselines: -BLOSUM62 (substitution matrix) and AlphaMissense
(precomputed scores, 100% coverage on this set). Metric: AUROC with bootstrap 95% CIs.

## Result 1: model scale flips the conclusion

The same method gives opposite verdicts at different model sizes.

| model      | TP53 AUROC | pooled AUROC |
|------------|-----------|--------------|
| ESM-2 35M  | 0.475     | 0.698 (~= BLOSUM 0.690) |
| ESM-2 650M | 0.916     | 0.889        |

At 35M, delta-surprisal is indistinguishable from a coin on TP53 and does not beat BLOSUM overall.
At 650M it is a strong predictor. The small-model negative is a false negative. An earlier mechanistic
story (that delta-surprisal is blind to DNA-contact hotspots like TP53's) was an artifact of model
size, not a real limitation.

## Result 2: strong signal, but dominated by AlphaMissense

ESM-2 650M, per gene, delta-surprisal vs the two references:

| gene   | n    | path/benign | delta-surprisal | AlphaMissense | -BLOSUM62 |
|--------|------|-------------|-----------------|---------------|-----------|
| POOLED | 2145 | 1269/876    | 0.889           | **0.955**     | 0.694 |
| BRCA1  | 400  | 200/200     | 0.946           | **0.954**     | 0.671 |
| TP53   | 353  | 200/153     | 0.916           | **0.980**     | 0.689 |
| MSH2   | 341  | 141/200     | 0.897           | **0.919**     | 0.827 |
| BRCA2  | 313  | 113/200     | 0.889           | **0.947**     | 0.691 |
| MLH1   | 221  | 185/36      | 0.888           | **0.927**     | 0.690 |
| VHL    | 201  | 166/35      | 0.935           | **0.959**     | 0.700 |
| RET    | 109  | 64/45       | 0.864           | **0.984**     | 0.768 |
| PTEN   | 207  | 200/7*      | 1.000*          | 0.997         | 0.706 |

*PTEN benign n = 7; its AUROC is not interpretable and is shown only for completeness.

Rank-average ensemble (delta-surprisal + AlphaMissense) = **0.940**, which is *below* AlphaMissense
alone (0.955). Combining the two does not help, which is the signature of redundancy: delta-surprisal
carries no signal that AlphaMissense is missing.

## Conclusion

Delta-surprisal is a real pathogenicity signal, unlike wild-type per-residue surprisal, and it
comfortably beats a substitution-matrix baseline at adequate model scale. But it reproduces, rather
than exceeds, known ESM zero-shot variant-effect performance, and it is dominated by AlphaMissense
with no orthogonal contribution. It is reported here as validation of the signal and as an honest
negative on novelty, not as a competitive variant-effect predictor.

Practical notes for anyone reproducing: use the 650M model or larger (35M is underpowered for this
task), always compare against BLOSUM as a floor, and check class balance per gene before trusting a
per-gene AUROC.
