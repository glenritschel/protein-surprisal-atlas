# Biological Association Analyses

## Overview
This report presents the association between biological constraints (gnomAD metrics and intrinsic disorder) and protein-language-model surprisal (bits per residue), after adjusting for sequence length, homolog abundance (log_family_size), and low-complexity fraction.

## Interpretation Caveats
- **Coverage Limitations**: Biological annotations have partial coverage across the proteome. The sample size ($n$) varies across tests, and associations are evaluated only on complete-case subsets.
- **Collinearity**: Disorder and low-complexity partly overlap (disorder regions are often low-complexity), so interpret their coefficients jointly, not in isolation.
- **Causality**: Associations are cross-sectional and do not establish biological mechanism or causality.

## Results

| Annotation | n | Effect Size (Coef) | Delta R² | Robust p-value | Adjusted p-value (FDR) | Spearman Rank Corr |
|---|---|---|---|---|---|---|
| gnomad_loeuf | 441 | -0.0256 | 0.0005 | 6.57e-01 | 6.57e-01 | -0.0170 |
| gnomad_pli | 441 | -0.0684 | 0.0034 | 2.22e-01 | 2.96e-01 | 0.0069 |
| disorder_fraction | 499 | 0.5171 | 0.0338 | 4.23e-05 | 1.69e-04 | 0.1729 |
| depmap_gene_effect | 431 | -0.0659 | 0.0033 | 1.61e-01 | 2.96e-01 | -0.0053 |
