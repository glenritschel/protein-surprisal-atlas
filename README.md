# protein-surprisal-atlas

Estimating protein language-model surprisal (pseudo-log-likelihood in bits) across the reviewed human proteome, and testing its biological signal beyond sequence length and homolog abundance.

## Scientific Warning

**This project estimates model-relative sequence surprisal (pseudo-log-likelihood in bits) under a trained protein language model. It does not calculate a true Shannon code length or achievable compression size, does not calculate Kolmogorov complexity, does not prove protein functionality, and does not provide clinical interpretation of variants. Surprisal is strongly influenced by homolog abundance in the model's training data.**

The quantity `S(x)` is a pseudo-log-likelihood: the per-position conditionals do not multiply into a proper joint distribution and `S(x)` is not an achievable code length. Treat it as a model surprisal score, useful primarily in relative and rank terms.

## Execution
```bash
make download
make build_pilot
make score_pilot
make validate_scoring
make run_analysis
```

See the provided `Makefile` and `config.yaml` for parameters.
