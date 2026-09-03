# TMPS reproducibility package

This folder provides a minimal, public-repository-ready implementation of the
Therapy-anchored Multi-omics Prioritization Score (TMPS) prepared for a Functional & Integrative Genomics methodology submission package.

TMPS is a transparent rule-based multi-omics prioritization algorithm for
reproducible candidate-mechanism ranking. It is not a machine-learning model,
diagnostic algorithm, treatment-selection model, drug-target discovery
algorithm, or causal validation framework.

## Files

- `data/tmps_input_evidence.csv`: candidate-level evidence table used for TMPS ranking.
- `data/tmps_scoring_parameters.csv`: prespecified evidence domains and interpretation boundaries.
- `data/example_input.csv`: small example input file for checking the scoring workflow.
- `tmps_scoring.py`: rebuilds total TMPS scores and evidence tiers.
- `tmps_ablation_robustness.py`: reruns leave-one-domain-out ablation and alternative weighting checks.
- `tmps_benchmark_rankings.py`: compares TMPS with simpler baseline ranking rules.
- `tmps_random_weight_robustness.py`: runs Monte Carlo random-weight perturbation robustness analysis.
- `outputs/TMPS_Table_S2_scores.csv`: ranked TMPS output mirrored in Supplementary Table S2.
- `outputs/TMPS_Table_S3_ablation_parameter_robustness.csv`: robustness output mirrored in Supplementary Table S3.
- `outputs/TMPS_benchmark_rankings.csv`: comparator output for full TMPS, simple evidence-sum, evidence-count, and genetic-evidence-only rankings.
- `outputs/TMPS_random_weight_summary.csv`: summary of the Monte Carlo random-weight robustness analysis.
- `outputs/TMPS_random_weight_detail.csv`: iteration-level Monte Carlo random-weight output.
- `outputs/example_output.csv`: expected output generated from `data/example_input.csv`.

## Quick start

```bash
python tmps_scoring.py --input data/tmps_input_evidence.csv --out outputs/TMPS_recomputed_scores.csv
python tmps_ablation_robustness.py --input outputs/TMPS_recomputed_scores.csv --out outputs/TMPS_recomputed_ablation_robustness.csv
python tmps_benchmark_rankings.py --input data/tmps_input_evidence.csv --out outputs/TMPS_benchmark_rankings.csv
python tmps_random_weight_robustness.py --input data/tmps_input_evidence.csv --summary-out outputs/TMPS_random_weight_summary.csv --detail-out outputs/TMPS_random_weight_detail.csv
python tmps_scoring.py --input data/example_input.csv --out outputs/example_output.csv
```

## Interpretation

Tier 1 indicates strong cross-resource support without major traceability
penalties. Tier 2 indicates multi-domain support with remaining boundaries.
Exploratory indicates incomplete, context-only, or hypothesis-generating
evidence requiring follow-up. The manuscript interprets APEX1 as the stable
Tier 1 candidate and retains PRDX2 and KDM1A as exploratory signals.

## Repository note for submission

Before final submission, this folder can be uploaded unchanged as supplementary
code or pushed to a public GitHub/Zenodo repository. If a repository is created,
replace the placeholder Code availability sentence in the manuscript with the
final URL or DOI.

## Portability

TMPS can be adapted to other immune-mediated diseases by replacing the
predefined clinically anchored search space, disease GWAS, molecular QTL
resources, and disease-context evidence while keeping the evidence domains,
weights, penalties, and interpretation boundaries prespecified before
outcome-aware interpretation.
