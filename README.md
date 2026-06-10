# Baseline Gut Multi-Omics Profiles Predict Donor-Specific Metabolomic Responses to Dietary Substrates

This repository contains the data and analysis code for the regression-based
machine-learning pipeline described in:

> Kim, G., Shin, H. *Baseline Gut Multi-Omics Profiles Predict Donor-Specific
> Metabolomic Responses to Dietary Substrates: A Regression-Based Machine
> Learning Approach.* (manuscript in preparation)

## Overview

Fecal samples from 19 donors were incubated *in vitro* for 24 h with three
dietary substrates (FOS, laminarin, onion) and a matched control. For each
donor and substrate, a metabolome "response score" was defined as the
treatment-induced shift along PC1 of the metabolome profile (treatment minus
matched control).

Using only each donor's **pre-treatment (baseline) gut microbiota (16S rRNA)
and metabolome profiles**, three regression models (ElasticNet, Random
Forest, and SVR with a linear kernel) were trained to predict this response
score under Leave-One-Donor-Out cross-validation (LODO-CV). Permutation
importance was then used to identify which baseline taxa/metabolites are
most predictive of the metabolomic response, both overall and per substrate.

## Repository structure

```
.
├── README.md
├── requirements.txt
├── data/
│   ├── Merge_metadata.txt                          # sample metadata (donor, group, etc.)
│   ├── feature-count_clr_transposed.tsv            # 16S rRNA ASV table (CLR-transformed)
│   ├── metabolite_common_log2_zscore_transposed.tsv# metabolome table (log2 + z-score)
│   ├── Microbiota_All_transposed.tsv               # full ASV table (raw counts)
│   └── alpha_div.txt                               # alpha diversity metrics
├── code/
│   └── regression_analysis.py                      # core LODO-CV regression pipeline
└── results/
    ├── Figure1_statistics.xlsx
    ├── Figure2_regression_statistics.xlsx
    ├── Figure3_feature_response_correlations.xlsx
    ├── FigureS2_Taxonomy_regression_results.xlsx
    ├── FigureS3_Metabolome_regression_results.xlsx
    ├── FigureS4_Taxonomy_feature_response_correlations.xlsx
    ├── FigureS5_Metabolome_feature_response_correlations.xlsx
    ├── MainFigure2_Metabolome_Differential_statistics_v6.xlsx
    └── ScreePlot_baronly_statistics.xlsx
```

`results/` contains the precomputed statistics underlying the manuscript's
figures and tables. `code/regression_analysis.py` reproduces the core
regression results (`Model_Performance`, `BestModel_Predictions`,
`Importance_BaselineOnly`, `TreatmentSpecific_Perf`, `TS_Importance_*`
sheets), saved to `results/regression_analysis_output.xlsx`.

## Setup

```bash
pip install -r requirements.txt
```

Requires Python 3.10+.

## Usage

```bash
cd code
python regression_analysis.py
```

The script automatically locates `../data` and writes its output to
`../results/regression_analysis_output.xlsx`. No path editing is required.

### Pipeline steps

1. Load raw metadata, 16S rRNA, and metabolome tables.
2. Keep only donors with complete paired Control/FOS/Laminarin/Onion samples.
3. Compute PC1-based metabolome response scores (treatment vs. matched control).
4. Build three baseline predictor sets from each donor's Control-condition
   profile: Taxonomy-only, Metabolome-only, and Combined.
5. Evaluate ElasticNet, Random Forest, and SVR (linear kernel) with
   Leave-One-Donor-Out cross-validation (19 donors x 3 substrates = 57 samples).
6. Select the best-performing input/model combination by Spearman correlation.
7. Compute permutation importance (n_repeats = 30, random_state = 42) for
   the overall baseline model and for each substrate separately.
8. Export all performance metrics, predictions, and importance scores to Excel.

With the settings used in the manuscript (500-tree Random Forest, 30
permutation repeats x 3 models x 3 inputs), the full run takes a few minutes.

## Data availability

Raw 16S rRNA sequencing data are not included in this repository due to file
size / privacy considerations. The processed feature tables provided here
(`feature-count_clr_transposed.tsv`, `metabolite_common_log2_zscore_transposed.tsv`)
are sufficient to reproduce all regression analyses and figures reported in
the manuscript.

## Citation

If you use this code or data, please cite the manuscript above (citation
details to be updated upon publication).
