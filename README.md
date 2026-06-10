# Gut Multi-Omics Response Prediction

Code accompanying the manuscript:
"Baseline Gut Multi-Omics Profiles Predict Donor-Specific Metabolomic
Responses to Dietary Substrates: A Regression-Based Machine Learning Approach"

## Contents
- `ML_code_v2.ipynb`: Regression analysis (ElasticNet, Random Forest, SVR)
  with leave-one-donor-out cross-validation, permutation importance, and
  figure generation (Figures 2-3, S4-S7)
- `data/`: processed feature tables used as model input
- `requirements.txt`: Python package versions

## Usage
1. Install dependencies: `pip install -r requirements.txt`
2. Update `BASE_DIR` in the notebook to point to the `data/` folder
3. Run all cells

## Raw data availability
16S rRNA sequencing data are available in Qiita under accession numbers
15816 (FOS), 15852 (laminarin), and 15658 (onion).
