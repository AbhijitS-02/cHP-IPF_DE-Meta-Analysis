# Differential Expression Meta-Analysis of NGS Datasets to Classify Chronic Hypersensitivity Pneumonitis and Idiopathic Pulmonary Fibrosis

## Abstract 
**Chronic hypersensitivity pneumonitis (cHP)** is an interstitial lung disease (ILD) driven by aberrant immune response that often progresses to pulmonary fibrosis and exhibits substantial clinical, radiological, and histopathological overlap with **idiopathic pulmonary fibrosis (IPF)**. This overlap frequently leads to misdiagnosis and delays in appropriate therapeutic intervention. To identify robust molecular signatures that discriminate between these overlapping ILD subtypes, we performed ***differential expression meta-analysis*** of next-generation sequencing (NGS) transcriptomic datasets (GSE184316 and GSE150910). Integration using the ***Fused Inverse Normal (FINM)*** framework identified **380 significantly dysregulated genes** between cHP and IPF. Gene ontology enrichment highlighted immune- and receptor-mediated signaling predominating in cHP, whereas IPF was enriched for developmental and fibrotic remodeling pathways. ***Minimum redundancy–maximum relevance (mRMR)*** feature ranking followed by ***LASSO regression*** yielded a parsimonious **three-gene panel (ZNF443, GBP4, RNF208)** with strong discriminatory accuracy **(AUC = 0.921; 95% CI: 0.851–0.991)**. Future validation in patient population is necessary to establish the candidate genes as a potential biomarker panel.

## Getting started with this repository
1. Clone the repo:
```bash
git clone https://github.com/AbhijitS-02/cHP-IPF_DE-Meta-Analysis.git
cd cHP-IPF_DE-Meta-Analysis
```

2. The notebook **03_mRMR.ipynb** was run on a conda based virtual environment:
```bash
conda create -n cHP-IPF_DE-Meta-Analysis python=3.13.
conda activate cHP-IPF_DE-Meta-Analysis
```
Install python dependencies-
```bash
conda install pandas numpy matplotlib scikit-learn
pip install mrmr_selection ipykernel

# Add conda env to jupyter notebook
python -m ipykernel install --user --name=cHP-IPF_DE-Meta-Analysis --display-name "cHP-IPF_DE-Meta-Analysis"
```


3. Install the required R packages:
```r
install.packages(c(
  "tidyverse",
  "here",
  "ggplot2",
  "reshape2",
  "sva",
  "glmnet",
  "pROC",
  "ggpubr"
))
```

