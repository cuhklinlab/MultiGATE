# MultiGATE Cross-Modality Attention Score Analysis for Cis-Regulation

## Overview

This repository contains the analysis pipeline for evaluating MultiGATE's cross-modality attention scores in predicting cis-regulatory relationships between genes and chromatin accessibility peaks. The analysis compares MultiGATE against established baseline methods using eQTL associations as ground truth.

## Objective

To assess the performance of MultiGATE attention scores for identifying gene-peak regulatory relationships in spatial multi-omics data compared to conventional computational approaches.

## Data

- **Multi-omics Data**: Human spatial RNA-seq and ATAC-seq datasets
- **Ground Truth**: eQTL (expression Quantitative Trait Loci) associations from publicly available databases
- **Distance Range**: Gene-peak pairs within 150kb genomic distance

## Methods

### Baseline Comparisons
1. **Spearman Correlation**: Expression-accessibility correlation using metacell aggregation
2. **LASSO Regression**: L1-regularized linear regression for feature selection
3. **Cicero**: Co-accessibility analysis for chromatin interaction prediction
4. **MultiGATE**: Cross-modality attention scores (our method)

### Analysis Pipeline
1. **Data Preprocessing**: Spatial multi-omics data integration and quality control
2. **Distance Graph Construction**: 150kb window-based gene-peak pairing
3. **Method Implementation**: Parallel computation of all baseline methods
4. **Performance Evaluation**: ROC-AUC and Precision-Recall analysis

## Key Results

Performance metrics for eQTL prediction (AUROC/AP scores):
- MultiGATE: Superior performance in both ROC-AUC and Average Precision
- Spearman Correlation: Traditional correlation-based approach
- LASSO Regression: Feature selection-based method  
- Cicero: Co-accessibility-based prediction

## Output Files

- `eqtl_roc.pdf`: ROC curves comparing all methods
- `eqtl_prc.pdf`: Precision-recall curves for method comparison
- `dist_binned_score_pchic.pdf`: Distance-stratified score distributions
- Results saved to `results_eqtl_only_20250613/` directory

## Dependencies

- Python 3.x
- AnnData, Scanpy (single-cell analysis)
- NetworkX (graph analysis)
- Scikit-learn (machine learning)
- Pandas, NumPy (data manipulation)
- Matplotlib, Seaborn (visualization)

## Usage

Execute the Jupyter notebook `eqtl_inference20250613.ipynb` to reproduce the complete analysis pipeline. The notebook contains sequential steps for data loading, method implementation, and comparative evaluation.
