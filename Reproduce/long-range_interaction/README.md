# Long-Range Gene-Peak Interaction Analysis

This directory contains the analysis of long-range chromatin interactions using MultiGATE's cross-modality attention mechanism, specifically designed to investigate interactions beyond the 150 kb genomic distance constraint.

## Overview

We developed this analysis to demonstrate that while MultiGATE initially employs a 150 kb genomic distance prior, its cross-modality attention framework can effectively identify and score distal regulatory relationships when provided with appropriate external evidence.

## Background

MultiGATE initially restricts candidate peak-gene edges to within 150 kb through its genomic distance prior. However, the attention mechanism itself places no hard limit on distance and can incorporate more distal links when supported by external data sources. To demonstrate this capability, we augmented the standard distance prior with enhancer-promoter contacts from the HiChIP database across 24 diverse human tissues.

## Key Findings

1. **Initial Short-range Analysis**: Under the 150 kb genomic distance prior alone, we identified 16,347 short-range peak-gene candidates.

2. **HiChIP Data Integration**: We augmented this prior by adding 205,627 enhancer-promoter contacts from HiChIP experiments across 24 human tissues, including 10,644 brain-specific interactions.

3. **Distance-stratified Analysis**: After retraining MultiGATE with augmented connections, we stratified all candidate edges into eight distance bins (0-150 kb, 150-300 kb, ..., >1.25 Mb) and compared attention scores.

4. **Tissue-specific Patterns**: Brain-specific HiChIP loops received significantly higher attention scores than other loops across all distance categories (Mann-Whitney U test, P<1×10⁻³).



## Files and Structure

### Main Analysis Notebook
- `longrange_detection_notebook.ipynb`: Complete analysis workflow implementing the long-range interaction detection pipeline

### Input Data
- `Results/MultiGATE_gene_peak_attention.npz`: Sparse matrix of gene-peak attention scores from MultiGATE model
- `Results/gencode.v25.chr_patch_hapl_scaff.annotation.gtf.gz`: Gene annotation file for coordinate mapping
- `Results/peak_gene_hic_df_alignchipdbintohumandata.csv`: HiChIP database with tissue-specific enhancer-promoter contacts

### Output Files
- `gene_peak_interactions_final.csv`: Complete processed dataset with all interactions and annotations
- `interaction_summary_stats.csv`: Statistical summary by distance category and tissue type
- `interaction_score_trends_with_stats_Final.pdf`: Publication-quality visualization

## Analysis Workflow

1. **Data Loading**: Load MultiGATE attention scores, RNA-seq and ATAC-seq data
2. **Processing**: Convert sparse matrix to DataFrame, filter meaningful interactions
3. **Annotation**: Integrate gene coordinates from GTF file
4. **HiChIP Integration**: Classify interactions by tissue specificity (brain vs. others)
5. **Distance Analysis**: Stratify into 8 distance bins (0-150 kb to >1.25 Mb)
6. **Statistics**: Mann-Whitney U tests, effect size calculations, confidence intervals
7. **Visualization**: Publication-quality plots with statistical annotations
8. **Export**: Save processed datasets and visualizations

## Technical Requirements

- Python 3.8+ with scientific computing libraries (scanpy, pandas, numpy, scipy, matplotlib, seaborn)
- GPU recommended for large matrix operations
- 16GB+ RAM for sparse matrix operations

## Statistical Methods

- **Distance bins**: 8 categories from 0-150 kb to >1.25 Mb
- **Statistical test**: Mann-Whitney U test with Cohen's d effect size
- **Tissue classification**: Brain-specific, other tissues, unclassified (based on HiChIP data)

## Key Results

- Brain-specific HiChIP loops show significantly higher attention scores across all distance categories (P<1×10⁻³)
- Long-range interactions demonstrate stronger tissue specificity than short-range interactions
- Identified interactions are supported by independent eQTL and GWAS evidence
- Target genes are enriched for neurological functions

## Reproducibility

To reproduce this analysis:

1. Ensure all input files are available in the `Results/` directory
2. Install required Python packages and dependencies
3. Run the `longrange_detection_notebook.ipynb` notebook sequentially
4. Output files will be generated in the current directory



