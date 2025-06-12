# TF_validation

This folder contains scripts and notebooks for validating transcription factor (TF) regulatory predictions using multi-omics data.

## Contents
- **SOX2_validation.ipynb**: Jupyter notebook for evaluating the performance of different methods in predicting SOX2 regulatory elements using RNA-seq, ATAC-seq, and ChIP-seq data.
- **Input files**: Required data files include processed RNA and ATAC AnnData objects, attention score CSVs, TF binding matrices, and ChIP-seq peak files.
- **tf_binding_all.pkl**: This file must be downloaded from [Google Drive](https://drive.google.com/file/d/1qraLX1k36OQMuAPUvY44ktrYBRuxmzbJ/view?usp=sharing) and placed in this folder.

## Main Steps
1. **Load Data**: Import RNA, ATAC, attention scores, and TF binding matrices.
2. **Preprocessing**: Filter and align peaks/genes across datasets.
3. **Scoring**: Compute cosine similarity, process attention scores, and prepare motif/TF binding matrices.
4. **Validation**: Compare predicted regulatory elements with ChIP-seq data using AUC and AUPR metrics.
5. **Visualization**: Plot and save comparison results as PDF.

## Usage
1. Place all required data files in the folder.
2. Download `tf_binding_all.pkl` from [Google Drive](https://drive.google.com/file/d/1qraLX1k36OQMuAPUvY44ktrYBRuxmzbJ/view?usp=sharing).
3. Open and run `SOX2_validation.ipynb` step by step.
4. Results and plots will be saved in the same directory. 