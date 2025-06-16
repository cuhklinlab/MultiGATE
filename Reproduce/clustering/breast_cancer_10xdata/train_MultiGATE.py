import os
import argparse
import warnings
import sys

# Add MultiGATE package path to system path
sys.path.insert(0, '/lustre/project/Stat/s1155202250/fastfolder/code/st/MultiGATE/MultiGATEgithub0607/MultiGATE/Reproduce/clustering')

import MultiGATE
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.decomposition import PCA
from scipy.sparse import issparse

# Set CUDA visible devices
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

# Suppress warnings
warnings.filterwarnings('ignore')

# Set output directory
base_path = "./outputbreast0307/"
os.makedirs(base_path, exist_ok=True)

# Load RNA data
file_name = "./Inputdata/adata_rna.h5ad" 
adata1 = sc.read_h5ad(file_name)
adata1.obsm["spatial"][:, 1] *= -1  # Adjust spatial coordinates

# Load protein data
file_name = "./Inputdata/adata_protein.h5ad"
adata2 = sc.read_h5ad(file_name)
adata2.obsm["spatial"][:, 1] *= -1  # Adjust spatial coordinates

# Filter RNA data to highly variable genes only
adata1 = adata1[:, adata1.var['highly_variable']]

# Calculate gene-protein network
MultiGATE.Cal_gene_protein_Net(adata1, adata2) 
adata1.uns['gene_peak_Net'] = adata2.uns['gene_peak_Net']

# Train MultiGATE model
adata1, adata2 = MultiGATE.train_MultiGATE(
    adata1, adata2,  
    type='protein', 
    n_epochs=1000, 
    save_attention=False, 
    temp=-12,
    protein_value=1,
)

# Set R environment variables for WNN integration
os.environ['R_HOME'] = "/lustre/project/Stat/s1155077016/condaenvs/Seurat4/lib/R" 
os.environ['R_USER'] = '/users/s1155077016/anaconda3/lib/python3.9/site-packages/rpy2'

# Perform WNN (Weighted Nearest Neighbor) integration using R
adata1, adata2 = MultiGATE.wnn_R(adata1, adata2, res=0.24)

# Create visualization plots
size = 20 
plt.rcParams["figure.figsize"] = (7, 3) 
fig, axs = plt.subplots(1, 2)

# Plot spatial and UMAP visualizations
sc.pl.embedding(adata1, basis="spatial", color="wnn", s=size, title='MultiGATE', ax=axs[0], show=False, legend_loc=None)
sc.pl.umap(adata1, color="wnn", title='MultiGATE', ax=axs[1], show=False)
plt.tight_layout()
plt.savefig(base_path + f'human_plot.png')
plt.close()

# Clean up adata1 object before saving
# Remove attention matrices if present
if 'MultiGATE_attention' in adata1.uns:
    del adata1.uns['MultiGATE_attention']
if 'MultiGATE_gene_peak_attention' in adata1.uns:
    del adata1.uns['MultiGATE_gene_peak_attention']

# Keep only necessary obs columns
obs_to_keep = ['wnn']
for col in list(adata1.obs.columns):
    if col not in obs_to_keep:
        del adata1.obs[col]

# Keep only necessary obsm matrices
obsm_to_keep = ['spatial', 'MultiGATE', 'X_umap']
for key in list(adata1.obsm.keys()):
    if key not in obsm_to_keep:
        del adata1.obsm[key]

# Store separate RNA and protein embeddings
adata1.obsm['MultiGATE_RNA'] = adata1.obsm['MultiGATE']
adata1.obsm['MultiGATE_protein'] = adata2.obsm['MultiGATE']

# Save the processed data
adata1.write_h5ad(base_path + f'breast_cancer_adata_MultiGATE.h5ad')
