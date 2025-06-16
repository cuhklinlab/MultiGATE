# %%
# Configure available GPUs
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

import warnings
warnings.filterwarnings('ignore')

# Add MultiGATE package path to system path
import sys
sys.path.insert(0, '/lustre/project/Stat/s1155202250/fastfolder/code/st/MultiGATE/MultiGATEgithub0607/MultiGATE/Reproduce/clustering')

# Import required packages
import MultiGATE
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd

# Training configuration
num_epoch = 3500

# %%
# Setup output directory
base_path = './outputSMST0209/'
os.makedirs(base_path, exist_ok=True)

# Load spatial transcriptomics (SPT) and spatial metabolomics (SPM) data
SPT = sc.read('./InputData/SPT.h5ad')
SPM = sc.read('./InputData/SPM.h5ad')

# Calculate and display spatial network for transcriptomics data
MultiGATE.Cal_Spatial_Net(SPT, rad_cutoff=1000)
MultiGATE.Stats_Spatial_Net(SPT)

# %%
# Calculate and display spatial network for metabolomics data
MultiGATE.Cal_Spatial_Net(SPM, rad_cutoff=1000)
MultiGATE.Stats_Spatial_Net(SPM)

# %%
# Load prior knowledge links between metabolites and genes
prior_links = pd.read_csv("./InputData/prior_links.csv")
print("Prior links data shape:", prior_links.shape)

# %%
# Display transcriptomics data info
print("Spatial transcriptomics data info:")
print(SPT)

# %%
# Calculate highly variable genes for transcriptomics data
sc.pp.highly_variable_genes(SPT, flavor="seurat_v3", n_top_genes=2000)

def Cal_mz_gene_Net(spm, spt, prior_links, verbose=True):
    """
    Construct the links between metabolites (m/z) and genes based on prior knowledge.

    Parameters
    ----------
    spm : AnnData
        AnnData object containing metabolomics data
    spt : AnnData
        AnnData object containing transcriptomics data
    prior_links : DataFrame
        DataFrame containing prior knowledge of m/z-gene links
    verbose : bool, optional
        Whether to print progress messages (default: True)

    Returns
    -------
    None
        The m/z-gene links are saved in both adata objects' uns['gene_peak_Net']
    """
    if verbose:
        print('------Calculating m/z-gene network...')
    
    # Extract metabolite and gene names
    mz_names = spm.var.values.flatten()
    gene_names = spt.var['gene_name'].values.flatten()
    gene_names = [gene.upper() for gene in gene_names]  # Convert to uppercase for matching
    
    # Initialize list for matched m/z-gene pairs
    matched_pairs = []
    
    # Process each row in prior_links to find valid m/z-gene connections
    for idx, row in prior_links.iterrows():
        # Check if RaMPdb entry contains gene information
        if (row['RaMPdb'] != 0 and 
            isinstance(row['RaMPdb'], str) and 
            row['RaMPdb'].startswith('{')):
            
            mz = row['mz_names']
            linked_genes = eval(row['RaMPdb'])  # Convert string to list
            
            # Check if m/z exists in metabolomics data
            if mz in mz_names:
                # Check each linked gene
                for gene in linked_genes:
                    gene = gene.upper()  # Convert to uppercase for matching
                    if gene in gene_names:
                        matched_pairs.append((mz, gene))
    
    # Create DataFrame of matched pairs and save to both datasets
    gene_peak_net = pd.DataFrame(matched_pairs, columns=["Peak", "Gene"])
    gene_peak_net = gene_peak_net[['Gene', 'Peak']]  # Reorder columns
    
    # Save network to both AnnData objects
    spm.uns['gene_peak_Net'] = gene_peak_net
    spt.uns['gene_peak_Net'] = gene_peak_net
    
    if verbose:
        print('------Network calculation completed')
        print(f'Found {len(matched_pairs)} m/z-gene links')

# %%
# Calculate m/z-gene network using prior knowledge
Cal_mz_gene_Net(SPM, SPT, prior_links)

# %%
# Import additional required modules for MultiGATE training
import numpy as np
import sklearn.neighbors
from sklearn.neighbors import NearestNeighbors
import scipy.sparse as sp
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import normalize
import networkx as nx
from networkx.algorithms.bipartite import biadjacency_matrix   
import MultiGATE.genomics as genomics

# Prepare data for MultiGATE training
# Use metabolomics data as adata1 and highly variable genes from transcriptomics as adata2
adata1 = SPM
adata2 = SPT[:, SPT.var['highly_variable']]

# Set proper indices for matching
adata1.var.index = adata1.var['mz']
adata2.var.index = adata2.var['gene_name'].str.upper()

# Recalculate m/z-gene network for prepared data
Cal_mz_gene_Net(adata1, adata2, prior_links)

# Train MultiGATE model
print("Starting MultiGATE training...")
adata1, adata2 = MultiGATE.train_MultiGATE(
    adata2, adata1,
    hidden_dims=[256, 10],
    type='protein',
    n_epochs=num_epoch,
    save_attention=False,
    temp=-3,
    protein_value=0.01
)

# %%
# Setup R environment for clustering
os.environ['R_HOME'] = "/lustre/project/Stat/s1155077016/condaenvs/Seurat4/lib/R"
os.environ['R_USER'] = '/users/s1155077016/anaconda3/lib/python3.9/site-packages/rpy2'

# Clustering configuration
size = 20
n_clusters = 11

# Perform clustering using mclust
print("Performing clustering...")
adata1 = MultiGATE.mclust_R(adata1, used_obsm='MultiGATE_clip_all', num_cluster=n_clusters)

# %%
# Evaluate clustering results and create visualizations
from sklearn.metrics import adjusted_rand_score


# Create visualization plots
plt.rcParams["figure.figsize"] = (12, 3)
fig, axs = plt.subplots(1, 3)

# Plot clustering results in spatial coordinates
sc.pl.embedding(adata1, basis='spatial', color='mclust', 
                title=f'MultiGATE', s=size, 
                show=False, ax=axs[0], legend_loc=None)

# Create UMAP visualization
sc.pp.neighbors(adata1, use_rep='MultiGATE_clip_all', key_added='avg')
sc.tl.umap(adata1, neighbors_key='avg')
sc.pl.umap(adata1, color="mclust", title='MultiGATE', ax=axs[2], show=False)

# Plot ground truth annotations
sc.pl.embedding(adata1, basis='spatial', color='annotations', 
                title='Ground Truth', s=size, ax=axs[1], legend_loc=None)

plt.tight_layout()
plt.savefig(f'{base_path}/clustering.pdf')
plt.close()

# %%
# Clean up and save results
print("Cleaning up and saving results...")

# Keep only necessary obsm matrices to reduce file size
obsm_to_keep = ['spatial', 'MultiGATE', 'X_umap']
for key in list(adata1.obsm.keys()):
    if key not in obsm_to_keep:
        del adata1.obsm[key]

# Remove attention matrices if present to save memory
if 'MultiGATE_attention' in adata1.uns:
    del adata1.uns['MultiGATE_attention']
if 'MultiGATE_gene_peak_attention' in adata1.uns:
    del adata1.uns['MultiGATE_gene_peak_attention']

# Save processed data
output_file = base_path + 'STSM_MultiGATE.h5ad'
adata1.write_h5ad(output_file)
print(f"Results saved to: {output_file}")