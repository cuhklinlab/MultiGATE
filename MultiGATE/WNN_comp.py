import pandas as pd;import numpy as np
import scanpy as sc;import matplotlib.pyplot as plt;import os;import sys
from sklearn.metrics.cluster import adjusted_rand_score

os.environ['R_HOME'] = "/lustre/project/Stat/s1155077016/condaenvs/Seurat4/lib/R"  #/users/s1155077016/anaconda3/envs/seurate4
os.environ['R_USER'] = '/users/s1155077016/anaconda3/lib/python3.9/site-packages/rpy2'


base_path ='/lustre/project/Stat/s1155077016/Spatial_RNA_Protein/mouse_spleen/output'


#'/lustre/project/Stat/s1155077016/Spatial_multi-omics_data/P22/output'
#'/lustre/project/Stat/s1155077016/Spatial_multi-omics_data/Human/output'
#'/lustre/project/Stat/s1155077016/Spatial_RNA_Protein/mouse_spleen/output'

file_name="/lustre/project/Stat/s1155077016/Spatial_RNA_Protein/mouse_spleen/spleen_RNA_processed.h5ad"

#"/lustre/project/Stat/s1155077016/Spatial_multi-omics_data/P22/P22_RNA.h5ad" #RNA
#"/lustre/project/Stat/s1155077016/Spatial_multi-omics_data/Human/Human_RNA.h5ad"
#"/lustre/project/Stat/s1155077016/Spatial_RNA_Protein/mouse_spleen/spleen_RNA_processed.h5ad"

stage = os.path.basename(file_name).split('_')[0]
adata1 = sc.read_h5ad(file_name)

#adata1.obsm["spatial"] = adata1.obsm["spatial"][:, [1, 0]]
#adata1.obsm["spatial"] = adata1.obsm["spatial"]*(-1)
#P22

#adata1.obsm["spatial"][:, 1]= adata1.obsm["spatial"][:, 1]*(-1)
#Human

adata1.obsm["spatial"] = adata1.obsm["spatial"][:, [1, 0]]
adata1.obsm["spatial"][:, 1]= adata1.obsm["spatial"][:, 1]*(-1)
#Protein
#adata1.obsm["spatial"] = adata1.obsm["spatial"]*(-1)
#slide-tags
file_name= "/lustre/project/Stat/s1155077016/Spatial_RNA_Protein/mouse_spleen/spleen_Pro_processed.h5ad"

#"/lustre/project/Stat/s1155077016/Spatial_multi-omics_data/P22/P22_ATAC_lsi.h5ad" #ATAC  _lsi
#"/lustre/project/Stat/s1155077016/Spatial_multi-omics_data/Human/Human_ATAC_lsi.h5ad" 
#"/lustre/project/Stat/s1155077016/Spatial_RNA_Protein/mouse_spleen/spleen_Pro_processed.h5ad"

#stage = os.path.basename(file_name).split('_')[0]
adata2 = sc.read_h5ad(file_name)
#adata2.obsm["spatial"] = adata2.obsm["spatial"][:, [1, 0]]
#adata2.obsm["spatial"] = adata2.obsm["spatial"]*(-1)
#P22

#adata2.obsm["spatial"][:, 1]= adata2.obsm["spatial"][:, 1]*(-1) 
#Human

adata2.obsm["spatial"] = adata2.obsm["spatial"][:, [1, 0]]
adata2.obsm["spatial"][:, 1]= adata2.obsm["spatial"][:, 1]*(-1)
#Protein
#adata1.obsm["spatial"] = adata1.obsm["spatial"]*(-1)
#slide-tags

from sklearn.decomposition import PCA
from scipy.sparse import issparse
def wnn(adata1, adata2, res=0.5, algo =  1):
    import rpy2.robjects as ro
    ro.r.library("Seurat")
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    import rpy2.robjects.pandas2ri
    rpy2.robjects.pandas2ri.activate()
    pca = PCA(n_components=50, svd_solver='arpack',random_state=2024)
    if issparse(adata1.X):
        data1 = adata1.X.A
    else:
        data1 = adata1.X        
    pca.fit(data1)
    data1 = pca.transform(data1)
    if issparse(adata2.X):
        data2 = adata2.X.A
    else:
        data2 = adata2.X       
    pca.fit(data2)
    data2 = pca.transform(data2)
    rna_matrix =  rpy2.robjects.numpy2ri.numpy2rpy(data1)#robjects.r['as.matrix'](rna_df)
    atac_matrix = rpy2.robjects.numpy2ri.numpy2rpy(data2)
    barcodes = rpy2.robjects.pandas2ri.py2rpy(adata1.obs.index)
    ro.r.assign("rna_matrix", rna_matrix)
    ro.r.assign("atac_matrix", atac_matrix)
    ro.r.assign("barcodes", barcodes)
    ro.r.assign("res", res)  #res
    ro.r.assign("algo", algo)  #algo
    ro.r(
    '''
    colnames(rna_matrix) <- paste0("pca_", 1:dim(rna_matrix)[2])
    colnames(atac_matrix) <- paste0("lsi_", 1:dim(rna_matrix)[2])
    rownames(rna_matrix) <- barcodes
    P22 <- CreateSeuratObject(counts = t(rna_matrix))  
    rownames(atac_matrix) <- barcodes
    atac_obj <- CreateAssayObject(counts = t(atac_matrix)) 
    P22[["ATAC"]] <- atac_obj
    P22[["pca"]] = CreateDimReducObject(embeddings = rna_matrix, key = "pca_", assay = "RNA")
    P22[["lsi"]] = CreateDimReducObject(embeddings = atac_matrix, key = "lsi_", assay = "ATAC")
    P22 <- FindMultiModalNeighbors(P22, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
    P22 <- RunUMAP(P22, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    P22 <- FindClusters(P22, graph.name = "wsnn", algorithm = algo, resolution = res,verbose = FALSE)
    clusters = P22@meta.data$seurat_clusters 
    wwnumap = P22[["wnn.umap"]]
    umap = Embeddings(object = wwnumap)#[1:nrow(rna_matrix), 1:2]##wwnumap[1:nrow(rna_matrix), 1:2]#P22[["wnn.umap"]]  wwnumap[[1:nrow(rna_matrix), 1:2]] obj@reductions$umap
    '''
    )
    wnn_UMAP = ro.r('umap')
    wnn = ro.r('clusters')
    adata1.obs['wnn'] = wnn
    adata1.obs['wnn'] = adata1.obs['wnn'].astype('int')
    adata1.obs['wnn'] = adata1.obs['wnn'].astype('category')
    adata2.obs['wnn'] = wnn
    adata2.obs['wnn'] = adata2.obs['wnn'].astype('int')
    adata2.obs['wnn'] = adata2.obs['wnn'].astype('category')
    adata1.obsm['X_umap'] = wnn_UMAP
    adata2.obsm['X_umap'] = wnn_UMAP
    return adata1

def wnn_protein(adata1, adata2, res=0.5, algo =  1):
    import rpy2.robjects as ro
    ro.r.library("Seurat")
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    import rpy2.robjects.pandas2ri
    rpy2.robjects.pandas2ri.activate()
    pca = PCA(n_components=50, svd_solver='arpack',random_state=2024)
    if issparse(adata1.X):
        data1 = adata1.X.A
    else:
        data1 = adata1.X        
    pca.fit(data1)
    data1 = pca.transform(data1)
    data2 = adata2.X
    rna_matrix =  rpy2.robjects.numpy2ri.numpy2rpy(data1)#robjects.r['as.matrix'](rna_df)
    atac_matrix = rpy2.robjects.numpy2ri.numpy2rpy(data2)
    barcodes = rpy2.robjects.pandas2ri.py2rpy(adata1.obs.index)
    ro.r.assign("rna_matrix", rna_matrix)
    ro.r.assign("atac_matrix", atac_matrix)
    ro.r.assign("barcodes", barcodes)
    ro.r.assign("res", res)  #res
    ro.r.assign("algo", algo)  #algo
    ro.r(
    '''
    colnames(rna_matrix) <- paste0("pca_", 1:dim(rna_matrix)[2])
    colnames(atac_matrix) <- paste0("lsi_", 1:dim(atac_matrix)[2])
    rownames(rna_matrix) <- barcodes
    P22 <- CreateSeuratObject(counts = t(rna_matrix))  
    rownames(atac_matrix) <- barcodes
    atac_obj <- CreateAssayObject(counts = t(atac_matrix)) 
    P22[["ATAC"]] <- atac_obj
    P22[["pca"]] = CreateDimReducObject(embeddings = rna_matrix, key = "pca_", assay = "RNA")
    P22[["lsi"]] = CreateDimReducObject(embeddings = atac_matrix, key = "lsi_", assay = "ATAC")
    P22 <- FindMultiModalNeighbors(P22, reduction.list = list("pca", "lsi"), dims.list = list(1:dim(rna_matrix)[2], 1:dim(atac_matrix)[2]))
    P22 <- RunUMAP(P22, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    P22 <- FindClusters(P22, graph.name = "wsnn", algorithm = algo, resolution = res,verbose = FALSE)
    clusters = P22@meta.data$seurat_clusters 
    wwnumap = P22[["wnn.umap"]]
    umap = Embeddings(object = wwnumap)#[1:nrow(rna_matrix), 1:2]##wwnumap[1:nrow(rna_matrix), 1:2]#P22[["wnn.umap"]]  wwnumap[[1:nrow(rna_matrix), 1:2]] obj@reductions$umap
    '''
    )
    wnn_UMAP = ro.r('umap')
    wnn = ro.r('clusters')
    adata1.obs['wnn'] = wnn
    adata1.obs['wnn'] = adata1.obs['wnn'].astype('int')
    adata1.obs['wnn'] = adata1.obs['wnn'].astype('category')
    adata2.obs['wnn'] = wnn
    adata2.obs['wnn'] = adata2.obs['wnn'].astype('int')
    adata2.obs['wnn'] = adata2.obs['wnn'].astype('category')
    adata1.obsm['X_umap'] = wnn_UMAP
    adata2.obsm['X_umap'] = wnn_UMAP
    return adata1

#spatialglue adata
#file_name = '/lustre/project/Stat/s1155077016/SpatialGlue_data/human/human.h5ad'
#adata = sc.read_h5ad(file_name)

#
#adata1 = wnn(adata1, adata2, res=0.5)

#df = pd.read_csv('/lustre/project/Stat/s1155077016/Spatial_multi-omics_data/Human/true_label_0202.csv',header = 0, index_col = 0)
#shared_index = df.index.intersection(adata1.obs_names)  
# adata spatial_glue results; adata1 other methods

# Update 'adata.obs['true_label']' with the values from the matched part of the first column in 'df'
#adata1.obs['true_label'] = pd.Series(df.loc[shared_index].iloc[:, 0], index=shared_index)
#adata1.obs['true_label'] = adata1.obs['true_label'].astype('category')
#obs_df = adata1.obs.dropna()

#shared_index = df.index.intersection(adata.obs_names)
#obs_df= obs_df.iloc[df.loc[shared_index].iloc[:, 0]] 
#ARI = adjusted_rand_score(obs_df['wnn'], obs_df["true_label"])
#print('ARI: ', ARI)

adata1 = wnn_protein(adata1, adata2, res=1.0)
# P22 0.5 Human 1.0 Protein 1.0 slide-tags 0.15
size = 15 # P22 15 Human 25 Protein 20 slide-tags 15
plt.rcParams["figure.figsize"] = (3, 3);
sc.pl.embedding(adata1, basis="spatial", color="wnn",s=size, show=False, title='Seurat')
plt.savefig(os.path.join(base_path,"0329_Seurat_"+stage+".png"), bbox_inches="tight",dpi=600) #pred_bp1k.png avremb
plt.axis('off')
plt.close()

plt.rcParams["figure.figsize"] = (3, 3);
sc.pl.umap(adata1, color="wnn",title='WNN');
plt.savefig(os.path.join(base_path,"0329_Seurat_umap_graph_"+stage+".png"), bbox_inches="tight",dpi=600) ;plt.axis('off');plt.close()

adata1.write(os.path.join(base_path,'adata_Seurat_'+stage+'.h5ad'))
print('finished')