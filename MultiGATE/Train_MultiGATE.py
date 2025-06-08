from sklearn.preprocessing import Normalizer
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
from .MultiGATE import MultiGATE
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

def train_MultiGATE(adata1, adata2, hidden_dims=[512, 30], bp_width=450, temp=1.0, type='ATAC', n_epochs=500, lr=0.0001, key_added='MultiGATE',
                    gradient_clipping=5, nonlinear=True, weight_decay=0.0001, verbose=False,
                    random_seed=2020, pre_labels=None, pre_resolution=0.2,
                    save_attention=False, save_loss=False, save_reconstrction=False,protein_value=0.001):
    """\
    Training graph attention auto-encoder.

    Parameters
    ----------
    adata
        AnnData object of scanpy package.
    hidden_dims
        The dimension of the encoder.

    n_epochs
        Number of total epochs in training.
    lr
        Learning rate for AdamOptimizer.
    key_added
        The latent embeddings are saved in adata.obsm[key_added].
    gradient_clipping
        Gradient Clipping.
    nonlinear
        If True, the nonlinear avtivation is performed.
    weight_decay
        Weight decay for AdamOptimizer.
    pre_labels
        The key in adata.obs for the manually designate the pre-clustering results. Only used when alpha>0.
    pre_resolution
        The resolution parameter of sc.tl.louvain for the pre-clustering. Only used when alpha>0 and per_labels==None.
    save_attention
        If True, the weights of the attention layers are saved in adata.uns['STAGATE_attention']
    save_loss
        If True, the training loss is saved in adata.uns['STAGATE_loss'].
    save_reconstrction
        If True, the reconstructed expression profiles are saved in adata.layers['STAGATE_ReX'].

    Returns
    -------
    AnnData
    """

    tf.reset_default_graph()
    np.random.seed(random_seed)
    tf.set_random_seed(random_seed)
    if 'highly_variable' in adata1.var.columns and 'highly_variable' in adata2.var.columns:
        adata_Vars1 = adata1[:, adata1.var['highly_variable']]
        adata_Vars2 = adata2[:, adata2.var['highly_variable']]
    else:
        adata_Vars1 = adata1
        adata_Vars2 = adata2
    # import pandas as pd
    if isinstance(adata_Vars1.X, np.ndarray):
        X1 = pd.DataFrame(
            adata_Vars1.X[:, ], index=adata_Vars1.obs.index, columns=adata_Vars1.var.index)
    else:
        X1 = pd.DataFrame(adata_Vars1.X.toarray()[
                          :, ], index=adata_Vars1.obs.index, columns=adata_Vars1.var.index)
    if isinstance(adata_Vars2.X, np.ndarray):
        X2 = pd.DataFrame(
            adata_Vars2.X[:, ], index=adata_Vars2.obs.index, columns=adata_Vars2.var.index)
    else:
        X2 = pd.DataFrame(adata_Vars2.X.toarray()[
                          :, ], index=adata_Vars2.obs.index, columns=adata_Vars2.var.index)
#    X1 = pd.DataFrame(adata_Vars1.X.toarray()[:, ], index=adata_Vars1.obs.index, columns=adata_Vars1.var.index)
#    X2 = pd.DataFrame(adata_Vars2.X.toarray()[:, ], index=adata_Vars2.obs.index, columns=adata_Vars2.var.index)
    if verbose:
        print('Size of Input of rna: ', adata_Vars1.shape)
        print('Size of Input of atac: ', adata_Vars2.shape)
    cells = np.array(X1.index)
    cells_id_tran = dict(zip(cells, range(cells.shape[0])))
    genes = np.array(X1.columns)
    peaks = np.array(X2.columns)
    genes_id_tran = dict(zip(genes, range(genes.shape[0])))
    peaks_id_tran = dict(zip(peaks, range(peaks.shape[0])))
    if 'Spatial_Net' not in adata1.uns.keys():
        raise ValueError(
            "Spatial_Net is not existed! Run Cal_Spatial_Net first!")
    Spatial_Net = adata_Vars1.uns['Spatial_Net']
    G_df = Spatial_Net.copy()
    G_df['Cell1'] = G_df['Cell1'].map(cells_id_tran)
    G_df['Cell2'] = G_df['Cell2'].map(cells_id_tran)
    G = sp.coo_matrix((np.ones(G_df.shape[0]), (G_df['Cell1'], G_df['Cell2'])), shape=(
        adata_Vars1.n_obs, adata_Vars1.n_obs))
    G_tf = prepare_graph_data(G)

    if 'gene_peak_Net' not in adata1.uns.keys():
        raise ValueError(
            "gene_peak_Net is not existed! Run Cal_gene_peak_Net first!")
    gene_peak_Net = adata_Vars1.uns['gene_peak_Net']
    # gene_peak_Net.columns = ['Gene','Peak']
    if type is 'protein':
        gene_peak_Net.columns = ['Gene','Peak']
    G_gp_df = gene_peak_Net.copy()
    G_gp_df['Gene'] = G_gp_df['Gene'].map(genes_id_tran)
    G_gp_df['Peak'] = G_gp_df['Peak'].map(peaks_id_tran) + adata_Vars1.n_vars
    if type is 'ATAC':
        dist = G_gp_df['Distance']
        weights = np.concatenate((((dist + bp_width) / bp_width) **
                                 (-0.75), ((dist + bp_width) / bp_width) ** (-0.75)), axis=0)
        temp=-10
        G_gp = sp.coo_matrix((weights, (np.concatenate((G_gp_df['Gene'], G_gp_df['Peak']), axis=0), np.concatenate((G_gp_df['Peak'], G_gp_df['Gene']), axis=0))), shape=(adata_Vars1.n_vars+adata_Vars2.n_vars,
                                                                                                                                                                          adata_Vars1.n_vars+adata_Vars2.n_vars))
    elif type is 'ATAC_RNA':
        dist = G_gp_df['Distance']
        bp_width=2000
        weights = np.concatenate((((dist + bp_width) / bp_width) **
                                    (-0.75), ((dist + bp_width) / bp_width) ** (-0.75)), axis=0)
        G_gp = sp.coo_matrix((weights, (np.concatenate((G_gp_df['Gene'], G_gp_df['Peak']), axis=0), np.concatenate((G_gp_df['Peak'], G_gp_df['Gene']), axis=0))), shape=(adata_Vars1.n_vars+adata_Vars2.n_vars,
                                                                                                                                                                            adata_Vars1.n_vars+adata_Vars2.n_vars))
    else:
        G_gp = sp.coo_matrix((np.ones(G_gp_df.shape[0]*2)*protein_value, (np.concatenate((G_gp_df['Gene'], G_gp_df['Peak']), axis=0), np.concatenate((G_gp_df['Peak'], G_gp_df['Gene']), axis=0))), shape=(adata_Vars1.n_vars+adata_Vars2.n_vars,
                                                                                                                                                                                             adata_Vars1.n_vars+adata_Vars2.n_vars))
   
    G_gp_tf = prepare_graph_data(G_gp)


    trainer = MultiGATE(hidden_dims1=[X1.shape[1]] + hidden_dims, hidden_dims2=[X2.shape[1]] + hidden_dims, spot_num=X1.shape[0],
                        temp=temp, n_epochs=n_epochs, lr=lr, gradient_clipping=gradient_clipping,
                        nonlinear=nonlinear, weight_decay=weight_decay, verbose=verbose,
                        random_seed=random_seed)

    trainer(G_tf, G_tf, G_gp_tf, X1, X2)

    loss_df = pd.DataFrame({
        'loss': trainer.loss_list,
        'loss_atac': trainer.loss_list_atac,
        'loss_rna': trainer.loss_list_rna,
        'loss_clip': trainer.loss_list_clip,
        'weight_decay_loss': trainer.weight_decay_loss_list,
    })

    loss_df.to_csv('losses.csv', index=False)
    embeddings_RNA, embeddings_ATAC, attentions1, attentions2, attentions_gp, loss, ReX_RNA, ReX_ATAC = trainer.infer(G_tf, G_tf,
                                                                                                                      G_gp_tf, X1, X2)

    # embeddings, attentions, loss, ReX
    cell_reps1 = pd.DataFrame(embeddings_RNA)
    cell_reps1.index = cells

    adata1.obsm[key_added] = cell_reps1.loc[adata1.obs_names, ].values

    cell_reps2 = pd.DataFrame(embeddings_ATAC)
    cell_reps2.index = cells

    adata2.obsm[key_added] = cell_reps2.loc[adata2.obs_names, ].values
    norm2 = Normalizer(norm='l2')

    adata1.obsm[key_added + '_clip_all'] = (norm2.fit_transform(cell_reps1.loc[adata1.obs_names,].values) +
                                            norm2.fit_transform(cell_reps2.loc[adata2.obs_names, ].values)) / 2.0
    adata2.obsm[key_added + '_clip_all'] = (norm2.fit_transform(cell_reps1.loc[adata1.obs_names,].values) +
                                            norm2.fit_transform(cell_reps2.loc[adata2.obs_names, ].values)) / 2.0
    if save_attention:
        adata1.uns['MultiGATE_attention'] = attentions1
        adata2.uns['MultiGATE_attention'] = attentions2
        adata1.uns['MultiGATE_gene_peak_attention'] = attentions_gp
        adata2.uns['MultiGATE_gene_peak_attention'] = attentions_gp
    if save_loss:
        adata1.uns['MultiGATE_loss_united'] = loss
        adata2.uns['MultiGATE_loss_united'] = loss
        # adata1.uns['MultiGATE_loss'] = loss
    if save_reconstrction:
        ReX_RNA = pd.DataFrame(ReX_RNA, index=X1.index, columns=X1.columns)
        ReX_RNA[ReX_RNA < 0] = 0
        adata1.layers['MultiGATE_ReX'] = ReX_RNA.values

        ReX_ATAC = pd.DataFrame(ReX_ATAC, index=X2.index, columns=X2.columns)
        ReX_ATAC[ReX_ATAC < 0] = 0
        adata2.layers['MultiGATE_ReX'] = ReX_ATAC.values
    return adata1, adata2


def prune_spatial_Net(Graph_df, label):
    print('------Pruning the graph...')
    print('%d edges before pruning.' % Graph_df.shape[0])
    pro_labels_dict = dict(zip(list(label.index), label))
    Graph_df['Cell1_label'] = Graph_df['Cell1'].map(pro_labels_dict)
    Graph_df['Cell2_label'] = Graph_df['Cell2'].map(pro_labels_dict)
    Graph_df = Graph_df.loc[Graph_df['Cell1_label']
                            == Graph_df['Cell2_label'],]
    print('%d edges after pruning.' % Graph_df.shape[0])
    return Graph_df


def prepare_graph_data(adj):
    # adapted from preprocess_adj_bias
    num_nodes = adj.shape[0]
    adj = adj + sp.eye(num_nodes)  # self-loop
    # data =  adj.tocoo().data
    # adj[adj > 0.0] = 1.0
    if not sp.isspmatrix_coo(adj):
        adj = adj.tocoo()
    adj = adj.astype(np.float32)
    indices = np.vstack((adj.col, adj.row)).transpose()
    return (indices, adj.data, adj.shape)


def recovery_Imputed_Count(adata, size_factor):
    assert ('ReX' in adata.uns)
    temp_df = adata.uns['ReX'].copy()
    sf = size_factor.loc[temp_df.index]
    temp_df = np.expm1(temp_df)
    temp_df = (temp_df.T * sf).T
    adata.uns['ReX_Count'] = temp_df
    return adata
