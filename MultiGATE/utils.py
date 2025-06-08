import pandas as pd
import numpy as np
import re
import sklearn.neighbors
from sklearn.neighbors import NearestNeighbors
import scipy.sparse as sp
from sklearn.metrics import pairwise_distances
from sklearn.preprocessing import normalize
import networkx as nx
import re
from networkx.algorithms.bipartite import biadjacency_matrix   
import MultiGATE.genomics as genomics

def calculate_weight(rna_ratio, atac_ratio):
    try:
        with np.errstate(divide='raise', invalid='raise'):
            rna_exp = np.exp(rna_ratio)
            atac_exp = np.exp(atac_ratio)
            rna_weight = rna_exp / (rna_exp + atac_exp)
    except (ZeroDivisionError, FloatingPointError, RuntimeError):
        rna_weight = np.where(np.isinf(rna_exp), 1.0, 0.0)
    return rna_weight
def get_lowest_jaccard_indices(jaccard_sim, embedding, n_neighbors=15):
    num_obs = jaccard_sim.shape[0]
    lowest_indices = []

    for i in range(num_obs):
        row = jaccard_sim[i]
        nonzero_indices = np.nonzero(row)[0]

        if len(nonzero_indices) > n_neighbors:
            rna_diff = embedding[i] - embedding[nonzero_indices]
            norms = np.linalg.norm(rna_diff, axis=1)
            sort_indices = np.lexsort((norms, row[nonzero_indices]))  # Sort by norms, then row values
            lowest_indices.append(nonzero_indices[sort_indices[:n_neighbors]])
        else:
            lowest_indices.append(nonzero_indices)

    return lowest_indices
def calculate_kernel_bandwidth(embedding, lowest_jaccard_indices):
    num_obs = embedding.shape[0]
    euclidean_distances = pairwise_distances(embedding)
    kernel_bandwidths = []
    for i in range(num_obs):
        indices = lowest_jaccard_indices[i]
        avg_distance = np.mean(euclidean_distances[i, indices])
        kernel_bandwidths.append(avg_distance)
    return kernel_bandwidths
def run_wnn_analysis(rna, atac, n_neighbors = 15):
    rna_embedding = rna.obsm['STAGATE']  # Replace 'X_emb' with the actual embedding key in your RNA object
    atac_embedding = atac.obsm['STAGATE']
    rna_embedding = normalize(rna_embedding, norm='l2')  # Replace 'X_emb' with the actual embedding key in your RNA object
    atac_embedding = normalize(atac_embedding, norm='l2')
    # rna_avg_embedding = rna_neighbors.dot(rna_embedding) / rna_neighbors.sum(axis=1)
    nn = NearestNeighbors(n_neighbors=n_neighbors+1, algorithm='auto')
    nn.fit(rna_embedding)
    # Find the 15 nearest neighbors for each observation
    _, indices = nn.kneighbors(rna_embedding, n_neighbors=n_neighbors+1)
    # Exclude the observation itself from the neighbors
    neighbors_indices = indices[:, 1:]
    # Calculate the average RNA embedding using the neighbor indices
    rna_avg_embedding = np.mean(rna_embedding[neighbors_indices], axis=1)
    atac_cross_modality_avg_embedding = np.mean(atac_embedding[neighbors_indices], axis=1)
    n_observations = len(rna_embedding)
    rna_neighbors = sp.lil_matrix((n_observations, n_observations), dtype=int)

    # Set the entries in the neighbor matrix to 1 for the neighbors
    for i, neighbors in enumerate(neighbors_indices):
        rna_neighbors[i, neighbors] = 1

    _, indices = nn.kneighbors(atac_embedding, n_neighbors=n_neighbors+1)
    neighbors_indices = indices[:, 1:]
    atac_avg_embedding = np.mean(atac_embedding[neighbors_indices], axis=1)
    rna_cross_modality_avg_embedding  = np.mean(rna_embedding[neighbors_indices], axis=1)
    # Calculate cross-modality average embeddings
    n_observations = len(rna_embedding)
    atac_neighbors = sp.lil_matrix((n_observations, n_observations), dtype=int)

    # Set the entries in the neighbor matrix to 1 for the neighbors
    for i, neighbors in enumerate(neighbors_indices):
        atac_neighbors[i, neighbors] = 1

    # Assuming rna_neighbors is a sparse matrix
    rna_jaccard_sim = pairwise_distances(rna_neighbors.toarray(), metric='jaccard')

    # Calculate the Jaccard similarity matrix for ATAC
    atac_jaccard_sim = pairwise_distances(atac_neighbors.toarray(), metric='jaccard')
    # Get the indices of the 15 observations with the lowest non-zero Jaccard similarity for RNA
    rna_lowest_jaccard_indices = get_lowest_jaccard_indices(rna_avg_embedding,rna_jaccard_sim,n_neighbors=n_neighbors)

    # Get the indices of the 15 observations with the lowest non-zero Jaccard similarity for ATAC
    atac_lowest_jaccard_indices = get_lowest_jaccard_indices(atac_avg_embedding,atac_jaccard_sim,n_neighbors=n_neighbors)
    # Calculate the observation-specific kernel bandwidths for RNA
    rna_kernel_bandwidths = calculate_kernel_bandwidth(rna_embedding, rna_lowest_jaccard_indices)
    # Calculate the observation-specific kernel bandwidths for ATAC
    atac_kernel_bandwidths = calculate_kernel_bandwidth(atac_embedding, atac_lowest_jaccard_indices)

    # Assign the observation-specific kernel bandwidths to the respective AnnData objects
    rna.obs['kernel_bandwidth'] = rna_kernel_bandwidths
    atac.obs['kernel_bandwidth'] = atac_kernel_bandwidths
    rna_distances = pairwise_distances(rna_embedding)
    np.fill_diagonal(rna_distances, np.inf)
    closest_indices_rna = np.argmin(rna_distances, axis=1)
    closest_distances_rna = np.linalg.norm(rna_embedding - rna_embedding[closest_indices_rna], axis=1)
    rna_theta_within = np.exp(-np.maximum(np.linalg.norm(rna_embedding - rna_avg_embedding, axis=1)-closest_distances_rna,0)
                        /(rna_kernel_bandwidths - closest_distances_rna))

    rna_theta_cross = np.exp(-np.maximum(np.linalg.norm(rna_embedding - rna_cross_modality_avg_embedding, axis=1)-closest_distances_rna,0)
                        /(rna_kernel_bandwidths -  closest_distances_rna))


    atac_distances = pairwise_distances(atac_embedding)
    np.fill_diagonal(atac_distances, np.inf)
    closest_indices_atac = np.argmin(atac_distances, axis=1)
    closest_distances_atac = np.linalg.norm(atac_embedding - atac_embedding[closest_indices_atac], axis=1)
    atac_theta_within = np.exp(-np.maximum(np.linalg.norm(atac_embedding - atac_avg_embedding, axis=1)-closest_distances_atac,0)
                        /(atac_kernel_bandwidths  - closest_distances_atac))


    atac_theta_cross = np.exp(-np.maximum(np.linalg.norm(atac_embedding - atac_cross_modality_avg_embedding, axis=1)-closest_distances_atac,0)
                        /(atac_kernel_bandwidths -  closest_distances_atac))

    epsilon = 1e-4

    rna_ratio = rna_theta_within /(rna_theta_cross+epsilon)
    atac_ratio = atac_theta_within / (atac_theta_cross + epsilon)

    rna_weight = calculate_weight(rna_ratio, atac_ratio)
    atac_weight = 1- rna_weight

    num_rna = len(rna_embedding)
    num_atac = len(atac_embedding)

    # Initialize the matrix with zeros
    wnn_mat = np.zeros((num_rna, num_atac))
    # Iterate over the indices i and j
    for i in range(num_rna):
        for j in range(num_atac):
            # Calculate theta_rna and theta_atac based on rna_embedding and atac_embedding
            theta_rna = np.exp(-np.maximum(np.linalg.norm(rna_embedding[i] -rna_embedding[j])-closest_distances_rna[i],0)
                        /(rna_kernel_bandwidths[i] - closest_distances_rna[i]))
            theta_atac = np.exp(-np.maximum(np.linalg.norm(atac_embedding[i] -atac_embedding[j])-closest_distances_atac[i],0)
                        /(atac_kernel_bandwidths[i] - closest_distances_atac[i]))
            
            # Calculate the entry (i, j) based on the given formula
            wnn_mat[i, j] = rna_weight[i] * theta_rna + atac_weight[i] * theta_atac

    wnn_mat = np.nan_to_num(wnn_mat, nan=1)
    rna.obsp['wnn'] = wnn_mat
def Cal_gene_peak_Net_new(rna, atac, range=2000, 
                          file = "/lustre/project/Stat/s1155077016/Spatial_multi-omics_data/P22/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz", verbose=True):
    var_by=rna.var_names
    gtf = file
    gtf_by = "gene_name"
    COLUMNS = pd.Index(["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])

    loaded = pd.read_csv(gtf, sep="\t", header=None, comment="#")
    loaded.columns = COLUMNS[:loaded.shape[1]]
    loaded = loaded.query("feature == 'gene'")

    # ?????????
    pattern = re.compile(r'([^\s]+) "([^"]+)";')
    splitted = pd.DataFrame.from_records(np.vectorize(lambda x: {
        key: val for key, val in pattern.findall(x)
    })(loaded["attribute"]), index=loaded.index)
    loaded = loaded.assign(**splitted)
    loaded = loaded.sort_values("seqname").drop_duplicates(subset=[gtf_by], keep="last")

    # ????? BED ?? DataFrame
    bed_df = pd.DataFrame(loaded, copy=True).loc[:, ("seqname", "start", "end", "score", "strand")]
    bed_df.insert(3, "name", np.repeat(".", len(bed_df)) if gtf_by is None else loaded[gtf_by])
    bed_df["start"] -= 1  # ??? 0-based
    bed_df.columns = ("chrom", "chromStart", "chromEnd", "name", "score", "strand")

    # ??????
    for item in bed_df.columns:
        if item in ("chromStart", "chromEnd"):
            bed_df[item] = bed_df[item].astype(int)
        else:
            bed_df[item] = bed_df[item].astype(str)

    # ??????? strand ??
    gene_df = pd.concat([
        pd.DataFrame(bed_df),
        pd.DataFrame(loaded).drop(columns=COLUMNS)  # ?????????
    ], axis=1).set_index(gtf_by).reindex(var_by)

    # ???? strand ?
    gene_df = gene_df[['chrom', 'chromStart', 'chromEnd', 'strand']]
    gene_df.dropna(inplace=True)
    rna=rna[:,gene_df.index]
    rna.var['chrom'] = gene_df.loc[rna.var_names, 'chrom']
    rna.var['chromStart'] = gene_df.loc[rna.var_names, 'chromStart']
    rna.var['chromEnd'] = gene_df.loc[rna.var_names, 'chromEnd']
    rna.var['strand'] = gene_df.loc[rna.var_names, 'strand']
    # ?????????????????,????DataFrame
    chrom_info = atac.var_names.str.split(r"[:\-]", expand=True).to_frame(index=False)
    chrom_info.columns = ['chrom', 'chromStart', 'chromEnd']

    # ???????????????????
    chrom_info['chromStart'] = chrom_info['chromStart'].astype(int)
    chrom_info['chromEnd'] = chrom_info['chromEnd'].astype(int)
    # ? atac.var_names ??? chrom_info ???
    chrom_info.index = atac.var_names

    # ?? chrom_info ???? atac.var ?????
    if chrom_info.index.equals(atac.var.index):
        # ? chrom_info ?????? atac.var ?
        atac.var[['chrom', 'chromStart', 'chromEnd']] = chrom_info
    else:
        print("?????,?????????!")
    rna=rna[:,rna.var['highly_variable']]
    # genes = Bed(rna.var.assign(name=rna.var_names))
    # peaks = Bed(atac.var.assign(name=atac.var_names))
    # tss = genes.strand_specific_start_site()
    # promoters = tss.expand(2000, 0)
    genes = genomics.Bed(rna.var.assign(name=rna.var_names))
    peaks = genomics.Bed(atac.var.assign(name=atac.var_names))
    tss = genes.strand_specific_start_site()
    promoters = tss.expand(2000, 0)
    dist_graph = genomics.window_graph(
        promoters, peaks, range,
        attr_fn=lambda l, r, d: { 
            "dist": abs(d),
            "weight": genomics.dist_power_decay(abs(d)),
            "type": "dist"
        }
    )
    dist_graph = nx.DiGraph(dist_graph)
    dist_graph.number_of_edges()
    dist = biadjacency_matrix(dist_graph, genes.index, peaks.index, weight="dist", dtype=np.float32).tocoo()
    df=pd.DataFrame({
        # 'gene':dist.row,
        # 'peak':dist.col,
        'Distance':dist.data
    })
    non_zero_rows = dist.row
    non_zero_cols = dist.col
    df['Gene'] = [genes.index[row] for row in non_zero_rows]
    df['Peak'] = [peaks.index[col] for col in non_zero_cols]
    if verbose:
        print('The graph contains %d edges, %d genes.' % (dist_graph.number_of_edges(), len(set(df['Gene']))))
        print('%.4f peaks per gene on average.' % (dist_graph.number_of_edges() / len(set(df['Gene']))))
    atac.uns['gene_peak_Net'] = df
    rna.uns['gene_peak_Net'] = df

def Cal_gene_peak_Net(rna, atac, range=2000, states = ['Top'], # 'Bottom' 'Overlap'
                      file = "/lustre/project/Stat/s1155077016/Spatial_multi-omics_data/P22/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz", verbose=True):
    """\
    Construct the spatial neighbor networks.

    Parameters
    ----------
    adata
        AnnData object of scanpy package.
    rad_cutoff
        radius cutoff when model='Radius'
    k_cutoff
        The number of nearest neighbors when model='KNN'
    model
        The network construction model. When model=='Radius', the spot is connected to spots whose distance is less than rad_cutoff. When model=='KNN', the spot is connected to its first k_cutoff nearest neighbors.

    Returns
    -------
    The spatial networks are saved in adata.uns['Spatial_Net']
    """
    if verbose:
        print('------Calculating gene-peak graph...')
    split = atac.var_names.str.split(r"[:\-]")
    atac.var["chrom"] = split.map(lambda x: x[0])
    atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
    atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)
    atac_Vars = atac[:, atac.var['highly_variable']]

    var_by = rna.var_names
    gtf = file
    gtf_by = "gene_name"
    COLUMNS = pd.Index(["seqname", "source", "feature", "start", "end", "score", "strand",
                        "frame", "attribute"])

    loaded = pd.read_csv(gtf, sep="\t", header=None, comment="#")
    loaded.columns = COLUMNS[:loaded.shape[1]]
    loaded = loaded.query("feature == 'gene'")

    pattern = re.compile(r'([^\s]+) "([^"]+)";')
    splitted = pd.DataFrame.from_records(np.vectorize(lambda x: {
        key: val for key, val in pattern.findall(x)
    })(loaded["attribute"]), index=loaded.index)
    loaded = loaded.assign(**splitted)
    loaded = loaded.sort_values("seqname").drop_duplicates(subset=[gtf_by], keep="last")
    bed_df = pd.DataFrame(loaded, copy=True).loc[
             :, ("seqname", "start", "end", "score", "strand")
             ]
    bed_df.insert(3, "name", np.repeat(
        ".", len(bed_df)
    ) if gtf_by is None else loaded[gtf_by])
    bed_df["start"] -= 1  # Convert to zero-based
    bed_df.columns = ("chrom", "chromStart", "chromEnd", "name", "score", "strand")
    for item in bed_df.columns:
        if item in bed_df:
            if item in ("chromStart", "chromEnd"):
                bed_df[item] = bed_df[item].astype(int)
            else:
                bed_df[item] = bed_df[item].astype(str)
        elif item not in ("chrom", "chromStart", "chromEnd"):
            bed_df[item] = "."

    merge_df = pd.concat([
        pd.DataFrame(bed_df),
        pd.DataFrame(loaded).drop(columns=COLUMNS)  # Only use the splitted attributes
    ], axis=1).set_index(gtf_by).reindex(var_by).set_index(rna.var.index)
    rna.var = rna.var.assign(**merge_df)

    chrom_mask = ~(rna.var.loc[:, ["chromStart", "chromEnd"]].isin([np.nan, np.inf]).all(axis=1))
    rna = rna[:, chrom_mask]

    rna_copy = rna.copy()
    # Modify the 'chromStart' column of the 'var' attribute
    rna_copy.var["chromStart"] = rna_copy.var["chromStart"].astype(int)
    rna_copy.var["chromEnd"] = rna_copy.var["chromEnd"].astype(int)
    # rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].head()
    rna_Vars = rna_copy[:, rna_copy.var['highly_variable']]


    gene_ranges = rna_Vars.var.loc[:, ["chrom", "chromStart", "chromEnd"]]
    gene_ranges = gene_ranges.rename_axis('Gene')
    gene_ranges = gene_ranges.reset_index()
    peak_ranges = atac_Vars.var.loc[:, ["chrom", "chromStart", "chromEnd"]]
    peak_ranges = peak_ranges.rename_axis('Peak')
    peak_ranges = peak_ranges.reset_index()
    merged_df = gene_ranges.merge(peak_ranges, on='chrom', suffixes=('_gene', '_peak'))

    # Calculate the distance between gene and peak using vectorized operations
    gene_start = merged_df['chromStart_gene'].values
    gene_end = merged_df['chromEnd_gene'].values
    peak_start = merged_df['chromStart_peak'].values
    peak_end = merged_df['chromEnd_peak'].values

#    distances = np.minimum(
#        np.abs(gene_start - peak_end),
#        np.abs(peak_start - gene_end)
#    )
    distances = np.abs(gene_start - peak_start)
    merged_df['Distance'] = distances

    # Determine the strand using vectorized operations
    overlap_condition = (
        (peak_start <= gene_end) & (peak_start>= gene_start) |
        (peak_end >= gene_start) & (peak_end <= gene_end)
    )

    top_condition = (peak_start < gene_start) & (peak_end < gene_start)

    merged_df['strand'] = np.select(
        [overlap_condition, top_condition],
        ['Overlap', 'Top'],
        default='Bottom'
    )

    # Filter the results based on distance and strand conditions
    results_df = merged_df[(merged_df['Distance'] <= range) & (merged_df['strand'].isin(states))]
    results_df = results_df[['Gene', 'Peak', 'Distance']].reset_index(drop=True)

    gene_peak_Net = results_df.copy()
    if verbose:
        print('The graph contains %d edges, %d genes.' % (gene_peak_Net.shape[0], rna.n_vars))
        print('%.4f peaks per gene on average.' % (gene_peak_Net.shape[0] / rna.n_vars))
    rna.uns['gene_peak_Net'] = gene_peak_Net
    atac.uns['gene_peak_Net'] = gene_peak_Net

def Cal_gene_protein_Net(rna, pro, verbose=True):
    """\
    Construct the spatial neighbor networks.

    Parameters
    ----------
    adata
        AnnData object of scanpy package.
    rad_cutoff
        radius cutoff when model='Radius'
    k_cutoff
        The number of nearest neighbors when model='KNN'
    model
        The network construction model. When model=='Radius', the spot is connected to spots whose distance is less than rad_cutoff. When model=='KNN', the spot is connected to its first k_cutoff nearest neighbors.

    Returns
    -------
    The spatial networks are saved in adata.uns['Spatial_Net']
    """
    if verbose:
        print('------Calculating gene-protein graph...')
    adata1_var_names = rna.var_names
    adata2_var_names = pro.var_names
    matching_pairs = []
    # Find the matched pairs of var_names
    matched_pairs = [(name1, name2) for name1 in adata1_var_names for name2 in adata2_var_names if name1.capitalize() == name2.capitalize()]
    df = pd.DataFrame(matched_pairs, columns=["Gene", "Peak"])
    rna.uns['gene_peak_Net'] = df
    pro.uns['gene_peak_Net'] = df
    if verbose:
        print('------Calculation finished...')

def Cal_Spatial_Net(adata, rad_cutoff=None, k_cutoff=None, model='Radius', verbose=True):
    """\
    Construct the spatial neighbor networks.

    Parameters
    ----------
    adata
        AnnData object of scanpy package.
    rad_cutoff
        radius cutoff when model='Radius'
    k_cutoff
        The number of nearest neighbors when model='KNN'
    model
        The network construction model. When model=='Radius', the spot is connected to spots whose distance is less than rad_cutoff. When model=='KNN', the spot is connected to its first k_cutoff nearest neighbors.
    
    Returns
    -------
    The spatial networks are saved in adata.uns['Spatial_Net']
    """

    assert(model in ['Radius', 'KNN'])
    if verbose:
        print('------Calculating spatial graph...')
    coor = pd.DataFrame(adata.obsm['spatial'])
    coor.index = adata.obs.index
    coor.columns = ['imagerow', 'imagecol']

    if model == 'Radius':
        nbrs = sklearn.neighbors.NearestNeighbors(radius=rad_cutoff).fit(coor)
        distances, indices = nbrs.radius_neighbors(coor, return_distance=True)
        KNN_list = []
        for it in range(indices.shape[0]):
            KNN_list.append(pd.DataFrame(zip([it]*indices[it].shape[0], indices[it], distances[it])))
    
    if model == 'KNN':
        nbrs = sklearn.neighbors.NearestNeighbors(n_neighbors=k_cutoff+1).fit(coor)
        distances, indices = nbrs.kneighbors(coor)
        KNN_list = []
        for it in range(indices.shape[0]):
            KNN_list.append(pd.DataFrame(zip([it]*indices.shape[1],indices[it,:], distances[it,:])))

    KNN_df = pd.concat(KNN_list)
    KNN_df.columns = ['Cell1', 'Cell2', 'Distance']

    Spatial_Net = KNN_df.copy()
    Spatial_Net = Spatial_Net.loc[Spatial_Net['Distance']>0,]
    id_cell_trans = dict(zip(range(coor.shape[0]), np.array(coor.index), ))
    Spatial_Net['Cell1'] = Spatial_Net['Cell1'].map(id_cell_trans)
    Spatial_Net['Cell2'] = Spatial_Net['Cell2'].map(id_cell_trans)
    if verbose:
        print('The graph contains %d edges, %d cells.' %(Spatial_Net.shape[0], adata.n_obs))
        print('%.4f neighbors per cell on average.' %(Spatial_Net.shape[0]/adata.n_obs))

    adata.uns['Spatial_Net'] = Spatial_Net


def Cal_Spatial_Net_3D(adata, rad_cutoff_2D, rad_cutoff_Zaxis,
                       key_section='Section_id', section_order=None, verbose=True):
    """\
    Construct the spatial neighbor networks.

    Parameters
    ----------
    adata
        AnnData object of scanpy package.
    rad_cutoff_2D
        radius cutoff for 2D SNN construction.
    rad_cutoff_Zaxis
        radius cutoff for 2D SNN construction for consturcting SNNs between adjacent sections.
    key_section
        The columns names of section_ID in adata.obs.
    section_order
        The order of sections. The SNNs between adjacent sections are constructed according to this order.
    
    Returns
    -------
    The 3D spatial networks are saved in adata.uns['Spatial_Net'].
    """
    adata.uns['Spatial_Net_2D'] = pd.DataFrame()
    adata.uns['Spatial_Net_Zaxis'] = pd.DataFrame()
    num_section = np.unique(adata.obs[key_section]).shape[0]
    if verbose:
        print('Radius used for 2D SNN:', rad_cutoff_2D)
        print('Radius used for SNN between sections:', rad_cutoff_Zaxis)
    for temp_section in np.unique(adata.obs[key_section]):
        if verbose:
            print('------Calculating 2D SNN of section ', temp_section)
        temp_adata = adata[adata.obs[key_section] == temp_section, ]
        Cal_Spatial_Net(
            temp_adata, rad_cutoff=rad_cutoff_2D, verbose=False)
        temp_adata.uns['Spatial_Net']['SNN'] = temp_section
        if verbose:
            print('This graph contains %d edges, %d cells.' %
                  (temp_adata.uns['Spatial_Net'].shape[0], temp_adata.n_obs))
            print('%.4f neighbors per cell on average.' %
                  (temp_adata.uns['Spatial_Net'].shape[0]/temp_adata.n_obs))
        adata.uns['Spatial_Net_2D'] = pd.concat(
            [adata.uns['Spatial_Net_2D'], temp_adata.uns['Spatial_Net']])
    for it in range(num_section-1):
        section_1 = section_order[it]
        section_2 = section_order[it+1]
        if verbose:
            print('------Calculating SNN between adjacent section %s and %s.' %
                  (section_1, section_2))
        Z_Net_ID = section_1+'-'+section_2
        temp_adata = adata[adata.obs[key_section].isin(
            [section_1, section_2]), ]
        Cal_Spatial_Net(
            temp_adata, rad_cutoff=rad_cutoff_Zaxis, verbose=False)
        spot_section_trans = dict(
            zip(temp_adata.obs.index, temp_adata.obs[key_section]))
        temp_adata.uns['Spatial_Net']['Section_id_1'] = temp_adata.uns['Spatial_Net']['Cell1'].map(
            spot_section_trans)
        temp_adata.uns['Spatial_Net']['Section_id_2'] = temp_adata.uns['Spatial_Net']['Cell2'].map(
            spot_section_trans)
        used_edge = temp_adata.uns['Spatial_Net'].apply(
            lambda x: x['Section_id_1'] != x['Section_id_2'], axis=1)
        temp_adata.uns['Spatial_Net'] = temp_adata.uns['Spatial_Net'].loc[used_edge, ]
        temp_adata.uns['Spatial_Net'] = temp_adata.uns['Spatial_Net'].loc[:, [
            'Cell1', 'Cell2', 'Distance']]
        temp_adata.uns['Spatial_Net']['SNN'] = Z_Net_ID
        if verbose:
            print('This graph contains %d edges, %d cells.' %
                  (temp_adata.uns['Spatial_Net'].shape[0], temp_adata.n_obs))
            print('%.4f neighbors per cell on average.' %
                  (temp_adata.uns['Spatial_Net'].shape[0]/temp_adata.n_obs))
        adata.uns['Spatial_Net_Zaxis'] = pd.concat(
            [adata.uns['Spatial_Net_Zaxis'], temp_adata.uns['Spatial_Net']])
    adata.uns['Spatial_Net'] = pd.concat(
        [adata.uns['Spatial_Net_2D'], adata.uns['Spatial_Net_Zaxis']])
    if verbose:
        print('3D SNN contains %d edges, %d cells.' %
            (adata.uns['Spatial_Net'].shape[0], adata.n_obs))
        print('%.4f neighbors per cell on average.' %
            (adata.uns['Spatial_Net'].shape[0]/adata.n_obs))

def Stats_Spatial_Net(adata):
    import matplotlib.pyplot as plt
    Num_edge = adata.uns['Spatial_Net']['Cell1'].shape[0]
    Mean_edge = Num_edge/adata.shape[0]
    plot_df = pd.value_counts(pd.value_counts(adata.uns['Spatial_Net']['Cell1']))
    plot_df = plot_df/adata.shape[0]
    fig, ax = plt.subplots(figsize=[3,2])
    plt.ylabel('Percentage')
    plt.xlabel('')
    plt.title('Number of Neighbors (Mean=%.2f)'%Mean_edge)
    ax.bar(plot_df.index, plot_df)
    return Mean_edge

def mclust_R(adata, num_cluster, modelNames='EEE', used_obsm='STAGATE', random_seed=2020):
    """\
    Clustering using the mclust algorithm.
    The parameters are the same as those in the R package mclust.
    """
    
    np.random.seed(random_seed)
    import rpy2.robjects as robjects
    robjects.r.library("mclust")

    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    r_random_seed = robjects.r['set.seed']
    r_random_seed(random_seed)
    rmclust = robjects.r['Mclust']

    res = rmclust(rpy2.robjects.numpy2ri.numpy2rpy(adata.obsm[used_obsm]), num_cluster, modelNames)
    mclust_res = np.array(res[-2])

    adata.obs['mclust'] = mclust_res
    adata.obs['mclust'] = adata.obs['mclust'].astype('int')
    adata.obs['mclust'] = adata.obs['mclust'].astype('category')
    return adata
def wnn_R(adata1, adata2, res=0.5, algo =  1, used_obsm='MultiGATE'):
    """
    Clustering using the mclust algorithm.
    The parameters are the same as those in the R package mclust.
    """
    import rpy2.robjects as ro
    ro.r.library("Seurat")
    ro.r.library("ggplot2")
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    import rpy2.robjects.pandas2ri
    rpy2.robjects.pandas2ri.activate()
    data1 = adata1.obsm[used_obsm]
    # row_names1 = adata1.obs.index

    data2 = adata2.obsm[used_obsm]



    rna_matrix =  rpy2.robjects.numpy2ri.numpy2rpy(data1)#robjects.r['as.matrix'](rna_df)
    atac_matrix = rpy2.robjects.numpy2ri.numpy2rpy(data2)
    #pro_matrix = rpy2.robjects.numpy2ri.numpy2rpy(data3)
    barcodes = rpy2.robjects.pandas2ri.py2rpy(adata1.obs.index)


    ro.r.assign("rna_matrix", rna_matrix)
    ro.r.assign("atac_matrix", atac_matrix)
    ro.r.assign("barcodes", barcodes)
    ro.r.assign("res", res)
    ro.r.assign("algo", algo)

    ro.r(
    '''
    colnames(rna_matrix) <- paste0("pca_", seq_len(ncol(rna_matrix)))
    colnames(atac_matrix) <- paste0("lsi_", seq_len(ncol(atac_matrix)))
    rownames(rna_matrix) <- barcodes

    #colnames(pro_matrix) <- pronames
    #rownames(pro_matrix) <- barcodes
    #Pro <- CreateSeuratObject(counts = t(pro_matrix))

    P22 <- CreateSeuratObject(counts = t(rna_matrix))
    rownames(atac_matrix) <- barcodes
    atac_obj <- CreateAssayObject(counts = t(atac_matrix))
    P22[["ATAC"]] <- atac_obj
    P22[["pca"]] = CreateDimReducObject(embeddings = rna_matrix, key = "pca_", assay = "RNA")
    P22[["lsi"]] = CreateDimReducObject(embeddings = atac_matrix, key = "lsi_", assay = "ATAC")
    P22 <- FindMultiModalNeighbors(
      P22,
      reduction.list = list("pca", "lsi"),
      dims.list      = list(seq_len(ncol(rna_matrix)), seq_len(ncol(atac_matrix)))
    )
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
    return adata1, adata2
def wnn_protein(adata1, adata2, res=0.5, algo =  1, used_obsm='MultiGATE'):
    """
    Clustering using the mclust algorithm.
    The parameters are the same as those in the R package mclust.
    """
    import rpy2.robjects as ro
    ro.r.library("Seurat")
    ro.r.library("ggplot2")
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    import rpy2.robjects.pandas2ri
    rpy2.robjects.pandas2ri.activate()
    data1 = adata1.obsm[used_obsm]
    # row_names1 = adata1.obs.index

    data2 = adata2.obsm[used_obsm]

    data3 = adata2.X
    # row_names2 = adata2.obs.index
    avg_label = adata1.obs['mclust'].astype('int')
    avg_label =  rpy2.robjects.numpy2ri.numpy2rpy(avg_label)

    rna_matrix =  rpy2.robjects.numpy2ri.numpy2rpy(data1)#robjects.r['as.matrix'](rna_df)
    atac_matrix = rpy2.robjects.numpy2ri.numpy2rpy(data2)
    pro_matrix = rpy2.robjects.numpy2ri.numpy2rpy(data3)
    barcodes = rpy2.robjects.pandas2ri.py2rpy(adata1.obs.index)

    #pro_matrix = rpy2.robjects.numpy2ri.numpy2rpy(data3)
    pronames = rpy2.robjects.pandas2ri.py2rpy(adata2.var_names)

    ro.r.assign("rna_matrix", rna_matrix)
    ro.r.assign("atac_matrix", atac_matrix)
    ro.r.assign("barcodes", barcodes)
    ro.r.assign("res", res)
    ro.r.assign("algo", algo)

    ro.r.assign("pronames", pronames)
    ro.r.assign("pro_matrix", pro_matrix)
    ro.r.assign("avg_label", avg_label)

    ro.r(
    '''
    colnames(rna_matrix) <- paste0("pca_", 1:30)
    colnames(atac_matrix) <- paste0("lsi_", 1:30)
    rownames(rna_matrix) <- barcodes

    colnames(pro_matrix) <- pronames
    rownames(pro_matrix) <- barcodes
    Pro <- CreateSeuratObject(counts = t(pro_matrix))

    P22 <- CreateSeuratObject(counts = t(rna_matrix))
    rownames(atac_matrix) <- barcodes
    atac_obj <- CreateAssayObject(counts = t(atac_matrix))
    P22[["ATAC"]] <- atac_obj
    P22[["pca"]] = CreateDimReducObject(embeddings = rna_matrix, key = "pca_", assay = "RNA")
    P22[["lsi"]] = CreateDimReducObject(embeddings = atac_matrix, key = "lsi_", assay = "ATAC")
    P22 <- FindMultiModalNeighbors(P22, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 1:30))
    P22 <- RunUMAP(P22, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    P22 <- FindClusters(P22, graph.name = "wsnn", algorithm = algo, resolution = res,verbose = FALSE)
    clusters = P22@meta.data$seurat_clusters 
    wwnumap = P22[["wnn.umap"]]
    umap = Embeddings(object = wwnumap)#[1:nrow(rna_matrix), 1:2]##wwnumap[1:nrow(rna_matrix), 1:2]#P22[["wnn.umap"]]  wwnumap[[1:nrow(rna_matrix), 1:2]] obj@reductions$umap
    VariableFeatures(Pro) = rownames(Pro@assays$RNA@counts)    
    #Pro@active.ident = as.factor(clusters)
    #heatmap = DoHeatmap(Pro, slot = "data", disp.max=1.5, size = 3) 
    #ggsave("/lustre/project/Stat/s1155077016/Spatial_RNA_Protein/mouse_spleen/heatmap_wnn.png", heatmap)
    Pro@active.ident = as.factor(avg_label)
    heatmap = DoHeatmap(Pro, slot = "data", disp.max=1.5, size = 3) 
    ggsave("/lustre/project/Stat/s1155077016/Spatial_RNA_Protein/mouse_spleen/heatmap_mclust_0311.png", heatmap)
    saveRDS(Pro, file = "/lustre/project/Stat/s1155077016/Spatial_RNA_Protein/mouse_spleen/pro_seurat_0311.rds")   
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
    return adata1, adata2