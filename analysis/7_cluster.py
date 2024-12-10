########################## file information #############################
# Author: Pengfei Ren
# Date: 2024-12-04 21:43:44
# LastEditTime: 2024-12-04 23:37:41
# Description: 
# FilePath: /undefined/Users/morsouron/Desktop/Project/benchmark/code/7_cluster.py
#########################################################################

def filter_low_quality(adata,unit_type):
    sc.pp.filter_genes(adata, min_cells=100)
    adata = qc(adata)
    if unit_type == 'cell': 
        lower = adata.obs['total_counts'].quantile(0.10)
    else:
        lower = adata.obs['total_counts'].quantile(0.20)
    filtered_adata = adata[adata.obs['total_counts'] >= lower]
    filtered_adata = qc(filtered_adata)
    return filtered_adata

def cluster(adata):
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=int(adata.shape[1]/10),flavor='seurat')
    sc.tl.pca(adata, n_comps = 30)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    silhouette_scores = silhouette_score(adata.obsm['X_pca'], adata.obs['cluster'])
    return(adata)