########################## file information #############################
# Author: Pengfei Ren
# Date: 2024-12-04 21:53:08
# LastEditTime: 2024-12-04 23:34:29
# Description: 
# FilePath: /undefined/Users/morsouron/Desktop/Project/benchmark/code/9_spatial_cluster.py
#########################################################################

# CODEX data
def spatial_cluster(adata,adata_path,coord_system):
    adata = adata[adata.obs['codex_common'] == 1,:]
    adata.obs['sample'] = 'codex'
    trvae = sca.models.TRVAE(
        adata=adata,
        condition_key='sample',
        conditions=['codex'],
        hidden_layer_sizes=[128, 128],
    )
    adata.X = adata.X.astype(np.float32)
    trvae.train(
        n_epochs=trvae_epochs,
        alpha_epoch_anneal=200,
        early_stopping_kwargs=early_stopping_kwargs
    )
    adata.obsm['X_trVAE'] = trvae.get_latent().astype(np.float32)
    sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True,spatial_key=coord_system)
    cc.gr.remove_long_links(adata)
    cc.gr.aggregate_neighbors(adata, n_layers=5, use_rep='X_trVAE')
    autok = cc.tl.ClusterAutoK(
        n_clusters=(5,10), 
        max_runs=10, 
        model_params=dict(
            random_state=12345,
            # If running on GPU
            trainer_params=dict(accelerator='gpu', devices=1)
        )
    )
    autok.fit(adata, use_rep='X_cellcharter')
    gmm = cc.tl.Cluster(
        n_clusters=autok.best_k, 
        random_state=12345,
        # If running on GPU
        #trainer_params=dict(accelerator='gpu', devices=1)
    )
    gmm.fit(adata, use_rep='X_cellcharter')
    adata.obs['spatial_cluster'] = gmm.predict(adata, use_rep='X_cellcharter')
    return([autok.best_k,tissue,platform])

# ST data
def spatial_cluster(adata,adata_path,coord_system):
    adata = adata[adata.obs['codex_common'] == 1,:]
    adata.obs['sample'] = 'transcriptome'
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    scvi.settings.seed = 12345
    scvi.model.SCVI.setup_anndata(
        adata, 
        layer="counts", 
        batch_key='sample'
    )
    model = scvi.model.SCVI(adata)
    model.train(early_stopping=True, enable_progress_bar=True)
    adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)
    sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True, spatial_key=coord_system)
    cc.gr.remove_long_links(adata)
    cc.gr.aggregate_neighbors(adata, n_layers=5, use_rep='X_scVI', out_key='X_cellcharter', sample_key='sample')
    
    platform = os.path.basename(os.path.dirname(adata_path))
    tissue = os.path.basename(adata_path).split('_')[0]
    unit = os.path.basename(adata_path).replace('_filter.h5ad','')
    
    k = best_ks.loc[(best_ks['tissue'] == tissue) & (best_ks['platform'] == platform),'k'].tolist()[0]
    
    gmm = cc.tl.Cluster(
        n_clusters=k, 
        random_state=12345,
        # If running on GPU
        trainer_params=dict(accelerator='gpu', devices=1)
    )
    gmm.fit(adata, use_rep='X_cellcharter')
    adata.obs['spatial_cluster'] = gmm.predict(adata, use_rep='X_cellcharter')

# correlation between CODEX and ST

# calculate the cell proportion of each cluster in each grid
def calculate_cluster_proportion(data, grid_coords):
    cluster_counts = {}
    for i, coords in enumerate(grid_coords):
        key = (coords[0], coords[1])  # grid index
        label = data[i, 2]  # cluster label
        if key not in cluster_counts:
            cluster_counts[key] = {}
        if label not in cluster_counts[key]:
            cluster_counts[key][label] = 0
        cluster_counts[key][label] += 1

    cluster_proportions = {}
    for key, clusters in cluster_counts.items():
        total_cells = sum(clusters.values())
        cluster_proportions[key] = {label: count / total_cells for label, count in clusters.items()}

    return cluster_proportions

# generate a cluster proportion vector and compare only the common grids in the two slices
def create_proportion_vectors(proportions, common_grids, num_clusters):
    vectors = {label: [] for label in range(num_clusters)}
    
    for grid in common_grids:
        for label in vectors.keys():
            if label in proportions.get(grid, {}):
                vectors[label].append(proportions[grid][label])
            else:
                vectors[label].append(0)
    
    return vectors

for tissue_idx, tissue in enumerate(tissues):
    count = 1
    for kind_idx, kind in enumerate(kinds):
        trans_np, codex_np = load_data_np(tissue, kind)

        data1 = trans_np
        data2 = codex_np

        # divide the coordinate by 100 and round it to get the grid index.
        grid_coords1 = np.floor(data1[:, :2] / 100).astype(int)
        grid_coords2 = np.floor(data2[:, :2] / 100).astype(int)
        
        proportions1 = calculate_cluster_proportion(data1, grid_coords1)
        proportions2 = calculate_cluster_proportion(data2, grid_coords2)
        
        
        # find common grids
        common_grids = set(proportions1.keys()).intersection(set(proportions2.keys()))
        
        num_clusters = int(max(np.max(data1[:, 2]), np.max(data2[:, 2]))) + 1  # 确保 num_clusters 是整数
        vectors1 = create_proportion_vectors(proportions1, common_grids, num_clusters)
        vectors2 = create_proportion_vectors(proportions2, common_grids, num_clusters)
        
        # calculate the correlation between clusters
        correlation_matrix = np.zeros((num_clusters, num_clusters))
        
        for i in range(num_clusters):
            for j in range(num_clusters):
                if len(vectors1[i]) > 0 and len(vectors2[j]) > 0:
                    corr, _ = pearsonr(vectors1[i], vectors2[j])
                    correlation_matrix[i, j] = corr
                else:
                    correlation_matrix[i, j] = 0  # if there is no data, set it to 0

        column_names = [f'codex_label_{i}' for i in range(correlation_matrix.shape[1])]
        index_names = [f'trans_label_{i}' for i in range(correlation_matrix.shape[0])]


