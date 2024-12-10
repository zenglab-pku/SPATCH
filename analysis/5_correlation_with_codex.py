########################## file information #############################
# Author: Pengfei Ren
# Date: 2024-12-04 21:38:41
# LastEditTime: 2024-12-04 23:34:04
# Description: 
# FilePath: /undefined/Users/morsouron/Desktop/Project/benchmark/code/5_correlation_with_codex.py
#########################################################################

def bin_adata(adata,coord_system,original_resolution,new_resolution):
    """
    Bin the spatial transcriptomics data in the `adata` object based on a specified resolution.
    The function aggregates the expression data of cells falling within the same grid cell after rescaling their spatial coordinates.

    Parameters:
        adata (AnnData): 
            The AnnData object containing the original spatial transcriptomics data

        coord_system (str): 
            The key in the `adata.obsm` dictionary that contains the spatial coordinates of the cells. These coordinates will be rescaled according to the new resolution.

        original_resolution (float): 
            The original resolution (in microns or other spatial units) of the spatial coordinates in `adata.obsm[coord_system]`.

        new_resolution (float): 
            The new resolution (in microns or other spatial units) to which the spatial coordinates should be rescaled. This will determine the size of the grid cells after binning.

    Returns:
        AnnData:
            A new AnnData object where the spatial coordinates are rescaled to the new resolution and the gene expression data is aggregated within each new grid cell.
    """
    spatial_coords = adata.obsm[coord_system]
    grid_indices = (spatial_coords * original_resolution // new_resolution).astype(int)
    unique_grid_indices, inverse_indices = np.unique(grid_indices, axis=0, return_inverse=True)
    old_X = adata.X.copy()
    new_X = np.zeros((len(unique_grid_indices), adata.X.shape[1]))

    for i in range(len(unique_grid_indices)):
        points_in_grid = (inverse_indices == i)
        new_X[i, :] = old_X[points_in_grid, :].sum(axis=0)

    new_adata = sc.AnnData(X=new_X)
    new_adata.obsm['spatial'] = unique_grid_indices
    new_adata.var_names = adata.var_names
    
    return(new_adata)

def tr_codex_cell_intensity(adata,adata_codex,marker_gene,tissue,platform):
    """
    This function calculates the cell intensity for marker genes and their corresponding proteins from both ST 
    and Codex data, respectively, by binning the spatial data and matching the gene and protein expression within the same bins.
    
    Parameters:
        adata (AnnData): 
            The AnnData object containing the single-cell RNA-seq data with spatial information.
        
        adata_codex (AnnData): 
            The AnnData object containing the Codex data with spatial information.
        
        marker_gene (dict): 
            A dictionary mapping gene symbols to corresponding protein markers.
        
        tissue (str): 
            The type of tissue for the analysis.
        
        platform (str): 
            The platform type used for spatial transcriptomics. 
            
    Returns:
        pd.DataFrame:
            A dataframe containing the binned spatial coordinates, gene expression, and protein expression for each bin, 
            with columns:
            - 'x_bin': X coordinate of the bin
            - 'y_bin': Y coordinate of the bin
            - 'gene_expr': Expression of the gene in the bin
            - 'protein_expr': Expression of the corresponding protein in the bin (with NaN values replaced by 0).
    """
    
    adata = adata[adata.obs['codex_common']==1,:]
    adata_codex = adata_codex[adata_codex.obs['codex_common']==1,:]
    
    bin_dfs = []
    for bin_size in [100,200,300,400,500]:
        
        adata_binned = bin_adata(adata,'codex_registered',1,bin_size)
        adata_codex_binned = bin_adata(adata_codex,'spatial',1,bin_size)
        
        for gene,protein in marker_gene.items():
            gene_binned = adata_binned[:, gene].X.toarray().flatten()
            protein_binned = adata_codex_binned[:, protein].X.toarray().flatten()
            
            bin_gene_df = pd.DataFrame({
                'x_bin': adata_binned.obsm['spatial'][:,0],
                'y_bin': adata_binned.obsm['spatial'][:,1],
                'gene_expr': gene_binned
            })
            bin_protein_df = pd.DataFrame({
                'x_bin': adata_codex_binned.obsm['spatial'][:,0],
                'y_bin': adata_codex_binned.obsm['spatial'][:,1],
                'protein_expr': protein_binned
            })
            
            bin_protein_gene = pd.merge(bin_gene_df,bin_protein_df,how='inner')
            bin_protein_gene.loc[bin_protein_gene['protein_expr'].isna(),'protein_expr'] = 0

    return(bin_protein_gene)