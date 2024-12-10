'''
Author: Pengfei Ren
Date: 2024-12-04 21:14:03
LastEditTime: 2024-12-04 21:14:04
Description: 
FilePath: /undefined/Users/morsouron/Desktop/Project/benchmark/code/load_data
'''

def load_stereo_cell_bin(mask, adata):
    """
    Aggregate bin-level expression data into cell-level data based on a given mask.

    Parameters:
    ----------
        mask : A 2D array where each pixel is assigned a cell ID, indicating the cell it belongs to.
        adata : A Scanpy AnnData object containing bin-level gene expression data 

    Returns:
    --------
        AnnData: A new AnnData object containing aggregated gene expression data and spatial information for each cell.
    """

    # Step 1: Group bins by their corresponding IDs in the mask
    cell_grids = {}  
    for idx, (x, y) in enumerate(zip(adata.obs['x'], adata.obs['y'])):
        cell_id = mask[y, x] 
        if cell_id not in cell_grids:
            cell_grids[cell_id] = []
        cell_grids[cell_id].append(idx)  

    # Step 2: Aggregate gene expression and calculate spatial centers for each bin
    cell_data = []  
    cell_obs = []   

    for cell_id, indices in cell_grids.items():
        # Aggregate gene expression data for all bins in the cell 
        cell_expression = np.sum(adata[indices].X, axis=0)

        # Calculate the spatial center of the bin by averaging the coordinates of its bins
        center_x = np.mean(adata.obs['x'][indices])
        center_y = np.mean(adata.obs['y'][indices])

        cell_data.append(cell_expression)  
        cell_obs.append([center_x, center_y]) 

    # Step 3: Create a new AnnData object for the aggregated data
    cell_X = csr_matrix(np.vstack(cell_data))  
    cell_obs = pd.DataFrame(cell_obs, columns=['x', 'y'])  

    # Create a new AnnData object with aggregated gene expression and spatial data
    cell_adata = sc.AnnData(X=cell_X, obs=cell_obs, var=adata.var)
    cell_adata.obsm['spatial'] = cell_adata.obs[["x", "y"]].values / 2

    return cell_adata

def load_hd_cell_bin(dir_base,polys):
    """
    Aggregate bin-level expression data into cell-level data based on a given mask.

    Parameters:
    ----------
        dir_base : Root data path.
        polys : A 2D array where each pixel is assigned a cell ID, indicating the cell it belongs to.

    Returns:
    --------
        AnnData: A new AnnData object containing aggregated gene expression data and spatial information for each cell.
    """
    # Creating a list to store Polygon geometries
    geometries = []

    # Iterating through each nuclei in the 'polys' DataFrame
    for nuclei in range(len(polys['coord'])):

        # Extracting coordinates for the current nuclei and converting them to (y, x) format
        coords = [(y, x) for x, y in zip(polys['coord'][nuclei][0], polys['coord'][nuclei][1])]

        # Creating a Polygon geometry from the coordinates
        geometries.append(Polygon(coords))

    # Creating a GeoDataFrame using the Polygon geometries
    gdf = gpd.GeoDataFrame(geometry=geometries)
    gdf['id'] = [f"ID_{i+1}" for i, _ in enumerate(gdf.index)]
    # Load Visium HD data
    raw_h5_file = dir_base+'filtered_feature_bc_matrix.h5'
    adata = sc.read_10x_h5(raw_h5_file)

    # Load the Spatial Coordinates
    tissue_position_file = dir_base+'spatial/tissue_positions.parquet'
    df_tissue_positions=pd.read_parquet(tissue_position_file)

    #Set the index of the dataframe to the barcodes
    df_tissue_positions = df_tissue_positions.set_index('barcode')

    # Create an index in the dataframe to check joins
    df_tissue_positions['index']=df_tissue_positions.index

    # Adding the tissue positions to the meta data
    adata.obs =  pd.merge(adata.obs, df_tissue_positions, left_index=True, right_index=True)

    # Create a GeoDataFrame from the DataFrame of coordinates
    geometry = [Point(xy) for xy in zip(df_tissue_positions['pxl_col_in_fullres'], df_tissue_positions['pxl_row_in_fullres'])]
    gdf_coordinates = gpd.GeoDataFrame(df_tissue_positions, geometry=geometry)
    # Perform a spatial join to check which coordinates are in a cell nucleus
    result_spatial_join = gpd.sjoin(gdf_coordinates, gdf, how='left', predicate='within')

    # Identify nuclei associated barcodes and find barcodes that are in more than one nucleus
    result_spatial_join['is_within_polygon'] = ~result_spatial_join['index_right'].isna()
    barcodes_in_overlaping_polygons = pd.unique(result_spatial_join[result_spatial_join.duplicated(subset=['index'])]['index'])
    result_spatial_join['is_not_in_an_polygon_overlap'] = ~result_spatial_join['index'].isin(barcodes_in_overlaping_polygons)

    # Remove barcodes in overlapping nuclei
    barcodes_in_one_polygon = result_spatial_join[result_spatial_join['is_within_polygon'] & result_spatial_join['is_not_in_an_polygon_overlap']]

    # The AnnData object is filtered to only contain the barcodes that are in non-overlapping polygon regions
    filtered_obs_mask = adata.obs_names.isin(barcodes_in_one_polygon['index'])
    filtered_adata = adata[filtered_obs_mask,:]

    # Add the results of the point spatial join to the Anndata object
    filtered_adata.obs =  pd.merge(filtered_adata.obs, barcodes_in_one_polygon[['index','geometry','id','is_within_polygon','is_not_in_an_polygon_overlap']], left_index=True, right_index=True)
    # Group the data by unique nucleous IDs
    groupby_object = filtered_adata.obs.groupby(['id'], observed=True)

    # Extract the gene expression counts from the AnnData object
    counts = filtered_adata.X

    N_groups = groupby_object.ngroups
    N_genes = counts.shape[1]

    # Initialize a sparse matrix to store the summed gene counts for each nucleus
    summed_counts = sparse.lil_matrix((N_groups, N_genes))

    polygon_id = []
    array_rows = []
    array_cols = []
    pxl_rows = []
    pxl_cols = []
    row = 0

    # Iterate over each unique polygon to calculate the sum of gene counts.
    for polygons, idx_ in groupby_object.indices.items():
        summed_counts[row] = counts[idx_].sum(0)
        barcodes = filtered_adata.obs.iloc[idx_,].index
        row += 1
        polygon_id.append(polygons)
        pxl_rows.append(np.mean(filtered_adata.obs.loc[barcodes,'pxl_row_in_fullres']))
        array_rows.append(np.mean(filtered_adata.obs.loc[barcodes,'array_row']))
        pxl_cols.append(np.mean(filtered_adata.obs.loc[barcodes,'pxl_col_in_fullres']))
        array_cols.append(np.mean(filtered_adata.obs.loc[barcodes,'array_col']))
    summed_counts = summed_counts.tocsr()
    grouped_filtered_adata = anndata.AnnData(X=summed_counts,obs=pd.DataFrame({'id':polygon_id,'pxl_row_in_fullres':pxl_rows,'pxl_col_in_fullres':pxl_cols,'array_row':array_rows,'array_col':array_cols},index=polygon_id),var=filtered_adata.var)
    grouped_filtered_adata.obsm['spatial'] =  grouped_filtered_adata.obs[["array_row", "array_col"]].values * 2
    return(grouped_filtered_adata)

def load_xenium(
    path : str
) -> AnnData :
    """
    Read *10x Genomics* Xenium formatted dataset.

    Parameters
    ----------
    path
        Path to the root directory containing *Xenium* files.

    Returns
    -------
    Annotated data object
    """ 
    # count matrix
    adata = sc.read_10x_h5(filename=os.path.join(path, "cell_feature_matrix.h5"))
    # meta data for cells
    df = pd.read_csv(os.path.join(path, "cells.csv.gz"))
    df.set_index(adata.obs_names, inplace=True)
    adata.obs = df.copy()
    # add spatial location (default um)
    adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()
    return adata

def load_hd(
    path: PathLike,
    *,
    counts_file: str = "filtered_feature_bc_matrix.h5", 
    library_id: str | None = None,
    load_images: bool = False,
    source_image_path: PathLike | None = None,
    **kwargs: Any,
) -> AnnData:
    """
    Read *10x Genomics* Visium formatted dataset.

    In addition to reading the regular *Visium* output, it looks for the *spatial* directory and loads the images,
    spatial coordinates and scale factors.

    Parameters
    ----------
    path
        Path to the root directory containing *Visium* files.
    counts_file
        Which file in the passed directory to use as the count file. Typically either *filtered_feature_bc_matrix.h5* or
        *raw_feature_bc_matrix.h5*.
    library_id
        Identifier for the *Visium* library. Useful when concatenating multiple :class:`anndata.AnnData` objects.
    kwargs
        Keyword arguments for :func:`scanpy.read_10x_h5`, :func:`anndata.read_mtx` or :func:`read_text`.

    Returns
    -------
    Annotated data object
    """ 
    
    if '002' in path:
        resolution = 2
    elif '008' in path:
        resolution = 8
    elif '016' in path:
        resolution = 16
    
    path = Path(path)
    adata, library_id = _read_counts(path, count_file=counts_file, library_id=library_id, **kwargs)

    tissue_positions_file = (
        path / "spatial/tissue_positions.parquet"
    )
    coords = pd.read_parquet(tissue_positions_file)
    coords.index = coords['barcode']
    coords = coords.drop(columns='barcode')
    coords.columns = ["in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]
    coords.set_index(coords.index.astype(adata.obs.index.dtype), inplace=True)

    adata.obs = pd.merge(adata.obs, coords, how="left", left_index=True, right_index=True)
    adata.obsm[Key.obsm.spatial] = adata.obs[["array_row", "array_col"]].values
    # convert to um
    adata.obsm[Key.obsm.spatial] = adata.obsm[Key.obsm.spatial]*resolution
    
    if not load_images:
        return adata

    adata.uns[Key.uns.spatial][library_id][Key.uns.image_key] = {
        res: _load_image(path / f"{Key.uns.spatial}/tissue_{res}_image.png") for res in ["hires", "lowres"]
    }
    adata.uns[Key.uns.spatial][library_id]["scalefactors"] = json.loads(
        (path / f"{Key.uns.spatial}/scalefactors_json.json").read_bytes()
    )

    if source_image_path is not None:
        source_image_path = Path(source_image_path).absolute()
        if not source_image_path.exists():
            logg.warning(f"Path to the high-resolution tissue image `{source_image_path}` does not exist")
        adata.uns["spatial"][library_id]["metadata"]["source_image_path"] = str(source_image_path)

    return adata

def load_nanostring(
    path: str
) -> AnnData:
    """
    Read *Nanostring* formatted dataset.

    In addition to reading the regular *Nanostring* output, it loads the metadata file, if present *CellComposite* and *CellLabels*
    directories containing the images and optionally the field of view file.

    Parameters
    ----------
    path
        Path to the root directory containing *Nanostring* files.

    Returns
    -------
    Annotated data object 
    """  # 
    counts_file = 'exprMat_file.csv.gz'
    meta_file = 'metadata_file.csv.gz'
    path, fov_key = Path(path), "fov"
    cell_id_key = "cell_ID"
    counts = pd.read_csv(path / counts_file, header=0, index_col=cell_id_key)
    counts.index = counts.index.astype(str).str.cat(counts.pop(fov_key).astype(str).values, sep="_")
    # counts = counts.drop('cell',axis = 1)

    obs = pd.read_csv(path / meta_file, header=0, index_col=cell_id_key)
    obs[fov_key] = pd.Categorical(obs[fov_key].astype(str))
    obs[cell_id_key] = obs.index.astype(np.int64)
    obs.rename_axis(None, inplace=True)
    obs.index = obs.index.astype(str).str.cat(obs[fov_key].values, sep="_")

    common_index = obs.index.intersection(counts.index)

    adata = AnnData(
        csr_matrix(counts.loc[common_index, :].values),
        dtype=counts.values.dtype,
        obs=obs.loc[common_index, :],
    )
    adata.var_names = counts.columns
    
    adata.obsm['spatial'] = adata.obs[['CenterX_global_px', 'CenterY_global_px']].values
    # convert to um
    adata.obsm['spatial'] = adata.obsm['spatial']*120.280945/1000
    return adata

def load_codex(
    path: str
    ):
    """
    Read *CODEX* data processed by Stardist and QuPath.

    Parameters
    ----------
    path
        Path to the root directory containing *CODEX* files.

    Returns
    -------
    Annotated data object 
    """  # 
    major_type = {'Pan-Cytokeratin':'Epithelial','CD34':'Endothelial','CD3e':'T','SMA':'Fibroblast','CD68':'Macrophage','CD20':'B','CD56':'NK','CD4':'CD4T','CD8':'CD8T','FOXP3':'Treg'}
    # load data from the output of qupath 
    data = pd.read_table(path + '/measurements_major.tsv')
    t =  pd.read_table(path + '/measurements_t.tsv')
    foxp3 = pd.read_table(path + '/measurements_foxp3.tsv')
    data['T'] = t['Classification']
    data['FOXP3'] = foxp3['Classification']
    data['major'] = data['Classification']
    data['minor'] = data['Classification']
    data.loc[(data['major'] == 'CD3e') & (data['T'] == 'CD4'),'minor'] = 'CD4'
    data.loc[(data['major'] == 'CD3e') & (data['T'] == 'CD8'),'minor'] = 'CD8'
    data.loc[(data['minor'] == 'CD4') & (data['FOXP3'] == 'FOXP3'),'minor'] = 'FOXP3'
    # remove DAPI channel
    data = data.loc[:, ~data.columns.str.contains("DAPI")]
    # extract counts matrix
    counts = data.filter(like="Cell: Mean")
    counts.columns = counts.columns.str.replace(": Cell: Mean", "", regex=False)
    adata = AnnData(counts)
    adata.obs_names = data['Object ID']
    # add spatial location
    adata.obsm['spatial'] = data.loc[:,['Centroid X µm','Centroid Y µm']].to_numpy()
    adata.obs['major'] = data['major'].tolist()
    adata.obs['major'] = adata.obs['major'].map(major_type)
    adata.obs['minor'] = data['minor'].tolist()
    adata.obs['minor'] = adata.obs['minor'].map(major_type)
    # remove low quality cells(based on stardist output)
    index = (data['Detection probability']>=np.quantile(data['Detection probability'],q = np.arange(0,1,0.1))[1]) & (~data['major'].isna())
    adata = adata[index,:]
    return(adata)