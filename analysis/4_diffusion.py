########################## file information #############################
# Author: Pengfei Ren
# Date: 2024-12-04 21:37:23
# LastEditTime: 2024-12-04 23:27:24
# Description: 
# FilePath: /undefined/Users/morsouron/Desktop/Project/benchmark/code/4_diffusion.py
#########################################################################

tissues = ['COAD','HCC','OV']

def qc(obj):
    obj.obs_names_make_unique()
    obj.var_names_make_unique()
    obj.var["mt"] = obj.var_names.str.startswith("MT-")
    # ribosomal genes
    obj.var["ribo"] = obj.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    obj.var["hb"] = obj.var_names.str.contains("^HB[^(P)]")
    sc.pp.calculate_qc_metrics(
        obj, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )
    return(obj)

for tissue in tissues:
    # hd
    hd = load_hd_all('/home/renpf/benchmark/data/raw/hd/' + tissue + '/outs/binned_outputs/square_008um/')
    hd_filter = sc.read_h5ad('/home/renpf/benchmark/data/data/hd/' + tissue + '_bin_raw.h5ad')
    hd = hd[:,hd.var_names.isin(hd_filter.var_names)]
    hd = qc(hd)

    #stereo
    bin_size = 16
    path = '/home/renpf/benchmark/data/raw/stereo/' + tissue + '/' 
    adata = stereo.io.read_gef(file_path= glob.glob('/home/renpf/benchmark/data/raw/stereo/' + tissue + '/*.gef')[2], bin_size=bin_size, gene_name_index=True)
    adata = stereo.io.stereo_to_anndata(adata,flavor='scanpy')

    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith('_tissue_cut.tif'):
                tifile = path + file

    # load tissue mask for Stereo-seq
    mask = io.imread(tifile)
    mask = np.transpose(mask)
    mask = np.where(mask == 1)
    mask = pd.DataFrame({'x':mask[0], 'y':mask[1]})
    mask = mask//bin_size
    mask = mask.drop_duplicates()

    adata = qc(adata)
    bin_loc = adata.obs.loc[:,['x','y']]//bin_size
    in_tissue = pd.merge(mask,bin_loc, on=['x', 'y'], how='inner')
    # extract spots within tissue region
    index = bin_loc.apply(tuple, 1).isin(in_tissue.apply(tuple, 1))
    
    #convert to um
    hd.obs[['x','y']] = hd.obs.loc[:,['array_row','array_col']]*8
    adata.obs[['x','y']] = adata.obs.loc[:,['x','y']]/2
    
    #common genes between Stereo-seq and Visium HD
    common_genes = reduce(np.intersect1d,[hd.var_names,adata.var_names])
    hd = hd[:,common_genes]
    hd = qc(hd)
    adata_common = adata[:,common_genes]
    adata_common = qc(adata_common)
    
    #save counts
    in_tissue = hd.obs.loc[hd.obs['in_tissue']==1,['total_counts','x','y']]
    out_tissue = hd.obs.loc[hd.obs['in_tissue']==0,['total_counts','x','y']]
    in_tissue.to_csv('/home/renpf/benchmark/res/gene/diffusion/tissue/' + tissue + '_hd_in.csv')
    out_tissue.to_csv('/home/renpf/benchmark/res/gene/diffusion/tissue/' + tissue + '_hd_out.csv')

    in_tissue = adata_common.obs.loc[index,['total_counts','x','y']]
    out_tissue = adata_common.obs.loc[~index,['total_counts','x','y']]
    in_tissue.to_csv('/home/renpf/benchmark/res/gene/diffusion/tissue/' + tissue + '_stereo_common_in.csv')
    out_tissue.to_csv('/home/renpf/benchmark/res/gene/diffusion/tissue/' + tissue + '_stereo_common_out.csv')
    
    #stereo-seq unique genes
    adata_uniq = adata[:,~adata.var_names.isin(common_genes)]
    adata_uniq = qc(adata_uniq)
    in_tissue = adata_uniq.obs.loc[index,['total_counts','x','y']]
    out_tissue = adata_uniq.obs.loc[~index,['total_counts','x','y']]
    in_tissue.to_csv('/home/renpf/benchmark/res/gene/diffusion/tissue/' + tissue + '_stereo_uniq_in.csv')
    out_tissue.to_csv('/home/renpf/benchmark/res/gene/diffusion/tissue/' + tissue + '_stereo_uniq_out.csv')
    
def calculate_distance(df1,df2,batch_size):
    """
    Calculate the minimum distance between points in two DataFrames (df1 and df2) based on their x and y coordinates.
    The distances are computed using the nearest neighbor algorithm in batches to handle large datasets efficiently.

    Parameters:
        df1 (pandas.DataFrame): 
            The first DataFrame containing x and y coordinates of the reference points. 

        df2 (pandas.DataFrame): 
            The second DataFrame containing x and y coordinates of the target points.

        batch_size (int): 
            The size of each batch for processing the distances. The function processes `df2` in batches of this 
            size to improve performance when handling large datasets.

    Returns:
        pandas.DataFrame: 
            The updated `df2` DataFrame, with an additional column 'distance' that contains the minimum distance 
            from each point in `df2` to the nearest point in `df1`.
    """
    X1 = df1[['x', 'y']].values
    X2 = df2[['x', 'y']].values
    nbrs = NearestNeighbors(n_neighbors=1).fit(X1)

    batch_size = batch_size
    min_distances = []

    for i in range(0, len(X2), batch_size):
        batch_X2 = X2[i:i+batch_size]
        distances, indices = nbrs.kneighbors(batch_X2)
        min_distances.extend(distances.flatten())  

    df2['distance'] = min_distances
    return df2
    
for tissue in tissues:
    for platform in ['hd','stereo_common','stereo_uniq']:
        df1 = pd.read_csv('/home/renpf/benchmark/res/gene/diffusion/tissue/' + tissue + '_' + platform + '_in.csv')
        df2 = pd.read_csv('/home/renpf/benchmark/res/gene/diffusion/tissue/' + tissue + '_' + platform + '_out.csv')
        df2 = calculate_distance(df1,df2,1000000)
        df2.to_csv('/home/renpf/benchmark/res/gene/diffusion/tissue/' + tissue + '_' + platform + '_out.csv', index = False)
print('All done')