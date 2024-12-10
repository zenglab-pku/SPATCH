########################## file information #############################
# Author: Pengfei Ren
# Date: 2024-12-04 21:16:02
# LastEditTime: 2024-12-04 23:24:19
# Description: 
# FilePath: /undefined/Users/morsouron/Desktop/Project/benchmark/code/2_8um_bin.py
#########################################################################

# bin xenium data at 8um resolution
def xenium_8um_bin(tr,tissue,resolution,platform):
    """
    Process transcriptomic data from Xenium into 8 µm resolution.

    Parameters:
        tr (str): Path to the transcript coordinate file. 
        tissue (str): Type of tissue being analyzed.
        resolution (float): Original resolution of the spatial data in µm. 
        platform (str): Sequencing platform type.

    Returns:
        AnnData with 8 µm resolution.
    """
    tr['binx'] = tr['x_location']*resolution//8
    tr['biny'] = tr['y_location']*resolution//8
    tr = tr.groupby(['binx', 'biny', 'feature_name']).size().unstack(fill_value=0)
    spatial = np.array(list(tr.index))
    tr.index = range(tr.shape[0])
    adata = sc.AnnData(csr_matrix(tr.values))
    adata.var_names = tr.columns
    adata.obs_names = tr.index
    adata.obsm['spatial'] = spatial * 8
    return(adata)

# bin cosmx data at 8um resolution    
def nano_8um_bin(tr,tissue,resolution,platform):
    """
    Process transcriptomic data from CoxMx into 8 µm resolution.

    Parameters:
        tr (str): Path to the transcript coordinate file. 
        tissue (str): Type of tissue being analyzed.
        resolution (float): Original resolution of the spatial data in µm. 
        platform (str): Sequencing platform type.

    Returns:
        AnnData with 8 µm resolution.
    """
    tr['binx'] = tr['x_global_px']*resolution//8
    tr['biny'] = tr['y_global_px']*resolution//8
    tr = tr.groupby(['binx', 'biny', 'target']).size().unstack(fill_value=0)
    spatial = np.array(list(tr.index))
    tr.index = range(tr.shape[0])
    adata = sc.AnnData(csr_matrix(tr.values))
    adata.var_names = tr.columns
    adata.obs_names = tr.index
    adata.obsm['spatial'] = spatial * 8
    return(adata)

# extract tissue mask based on the DAPI image
def dapi_mask(image,kernel_width):
    """
    Generate a binary mask for the largest contour in a DAPI-stained image.

    Parameters:
        image (numpy.ndarray): 
            A 2D array representing the grayscale image of a DAPI stain. 
            The image values are expected to be raw intensity values that will be normalized within the function.
        
        kernel_width (int): 
            The width of the square structuring element (kernel) used for morphological operations 
            (dilation and erosion). Determines the size of the neighborhood for processing.

    Returns:
        numpy.ndarray:
            A binary mask (2D array) with the same dimensions as the input image. 
            The mask highlights the largest connected component identified in the processed image.
    """
    image = (image/np.max(image)*255).astype(np.uint8)
    _, otsu_thresh = cv2.threshold(image, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)
    kernel = np.ones((kernel_width, kernel_width), np.uint8)
    dilated_image = cv2.dilate(otsu_thresh, kernel, iterations=2)
    eroded_image = cv2.erode(dilated_image, kernel, iterations=1)
    contours, _ = cv2.findContours(eroded_image, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    max_contour = max(contours, key=cv2.contourArea)
    mask = np.zeros_like(image)
    cv2.drawContours(mask, [max_contour], -1, 255, thickness=cv2.FILLED)
    return(mask)

# extract data within tissue mask
def tissue_extraction(adata,coord_system,mask):
    """
    Extract tissue data from a spatial transcriptomics dataset using a mask.

    Parameters:
        adata (AnnData): 
            The input AnnData object containing spatial transcriptomics data. 
        
        coord_system (str): 
            Specifies the coordinate system used in the dataset. 
        
        mask (numpy.ndarray): 
            A binary mask (2D array) representing the tissue region. 

    Returns:
        AnnData:
            A filtered AnnData object containing only the data corresponding to the tissue region defined by the mask.
    """
    loc = adata.obsm[coord_system].astype(int)
    common = []
    for i in range(adata.shape[0]):
        x = loc[i][0]
        y = loc[i][1]
        if y <= (mask.shape[0]-1) and x <= (mask.shape[1]-1) and mask[y,x] == 255:
            common.append(1)
        else:
            common.append(0)
    adata.obs['tissue_cut'] = common
    return(adata)

for tissue in tissues:
    nano = sc.read_h5ad('/home/renpf/benchmark/data/data/nanostring/' + tissue + '_bin_raw.h5ad')
    xenium = sc.read_h5ad('/home/renpf/benchmark/data/data/xenium/' + tissue + '_bin_raw.h5ad')
    nano_dapi = tifffile.imread('/home/renpf/benchmark/data/raw/nanostring/' + tissue + '/' + nano_files[tissues.index(tissue)])
    xenium_file = '/data200T/benchmark/xenium/' + '/' + xenium_files[tissues.index(tissue)]
    xenium_dapi = []
    with tifffile.TiffFile(xenium_file) as tif:
        for page in tif.pages:
            xenium_dapi.append(page.asarray())
    xenium_dapi = xenium_dapi[0]     
    nano_mask = dapi_mask(nano_dapi,1000)
    xenium_mask = dapi_mask(xenium_dapi,500)
    
    nano.obsm['tissue_cut'] = nano.obsm['spatial'].copy()
    nano.obsm['tissue_cut'][:,1] = (nano_dapi.shape[0] - nano.obsm['spatial'][:,1]/nano_res - fov_size)
    nano.obsm['tissue_cut'][:,0] = nano.obsm['spatial'][:,0]/nano_res
    
    nano = tissue_extraction(nano,'tissue_cut',nano_mask)
    nano.obsm['tissue_cut'] = nano.obsm['tissue_cut'] * nano_res
    nano = nano[nano.obs['tissue_cut']==1,:]
    nano.write('/home/renpf/benchmark/data/data/nanostring/' + tissue + '_bin_raw.h5ad')
    
    xenium.obsm['tissue_cut'] = xenium.obsm['spatial'].copy()
    xenium.obsm['tissue_cut'] = xenium.obsm['tissue_cut']/xenium_res
    xenium = tissue_extraction(xenium,'tissue_cut',xenium_mask)
    xenium.obsm['tissue_cut'] = xenium.obsm['tissue_cut'] * xenium_res
    xenium = xenium[xenium.obs['tissue_cut']==1,:]
    xenium.write('/home/renpf/benchmark/data/data/xenium/' + tissue + '_bin_raw.h5ad')

# bin sST data at 8um resolution
def bin_8um(adata,coord_system,original_resolution):
    """
    Process transcriptomic data for sST platforms into 8 µm resolution.

    Parameters:
        adata (AnnData): The input AnnData object containing spatial transcriptomics data. 
        coord_system (str): Specifies the coordinate system used in the dataset. 
        original_resolution (float): Original resolution of the spatial data in µm. 

    Returns:
        AnnData with 8 µm resolution.
    """
    spatial_coords = adata.obsm[coord_system]
    new_resolution = 8.0       

    grid_indices = (spatial_coords * original_resolution // new_resolution).astype(int)
    unique_grid_indices, inverse_indices = np.unique(grid_indices, axis=0, return_inverse=True)
    new_X = np.zeros((len(unique_grid_indices), adata.X.shape[1]))

    for i in range(len(unique_grid_indices)):
        points_in_grid = (inverse_indices == i)
        new_X[i, :] = adata.X[points_in_grid, :].sum(axis=0)

    new_adata = sc.AnnData(X=new_X)
    new_adata.obsm['spatial'] = unique_grid_indices
    new_adata.var_names = adata.var_names
    
    return(new_adata)