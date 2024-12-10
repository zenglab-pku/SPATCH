'''
Author: Pengfei Ren
Date: 2024-12-04 21:21:27
LastEditTime: 2024-12-04 21:23:29
Description: 
FilePath: /undefined/Users/morsouron/Desktop/Project/benchmark/code/registeration
'''
# resize images
def resample_image(image, fix_img_res, codex_res):
    """
    Resample an image to a new resolution using SimpleITK.

    Parameters:
        image (numpy.ndarray): 
            A 2D or 3D array representing the image to be resampled. 
        
        fix_img_res (float): 
            The spatial resolution (pixel size) of the input image. 

        codex_res (float): 
            The target spatial resolution (pixel size) for the resampled image. 

    Returns:
        numpy.ndarray: 
            A 2D or 3D array representing the resampled image with the target resolution. 
    """
    image = sitk.GetImageFromArray(image.astype(np.float32))
    image.SetSpacing([fix_img_res, fix_img_res])
    
    original_spacing = image.GetSpacing()
    original_size = image.GetSize()
    new_spacing = [codex_res,codex_res]
    new_size = [int(round(osz * ospc / nspc)) for osz, ospc, nspc in zip(original_size, original_spacing, new_spacing)]
    
    resample = sitk.ResampleImageFilter()
    resample.SetOutputSpacing(new_spacing)
    resample.SetSize(new_size)
    resample.SetOutputDirection(image.GetDirection())
    resample.SetOutputOrigin(image.GetOrigin())
    resample.SetInterpolator(sitk.sitkLinear)
    
    return sitk.GetArrayFromImage(resample.Execute(image))

# register 
def register(fix_img,mv_img,fix_img_res,codex_res,landmark_path,niter,platform):
    """
    Perform image registration between a fixed image (`fix_img`) and a moving image (`mv_img`)
    using SimpleITK. The registration uses landmarks to define the initial transform
    and refines it using a multi-resolution approach with Mattes Mutual Information.

    Parameters:
        fix_img (numpy.ndarray): 
            The fixed reference image to which the moving image will be aligned.
        
        mv_img (numpy.ndarray): 
            The moving image that will be transformed to align with the fixed image.
        
        fix_img_res (float): 
            The resolution (pixel spacing) of the fixed image.
        
        codex_res (float): 
            The resolution (pixel spacing) for the moving image and target resolution
            for registration.
        
        landmark_path (str): 
            Path to the directory containing landmark files for both the fixed and moving images.
        
        niter (int): 
            Number of iterations for the gradient descent optimizer in the registration process.
        
        platform (str): 
            Specifies the platform type. If `platform='nano'`, additional preprocessing is applied 
            (e.g., flipping and adjusting landmarks).

    Returns:
        tuple: 
            - `fix_img` (SimpleITK.Image): The fixed image resampled to `codex_res` resolution.
            - `final_transform` (SimpleITK.Transform): The final transformation that aligns `mv_img` to `fix_img`.
    """
    
    # Use landmarks to define the initial transform
    if platform=='nano':
        mv_img = np.flip(mv_img, axis=0)
    mv_landmarks = pd.read_csv(landmark_path + '/codex.tif-points.tsv',header=0,sep='\t')
    mv_landmarks = np.array((mv_landmarks.loc[:,['x','y']]*codex_res).astype(np.float32))
    if platform=='nano':
        mv_landmarks[:, 1] = mv_img.shape[0]*codex_res - mv_landmarks[:, 1]
        mv_landmarks = mv_landmarks.reshape(-1).tolist()
    mv_img = sitk.GetImageFromArray(mv_img.astype(np.float32))
    mv_img.SetSpacing([codex_res, codex_res])
    
    fix_landmarks = pd.read_csv(landmark_path + '/nano.tif-points.tsv',header=0,sep='\t')
    fix_landmarks = np.array((fix_landmarks.loc[:,['x','y']]*codex_res).astype(np.float32)).reshape(-1).tolist()
    fix_img = resample_image(fix_img, fix_img_res, codex_res)
    fix_img = sitk.GetImageFromArray(fix_img.astype(np.float32))
    fix_img.SetSpacing([codex_res, codex_res])
    
    initial_transform = sitk.LandmarkBasedTransformInitializer(sitk.Similarity2DTransform(),
                                                           fix_landmarks,
                                                           mv_landmarks)
    
    # RegistrationMethod
    registration_method = sitk.ImageRegistrationMethod()

    # Similarity metric settings.
    registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=50)
    registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
    registration_method.SetMetricSamplingPercentage(0.01)
    registration_method.SetInterpolator(sitk.sitkLinear)

    # Optimizer settings.
    registration_method.SetOptimizerAsGradientDescent(learningRate=1, numberOfIterations=niter, convergenceMinimumValue=1e-6, convergenceWindowSize=10)
    registration_method.SetOptimizerScalesFromPhysicalShift()

    # Setup for the multi-resolution framework.            
    registration_method.SetShrinkFactorsPerLevel(shrinkFactors = [4,2,1])
    registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas=[2,1,0])
    registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

    # Don't optimize in-place, we would possibly like to run this cell multiple times.
    registration_method.SetInitialTransform(initial_transform, inPlace=False)

    # Execute registration
    final_transform = registration_method.Execute(sitk.Cast(fix_img, sitk.sitkFloat32), 
                                                   sitk.Cast(mv_img, sitk.sitkFloat32))
    # registration
    mv_img_resampled = sitk.Resample(mv_img, fix_img, final_transform, sitk.sitkLinear, 0.0, mv_img.GetPixelID())
    print(sitk.GetArrayFromImage(fix_img).shape)
    plt.imshow(sitk.GetArrayFromImage(fix_img))
    plt.show()
    print(sitk.GetArrayFromImage(mv_img_resampled).shape)
    plt.imshow(sitk.GetArrayFromImage(mv_img_resampled))
    plt.show()
    return(fix_img,final_transform)

def apply_transform_to_codex(codex,codex_meta,fix_img,final_transform,path_out,platform):
    """
    Apply a spatial transformation to a set of CODEX imaging data channels, aligning them to a reference image.
    The transformed data is saved as a multi-dimensional TIFF file.

    Parameters:
        codex (list of numpy.ndarray): 
            A list of CODEX image channels to be registered. Each channel is a 2D NumPy array.
        
        codex_meta (list of metadata objects): 
            Metadata for the CODEX image channels, containing information such as resolution and channel names.
        
        fix_img (SimpleITK.Image): 
            The fixed reference image used in the registration process.
        
        final_transform (SimpleITK.Transform): 
            The transformation object obtained from registration, applied to align the CODEX channels.
        
        path_out (str): 
            Path to save the registered output as a multi-channel OME-TIFF file.
        
        platform (str): 
            Specifies the imaging platform type. If `'nano'`, the input images are preprocessed by flipping along the y-axis.

    Returns:
        numpy.ndarray: 
            A stacked 3D NumPy array of the registered CODEX image channels, with dimensions (C, Y, X).
    """
    # apply transform to each channel
    registered_codex = []
    resampler = sitk.ResampleImageFilter()
    resampler.SetReferenceImage(fix_img)
    resampler.SetTransform(final_transform)
    resampler.SetInterpolator(sitk.sitkLinear)
    codex_res = 10000/(codex_meta[0].values()[9].value[0]/codex_meta[0].values()[9].value[1])
    for channel in codex:
        if platform=='nano':
            channel = np.flip(channel, axis=0)
        channel = sitk.GetImageFromArray(channel)
        channel = sitk.Cast(channel, sitk.sitkFloat32)
        channel.SetSpacing([codex_res, codex_res])
        channel = resampler.Execute(channel)
        registered_codex.append(sitk.GetArrayFromImage(channel))
    plt.imshow(registered_codex[0])
    plt.show()
    
    # convert data type from float to int to save memory
    registered_codex = [np.round(array).astype(np.uint16) for array in registered_codex]
    
    # Convert list of images to a multi-dimensional numpy array (C, Y, X)
    registered_codex_stack = np.stack(registered_codex, axis=0)

    # Determine resolution from metadata (assuming it's in the first image metadata)
    resolution = codex_meta[0].get('XResolution', (1, 1))
    resolution = 10000/(resolution.value[0]/resolution.value[1])

    # obtain channel names from metadata 
    channel_names = []
    for meta in codex_meta:
        channel_names.append(xmltodict.parse(meta.values()[6].value)['PerkinElmer-QPI-ImageDescription']['Biomarker'])
    
    #save data
    output_filename = path_out + '/codex_registered.tif'
    tifffile.imwrite(output_filename, registered_codex_stack, 
                 metadata={'Pixels': {
                                'PhysicalSizeX': resolution,
                                'PhysicalSizeXUnit': 'µm',
                                'PhysicalSizeY': resolution,
                                'PhysicalSizeYUnit': 'µm'}, 
                            "axes": 'CYX', 
                            "Channel": {"Name": channel_names}}, ome=True, compression = 'lzw')
    return(registered_codex_stack)

def apply_transform_to_scanpy(adata, coord_system, transform, mv_landmark):
    """
    Apply a spatial transformation to the coordinates stored in a Scanpy AnnData object.
    
    Parameters:
        adata (AnnData): 
            The Scanpy AnnData object containing the spatial coordinates and other data.
        
        coord_system (str): 
            The key in `adata.obsm` where the spatial coordinates are stored. 
            This is typically a 2D array with x and y coordinates for each cell.
        
        transform (SimpleITK.Transform): 
            The transformation object that will be applied to the coordinates. 
            This is typically the result of an image registration, representing translation and rotation.
        
        mv_landmark (pandas.DataFrame): 
            A DataFrame containing the landmark coordinates of the moving image in the form of x and y values.
    
    Returns:
        numpy.ndarray: 
            The transformed coordinates as a NumPy array.
    """
    coordinates = adata.obsm[coord_system]/codex_resolution  
    
    transform = transform.GetNthTransform(0)
    # center = np.array(transform.GetCenter())
    matrix = np.array(transform.GetMatrix()).reshape((2, 2))  
    translation = np.array(transform.GetTranslation())
    center = np.array([np.mean(mv_landmark['x']),np.mean(mv_landmark['y'])])

    transformed_coords = (coordinates - center) @ matrix + center - translation
    return transformed_coords # codex pixel


# extract common region
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
        if mask[y,x] !=0:
            common.append(1)
        else:
            common.append(0)
    adata.obs['codex_common'] = common
    return adata

mask1 = dapi_mask(image1,200)
mask2 = dapi_mask(image2,200)
intersection_mask = cv2.bitwise_and(mask1, mask2)