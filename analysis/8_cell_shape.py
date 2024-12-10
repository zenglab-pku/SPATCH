########################## file information #############################
# Author: Pengfei Ren
# Date: 2024-12-04 23:35:08
# LastEditTime: 2024-12-04 23:35:09
# Description: 
# FilePath: /undefined/Users/morsouron/Desktop/Project/benchmark/code/cell_shape.py
#########################################################################

def calculate_contour_metrics(contour):
    """
    Calculate shape-based metrics for a given contour.
    
    Parameters:
    contour (numpy.ndarray): Contour points for the object.
    
    Returns:
    tuple: Contour solidity, circularity, and aspect ratio.
    """
    area = cv2.contourArea(contour)
    hull = cv2.convexHull(contour)
    hull_area = cv2.contourArea(hull)
    
    # Calculate Solidity
    solidity = float(area) / hull_area if hull_area > 0 else 0
    
    perimeter = cv2.arcLength(contour, True)
    # Calculate Circularity
    circularity = 4 * np.pi * (area / (perimeter ** 2)) if perimeter > 0 else 0
    
    # Get bounding rectangle for aspect ratio
    x, y, w, h = cv2.boundingRect(contour)
    aspect_ratio = float(w) / h if h > 0 else 0

    return solidity, circularity, aspect_ratio

def process_image(image, tissue_label, platform_label):
    """
    Process the binary mask image and extract geometric metrics for each detected contour.
    
    Parameters:
    image (numpy.ndarray): Binary mask image where non-zero pixels represent object regions.
    tissue_label (str): Label for the tissue type.
    platform_label (str): Label for the platform used in the experiment.
    
    Returns:
    tuple: Lists of extracted metrics (solidity, circularity, aspect ratio) and corresponding labels.
    """
    # Initialize lists to store metrics and labels
    solidity_lst, circularity_lst, aspect_ratio_lst = [], [], []
    tissue_lst, platform_lst = [], []

    # Ensure the mask image is binary
    image[image != 0] = 255

    # Contour detection using OpenCV
    contours, _ = cv2.findContours(image, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # Calculate metrics for each contour and store results
    for contour in contours:
        solidity, circularity, aspect_ratio = calculate_contour_metrics(contour)
        solidity_lst.append(solidity)
        circularity_lst.append(circularity)
        aspect_ratio_lst.append(aspect_ratio)
        tissue_lst.append(tissue_label)
        platform_lst.append(platform_label)

    return solidity_lst, circularity_lst, aspect_ratio_lst, tissue_lst, platform_lst

# File paths for the tissue mask images
mask_paths_list = [
    '../stereo/COAD/mask.tif',  # COAD
    '../stereo/HCC/mask.tif',   # HCC
    '../stereo/OV/mask.tif',    # OV
    '...',
    '...',
]

# Tissue and platform labels corresponding to each mask image
tissue_list = ['COAD', 'HCC', 'OV']
platform_list = ['Xenium', 'Nanostring', 'Stereo-seq']

# DataFrame to store the results for each metric
df_solidity = pd.DataFrame()
df_circularity = pd.DataFrame()
df_aspect_ratio = pd.DataFrame()

# Process each tissue mask and extract metrics
for i, mask_path in enumerate(mask_paths_list):
    platform = platform_list[i % len(platform_list)]  # Assume platform cycles through the list

    # Read the mask image using tifffile
    with tifffile.TiffFile(mask_path) as tif:
        mask_data = [page.asarray() for page in tif.pages]
        mask_orig = mask_data[0]  # Assuming the first page contains the relevant data
    
    # Extract metrics for the current image
    solidity_lst, circularity_lst, aspect_ratio_lst, tissue_lst, platform_lst = process_image(mask_orig, tissue_list[i], platform)
    
    # Append the results to the corresponding DataFrame
    df_solidity = pd.concat([df_solidity, pd.DataFrame({'tissue': tissue_lst, 'platform': platform_lst, 'value': solidity_lst})], ignore_index=True)
    df_circularity = pd.concat([df_circularity, pd.DataFrame({'tissue': tissue_lst, 'platform': platform_lst, 'value': circularity_lst})], ignore_index=True)
    df_aspect_ratio = pd.concat([df_aspect_ratio, pd.DataFrame({'tissue': tissue_lst, 'platform': platform_lst, 'value': aspect_ratio_lst})], ignore_index=True)