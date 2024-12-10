########################## file information #############################
# Author: Pengfei Ren
# Date: 2024-12-04 21:28:21
# LastEditTime: 2024-12-04 23:33:42
# Description: 
# FilePath: /undefined/Users/morsouron/Desktop/Project/benchmark/code/6_scrna.py
#########################################################################

rm_doublet = function(obj){
    double_percents = c(2.3,3.1,3.9,4.6,5.4,6.1,6.9,7.6)
    # doublelet 
    n_trans = round(ncol(obj)/1000)-2
    if (n_trans>=8){
        nExp <- round(ncol(obj) * double_percents[8]/100) 
    } else {
        nExp <- round(ncol(obj) * double_percents[n_trans]/100) 
    }
    obj = NormalizeData(obj, verbose = FALSE)
    obj = FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    obj <- ScaleData(obj, verbose = FALSE)
    obj <- RunPCA(obj, npcs = 50, verbose = FALSE, features = VariableFeatures(obj))

    sweep.res.list <- paramSweep_v3(obj, PCs = 1:30, sct = FALSE,num.cores = 1)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pk = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    obj <- doubletFinder_v3(obj, pN = 0.25, pK = pk, nExp = nExp, PCs = 1:30)
    colnames(obj@meta.data)[str_detect(colnames(obj@meta.data),'classifications')] = 'DF.classifications'
    obj = obj[,obj$DF.classifications!='Doublet']
    return(obj)
}

COAD = Read10X_h5('~/benchmark/data/raw/scrna/expr/COAD/outs/filtered_feature_bc_matrix.h5')
COAD = CreateSeuratObject(counts = COAD, project = "COAD",min.cells = 10, min.features = 20)
COAD[["percent.mt"]] <- PercentageFeatureSet(COAD, pattern = "^MT-")
COAD=COAD[,COAD$nFeature_RNA>=500 & COAD$nCount_RNA>=1000 & COAD$percent.mt<=10]
COAD=COAD[,COAD$nFeature_RNA<=5000 & COAD$nCount_RNA<=25000]
COAD = rm_doublet(COAD)
COAD <- NormalizeData(COAD, verbose = FALSE)
COAD <- FindVariableFeatures(COAD, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
COAD <- ScaleData(COAD, verbose = FALSE)
COAD <- RunPCA(COAD, features = VariableFeatures(COAD), npcs = 30)
COAD = RunUMAP(COAD, reduction = "pca", dims = 1:30) %>%  
        FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
        FindClusters(resolution = 0.6, verbose = FALSE)