# cell tyope number
res = list()
for (tissue in c('COAD','HCC','OV')){
    res_tissue = list()
    seurat_obj = readRDS(paste0('~/benchmark/data/data/scrna/seurat/', tolower(tissue),'_annotated.rds'))
    ref_cts = unique(seurat_obj$major)
    for (platform in platforms){
        
        path = file.path('~/benchmark/res/cell/annotation/',str_split(platform,'_')[[1]][1],tissue)
        
        # List all files in the folder
        if (str_detect(platform,'bin')){
            files <- list.files(path = path, pattern = "bin.txt")
        } else if (str_detect(platform,'cell')){
            files <- list.files(path = path, pattern = "cell.txt")
        } else {
            files <- list.files(path = path, pattern = "*.txt")
        }

        # Read all the remaining files
        data_list <- lapply(files, function(file) {
          df = read.table(file.path(path, file), header = FALSE, stringsAsFactors = FALSE)
        })

        # Combine the first column of all files into a data frame
        merged_data <- do.call(cbind, data_list)
        merged_data = as.data.frame(merged_data)
        
        # number of cell types
        ncts = c()
        for (i in 1:ncol(merged_data)){
            ncts = c(ncts,length(unique(merged_data[,i])))
        }
        res_tissue[[platform]] <- ncts
    }
    res_tissue = as.data.frame(res_tissue %>% reduce(cbind))
    colnames(res_tissue) = platforms
    res_tissue$tissue = tissue
    res[[tissue]] = res_tissue
}

# percentage of cells annotated with the same cell type
res = list()
for (tissue in c('COAD','HCC','OV')){
    res_tissue = list()
    for (platform in platforms){
        path = file.path('~/benchmark/res/cell/annotation/',str_split(platform,'_')[[1]][1],tissue)
        
        # List all files in the folder
        if (str_detect(platform,'bin')){
            files <- list.files(path = path, pattern = "bin.txt")
        } else if (str_detect(platform,'cell')){
            files <- list.files(path = path, pattern = "cell.txt")
        } else {
            files <- list.files(path = path, pattern = "*.txt")
        }

        # Read all the remaining files
        data_list <- lapply(files, function(file) {
          read.table(file.path(path, file), header = FALSE, stringsAsFactors = FALSE)
        })

        # Combine the first column of all files into a data frame
        merged_data <- do.call(cbind, data_list)

        # Calculate the maximum number of tools with the same annotation for each row
        max_counts <- apply(merged_data, 1, function(x) {
          max(table(x))
        })

        # Count the distribution of the maximum values across all cells
        res_tissue[[platform]] <- as.vector(table(max_counts))/length(max_counts)
    }
    res_tissue = as.data.frame(res_tissue %>% reduce(cbind))
    colnames(res_tissue) = platforms
    res_tissue$tissue = tissue
    res_tissue$Method = 1:5
    res[[tissue]] = res_tissue
}

# entropy of cell number
res = list()
for (tissue in c('COAD','HCC','OV')){
    res_tissue = list()
    seurat_obj = readRDS(paste0('~/benchmark/data/data/scrna/seurat/', tolower(tissue),'_annotated.rds'))
    ref_cts = unique(seurat_obj$major)
    for (platform in platforms){
        
        path = file.path('~/benchmark/res/cell/annotation/',str_split(platform,'_')[[1]][1],tissue)
        
        # List all files in the folder
        if (str_detect(platform,'bin')){
            files <- list.files(path = path, pattern = "bin.txt")
        } else if (str_detect(platform,'cell')){
            files <- list.files(path = path, pattern = "cell.txt")
        } else {
            files <- list.files(path = path, pattern = "*.txt")
        }

        # Read all the remaining files
        data_list <- lapply(files, function(file) {
          df = read.table(file.path(path, file), header = FALSE, stringsAsFactors = FALSE)
        })

        # Combine the first column of all files into a data frame
        merged_data <- do.call(cbind, data_list)
        merged_data = as.data.frame(merged_data)
        
        # entropy
        freqs = list()
        for (i in 1:ncol(merged_data)){
            merged_data[,i] = factor(merged_data[,i],levels = ref_cts)
            freqs[[i]] = table(merged_data[,i])
        }
        freqs = reduce(freqs,cbind)
        entropies = apply(freqs,1,entropy)
        res_tissue[[platform]] <- entropies
    }
    res_tissue = as.data.frame(res_tissue %>% reduce(cbind))
    colnames(res_tissue) = platforms
    res_tissue$ct = rownames(res_tissue)
    res_tissue$tissue = tissue
    res[[tissue]] = res_tissue
}