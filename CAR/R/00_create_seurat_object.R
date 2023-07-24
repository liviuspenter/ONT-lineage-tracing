# create Seurat object

library(ggplot2)
library(Seurat)

seurat.data = Read10X(data.dir = paste0('./data/CAR/Pool97_6/filtered_feature_bc_matrix/'))

rownames(x = seurat.data[['Antibody Capture']]) = gsub(pattern = '*_CITEseq', replacement = '', rownames(seurat.data[['Antibody Capture']])) 
so <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = '97_6', min.cells = 3, min.features = 200)
message(paste(ncol(so), ' barcodes'))
so[['percent.mt']] <- PercentageFeatureSet(so, pattern = '^MT-')
so[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[,colnames(x = so)])
so = NormalizeData(so, assay = 'ADT', normalization.method = 'CLR')
DefaultAssay(so) = 'RNA'
so = NormalizeData(so)
so = ScaleData(so)
so = FindVariableFeatures(so)
so = RunPCA(so)
so = FindNeighbors(so)
so = FindClusters(so)
so = RunUMAP(so, dims = 1:20)
saveRDS(file='./data/CAR/objects/20220819_CAR_so.rds', so)