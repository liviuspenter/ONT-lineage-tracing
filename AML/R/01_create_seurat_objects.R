# create individual Seurat objects of AML cases for easier handling

library(ComplexHeatmap)
library(ggplot2)
library(Seurat)
library(SummarizedExperiment)

bm <- readRDS(file = "./data/reference/bm.reference.rds")
bm[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(
  object = bm[["spca.annoy.neighbors"]],
  file = "./data/reference/reftmp.idx"
)

## AML1002.13
seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool104-2_10/filtered_feature_bc_matrix/"))
so <- CreateSeuratObject(counts = seurat.data, project = "AML1002.1", min.cells = 3, min.features = 200)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- RenameCells(so, add.cell.id = "AML1002.1")
seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool104-2_12/filtered_feature_bc_matrix/"))
so2 <- CreateSeuratObject(counts = seurat.data, project = "AML1002.3", min.cells = 3, min.features = 200)
so2[["percent.mt"]] <- PercentageFeatureSet(so2, pattern = "^MT-")
so2 <- RenameCells(so2, add.cell.id = "AML1002.3")
so <- merge(so, so2)
so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so)
so <- RunPCA(so)
so <- RunUMAP(so, dims = 1:30)
so <- FindNeighbors(so)
so <- FindClusters(so, resolution = 0.3)

anchors <- FindTransferAnchors(
  reference = bm, query = so, k.filter = NA, reference.reduction = "spca",
  reference.neighbors = "spca.annoy.neighbors", dims = 1:50
)
AML1002 <- MapQuery(
  anchorset = anchors, so, reference = bm,
  refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap"
)
saveRDS(file = "./data/AML/objects/20220911_AML1002.rds", AML1002)

## AML1002 all
seurat.data <- Read10X(data.dir = paste0("./data/AML/objects/Pool104-2_10/filtered_feature_bc_matrix/"))
so <- CreateSeuratObject(counts = seurat.data, project = "AML1002.1", min.cells = 3, min.features = 200)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- RenameCells(so, add.cell.id = "AML1002.1")
seurat.data <- Read10X(data.dir = paste0("./data/AML/objects/Pool104-2_11/filtered_feature_bc_matrix/"))
so2 <- CreateSeuratObject(counts = seurat.data, project = "AML1002.2", min.cells = 3, min.features = 200)
so2[["percent.mt"]] <- PercentageFeatureSet(so2, pattern = "^MT-")
so2 <- RenameCells(so2, add.cell.id = "AML1002.2")
seurat.data <- Read10X(data.dir = paste0("./data/AML/objects/Pool104-2_12/filtered_feature_bc_matrix/"))
so3 <- CreateSeuratObject(counts = seurat.data, project = "AML1002.3", min.cells = 3, min.features = 200)
so3[["percent.mt"]] <- PercentageFeatureSet(so3, pattern = "^MT-")
so3 <- RenameCells(so3, add.cell.id = "AML1002.3")
so <- merge(so, list(so2, so3))
so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so)
so <- RunPCA(so)
so <- RunUMAP(so, dims = 1:30)
so <- FindNeighbors(so)
so <- FindClusters(so, resolution = 0.3)

anchors <- FindTransferAnchors(
  reference = bm, query = so, k.filter = NA, reference.reduction = "spca",
  reference.neighbors = "spca.annoy.neighbors", dims = 1:50
)
AML1002 <- MapQuery(
  anchorset = anchors, so, reference = bm,
  refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap"
)
saveRDS(file = "./data/AML/objects/20220930_AML1002_all.rds", AML1002)

# AML1007.135
seurat.data <- Read10X(data.dir = paste0("./data/AML/objects/Pool96_4/filtered_feature_bc_matrix/"))
rownames(x = seurat.data[["Antibody Capture"]]) <- gsub(pattern = "*_CITEseq", replacement = "", rownames(seurat.data[["Antibody Capture"]]))
so <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = "AML1007.1", min.cells = 3, min.features = 200)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[, colnames(x = so)])
so <- NormalizeData(so, assay = "ADT", normalization.method = "CLR")
so <- RenameCells(so, add.cell.id = "AML1007.1")
seurat.data <- Read10X(data.dir = paste0("./data/AML/objects/Pool96_6/filtered_feature_bc_matrix/"))
rownames(x = seurat.data[["Antibody Capture"]]) <- gsub(pattern = "*_CITEseq", replacement = "", rownames(seurat.data[["Antibody Capture"]]))
so2 <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = "AML1007.3", min.cells = 3, min.features = 200)
so2[["percent.mt"]] <- PercentageFeatureSet(so2, pattern = "^MT-")
so2[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[, colnames(x = so2)])
so2 <- NormalizeData(so2, assay = "ADT", normalization.method = "CLR")
so2 <- RenameCells(so2, add.cell.id = "AML1007.3")
seurat.data <- Read10X(data.dir = paste0("./data/AML/objects/Pool96_8/filtered_feature_bc_matrix/"))
rownames(x = seurat.data[["Antibody Capture"]]) <- gsub(pattern = "*_CITEseq", replacement = "", rownames(seurat.data[["Antibody Capture"]]))
so3 <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = "AML1007.5", min.cells = 3, min.features = 200)
so3[["percent.mt"]] <- PercentageFeatureSet(so3, pattern = "^MT-")
so3[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[, colnames(x = so3)])
so3 <- NormalizeData(so3, assay = "ADT", normalization.method = "CLR")
so3 <- RenameCells(so3, add.cell.id = "AML1007.5")
so <- merge(so, list(so2, so3))
so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so)
so <- RunPCA(so)
so <- RunUMAP(so, dims = 1:30)
so <- FindNeighbors(so)
so <- FindClusters(so, resolution = 0.3)

anchors <- FindTransferAnchors(
  reference = bm, query = so, k.filter = NA, reference.reduction = "spca",
  reference.neighbors = "spca.annoy.neighbors", dims = 1:50
)
AML1007 <- MapQuery(
  anchorset = anchors, so, reference = bm,
  refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap"
)
saveRDS(file = "./data/AML/objects/20220825_AML1007.rds", AML1007)

# AML1010.1
seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool93_35/filtered_feature_bc_matrix/"))
rownames(x = seurat.data[["Antibody Capture"]]) <- gsub(pattern = "*_CITEseq", replacement = "", rownames(seurat.data[["Antibody Capture"]]))
so <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = "AML1010.1", min.cells = 3, min.features = 200)
message(paste(ncol(so), " barcodes"))
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[, colnames(x = so)])
so <- NormalizeData(so, assay = "ADT", normalization.method = "CLR")
so <- RenameCells(so, add.cell.id = "AML1010.1")

anchors <- FindTransferAnchors(
  reference = bm, query = so, k.filter = NA, reference.reduction = "spca",
  reference.neighbors = "spca.annoy.neighbors", dims = 1:50
)
AML1010.1 <- MapQuery(
  anchorset = anchors, so, reference = bm,
  refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap"
)
saveRDS(file = "./data/AML/objects/20220711_AML1010_1.rds", AML1010.1)

# AML1019.1
seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool96_30/filtered_feature_bc_matrix/"))
so <- CreateSeuratObject(counts = seurat.data, project = "AML1019.1", min.cells = 3, min.features = 200)
message(paste(ncol(so), " barcodes"))
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- RenameCells(so, add.cell.id = "AML1019.1")
so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so)
so <- RunPCA(so)
so <- RunUMAP(so, dims = 1:30)
so <- FindNeighbors(so)
so <- FindClusters(so, resolution = 0.3)

anchors <- FindTransferAnchors(
  reference = bm, query = so, k.filter = NA, reference.reduction = "spca",
  reference.neighbors = "spca.annoy.neighbors", dims = 1:50
)
AML1019.1 <- MapQuery(
  anchorset = anchors, so, reference = bm,
  refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap"
)
saveRDS(file = "./data/AML/objects/20220822_AML1019_1.rds", AML1019.1)

# AML1022.1
seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool97_3/filtered_feature_bc_matrix/"))
rownames(x = seurat.data[["Antibody Capture"]]) <- gsub(pattern = "*_CITEseq", replacement = "", rownames(seurat.data[["Antibody Capture"]]))
so <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = "AML1022.1", min.cells = 3, min.features = 200)
message(paste(ncol(so), " barcodes"))
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[, colnames(x = so)])
so <- NormalizeData(so, assay = "ADT", normalization.method = "CLR")
so <- RenameCells(so, add.cell.id = "AML1022.1")

anchors <- FindTransferAnchors(
  reference = bm, query = so, k.filter = NA, reference.reduction = "spca",
  reference.neighbors = "spca.annoy.neighbors", dims = 1:50
)
AML1022.1 <- MapQuery(
  anchorset = anchors, so, reference = bm,
  refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap"
)
saveRDS(file = "./data/AML/objects/20220729_AML1022_1.rds", AML1022.1)


# AML1026.2 - for mtDNA analysis
seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool96_16/filtered_feature_bc_matrix/"))
rownames(x = seurat.data[["Antibody Capture"]]) <- gsub(pattern = "*_CITEseq", replacement = "", rownames(seurat.data[["Antibody Capture"]]))
so <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = "AML1026.2", min.cells = 3, min.features = 200)
message(paste(ncol(so), " barcodes"))
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[, colnames(x = so)])
so <- NormalizeData(so, assay = "ADT", normalization.method = "CLR")
so <- RenameCells(so, add.cell.id = "AML1026.2")

anchors <- FindTransferAnchors(
  reference = bm, query = so, k.filter = NA, reference.reduction = "spca",
  reference.neighbors = "spca.annoy.neighbors", dims = 1:50
)
AML1026.2 <- MapQuery(
  anchorset = anchors, so, reference = bm,
  refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap"
)
saveRDS(file = "./data/AML/objects/20220711_AML1026_2.rds", AML1026.2)


# AML3003.14
seurat.data <- Read10X(data.dir = paste0("/Users/shaka87/dfci/10026/10x/data/Pool104-2_29/filtered_feature_bc_matrix/"))
rownames(x = seurat.data[["Antibody Capture"]]) <- gsub(pattern = "*_CITEseq", replacement = "", rownames(seurat.data[["Antibody Capture"]]))
so <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = "AML3003.1", min.cells = 3, min.features = 200)
so[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[, colnames(x = so)])
so <- NormalizeData(so, assay = "ADT", normalization.method = "CLR")
so <- RenameCells(so, add.cell.id = "AML3003.1")

seurat.data <- Read10X(data.dir = paste0("/Users/shaka87/dfci/10026/10x/data/Pool104-2_32/filtered_feature_bc_matrix/"))
rownames(x = seurat.data[["Antibody Capture"]]) <- gsub(pattern = "*_CITEseq", replacement = "", rownames(seurat.data[["Antibody Capture"]]))
so2 <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = "AML3003.4", min.cells = 3, min.features = 200)
so2[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[, colnames(x = so2)])
so2 <- NormalizeData(so2, assay = "ADT", normalization.method = "CLR")
so2 <- RenameCells(so2, add.cell.id = "AML3003.4")
so <- merge(so, so2)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

anchors <- FindTransferAnchors(
  reference = bm, query = so, k.filter = NA, reference.reduction = "spca",
  reference.neighbors = "spca.annoy.neighbors", dims = 1:50
)
AML3003.14 <- MapQuery(
  anchorset = anchors, so, reference = bm,
  refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap"
)
saveRDS(file = "./data/AML/objects/202208_AML3003_14.rds", AML3003.14)

# AML3005.13
seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool106-1_8/filtered_feature_bc_matrix/"))
so <- CreateSeuratObject(counts = seurat.data, project = "AML3005.1", min.cells = 3, min.features = 200)
so <- RenameCells(so, add.cell.id = "AML3005.1")
seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool106-1_10/filtered_feature_bc_matrix/"))
so2 <- CreateSeuratObject(counts = seurat.data, project = "AML3005.3", min.cells = 3, min.features = 200)
so2 <- RenameCells(so2, add.cell.id = "AML3005.3")
so <- merge(so, so2)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

anchors <- FindTransferAnchors(
  reference = bm, query = so, k.filter = NA, reference.reduction = "spca",
  reference.neighbors = "spca.annoy.neighbors", dims = 1:50
)
AML3005.13 <- MapQuery(
  anchorset = anchors, so, reference = bm,
  refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap"
)

AML3005.13$manual.cluster <- "CD4 T cell"
AML3005.13$manual.cluster[which(AML3005.13$seurat_clusters %in% c(6, 8, 9, 10))] <- "Erythropoiesis"
AML3005.13$manual.cluster[which(AML3005.13$seurat_clusters %in% c(5, 11))] <- "Progenitor"
AML3005.13$manual.cluster[which(AML3005.13$seurat_clusters %in% c(12))] <- "B lymphopoiesis"
AML3005.13$manual.cluster[which(AML3005.13$seurat_clusters %in% c(3))] <- "NK cell"
AML3005.13$manual.cluster[which(AML3005.13$seurat_clusters %in% c(7))] <- "CD8 T cell"

saveRDS(file = "./data/AML/objects/20220718_AML3005_13.rds", AML3005.13)

### AML8007
seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool106-1_15/filtered_feature_bc_matrix/"))
rownames(x = seurat.data[["Antibody Capture"]]) <- gsub(pattern = "*_CITEseq", replacement = "", rownames(seurat.data[["Antibody Capture"]]))
so <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = "AML8007.1", min.cells = 3, min.features = 200)
so[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[, colnames(x = so)])
so <- NormalizeData(so, assay = "ADT", normalization.method = "CLR")
so <- RenameCells(so, add.cell.id = "AML8007.1")
seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool106-1_18/filtered_feature_bc_matrix/"))
rownames(x = seurat.data[["Antibody Capture"]]) <- gsub(pattern = "*_CITEseq", replacement = "", rownames(seurat.data[["Antibody Capture"]]))
so2 <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = "AML8007.4", min.cells = 3, min.features = 200)
so2[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[, colnames(x = so2)])
so2 <- NormalizeData(so2, assay = "ADT", normalization.method = "CLR")
so2 <- RenameCells(so2, add.cell.id = "AML8007.4")
so <- merge(so, so2)
so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so)
so <- RunPCA(so)
so <- RunUMAP(so, dims = 1:30)
so <- FindNeighbors(so)
so <- FindClusters(so, resolution = 0.3)

anchors <- FindTransferAnchors(
  reference = bm, query = so, k.filter = NA, reference.reduction = "spca",
  reference.neighbors = "spca.annoy.neighbors", dims = 1:50
)
AML8007 <- MapQuery(
  anchorset = anchors, so, reference = bm,
  refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap"
)
saveRDS(file = "./data/AML/objects/20220911_AML8007.rds", AML8007)

### AML8007 all
seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool106-1_15/filtered_feature_bc_matrix/"))
rownames(x = seurat.data[["Antibody Capture"]]) <- gsub(pattern = "*_CITEseq", replacement = "", rownames(seurat.data[["Antibody Capture"]]))
so <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = "AML8007.1", min.cells = 3, min.features = 200)
so[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[, colnames(x = so)])
so <- NormalizeData(so, assay = "ADT", normalization.method = "CLR")
so <- RenameCells(so, add.cell.id = "AML8007.1")
seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool106-1_16/filtered_feature_bc_matrix/"))
rownames(x = seurat.data[["Antibody Capture"]]) <- gsub(pattern = "*_CITEseq", replacement = "", rownames(seurat.data[["Antibody Capture"]]))
so.2 <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = "AML8007.2", min.cells = 3, min.features = 200)
so.2[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[, colnames(x = so.2)])
so.2 <- NormalizeData(so.2, assay = "ADT", normalization.method = "CLR")
so.2 <- RenameCells(so.2, add.cell.id = "AML8007.2")
seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool106-1_17/filtered_feature_bc_matrix/"))
rownames(x = seurat.data[["Antibody Capture"]]) <- gsub(pattern = "*_CITEseq", replacement = "", rownames(seurat.data[["Antibody Capture"]]))
so.3 <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = "AML8007.3", min.cells = 3, min.features = 200)
so.3[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[, colnames(x = so.3)])
so.3 <- NormalizeData(so.3, assay = "ADT", normalization.method = "CLR")
so.3 <- RenameCells(so.3, add.cell.id = "AML8007.3")
seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool106-1_18/filtered_feature_bc_matrix/"))
rownames(x = seurat.data[["Antibody Capture"]]) <- gsub(pattern = "*_CITEseq", replacement = "", rownames(seurat.data[["Antibody Capture"]]))
so.4 <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = "AML8007.4", min.cells = 3, min.features = 200)
so.4[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[, colnames(x = so.4)])
so.4 <- NormalizeData(so.4, assay = "ADT", normalization.method = "CLR")
so.4 <- RenameCells(so.4, add.cell.id = "AML8007.4")

so <- merge(so, list(so.2, so.3, so.4))
so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so)
so <- RunPCA(so)
so <- RunUMAP(so, dims = 1:30)
so <- FindNeighbors(so)
so <- FindClusters(so, resolution = 0.3)

anchors <- FindTransferAnchors(
  reference = bm, query = so, k.filter = NA, reference.reduction = "spca",
  reference.neighbors = "spca.annoy.neighbors", dims = 1:50
)
AML8007 <- MapQuery(
  anchorset = anchors, so, reference = bm,
  refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap"
)
saveRDS(file = "./data/AML/objects/20220929_AML8007_all.rds", AML8007)

### de-novo AML (aka FLT3-ITD cases because they all have FLT3-ITD)
# de-novo AML1: Pool129-1_3
# de-novo AML2: Pool129-1_5
# de-novo AML3: Pool129_1_4

seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool129-1_3/filtered_feature_bc_matrix/"))
so <- CreateSeuratObject(counts = seurat.data, project = "FLT3-ITD1", min.cells = 3, min.features = 200)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- RenameCells(so, add.cell.id = "FLT3-ITD1")
seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool129-1_4/filtered_feature_bc_matrix/"))
so2 <- CreateSeuratObject(counts = seurat.data, project = "FLT3-ITD2", min.cells = 3, min.features = 200)
so2[["percent.mt"]] <- PercentageFeatureSet(so2, pattern = "^MT-")
so2 <- RenameCells(so2, add.cell.id = "FLT3-ITD2")
seurat.data <- Read10X(data.dir = paste0("./data/AML/Pool129-1_5/filtered_feature_bc_matrix/"))
so3 <- CreateSeuratObject(counts = seurat.data, project = "FLT3-ITD3", min.cells = 3, min.features = 200)
so3[["percent.mt"]] <- PercentageFeatureSet(so3, pattern = "^MT-")
so3 <- RenameCells(so3, add.cell.id = "FLT3-ITD3")

so <- merge(so, list(so2, so3))
so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so)
so <- RunPCA(so)
so <- RunUMAP(so, dims = 1:30)
so <- FindNeighbors(so)
so <- FindClusters(so, resolution = 0.3)

anchors <- FindTransferAnchors(
  reference = bm, query = so, k.filter = NA, reference.reduction = "spca",
  reference.neighbors = "spca.annoy.neighbors", dims = 1:50
)
FLT3_ITD <- MapQuery(
  anchorset = anchors, so, reference = bm,
  refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap"
)
saveRDS(file = "./data/AML/objects/20230523_FLT3-ITD.rds", FLT3_ITD)
