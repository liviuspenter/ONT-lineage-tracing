# create seurat object for comparison GoT versus nanoranger

# AML1022.1: standard 10x cDNA
# AML1022.G: 10x cDNA with RT spike-in primers for DNMT3A, RUNX1 and SF3B1 loci

library(ComplexHeatmap)
library(ggplot2)
library(Seurat)
library(SummarizedExperiment)

bm <- readRDS(file = "./data/reference/bm.reference.rds")
bm[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(
  object = bm[["spca.annoy.neighbors"]],
  file = "./data/reference/reftmp.idx"
)

seurat.data <- Read10X(data.dir = paste0("./data/benchmarking/GoT_nanoranger/AML1022.1/filtered_feature_bc_matrix/"))
so <- CreateSeuratObject(counts = seurat.data, project = "AML1022.1", min.cells = 3, min.features = 200)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- RenameCells(so, add.cell.id = "AML1022.1")
seurat.data <- Read10X(data.dir = paste0("./data/benchmarking/GoT_nanoranger/AML1022.1_GoT_with_amplicon/filtered_feature_bc_matrix/"))
so2 <- CreateSeuratObject(counts = seurat.data, project = "AML1022.1G", min.cells = 3, min.features = 200)
so2[["percent.mt"]] <- PercentageFeatureSet(so2, pattern = "^MT-")
so2 <- RenameCells(so2, add.cell.id = "AML1022.1G")
so <- merge(so, list(so2))

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
AML1022 <- MapQuery(
  anchorset = anchors, so, reference = bm,
  refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap"
)
saveRDS(file = "./data/benchmarking/objects/20230812_AML1022_nanoranger_GoT.rds", AML1022)
