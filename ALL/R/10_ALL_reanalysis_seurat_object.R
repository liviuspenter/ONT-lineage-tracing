# reanalysis of data from Witkowski et al., Cancer Cell 2020 and Caron et al., Scientific Reports 2020

# create Seurat object and annotate using reference bone marrow dataset

library(dplyr)
library(ggplot2)
library(numbat)
library(Seurat)

samples <- c(
  "ETV001", "ETV002", "ETV003", "ETV004", "ETV005", "PH001", "PH002", "ETV6-RUNX1_1", "ETV6-RUNX1_2", "ETV6-RUNX1_3", "ETV6-RUNX1_4",
  "HHD_1", "HHD_2", "PRE-T_1", "PRE-T_2", "PH001_relapse", "PH002_relapse"
)

# create Seurat data
for (s in samples) {
  message(s)
  seurat.data <- Read10X(data.dir = paste0("./data/ALL/reanalysis/", s, "/filtered_feature_bc_matrix/"))

  so <- CreateSeuratObject(counts = seurat.data, project = s, min.cells = 3, min.features = 200)
  message(paste(s, " with ", ncol(so), " barcodes"))
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

  so <- subset(so, subset = nCount_RNA < 20000 & nFeature_RNA < 4000 & percent.mt < 20)

  so <- RenameCells(so, add.cell.id = s)

  if (s == samples[1]) {
    so.combined <- list(so)
  } else {
    so.combined <- c(so.combined, so)
  }
}

ALL.validation <- merge(x = so.combined[[1]], y = so.combined[2:length(samples)], project = "ALL.validation")
ALL.validation <- FindVariableFeatures(ALL.validation)
ALL.validation <- NormalizeData(ALL.validation)
ALL.validation <- ScaleData(ALL.validation)
ALL.validation <- RunPCA(ALL.validation)
set.seed(1987)
ALL.validation <- RunUMAP(ALL.validation, dims = 1:30)
ALL.validation <- FindNeighbors(ALL.validation)
ALL.validation <- FindClusters(ALL.validation, resolution = 0.3)

bm <- readRDS(file = "./data/reference/bm.reference.rds")
bm[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(
  object = bm[["spca.annoy.neighbors"]],
  file = "./data/reference/reftmp.idx"
)
Misc(bm[["wnn.umap"]], slot = "model")$num_precomputed_nns <- 1

anchors <- FindTransferAnchors(
  reference = bm, query = ALL.validation, k.filter = NA, reference.reduction = "spca",
  reference.neighbors = "spca.annoy.neighbors", dims = 1:50
)
ALL.validation <- MapQuery(
  anchorset = anchors, ALL.validation, reference = bm,
  refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap"
)

saveRDS(file = "./data/ALL/objects/20221027_ALL_validation.rds", ALL.validation)
