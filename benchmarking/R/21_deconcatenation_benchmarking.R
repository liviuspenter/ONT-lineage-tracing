# benchmark deconcatenation of nanoranger versus skera

library(dplyr)
library(ggplot2)
library(Seurat)
library(SeuratDisk)

# deconcatenation benchmarking provided by Mehdi
deconcat.data <- data.table::fread("./data/benchmarking/deconcatenation/nanoranger_skera.csv")

# generate plots
deconcat.data.filtered <- deconcat.data[which(deconcat.data$reads > 1000), ]
ggplot(deconcat.data.filtered, aes(x = skera_segments, y = nanoranger_segments, fill = log10(reads))) +
  geom_abline(slope = 1) +
  geom_point(shape = 21, size = 1.5) +
  scale_x_continuous("segments skera", limits = c(0, 17)) +
  scale_y_continuous("segments nanoranger", limits = c(0, 17)) +
  scale_fill_gradientn(colours = BuenColors::jdb_palette(name = "brewer_red")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/deconcatenation/20230814_deconcatenation_nanoranger_skera_segments.svg", width = 2, height = 2)

ggplot() +
  geom_col(
    data = deconcat.data %>% group_by(skera_segments) %>% summarize(reads = sum(reads)),
    aes(x = skera_segments, y = reads / 1000000),
    fill = "grey", color = "black"
  ) +
  geom_col(
    data = deconcat.data %>% group_by(nanoranger_segments) %>% summarize(reads = sum(reads)),
    aes(x = nanoranger_segments, y = reads / 1000000),
    fill = "orange", alpha = 0.5, color = "black"
  ) +
  scale_x_continuous("segments", limits = c(0, 17)) +
  scale_y_continuous("million reads", breaks = c(0, 2, 4, 6)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/deconcatenation/20230814_deconcatenation_nanoranger_skera_segments_histogram.svg",
  width = 2, height = 2
)

# gene counts 
gene.counts <- data.table::fread("./data/benchmarking/deconcatenation/gene_counts_PacBio_nanoranger.csv.gz")

ggplot(gene.counts, aes(x = log10_umis_PacBio, y = log10_umis_nanoranger)) +
  geom_abline(slope = 1) +
  ggrastr::rasterize(ggpointdensity::geom_pointdensity(size = 0.5), dpi = 600) +
  scale_x_continuous("log10(UMIs per gene) skera", limits = c(0, 5)) +
  scale_y_continuous("log10(UMIs per gene) nanoranger", limits = c(0, 5)) +
  scale_color_viridis_c() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/deconcatenation/20230814_deconcatenation_nanoranger_skera_UMIs_per_gene.svg",
  width = 2, height = 2
)

ggplot(gene.counts, aes(x = log10_cells_PacBio, y = log10_cells_nanoranger)) +
  geom_abline(slope = 1) +
  ggrastr::rasterize(ggpointdensity::geom_pointdensity(size = 0.5), dpi = 600) +
  scale_x_continuous("log10(cells per gene) skera", limits = c(0, 4)) +
  scale_y_continuous("log10(cells per gene) nanoranger", limits = c(0, 4)) +
  scale_color_viridis_c() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/deconcatenation/20230814_deconcatenation_nanoranger_skera_cells_per_gene.svg",
  width = 2, height = 2
)

# barcode counts
bc.counts <- data.table::fread("./data/benchmarking/deconcatenation/barcode_counts_PacBio_nanoranger.csv.gz")

ggplot(bc.counts, aes(x = log10_umis_PacBio, y = log10_umis_nanoranger)) +
  geom_abline(slope = 1) +
  ggrastr::rasterize(ggpointdensity::geom_pointdensity(size = 0.5), dpi = 600) +
  scale_x_continuous("log10(UMIs per cell) skera", limits = c(2, 5)) +
  scale_y_continuous("log10(UMIs per cell) nanoranger", limits = c(2, 5)) +
  scale_color_viridis_c() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/deconcatenation/20230814_deconcatenation_nanoranger_skera_UMIs_per_cell.svg",
  width = 2, height = 2
)

ggplot(bc.counts, aes(x = log10_genes_PacBio, y = log10_genes_nanoranger)) +
  geom_abline(slope = 1) +
  ggrastr::rasterize(ggpointdensity::geom_pointdensity(size = 0.5), dpi = 600) +
  scale_x_continuous("log10(genes per cell) skera", limits = c(2, 4)) +
  scale_y_continuous("log10(genes per cell) nanoranger", limits = c(2, 4)) +
  scale_color_viridis_c() +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/deconcatenation/20230814_deconcatenation_nanoranger_skera_genes_per_cell.svg",
  width = 2, height = 2
)

# modify feature matrix so it can be processed with Seurat
boo <- data.table::fread("./data/benchmarking/deconcatenation/nanoranger/features.original.tsv.gz", header = F) %>%
  as.data.frame()
boo$V2 <- "Gene Expression"
boo$V3 <- seq(1, nrow(boo))
write.table(boo[, c("V3", "V1", "V2")],
  file = "./data/benchmarking/deconcatenation/nanoranger/features.tsv",
  quote = F, sep = "\t", row.names = F, col.names = F
)

boo <- data.table::fread("./data/benchmarking/deconcatenation/pacbio/features.tsv.original.gz", header = F) %>%
  as.data.frame()
boo$V2 <- "Gene Expression"
boo$V3 <- seq(1, nrow(boo))
write.table(boo[, c("V3", "V1", "V2")],
  file = "./data/Mehdi/revision/deconcatenation/pacbio/features.tsv",
  quote = F, sep = "\t", row.names = F, col.names = F
)

# create Seurat object to compare gene counts between both deconcatenation approaches
seurat.data <- Read10X(data.dir = paste0("./data/benchmarking/deconcatenation/nanoranger/"))
so <- CreateSeuratObject(counts = seurat.data, project = "nanoranger", min.cells = 3, min.features = 200)
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- RenameCells(so, add.cell.id = "nanoranger")
seurat.data <- Read10X(data.dir = paste0("./data/benchmarking/deconcatenation/pacbio/"))
so2 <- CreateSeuratObject(counts = seurat.data, project = "skera", min.cells = 3, min.features = 200)
so2[["percent.mt"]] <- PercentageFeatureSet(so2, pattern = "^MT-")
so2 <- RenameCells(so2, add.cell.id = "skera")

so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so)
so <- RunPCA(so)
so <- RunUMAP(so, dims = 1:30)
so <- FindNeighbors(so)
so <- FindClusters(so, resolution = 0.3)

so2 <- NormalizeData(so2)
so2 <- FindVariableFeatures(so2)
so2 <- ScaleData(so2)
so2 <- RunPCA(so2)
so2 <- RunUMAP(so2, dims = 1:30)
so2 <- FindNeighbors(so2)
so2 <- FindClusters(so2, reso2lution = 0.3)


# map to reference
reference <- LoadH5Seurat("./data/pbmc_multimodal.h5seurat")
so <- SCTransform(so, verbose = FALSE)
so <- CellCycleScoring(so, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
so2 <- SCTransform(so2, verbose = FALSE)
so2 <- CellCycleScoring(so2, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

anchors <- FindTransferAnchors(
  reference = reference,
  query = so,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)
so <- MapQuery(
  anchorset = anchors,
  query = so,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)


anchors <- FindTransferAnchors(
  reference = reference,
  query = so2,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)
so2 <- MapQuery(
  anchorset = anchors,
  query = so2,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

cols <- c(
  "B" = "orange", "CD4 T" = "lightblue", "CD8 T" = "blue", "DC" = "green",
  "Mono" = "darkgreen", "NK" = "purple", "other" = "grey", "other T" = "cyan"
)

p <- DimPlot(so, group.by = "predicted.celltype.l1", reduction = "umap", cols = cols) +
  theme(plot.title = element_blank()) +
  NoAxes() +
  NoLegend()
q <- DimPlot(so2, group.by = "predicted.celltype.l1", reduction = "umap", cols = cols) +
  theme(plot.title = element_blank()) +
  NoAxes() +
  NoLegend()
ggsave("./benchmarking/figures/deconcatenation/20230814_deconcatenation_nanoranger.png",
  width = 4, height = 4,
  dpi = 600, plot = p
)
ggsave("./benchmarking/figures/deconcatenation/20230814_deconcatenation_skera.png",
  width = 4, height = 4,
  dpi = 600, plot = q
)

saveRDS(so, file = "./data/benchmarking/objects//20230814_deconcatenation_nanoranger.rds")
saveRDS(so2, file = "./data/benchmarking/objects/20230814_deconcatenation_pacbio.rds")

so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
so.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

### find differential gene expression between nanoranger and pacbio
so <- readRDS(file = "./data/benchmarking/objects//20230814_deconcatenation_nanoranger.rds")
so2 <- readRDS(file = "./data/benchmarking/objects/20230814_deconcatenation_pacbio.rds")

boo <- merge(so, so2)

boo.sub <- subset(boo, predicted.celltype.l1 == "NK")
boo.sub <- PrepSCTFindMarkers(boo.sub)
NK.markers <- FindMarkers(boo.sub, group.by = "orig.ident", ident.1 = "nanoranger", ident.2 = "skera")

boo.sub <- subset(boo, predicted.celltype.l1 == "Mono")
boo.sub <- PrepSCTFindMarkers(boo.sub)
Mono.markers <- FindMarkers(boo.sub, group.by = "orig.ident", ident.1 = "nanoranger", ident.2 = "skera")

boo.sub <- subset(boo, predicted.celltype.l1 == "B")
boo.sub <- PrepSCTFindMarkers(boo.sub)
B.markers <- FindMarkers(boo.sub, group.by = "orig.ident", ident.1 = "nanoranger", ident.2 = "skera")

# plot selected genes: HLA-C, HLA-E, ITGAL, IGHM, IL17RA
for (feature in c("HLA-C", "HLA-E", "ITGAL", "IGHM", "IL17RA")) {
  p <- FeaturePlot(so, reduction = "umap", features = feature, min.cutoff = "q05", max.cutoff = "q95") +
    scale_color_gradientn(colours = c("grey90", BuenColors::jdb_palette(name = "brewer_green"))) +
    NoLegend() + NoAxes() + theme(plot.title = element_blank())
  ggsave(paste0("./benchmarking/figures/deconcatenation/20230824_UMAP_nanoranger_", feature, ".png"), width = 4, height = 4, plot = p)
  q <- FeaturePlot(so2, reduction = "umap", features = feature, min.cutoff = "q05", max.cutoff = "q95") +
    scale_color_gradientn(colours = c("grey90", BuenColors::jdb_palette(name = "brewer_green"))) +
    NoLegend() + NoAxes() + theme(plot.title = element_blank())
  ggsave(paste0("./benchmarking/figures/deconcatenation/20230824_UMAP_skera_", feature, ".png"), width = 4, height = 4, plot = q)
}

svglite::svglite("./benchmarking/figures/deconcatenation/brewer_green.svg", width = 2, height = 0.5)
BuenColors::jdb_palette(name = "brewer_green", type = "continuous")
dev.off()
