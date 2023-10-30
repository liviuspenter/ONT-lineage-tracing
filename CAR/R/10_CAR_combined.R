# systematic analysis of 8 CAR amplicon libraries

library(ggplot2)
library(ggsignif)
library(Seurat)
library(SeuratDisk)

source("./R/20210213_read_tcr_bcr.R")

T.cell.subsets <- c(
  "CD4 TCM", "CD4 TEM", "CD4 Proliferating", "Treg",
  "CD8 Naive", "CD8 TCM", "CD8 TEM", "CD8 Proliferating"
)

T.cell.subset.colors <- c(
  "CD4 TCM" = RColorBrewer::brewer.pal(name = "Blues", n = 6)[3],
  "CD4 TEM" = RColorBrewer::brewer.pal(name = "Blues", n = 6)[4],
  "CD4 Proliferating" = RColorBrewer::brewer.pal(name = "Blues", n = 6)[5],
  "Treg" = RColorBrewer::brewer.pal(name = "Blues", n = 6)[6],
  "CD8 Naive" = RColorBrewer::brewer.pal(name = "Reds", n = 6)[1],
  "CD8 TCM" = RColorBrewer::brewer.pal(name = "Reds", n = 6)[2],
  "CD8 TEM" = RColorBrewer::brewer.pal(name = "Reds", n = 6)[3],
  "CD8 Proliferating" = RColorBrewer::brewer.pal(name = "Reds", n = 6)[4],
  "NK" = RColorBrewer::brewer.pal(name = "Reds", n = 6)[5],
  "NK Proliferating" = RColorBrewer::brewer.pal(name = "Reds", n = 6)[6]
)

libraries <- c("Pool97_6", "Pool104-2_4", "Pool104-2_6", "Pool104-2_7", "Pool106-2_3", "Pool106-2_5", "Pool106-2_16", "Pool106-2_18")

# create combined Seurat object
so.list <- list()
for (s in libraries) {
  message(s)
  seurat.data <- Read10X(paste0("./data/CAR/", s, "/filtered_feature_bc_matrix/"))
  if (length(seurat.data) == 2) {
    CAR.data <- CreateSeuratObject(seurat.data$`Gene Expression`, project = s, min.cells = 3, min.features = 200)
    CAR.data[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[, colnames(x = CAR.data)])
    CAR.data <- NormalizeData(CAR.data, assay = "ADT", normalization.method = "CLR")
  } else {
    CAR.data <- CreateSeuratObject(seurat.data, project = s, min.cells = 3, min.features = 200)
  }

  CAR.data[["percent.mt"]] <- PercentageFeatureSet(CAR.data, pattern = "^MT-")
  CAR.data <- RenameCells(CAR.data, add.cell.id = s)

  if (s == libraries[1]) {
    so.list <- CAR.data
  } else {
    so.list <- c(so.list, CAR.data)
  }
}
CAR.data <- merge(x = so.list[[1]], y = so.list[2:length(so.list)], project = "CAR.data")
CAR.data$timepoint <- ifelse(CAR.data$orig.ident %in% c("Pool97_6", "Pool104-2_6", "Pool106-2_5", "Pool106-2_16"), "IP", "D7")

CAR.data <- NormalizeData(CAR.data)
CAR.data <- FindVariableFeatures(CAR.data)
CAR.data <- ScaleData(CAR.data)
CAR.data <- RunPCA(CAR.data)
CAR.data <- RunUMAP(CAR.data, dims = 1:30)
CAR.data <- FindNeighbors(CAR.data)
CAR.data <- FindClusters(CAR.data, resolution = 0.4)

# read combined native TCR data
tcr.data <- data.frame()
tcr.data <- rbind(tcr.data, read_tcr_bcr("./data/CAR/Pool115-2_29/all_contig_annotations.csv",
  "./data/CAR/Pool115-2_29/clonotypes.csv",
  prefix = "Pool104-2_6", bcrtcr = "TCR"
))
tcr.data <- rbind(tcr.data, read_tcr_bcr("./data/CAR/Pool115-1_33/all_contig_annotations.csv",
  "./data/CAR/Pool115-1_33/clonotypes.csv",
  prefix = "Pool106-2_5", bcrtcr = "TCR"
))
tcr.data <- rbind(tcr.data, read_tcr_bcr("./data/CAR/Pool115-2_36/all_contig_annotations.csv",
  "./data/CAR/Pool115-2_36/clonotypes.csv",
  prefix = "Pool106-2_16", bcrtcr = "TCR"
))
tcr.data <- rbind(tcr.data, read_tcr_bcr("./data/CAR/Pool115-2_27/all_contig_annotations.csv",
  "./data/CAR/Pool115-2_27/clonotypes.csv",
  prefix = "Pool104-2_4", bcrtcr = "TCR"
))
tcr.data <- rbind(tcr.data, read_tcr_bcr("./data/CAR/Pool115-2_30/all_contig_annotations.csv",
  "./data/CAR/Pool115-2_30/clonotypes.csv",
  prefix = "Pool104-2_7", bcrtcr = "TCR"
))
tcr.data <- rbind(tcr.data, read_tcr_bcr("./data/CAR/Pool115-1_31/all_contig_annotations.csv",
  "./data/CAR/Pool115-1_31/clonotypes.csv",
  prefix = "Pool106-2_3", bcrtcr = "TCR"
))
tcr.data <- rbind(tcr.data, read_tcr_bcr("./data/CAR/Pool115-2_38/all_contig_annotations.csv",
  "./data/CAR/Pool115-2_38/clonotypes.csv",
  prefix = "Pool106-2_18", bcrtcr = "TCR"
))
tcr.data$sample <- stringr::str_split_fixed(tcr.data$TCR_clonotype_id, pattern = "_clonotype", n = 2)[, 1]
tcr.data$bc <- rownames(tcr.data)

tcr.data <- tcr.data %>%
  group_by(sample) %>%
  mutate(TCR.sample = length(sample))
tcr.data <- tcr.data %>%
  group_by(TCR_clonotype_id) %>%
  mutate(TCR.cells = length(TCR_clonotype_id))
tcr.data$freq <- tcr.data$TCR.cells / tcr.data$TCR.sample
tcr.data <- as.data.frame(tcr.data)
rownames(tcr.data) <- tcr.data$bc

CAR.data <- AddMetaData(CAR.data, metadata = tcr.data)

# map to reference
reference <- LoadH5Seurat("./data/pbmc_multimodal.h5seurat")
CAR.data <- SCTransform(CAR.data, verbose = FALSE)
CAR.data <- CellCycleScoring(CAR.data, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

anchors <- FindTransferAnchors(
  reference = reference,
  query = CAR.data,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50
)

CAR.data <- MapQuery(
  anchorset = anchors,
  query = CAR.data,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

# plot predicted cell types
p <- DimPlot(CAR.data,
  cols = "grey90", cells.highlight = list(
    "CD4" = colnames(CAR.data)[which(CAR.data$predicted.celltype.l1 == "CD4 T")],
    "CD8" = colnames(CAR.data)[which(CAR.data$predicted.celltype.l1 == "CD8 T")]
  ),
  cols.highlight = c("CD8" = "red", "CD4" = "blue"), sizes.highlight = 0.5
) +
  NoLegend() +
  NoAxes()
ggsave("./CAR/figures/UMAP/20230210_CAR_CD4_CD8.png", width = 4, height = 4, dpi = 600, plot = p)

p <- DimPlot(CAR.data, group.by = "predicted.celltype.l2", cols = T.cell.subset.colors, na.value = "grey90") +
  NoLegend() +
  NoAxes() +
  theme(plot.title = element_blank())
ggsave("./CAR/figures/UMAP/20230210_predicted_celltype.png", width = 4, height = 4, dpi = 600, plot = p)

# read CAR amplicons
CAR.df <- data.frame()
for (s in libraries) {
  CAR.ONT <- nanoranger.R::extract_fusion_gene(paste0("./data/CAR/amplicons/", s, "_CD28_CARTmod.csv.gz"), WILDTYPE = "CD28", FUSION = "CARTmod", filter = 10)
  CAR.ONT$bc <- paste0(s, "_", CAR.ONT$bc, "-1")
  CAR.ONT <- CAR.ONT[which(CAR.ONT$bc %in% colnames(CAR.data)), ]
  rownames(CAR.ONT) <- CAR.ONT$bc
  CAR.df <- rbind(CAR.ONT, CAR.df)
}
write.csv2(file = "./data/CAR/amplicons/20230209_CAR_processed.csv", CAR.df)
CAR.df <- read.table("./data/CAR/amplicons//20230209_CAR_processed.csv", row.names = 1, sep = ";", header = 1, dec = ",")

CAR.df$timepoint <- CAR.data$timepoint[CAR.df$bc]
CAR.df$predicted.celltype.l1 <- CAR.data$predicted.celltype.l1[CAR.df$bc]
CAR.df$predicted.celltype.l2 <- CAR.data$predicted.celltype.l2[CAR.df$bc]
CAR.df$S <- CAR.data$S.Score[CAR.df$bc]
CAR.df$G2M <- CAR.data$G2M.Score[CAR.df$bc]
CAR.df$B2M <- GetAssayData(CAR.data)["B2M", CAR.df$bc]
CAR.df$counts <- CAR.data$nCount_RNA[CAR.df$bc]
CAR.df$sample <- CAR.data$orig.ident[CAR.df$bc]
CAR.df$UMAP1 <- Embeddings(CAR.data, reduction = "ref.umap")[CAR.df$bc, "refUMAP_1"]
CAR.df$UMAP2 <- Embeddings(CAR.data, reduction = "ref.umap")[CAR.df$bc, "refUMAP_2"]
CAR.df$ratio <- CAR.df$CARTmod / CAR.df$CD28
CAR.df$CloneID <- CAR.data$TCR_clonotype_id[CAR.df$bc]
CAR.df$freq <- CAR.data$freq[CAR.df$bc]

# CAR expression
p <- ggplot() +
  geom_point(data = as.data.frame(Embeddings(CAR.data, reduction = "ref.umap")), aes(x = refUMAP_1, y = refUMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = CAR.df[order(CAR.df$CARTmod), ], aes(x = UMAP1, y = UMAP2, color = log10(CARTmod)), size = 0.5) +
  scale_color_gradientn(colors = BuenColors::jdb_palette(name = "solar_glare"), na.value = "grey90") +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./CAR/figures/UMAP/20230210_CAR_expression.png", width = 4, height = 4, dpi = 600, plot = p)

# CAR expression with legend
p <- ggplot() +
  geom_point(data = as.data.frame(Embeddings(CAR.data, reduction = "ref.umap")), aes(x = refUMAP_1, y = refUMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = CAR.df[order(CAR.df$CARTmod), ], aes(x = UMAP1, y = UMAP2, color = log10(CARTmod)), size = 0.5) +
  scale_color_gradientn(colors = BuenColors::jdb_palette(name = "solar_glare"), na.value = "grey90") +
  theme_classic() +
  NoAxes()
ggsave("./CAR/figures/UMAP/20230210_CAR_expression.svg", width = 4, height = 4, dpi = 600, plot = p)

# T cell expansion
boo <- CAR.df[which(!is.na(CAR.df$freq)), ]
p <- ggplot() +
  geom_point(data = as.data.frame(Embeddings(CAR.data, reduction = "ref.umap")), aes(x = refUMAP_1, y = refUMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = boo[order(boo$freq), ], aes(x = UMAP1, y = UMAP2, color = log10(freq)), size = 0.5) +
  scale_color_gradientn(colors = BuenColors::jdb_palette(name = "solar_glare"), na.value = "grey90") +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./CAR/figures/UMAP/20230210_Tcell_expansion.png", width = 4, height = 4, dpi = 600, plot = p)

# T cell expansion with legend
p <- ggplot() +
  geom_point(data = as.data.frame(Embeddings(CAR.data, reduction = "ref.umap")), aes(x = refUMAP_1, y = refUMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = boo[order(boo$freq), ], aes(x = UMAP1, y = UMAP2, color = log10(freq)), size = 0.5) +
  scale_color_gradientn(colors = BuenColors::jdb_palette(name = "solar_glare"), na.value = "grey90") +
  theme_classic() +
  NoAxes()
ggsave("./CAR/figures/UMAP/20230210_Tcell_expansion.svg", width = 4, height = 4, dpi = 600, plot = p)

CAR.data$freq.log <- log10(CAR.data$freq)
p <- FeaturePlot(CAR.data,
  cells = CAR.df$bc,
  features = "freq.log", order = T
) +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = "solar_glare"), na.value = "grey90") +
  NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave("./CAR/figures/UMAP/20230210_Tcell_expansion.png", width = 4, height = 4, dpi = 600, plot = p)

CAR.df.plot <- CAR.df %>%
  filter(CAR.df$predicted.celltype.l2 %in% c(T.cell.subsets, "NK", "NK Proliferating") & CAR.df$CARTmod > 0) %>%
  group_by(predicted.celltype.l2) %>%
  arrange(desc(-CARTmod)) %>%
  mutate(pos = seq(1, length(predicted.celltype.l2)) / length(predicted.celltype.l2))

stats <- as.data.frame(CAR.df.plot %>%
  group_by(predicted.celltype.l2) %>%
  summarize(mean = mean(CARTmod)) %>%
  arrange(desc(-mean)))
boo <- seq(1, 19, 2)
names(boo) <- stats$predicted.celltype.l2
stats$pos <- seq(1.8, 19.8, 2)

CAR.df.plot$abs.pos <- boo[CAR.df.plot$predicted.celltype.l2] + CAR.df.plot$pos
CAR.df.plot$CARTmod[which(CAR.df.plot$CARTmod == 0)] <- 0.1

# CAR expression by cell type
ggplot(CAR.df.plot, aes(x = abs.pos, y = CARTmod)) +
  geom_point(size = 0.5, aes(color = predicted.celltype.l2)) +
  scale_y_log10("CAR expression") +
  scale_color_manual(values = T.cell.subset.colors) +
  stat_summary(data = stats, aes(x = pos, y = mean), fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 1, size = 0.5, color = "black") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )
ggsave("./CAR/figures/plots/20230210_CART_expression_by_celltype.svg", width = 2.5, height = 2)

View(CAR.df.plot %>% group_by(sample, predicted.celltype.l2) %>% summarize(high = length(which(CARTmod > 10)) / length(predicted.celltype.l2)))

boo <- CAR.df %>%
  filter(CARTmod != 0) %>%
  group_by(sample, timepoint, predicted.celltype.l1) %>%
  mutate(CART.zscore = (CARTmod - mean(CARTmod) / sd(CARTmod)))

# CAR expression by CD4 or CD8 T cell
p <- ggplot(boo[which(boo$predicted.celltype.l1 %in% c("CD4 T", "CD8 T")), ], aes(x = predicted.celltype.l1, y = CART.zscore, color = predicted.celltype.l1)) +
  ggrastr::rasterize(geom_jitter(size = 0.5), dpi = 600) +
  geom_violin(color = "black") +
  scale_y_continuous("Expression CAR", limits = c(-2, 120)) +
  scale_x_discrete(labels = c("CD4+", "CD8+")) +
  scale_color_manual(values = c("CD4 T" = "blue", "CD8 T" = "red")) +
  geom_signif(comparisons = list(c("CD4 T", "CD8 T")), y_position = 100, color = "black", test = "t.test") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank()
  )
ggsave("./CAR/figures/plots/20230210_CART_expression_CD4_CD8.svg", width = 1.5, height = 2, plot = p)

boo <- boo %>%
  summarize(CART.mean = mean(CART.zscore))

# CAR expression by CD4 or CD8 T cell
ggplot(boo[which(boo$predicted.celltype.l1 %in% c("CD4 T", "CD8 T")), ], aes(x = predicted.celltype.l1, y = CART.mean)) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, width = 0.5, size = 0.5, color = "orange") +
  geom_line(aes(group = sample)) +
  geom_point(aes(color = timepoint), size = 0.5) +
  geom_signif(comparisons = list(c("CD4 T", "CD8 T")), test = "t.test", test.args = c("paired" = TRUE, alternative = "greater")) +
  scale_y_continuous("Mean expression CAR") +
  scale_x_discrete(labels = c("CD4+", "CD8+")) +
  scale_color_manual(values = c("IP" = "firebrick", "D7" = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank()
  )
ggsave("./CAR/figures/plots/20230210_CART_expression_mean_CD4_CD8.svg", width = 1.5, height = 2)

# CAR versus wildtype CD28 expression
p <- ggplot() +
  geom_abline(slope = 1) +
  ggrastr::rasterize(geom_point(data = CAR.df[which(CAR.df$timepoint == "IP"), ], aes(x = CD28, y = CARTmod), color = "firebrick", size = 0.5), dpi = 600) +
  ggrastr::rasterize(geom_point(data = CAR.df[which(CAR.df$timepoint == "D7"), ], aes(x = CD28, y = CARTmod), color = "grey", size = 0.5), dpi = 600) +
  # geom_point(aes(color=timepoint)) +
  scale_x_sqrt("CD28 expression", limits = c(0, 300), breaks = c(1, 10, 50, 100, 200, 300)) +
  scale_y_sqrt("CAR expression", limits = c(0, 300), breaks = c(1, 10, 50, 100, 200, 300)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./CAR/figures/plots/20230210_CART_CD28_expression.svg", width = 2, height = 2, plot = p)

# CAR expression versus native TCR expansion
p <- ggplot(data = CAR.df[which(!is.na(CAR.df$freq)), ], aes(x = freq, y = CARTmod, color = log10(freq))) +
  ggrastr::rasterize(geom_point(size = 0.5), dpi = 600) +
  scale_x_log10("%frequency native TCR", breaks = c(10^-4, 10^-3, 10^-2, 10^-1), labels = c(0.01, 0.1, 1, 10)) +
  scale_y_continuous("CAR expression") +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = "solar_glare")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./CAR/figures/plots/20230210_CART_expression_native_TCR_expansion.svg", width = 2, height = 2, plot = p)

p <- ggplot(
  data = CAR.df[which(!is.na(CAR.df$freq) & CAR.df$predicted.celltype.l2 %in% names(T.cell.subset.colors)), ],
  aes(x = freq, y = CARTmod, color = predicted.celltype.l2)
) +
  ggrastr::rasterize(geom_point(size = 0.5), dpi = 600) +
  scale_x_log10("%frequency native TCR", breaks = c(10^-4, 10^-3, 10^-2, 10^-1), labels = c(0.01, 0.1, 1, 10)) +
  scale_y_continuous("CAR expression") +
  scale_color_manual(values = T.cell.subset.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./CAR/figures/plots/20230210_CART_expression_native_TCR_expansion_celltype.svg", width = 2, height = 2, plot = p)

# CAR expression versus native TCR expansion in CD4+ T cells
p <- ggplot(data = CAR.df[which(!is.na(CAR.df$freq) & CAR.df$predicted.celltype.l1 == "CD4 T"), ], aes(x = freq, y = CARTmod)) +
  ggrastr::rasterize(geom_point(size = 0.5, aes(color = log10(freq))), dpi = 600) +
  scale_x_log10("%frequency native TCR", breaks = c(10^-4, 10^-3, 10^-2, 10^-1), labels = c(0.01, 0.1, 1, 10)) +
  scale_y_continuous("CAR expression") +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = "solar_glare")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./CAR/figures/plots/20230210_CART_expression_native_TCR_expansion_CD4.svg", width = 2, height = 2, plot = p)

# CAR expression versus native TCR expansion in CD8+ T cells
p <- ggplot(data = CAR.df[which(!is.na(CAR.df$freq) & CAR.df$predicted.celltype.l1 == "CD8 T"), ], aes(x = freq, y = CARTmod, color = log10(freq))) +
  ggrastr::rasterize(geom_point(size = 0.5), dpi = 600) +
  scale_x_log10("%frequency native TCR", breaks = c(10^-4, 10^-3, 10^-2, 10^-1), labels = c(0.01, 0.1, 1, 10)) +
  scale_y_continuous("CAR expression") +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = "solar_glare")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./CAR/figures/plots/20230210_CART_expression_native_TCR_expansion_CD8.svg", width = 2, height = 2, plot = p)

saveRDS(file = "./data/CAR/objects/20230210_CAR_data.rds", CAR.data)
CAR.data <- readRDS(file = "./data/CAR/objects/20230210_CAR_data.rds")

CAR.data <- AddModuleScore(CAR.data,
  features = list("Oliveira.tumor.score" = c("PDCD1", "CTLA4", "TIGIT", "HAVCR2", "TOX", "LAG3", "ENTPD1")),
  ctrl = 5, name = "Oliveira.tumor.score"
)
CAR.data <- AddModuleScore(CAR.data,
  features = list("Oliveira.memory.score" = c("TCF7", "IL7R", "SELL", "CCR7", "CD28")),
  ctrl = 5, name = "Oliveira.memory.score"
)

CAR.df$Oliveira.tumor.score <- CAR.data$Oliveira.tumor.score1[CAR.df$bc]
CAR.df$Oliveira.memory.score <- CAR.data$Oliveira.memory.score1[CAR.df$bc]

CAR.df <- CAR.df[CAR.df$CARTmod > 0, ]

CAR.df$exhausted <- ifelse(CAR.df$Oliveira.memory.score < 0.5 & CAR.df$Oliveira.tumor.score > 0.5, "exhausted", "non.exhausted")
CAR.df$high.expression <- ifelse(CAR.df$CARTmod > 20, "high", "low")

boo <- CAR.df %>%
  group_by(exhausted, high.expression) %>%
  tally() %>%
  tidyr::pivot_wider(names_from = c("exhausted", "high.expression"), values_from = "n")

# CAR T expression versus T cell exhaustion
p <- ggplot(data = CAR.df[order(CAR.df$CARTmod), ], aes(x = Oliveira.tumor.score, y = CARTmod, color = log10(CARTmod))) +
  ggrastr::rasterize(geom_point(size = 0.5), dpi = 600) +
  geom_hline(yintercept = 20) +
  geom_vline(xintercept = 0.5) +
  scale_x_continuous("exhaustion score", limits = c(-0.5, 1.5)) +
  scale_y_continuous("CAR expression", limits = c(0, 120), breaks = c(0, 30, 60, 90, 120)) +
  scale_color_gradientn(colours = BuenColors::jdb_palette(name = "solar_glare")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./CAR/figures/plots/20230223_CART_expression_Tcell_exhaustion.svg", width = 2, height = 2, plot = p)
