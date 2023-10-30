# analysis of cells with detectable wildtype CD28 or CAR sequence

library(dplyr)
library(ggplot2)
library(nanoranger.R)
library(Seurat)

technology.colors <- c("Both" = "#283350", "ONT" = "#f93800", "Illumina" = "#ffb500")

so <- readRDS(file = "./data/CAR/objects/20220819_CAR_so.rds")

### look into CAR data from ONT
CAR.ONT <- extract_fusion_gene("./data/CAR/97_6_CART_CD28.csv.gz", WILDTYPE = "CD28", FUSION = "CARTmod")
CAR.ONT$bc <- paste0(CAR.ONT$bc, "-1")
CAR.ONT <- CAR.ONT[which(CAR.ONT$bc %in% colnames(so)), ]
rownames(CAR.ONT) <- CAR.ONT$bc

so.subset <- subset(so, cells = which(colnames(so) %in% CAR.ONT$bc))
so.subset <- subset(so.subset, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20)
so.subset <- NormalizeData(so.subset)
so.subset <- FindVariableFeatures(so.subset)
so.subset <- ScaleData(so.subset)
so.subset <- RunPCA(so.subset)
so.subset <- FindNeighbors(so.subset)
so.subset <- FindClusters(so.subset, resolution = 0.3)
so.subset <- RunUMAP(so.subset, dims = 1:20)
so.subset$mutated <- CAR.ONT[colnames(so.subset), "mutated"]
saveRDS(file = "./data/CAR/objects/20220819_CAR_so_subset.rds", so.subset)
so.subset <- readRDS(file = "./data/CAR/objects/20220819_CAR_so_subset.rds")

p <- DimPlot(so.subset, group.by = "mutated", cols = c("mutated" = "red", "wildtype" = "#00aeef")) +
  NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave("./CAR/figures/UMAP/20220819_CAR_CD28.png", width = 3, height = 3, dpi = 600, plot = p)

p <- FeaturePlot(so.subset, features = "CD4", min.cutoff = "q05", max.cutoff = "q95") + NoLegend() + NoAxes() + theme(plot.title = element_blank()) +
  scale_color_gradient(low = "grey", high = "purple")
ggsave("./CAR/figures/UMAP/20220819_CD4.png", width = 3, height = 3, dpi = 600, plot = p)

p <- FeaturePlot(so.subset, features = "CD8A", min.cutoff = "q05", max.cutoff = "q95") + NoLegend() + NoAxes() + theme(plot.title = element_blank()) +
  scale_color_gradient(low = "grey", high = "purple")
ggsave("./CAR/figures/UMAP/20220819_CD8A.png", width = 3, height = 3, dpi = 600, plot = p)

p <- FeaturePlot(so.subset, features = "MKI67", min.cutoff = "q05", max.cutoff = "q95") + NoLegend() + NoAxes() + theme(plot.title = element_blank()) +
  scale_color_gradient(low = "grey", high = "purple")
ggsave("./CAR/figures/UMAP/20220819_MKI67.png", width = 3, height = 3, dpi = 600, plot = p)

p <- FeaturePlot(so.subset, features = "IL7R", min.cutoff = "q05", max.cutoff = "q95") + NoLegend() + NoAxes() + theme(plot.title = element_blank()) +
  scale_color_gradient(low = "grey", high = "purple")
ggsave("./CAR/figures/UMAP/20220819_IL7R.png", width = 3, height = 3, dpi = 600, plot = p)

markers <- FindMarkers(so.subset, group.by = "mutated", ident.1 = "mutated", ident.2 = "wildtype")
markers$gene <- rownames(markers)
markers$highlight <- ifelse(-log10(markers$p_val_adj) > 4 & abs(markers$avg_log2FC) > 0.5, "yes", "no")

ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = highlight)) +
  geom_point() +
  ggrepel::geom_label_repel(
    data = markers[which(markers$highlight == "yes"), ], aes(x = avg_log2FC, y = -log10(p_val_adj), color = highlight, label = gene),
    size = 3, label.size = 0
  ) +
  geom_hline(yintercept = 4) +
  geom_vline(xintercept = c(0.5, -0.5)) +
  scale_x_continuous("Log2FC", limits = c(-1.2, 1.2)) +
  scale_y_continuous("-log10(FDR)") +
  scale_color_manual(values = c("yes" = "black", "no" = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./CAR/figures/plots/20220819_volcano_CAR_CD28.svg", width = 2, height = 2)

so.subset$manual.cluster <- "CD4"
so.subset$manual.cluster[which(so.subset$seurat_clusters %in% c(2, 4))] <- "CD8"

VlnPlot(so.subset, group.by = "manual.cluster", split.by = "mutated", features = "IL7R", cols = c("mutated" = "red", "wildtype" = "green")) + NoLegend() +
  scale_x_discrete(labels = c("mutated" = "CAR", "wildtype" = "no CAR")) +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.title.x = element_blank(),
    plot.title = element_text("Arial", size = 10, color = "black", face = "bold.italic")
  )
ggsave("./CAR/figures/plots/20220819_CAR_IL7R.svg", width = 2, height = 2)

### clusters
pbmc.markers <- FindAllMarkers(so.subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

p <- DimPlot(so.subset, group.by = "manual.cluster", cols = c("CD4" = "grey", "CD8" = "black")) +
  NoLegend() + NoAxes() + theme(plot.title = element_blank())
ggsave("./CAR/figures/UMAP/20220819_CAR_CD4_CD8.png", width = 3, height = 3, dpi = 600, plot = p)

Idents(so.subset) <- "mutated"
ggplot(as.data.frame(prop.table(table(Idents(so.subset), so.subset$manual.cluster), margin = 2)), aes(x = Var2, y = 100 * Freq, fill = Var1)) +
  geom_col() +
  scale_fill_manual(values = c("mutated" = "red", "wildtype" = "green")) +
  scale_x_discrete(labels = c("mutated" = "CAR", "wildtype" = "no CAR")) +
  scale_y_continuous("% cells") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank()
  )
ggsave("./CAR/figures/plots/20220819_CAR_CD4_CD8.svg", width = 1.3, height = 2)
