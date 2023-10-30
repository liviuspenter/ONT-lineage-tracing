# compare IL7R isoform expression across libraries

library(dplyr)
library(nanoranger.R)
library(ggplot2)
library(Seurat)

# read libraries
CAR.so <- readRDS("./data/CAR/objects/20220819_CAR_so.rds")
AML3005 <- readRDS("./data/AML/objects/20220718_AML3005_13.rds")
AML1007 <- readRDS("./data/AML/objects/20220825_AML1007.rds")
AML1010 <- readRDS("./data/AML/objects/20220711_AML1010_1.rds")
AML1026 <- readRDS("./data/AML/objects/20220711_AML1026_2.rds")
TIL3 <- readRDS("./data/TCR/objects/20220820_TCR3.rds")
TIL4 <- readRDS("./data/TCR/objects/20230104_TIL4.rds")
PBMC7 <- readRDS("./data/TCR/objects/20230104_PBMC7.rds")
PBMC8 <- readRDS("./data/TCR/objects/20220820_TCR8.rds")

IL7R.97_6 <- extract_isoforms("./data/isoforms/97_6/97_6_exon_list.gz", GENE = "IL7R", EXON = "exon6", filter = 10)
IL7R.97_6$bc <- paste0(IL7R.97_6$bc, "-1")
IL7R.97_6 <- IL7R.97_6[which(IL7R.97_6$bc %in% colnames(CAR.so)), ]
IL7R.97_6$ratio <- IL7R.97_6$detected / (IL7R.97_6$detected + IL7R.97_6$not.detected)
IL7R.97_6$UMAP1 <- CAR.so@reductions$umap@cell.embeddings[IL7R.97_6$bc, "UMAP_1"]
IL7R.97_6$UMAP2 <- CAR.so@reductions$umap@cell.embeddings[IL7R.97_6$bc, "UMAP_2"]
IL7R.97_6$sample <- "97_6"

IL7R.1007_1 <- extract_isoforms("./data/isoforms/AML1007.1/AML1007_1_exons.gz", GENE = "IL7R", EXON = "exon6", filter = 10)
IL7R.1007_1$bc <- paste0("AML1007.1_", IL7R.1007_1$bc, "-1")
IL7R.1007_1 <- IL7R.1007_1[which(IL7R.1007_1$bc %in% colnames(AML1007)), ]
IL7R.1007_1$ratio <- IL7R.1007_1$detected / (IL7R.1007_1$detected + IL7R.1007_1$not.detected)
IL7R.1007_1$UMAP1 <- AML1007@reductions$umap@cell.embeddings[IL7R.1007_1$bc, "UMAP_1"]
IL7R.1007_1$UMAP2 <- AML1007@reductions$umap@cell.embeddings[IL7R.1007_1$bc, "UMAP_2"]
IL7R.1007_1$sample <- "AML1007_1"

IL7R.1007_3 <- extract_isoforms("./data/isoforms/AML1007.3/AML1007_3_exons.gz", GENE = "IL7R", EXON = "exon6", filter = 10)
IL7R.1007_3$bc <- paste0("AML1007.3_", IL7R.1007_3$bc, "-1")
IL7R.1007_3 <- IL7R.1007_3[which(IL7R.1007_3$bc %in% colnames(AML1007)), ]
IL7R.1007_3$ratio <- IL7R.1007_3$detected / (IL7R.1007_3$detected + IL7R.1007_3$not.detected)
IL7R.1007_3$UMAP1 <- AML1007@reductions$umap@cell.embeddings[IL7R.1007_3$bc, "UMAP_1"]
IL7R.1007_3$UMAP2 <- AML1007@reductions$umap@cell.embeddings[IL7R.1007_3$bc, "UMAP_2"]
IL7R.1007_3$sample <- "AML1007_3"

IL7R.1007_5 <- extract_isoforms("./data/isoforms/AML1007.5/AML1007_5_exons.gz", GENE = "IL7R", EXON = "exon6", filter = 10)
IL7R.1007_5$bc <- paste0("AML1007.5_", IL7R.1007_5$bc, "-1")
IL7R.1007_5 <- IL7R.1007_5[which(IL7R.1007_5$bc %in% colnames(AML1007)), ]
IL7R.1007_5$ratio <- IL7R.1007_5$detected / (IL7R.1007_5$detected + IL7R.1007_5$not.detected)
IL7R.1007_5$UMAP1 <- AML1007@reductions$umap@cell.embeddings[IL7R.1007_5$bc, "UMAP_1"]
IL7R.1007_5$UMAP2 <- AML1007@reductions$umap@cell.embeddings[IL7R.1007_5$bc, "UMAP_2"]
IL7R.1007_5$sample <- "AML1007_5"

IL7R.1010_1 <- extract_isoforms("./data/isoforms/AML1010.1/AML1010_1_exons.gz", GENE = "IL7R", EXON = "exon6", filter = 10)
IL7R.1010_1$bc <- paste0("AML1010.1_", IL7R.1010_1$bc, "-1")
IL7R.1010_1 <- IL7R.1010_1[which(IL7R.1010_1$bc %in% colnames(AML1010)), ]
IL7R.1010_1$ratio <- IL7R.1010_1$detected / (IL7R.1010_1$detected + IL7R.1010_1$not.detected)
IL7R.1010_1$UMAP1 <- AML1010@reductions$ref.umap@cell.embeddings[IL7R.1010_1$bc, "refUMAP_1"]
IL7R.1010_1$UMAP2 <- AML1010@reductions$ref.umap@cell.embeddings[IL7R.1010_1$bc, "refUMAP_2"]
IL7R.1010_1$sample <- "AML1010_1"

IL7R.1026_2 <- extract_isoforms("./data/isoforms/AML1026.2/AML1026_2_exons.gz", GENE = "IL7R", EXON = "exon6", filter = 10)
IL7R.1026_2$bc <- paste0("AML1026.2_", IL7R.1026_2$bc, "-1")
IL7R.1026_2 <- IL7R.1026_2[which(IL7R.1026_2$bc %in% colnames(AML1026)), ]
IL7R.1026_2$ratio <- IL7R.1026_2$detected / (IL7R.1026_2$detected + IL7R.1026_2$not.detected)
IL7R.1026_2$UMAP1 <- AML1026@reductions$umap@cell.embeddings[IL7R.1026_2$bc, "UMAP_1"]
IL7R.1026_2$UMAP2 <- AML1026@reductions$umap@cell.embeddings[IL7R.1026_2$bc, "UMAP_2"]
IL7R.1026_2$sample <- "AML1026_2"

IL7R.3005_1 <- extract_isoforms("./data/isoforms/AML3005.1/AML3005_1_exons.gz", GENE = "IL7R", EXON = "exon6", filter = 10)
IL7R.3005_1$bc <- paste0("AML3005.1_", IL7R.3005_1$bc, "-1")
IL7R.3005_1 <- IL7R.3005_1[which(IL7R.3005_1$bc %in% colnames(AML3005)), ]
IL7R.3005_1$ratio <- IL7R.3005_1$detected / (IL7R.3005_1$detected + IL7R.3005_1$not.detected)
IL7R.3005_1$UMAP1 <- AML3005@reductions$umap@cell.embeddings[IL7R.3005_1$bc, "UMAP_1"]
IL7R.3005_1$UMAP2 <- AML3005@reductions$umap@cell.embeddings[IL7R.3005_1$bc, "UMAP_2"]
IL7R.3005_1$sample <- "AML3005_1"

IL7R.3005_3 <- extract_isoforms("./data/isoforms/AML3005.3/AML3005_3_exons.gz", GENE = "IL7R", EXON = "exon6", filter = 10)
IL7R.3005_3$bc <- paste0("AML3005.3_", IL7R.3005_3$bc, "-1")
IL7R.3005_3 <- IL7R.3005_3[which(IL7R.3005_3$bc %in% colnames(AML3005)), ]
IL7R.3005_3$ratio <- IL7R.3005_3$detected / (IL7R.3005_3$detected + IL7R.3005_3$not.detected)
IL7R.3005_3$UMAP1 <- AML3005@reductions$umap@cell.embeddings[IL7R.3005_3$bc, "UMAP_1"]
IL7R.3005_3$UMAP2 <- AML3005@reductions$umap@cell.embeddings[IL7R.3005_3$bc, "UMAP_2"]
IL7R.3005_3$sample <- "AML3005_3"

IL7R.TIL3 <- extract_isoforms("./data/isoforms/TIL3/TIL3_exons.gz", GENE = "IL7R", EXON = "exon6", filter = 10)
IL7R.TIL3 <- IL7R.TIL3[which(IL7R.TIL3$bc %in% colnames(TIL3)), ]
IL7R.TIL3$ratio <- IL7R.TIL3$detected / (IL7R.TIL3$detected + IL7R.TIL3$not.detected)
IL7R.TIL3$UMAP1 <- TIL3@reductions$umap@cell.embeddings[IL7R.TIL3$bc, "UMAP_1"]
IL7R.TIL3$UMAP2 <- TIL3@reductions$umap@cell.embeddings[IL7R.TIL3$bc, "UMAP_2"]
IL7R.TIL3$sample <- "TIL3"

IL7R.TIL4 <- extract_isoforms("./data/isoforms/TIL4/TIL4_exons.gz", GENE = "IL7R", EXON = "exon6", filter = 10)
IL7R.TIL4 <- IL7R.TIL4[which(IL7R.TIL4$bc %in% colnames(TIL4)), ]
IL7R.TIL4$ratio <- IL7R.TIL4$detected / (IL7R.TIL4$detected + IL7R.TIL4$not.detected)
IL7R.TIL4$UMAP1 <- TIL4@reductions$umap@cell.embeddings[IL7R.TIL4$bc, "UMAP_1"]
IL7R.TIL4$UMAP2 <- TIL4@reductions$umap@cell.embeddings[IL7R.TIL4$bc, "UMAP_2"]
IL7R.TIL4$sample <- "TIL4"

IL7R.PBMC7 <- extract_isoforms("./data/isoforms/PBMC7/PBMC7_exons.gz", GENE = "IL7R", EXON = "exon6", filter = 10)
IL7R.PBMC7 <- IL7R.PBMC7[which(IL7R.PBMC7$bc %in% colnames(PBMC7)), ]
IL7R.PBMC7$ratio <- IL7R.PBMC7$detected / (IL7R.PBMC7$detected + IL7R.PBMC7$not.detected)
IL7R.PBMC7$UMAP1 <- PBMC7@reductions$umap@cell.embeddings[IL7R.PBMC7$bc, "UMAP_1"]
IL7R.PBMC7$UMAP2 <- PBMC7@reductions$umap@cell.embeddings[IL7R.PBMC7$bc, "UMAP_2"]
IL7R.PBMC7$sample <- "PBMC7"

IL7R.PBMC8 <- extract_isoforms("./data/isoforms/PBMC8/PBMC8_exons.gz", GENE = "IL7R", EXON = "exon6", filter = 10)
IL7R.PBMC8 <- IL7R.PBMC8[which(IL7R.PBMC8$bc %in% colnames(PBMC8)), ]
IL7R.PBMC8$ratio <- IL7R.PBMC8$detected / (IL7R.PBMC8$detected + IL7R.PBMC8$not.detected)
IL7R.PBMC8$UMAP1 <- PBMC8@reductions$umap@cell.embeddings[IL7R.PBMC8$bc, "UMAP_1"]
IL7R.PBMC8$UMAP2 <- PBMC8@reductions$umap@cell.embeddings[IL7R.PBMC8$bc, "UMAP_2"]
IL7R.PBMC8$sample <- "PBMC8"

IL7R <- dplyr::bind_rows(
  IL7R.1007_1[, c("sample", "bc", "detected", "not.detected", "ratio")],
  IL7R.1007_3[, c("sample", "bc", "detected", "not.detected", "ratio")],
  IL7R.1007_5[, c("sample", "bc", "detected", "not.detected", "ratio")],
  IL7R.3005_1[, c("sample", "bc", "detected", "not.detected", "ratio")],
  IL7R.3005_3[, c("sample", "bc", "detected", "not.detected", "ratio")],
  # IL7R.1010_1[,c('sample', 'bc', 'detected', 'not.detected', 'ratio')],
  IL7R.1026_2[, c("sample", "bc", "detected", "not.detected", "ratio")],
  IL7R.TIL3[, c("sample", "bc", "detected", "not.detected", "ratio")],
  IL7R.TIL4[, c("sample", "bc", "detected", "not.detected", "ratio")],
  IL7R.PBMC7[, c("sample", "bc", "detected", "not.detected", "ratio")],
  IL7R.PBMC8[, c("sample", "bc", "detected", "not.detected", "ratio")],
  IL7R.97_6[, c("sample", "bc", "detected", "not.detected", "ratio")]
)

IL7R$tissue <- "AML"
IL7R$tissue[which(IL7R$sample %in% c("TIL3", "TIL4"))] <- "TIL"
IL7R$tissue[which(IL7R$sample %in% c("PBMC7", "PBMC8"))] <- "TIL.PB"
IL7R$tissue[which(IL7R$sample %in% c("97_6"))] <- "CAR"

ggplot(IL7R, aes(x = reorder(sample, ratio), y = 100 * ratio, color = tissue)) +
  geom_jitter(size = 0.5) +
  stat_summary(geom = "crossbar", fun = mean, fun.min = mean, fun.max = mean, color = "black") +
  scale_x_discrete() +
  scale_y_continuous("% IL7R exon6") +
  scale_color_manual(values = c(
    "AML" = "darkgreen",
    "TIL.PB" = "darkblue",
    "TIL" = "orange",
    "CAR" = "firebrick"
  )) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./isoforms/figures/combined/plots/20230104_IL7R_exon6_overview.svg", width = 2.5, height = 2.5)

p <- ggplot() +
  geom_point(data = as.data.frame(CAR.so@reductions$umap@cell.embeddings), aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = IL7R.97_6[order(-IL7R.97_6$ratio), ], aes(x = UMAP1, y = UMAP2, color = ratio), size = 0.5) +
  scale_color_gradient(low = "black", high = "lightgreen") +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./isoforms/figures/CAR/UMAP/20230104_97_6_IL7R.png", width = 4, height = 4, dpi = 600, plot = p)

p <- ggplot() +
  geom_point(data = as.data.frame(TIL3@reductions$umap@cell.embeddings), aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = IL7R.TIL3[order(-IL7R.TIL3$ratio), ], aes(x = UMAP1, y = UMAP2, color = ratio), size = 0.5) +
  scale_color_gradient(low = "black", high = "lightgreen") +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./isoforms/figures/melanoma/UMAP/20230104_TIL3_IL7R.png", width = 4, height = 4, dpi = 600, plot = p)

p <- ggplot() +
  geom_point(data = as.data.frame(PBMC8@reductions$umap@cell.embeddings), aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = IL7R.PBMC8[order(-IL7R.PBMC8$ratio), ], aes(x = UMAP1, y = UMAP2, color = ratio), size = 0.5) +
  scale_color_gradient(low = "black", high = "lightgreen") +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./isoforms/figures/melanoma/UMAP/20230104_PBMC8_IL7R.png", width = 4, height = 4, dpi = 600, plot = p)

p <- ggplot() +
  geom_point(data = as.data.frame(TIL4@reductions$umap@cell.embeddings), aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = IL7R.TIL4[order(-IL7R.TIL4$ratio), ], aes(x = UMAP1, y = UMAP2, color = ratio), size = 0.5) +
  scale_color_gradient(low = "black", high = "lightgreen") +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./isoforms/figures/melanoma/UMAP/20230104_TIL4_IL7R.png", width = 4, height = 4, dpi = 600, plot = p)

p <- ggplot() +
  geom_point(data = as.data.frame(PBMC7@reductions$umap@cell.embeddings), aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = IL7R.PBMC7[order(-IL7R.PBMC7$ratio), ], aes(x = UMAP1, y = UMAP2, color = ratio), size = 0.5) +
  scale_color_gradient(low = "black", high = "lightgreen") +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./isoforms/figures/melanoma/UMAP/20230104_PBMC7_IL7R.png", width = 4, height = 4, dpi = 600, plot = p)

p <- ggplot() +
  geom_point(data = as.data.frame(AML1007@reductions$umap@cell.embeddings), aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = IL7R.1007_1[order(-IL7R.1007_1$ratio), ], aes(x = UMAP1, y = UMAP2, color = ratio), size = 0.5) +
  scale_color_gradient(low = "black", high = "lightgreen") +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./isoforms/figures/AML1007/UMAP/20230104_AML1007_1_IL7R.png", width = 4, height = 4, dpi = 600, plot = p)

p <- ggplot() +
  geom_point(data = as.data.frame(AML1007@reductions$umap@cell.embeddings), aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = IL7R.1007_3[order(-IL7R.1007_3$ratio), ], aes(x = UMAP1, y = UMAP2, color = ratio), size = 0.5) +
  scale_color_gradient(low = "black", high = "lightgreen") +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./isoforms/figures/AML1007/UMAP/20230104_AML1007_3_IL7R.png", width = 4, height = 4, dpi = 600, plot = p)

p <- ggplot() +
  geom_point(data = as.data.frame(AML1007@reductions$umap@cell.embeddings), aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = IL7R.1007_5[order(-IL7R.1007_5$ratio), ], aes(x = UMAP1, y = UMAP2, color = ratio), size = 0.5) +
  scale_color_gradient(low = "black", high = "lightgreen") +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./isoforms/figures/AML1007/UMAP/20230104_AML1007_5_IL7R.png", width = 4, height = 4, dpi = 600, plot = p)

p <- ggplot() +
  geom_point(data = as.data.frame(AML3005@reductions$umap@cell.embeddings), aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = IL7R.3005_1[order(-IL7R.3005_1$ratio), ], aes(x = UMAP1, y = UMAP2, color = ratio), size = 0.5) +
  scale_color_gradient(low = "black", high = "lightgreen") +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./isoforms/figures/AML3005/UMAP/20230104_AML3005_1_IL7R.png", width = 4, height = 4, dpi = 600, plot = p)

q <- ggplot() +
  geom_point(data = as.data.frame(AML3005@reductions$umap@cell.embeddings), aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = IL7R.3005_3[order(-IL7R.3005_3$ratio), ], aes(x = UMAP1, y = UMAP2, color = ratio), size = 0.5) +
  scale_color_gradient(low = "black", high = "lightgreen") +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./isoforms/figures/AML3005/UMAP/20230104_AML3005_3_IL7R.png", width = 4, height = 4, dpi = 600, plot = q)

p <- FeaturePlot(AML3005, features = "IL7R", order = T, reduction = "umap") + NoAxes() + NoLegend() +
  theme(plot.title = element_blank(), panel.spacing = unit(0, "pt")) + scale_color_gradientn(colours = c("grey90", BuenColors::jdb_palette(name = "brewer_green")))
ggsave("./isoforms/figures/AML3005/UMAP/20230104_AML3005_1_IL7R_original.png", width = 4, height = 4, dpi = 600, plot = p)

p <- FeaturePlot(AML3005, features = "IL7R", order = T, reduction = "umap") +
  theme(plot.title = element_blank(), panel.spacing = unit(0, "pt")) + scale_color_gradientn(colours = c("grey90", BuenColors::jdb_palette(name = "brewer_green")))
ggsave("./isoforms/figures/AML3005/UMAP/20230104_AML3005_1_IL7R_original.svg", width = 4, height = 4, dpi = 600, plot = p)
