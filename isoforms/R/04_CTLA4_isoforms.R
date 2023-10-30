# compare CTLA-4 isoform expression across libraries

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

CTLA4.97_6 <- extract_isoforms("./data/isoforms/97_6/97_6_exon_list.gz", GENE = "CTLA4", EXON = "exon2", filter = 10)
CTLA4.97_6$bc <- paste0(CTLA4.97_6$bc, "-1")
CTLA4.97_6 <- CTLA4.97_6[which(CTLA4.97_6$bc %in% colnames(CAR.so)), ]
CTLA4.97_6$ratio <- CTLA4.97_6$detected / (CTLA4.97_6$detected + CTLA4.97_6$not.detected)
CTLA4.97_6$UMAP1 <- CAR.so@reductions$umap@cell.embeddings[CTLA4.97_6$bc, "UMAP_1"]
CTLA4.97_6$UMAP2 <- CAR.so@reductions$umap@cell.embeddings[CTLA4.97_6$bc, "UMAP_2"]
CTLA4.97_6$sample <- "97_6"

CTLA4.1007_1 <- extract_isoforms("./data/isoforms/AML1007.1/AML1007_1_exons.gz", GENE = "CTLA4", EXON = "exon2", filter = 10)
CTLA4.1007_1$bc <- paste0("AML1007.1_", CTLA4.1007_1$bc, "-1")
CTLA4.1007_1 <- CTLA4.1007_1[which(CTLA4.1007_1$bc %in% colnames(AML1007)), ]
CTLA4.1007_1$ratio <- CTLA4.1007_1$detected / (CTLA4.1007_1$detected + CTLA4.1007_1$not.detected)
CTLA4.1007_1$UMAP1 <- AML1007@reductions$umap@cell.embeddings[CTLA4.1007_1$bc, "UMAP_1"]
CTLA4.1007_1$UMAP2 <- AML1007@reductions$umap@cell.embeddings[CTLA4.1007_1$bc, "UMAP_2"]
CTLA4.1007_1$sample <- "AML1007_1"

CTLA4.1007_3 <- extract_isoforms("./data/isoforms/AML1007.3/AML1007_3_exons.gz", GENE = "CTLA4", EXON = "exon2", filter = 10)
CTLA4.1007_3$bc <- paste0("AML1007.3_", CTLA4.1007_3$bc, "-1")
CTLA4.1007_3 <- CTLA4.1007_3[which(CTLA4.1007_3$bc %in% colnames(AML1007)), ]
CTLA4.1007_3$ratio <- CTLA4.1007_3$detected / (CTLA4.1007_3$detected + CTLA4.1007_3$not.detected)
CTLA4.1007_3$UMAP1 <- AML1007@reductions$umap@cell.embeddings[CTLA4.1007_3$bc, "UMAP_1"]
CTLA4.1007_3$UMAP2 <- AML1007@reductions$umap@cell.embeddings[CTLA4.1007_3$bc, "UMAP_2"]
CTLA4.1007_3$sample <- "AML1007_3"

CTLA4.1007_5 <- extract_isoforms("./data/isoforms/AML1007.5/AML1007_5_exons.gz", GENE = "CTLA4", EXON = "exon2", filter = 10)
CTLA4.1007_5$bc <- paste0("AML1007.5_", CTLA4.1007_5$bc, "-1")
CTLA4.1007_5 <- CTLA4.1007_5[which(CTLA4.1007_5$bc %in% colnames(AML1007)), ]
CTLA4.1007_5$ratio <- CTLA4.1007_5$detected / (CTLA4.1007_5$detected + CTLA4.1007_5$not.detected)
CTLA4.1007_5$UMAP1 <- AML1007@reductions$umap@cell.embeddings[CTLA4.1007_5$bc, "UMAP_1"]
CTLA4.1007_5$UMAP2 <- AML1007@reductions$umap@cell.embeddings[CTLA4.1007_5$bc, "UMAP_2"]
CTLA4.1007_5$sample <- "AML1007_5"

CTLA4.1010_1 <- extract_isoforms("./data/isoforms/AML1010.1/AML1010_1_exons.gz", GENE = "CTLA4", EXON = "exon2", filter = 10)
CTLA4.1010_1$bc <- paste0("AML1010.1_", CTLA4.1010_1$bc, "-1")
CTLA4.1010_1 <- CTLA4.1010_1[which(CTLA4.1010_1$bc %in% colnames(AML1010)), ]
CTLA4.1010_1$ratio <- CTLA4.1010_1$detected / (CTLA4.1010_1$detected + CTLA4.1010_1$not.detected)
CTLA4.1010_1$UMAP1 <- AML1010@reductions$ref.umap@cell.embeddings[CTLA4.1010_1$bc, "refUMAP_1"]
CTLA4.1010_1$UMAP2 <- AML1010@reductions$ref.umap@cell.embeddings[CTLA4.1010_1$bc, "refUMAP_2"]
CTLA4.1010_1$sample <- "AML1010_1"

CTLA4.1026_2 <- extract_isoforms("./data/isoforms/AML1026.2/AML1026_2_exons.gz", GENE = "CTLA4", EXON = "exon2", filter = 10)
CTLA4.1026_2$bc <- paste0("AML1026.2_", CTLA4.1026_2$bc, "-1")
CTLA4.1026_2 <- CTLA4.1026_2[which(CTLA4.1026_2$bc %in% colnames(AML1026)), ]
CTLA4.1026_2$ratio <- CTLA4.1026_2$detected / (CTLA4.1026_2$detected + CTLA4.1026_2$not.detected)
CTLA4.1026_2$UMAP1 <- AML1026@reductions$umap@cell.embeddings[CTLA4.1026_2$bc, "UMAP_1"]
CTLA4.1026_2$UMAP2 <- AML1026@reductions$umap@cell.embeddings[CTLA4.1026_2$bc, "UMAP_2"]
CTLA4.1026_2$sample <- "AML1026_2"

CTLA4.3005_1 <- extract_isoforms("./data/isoforms/AML3005.1/AML3005_1_exons.gz", GENE = "CTLA4", EXON = "exon2", filter = 10)
CTLA4.3005_1$bc <- paste0("AML3005.1_", CTLA4.3005_1$bc, "-1")
CTLA4.3005_1 <- CTLA4.3005_1[which(CTLA4.3005_1$bc %in% colnames(AML3005)), ]
CTLA4.3005_1$ratio <- CTLA4.3005_1$detected / (CTLA4.3005_1$detected + CTLA4.3005_1$not.detected)
CTLA4.3005_1$UMAP1 <- AML3005@reductions$umap@cell.embeddings[CTLA4.3005_1$bc, "UMAP_1"]
CTLA4.3005_1$UMAP2 <- AML3005@reductions$umap@cell.embeddings[CTLA4.3005_1$bc, "UMAP_2"]
CTLA4.3005_1$sample <- "AML3005_1"

CTLA4.3005_3 <- extract_isoforms("./data/isoforms/AML3005.3/AML3005_3_exons.gz", GENE = "CTLA4", EXON = "exon2", filter = 10)
CTLA4.3005_3$bc <- paste0("AML3005.3_", CTLA4.3005_3$bc, "-1")
CTLA4.3005_3 <- CTLA4.3005_3[which(CTLA4.3005_3$bc %in% colnames(AML3005)), ]
CTLA4.3005_3$ratio <- CTLA4.3005_3$detected / (CTLA4.3005_3$detected + CTLA4.3005_3$not.detected)
CTLA4.3005_3$UMAP1 <- AML3005@reductions$umap@cell.embeddings[CTLA4.3005_3$bc, "UMAP_1"]
CTLA4.3005_3$UMAP2 <- AML3005@reductions$umap@cell.embeddings[CTLA4.3005_3$bc, "UMAP_2"]
CTLA4.3005_3$sample <- "AML3005_3"

CTLA4.TIL3 <- extract_isoforms("./data/isoforms/TIL3/TIL3_exons.gz", GENE = "CTLA4", EXON = "exon2", filter = 10)
# CTLA4.TIL3$bc = paste0('AML3005.3_',CTLA4.TIL3$bc, '-1')
CTLA4.TIL3 <- CTLA4.TIL3[which(CTLA4.TIL3$bc %in% colnames(TIL3)), ]
CTLA4.TIL3$ratio <- CTLA4.TIL3$detected / (CTLA4.TIL3$detected + CTLA4.TIL3$not.detected)
CTLA4.TIL3$UMAP1 <- TIL3@reductions$umap@cell.embeddings[CTLA4.TIL3$bc, "UMAP_1"]
CTLA4.TIL3$UMAP2 <- TIL3@reductions$umap@cell.embeddings[CTLA4.TIL3$bc, "UMAP_2"]
CTLA4.TIL3$sample <- "TIL3"

CTLA4.TIL4 <- extract_isoforms("./data/isoforms/TIL4/TIL4_exons.gz", GENE = "CTLA4", EXON = "exon2", filter = 10)
CTLA4.TIL4 <- CTLA4.TIL4[which(CTLA4.TIL4$bc %in% colnames(TIL4)), ]
CTLA4.TIL4$ratio <- CTLA4.TIL4$detected / (CTLA4.TIL4$detected + CTLA4.TIL4$not.detected)
CTLA4.TIL4$UMAP1 <- TIL4@reductions$umap@cell.embeddings[CTLA4.TIL4$bc, "UMAP_1"]
CTLA4.TIL4$UMAP2 <- TIL4@reductions$umap@cell.embeddings[CTLA4.TIL4$bc, "UMAP_2"]
CTLA4.TIL4$sample <- "TIL4"

CTLA4.PBMC7 <- extract_isoforms("./data/isoforms/PBMC7/PBMC7_exons.gz", GENE = "CTLA4", EXON = "exon2", filter = 10)
CTLA4.PBMC7 <- CTLA4.PBMC7[which(CTLA4.PBMC7$bc %in% colnames(PBMC7)), ]
CTLA4.PBMC7$ratio <- CTLA4.PBMC7$detected / (CTLA4.PBMC7$detected + CTLA4.PBMC7$not.detected)
CTLA4.PBMC7$UMAP1 <- PBMC7@reductions$umap@cell.embeddings[CTLA4.PBMC7$bc, "UMAP_1"]
CTLA4.PBMC7$UMAP2 <- PBMC7@reductions$umap@cell.embeddings[CTLA4.PBMC7$bc, "UMAP_2"]
CTLA4.PBMC7$sample <- "PBMC7"

CTLA4.PBMC8 <- extract_isoforms("./data/isoforms/PBMC8/PBMC8_exons.gz", GENE = "CTLA4", EXON = "exon2", filter = 10)
# CTLA4.PBMC8$bc = paste0('AML3005.3_',CTLA4.PBMC8$bc, '-1')
CTLA4.PBMC8 <- CTLA4.PBMC8[which(CTLA4.PBMC8$bc %in% colnames(PBMC8)), ]
CTLA4.PBMC8$ratio <- CTLA4.PBMC8$detected / (CTLA4.PBMC8$detected + CTLA4.PBMC8$not.detected)
CTLA4.PBMC8$UMAP1 <- PBMC8@reductions$umap@cell.embeddings[CTLA4.PBMC8$bc, "UMAP_1"]
CTLA4.PBMC8$UMAP2 <- PBMC8@reductions$umap@cell.embeddings[CTLA4.PBMC8$bc, "UMAP_2"]
CTLA4.PBMC8$sample <- "PBMC8"

# removing libraries with too few evenst
CTLA4 <- dplyr::bind_rows( # CTLA4.1007_1[,c('sample', 'bc', 'detected', 'not.detected', 'ratio')],
  # CTLA4.1007_3[,c('sample', 'bc', 'detected', 'not.detected', 'ratio')],
  CTLA4.1007_5[, c("sample", "bc", "detected", "not.detected", "ratio")],
  CTLA4.3005_1[, c("sample", "bc", "detected", "not.detected", "ratio")],
  CTLA4.3005_3[, c("sample", "bc", "detected", "not.detected", "ratio")],
  # CTLA4.1010_1[,c('sample', 'bc', 'detected', 'not.detected', 'ratio')],
  # CTLA4.1026_2[,c('sample', 'bc', 'detected', 'not.detected', 'ratio')],
  CTLA4.TIL3[, c("sample", "bc", "detected", "not.detected", "ratio")],
  CTLA4.TIL4[, c("sample", "bc", "detected", "not.detected", "ratio")],
  CTLA4.PBMC7[, c("sample", "bc", "detected", "not.detected", "ratio")],
  CTLA4.PBMC8[, c("sample", "bc", "detected", "not.detected", "ratio")],
  CTLA4.97_6[, c("sample", "bc", "detected", "not.detected", "ratio")]
)

CTLA4$tissue <- "AML"
CTLA4$tissue[which(CTLA4$sample %in% c("TIL3", "TIL4"))] <- "TIL"
CTLA4$tissue[which(CTLA4$sample %in% c("PBMC7", "PBMC8"))] <- "TIL.PB"
CTLA4$tissue[which(CTLA4$sample %in% c("97_6"))] <- "CAR"

ggplot(CTLA4, aes(x = reorder(sample, ratio), y = 100 * ratio, color = tissue)) +
  geom_jitter(size = 0.5) +
  stat_summary(geom = "crossbar", fun = mean, fun.min = mean, fun.max = mean, color = "black") +
  scale_x_discrete() +
  scale_y_continuous("% CTLA-4 exon3") +
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
ggsave("./isoforms/figures/combined/plots/20230104_CTLA4_exon2_overview.svg", width = 2.5, height = 2.5)

p <- ggplot() +
  geom_point(data = as.data.frame(CAR.so@reductions$umap@cell.embeddings), aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = CTLA4.97_6, aes(x = UMAP1, y = UMAP2, color = ratio), size = 0.5) +
  scale_color_gradient(low = "black", high = "yellow") +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./isoforms/figures/CAR/UMAP/20230104_97_6_CTLA.png", width = 4, height = 4, dpi = 600, plot = p)

p <- FeaturePlot(CAR.so, features = "CTLA4", order = T) + NoAxes() + NoLegend() +
  theme(plot.title = element_blank(), panel.spacing = unit(0, "pt")) + scale_color_gradientn(colours = c("grey90", BuenColors::jdb_palette(name = "brewer_green")))
ggsave("./isoforms/figures/CAR/UMAP/20230104_97_6_CTLA4_original.png", width = 4, height = 4, dpi = 600, plot = p)

p <- ggplot() +
  geom_point(data = as.data.frame(TIL3@reductions$umap@cell.embeddings), aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = CTLA4.TIL3[order(-CTLA4.TIL3$ratio), ], aes(x = UMAP1, y = UMAP2, color = ratio), size = 0.5) +
  scale_color_gradient(low = "black", high = "yellow") +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./isoforms/figures/melanoma/UMAP/20230104_TIL3_CTLA.png", width = 4, height = 4, dpi = 600, plot = p)

p <- FeaturePlot(TIL3, features = "CTLA4", order = T) + NoAxes() + NoLegend() +
  theme(plot.title = element_blank(), panel.spacing = unit(0, "pt")) + scale_color_gradientn(colours = c("grey90", BuenColors::jdb_palette(name = "brewer_green")))
ggsave("./isoforms/figures/melanoma/UMAP/20230104_TIL3_CTLA4_original.png", width = 4, height = 4, dpi = 600, plot = p)

p <- ggplot() +
  geom_point(data = as.data.frame(AML3005@reductions$umap@cell.embeddings), aes(x = UMAP_1, y = UMAP_2), color = "grey90", size = 0.5) +
  geom_point(data = CTLA4.3005_3[order(-CTLA4.3005_3$ratio), ], aes(x = UMAP1, y = UMAP2, color = ratio), size = 0.5) +
  scale_color_gradient(low = "black", high = "yellow") +
  theme_classic() +
  NoAxes() +
  NoLegend()
ggsave("./isoforms/figures/AML3005/UMAP/20230104_AML3005_1_CTLA.png", width = 4, height = 4, dpi = 600, plot = p)

p <- FeaturePlot(AML3005, features = "CTLA4", order = T, reduction = "umap") + NoAxes() + NoLegend() +
  theme(plot.title = element_blank(), panel.spacing = unit(0, "pt")) + scale_color_gradientn(colours = c("grey90", BuenColors::jdb_palette(name = "brewer_green")))
ggsave("./isoforms/figures/AML3005/UMAP/20230104_AML3005_1_CTLA4_original.png", width = 4, height = 4, dpi = 600, plot = p)

p <- FeaturePlot(AML3005, features = "CTLA4", order = T, reduction = "umap") +
  theme(plot.title = element_blank(), panel.spacing = unit(0, "pt")) + scale_color_gradientn(colours = c("grey90", BuenColors::jdb_palette(name = "brewer_green")))
ggsave("./isoforms/figures/AML3005/UMAP/20230104_AML3005_1_CTLA4_original.svg", width = 4, height = 4, dpi = 600, plot = p)
