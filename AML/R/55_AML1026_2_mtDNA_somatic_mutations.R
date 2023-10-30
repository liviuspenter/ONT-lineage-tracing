### combined with NRAS and SF3B1 mutation

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(nanoranger.R)
library(Seurat)
library(SummarizedExperiment)


AML1026.2 <- readRDS("./data/AML/objects/20220711_AML1026_2.rds")

# functionality from MAESTER (https://github.com/petervangalen/MAESTER-2021)
source("./R/210215_FunctionsGeneral.R")

# get mtDNA variants
magtk.output <- readRDS("./data/AML/mtDNA/AML1026_2.maegatk/AML1026_2.rds")
af.dm <- data.matrix(computeAFMutMatrix(magtk.output)) * 100
rownames(af.dm) <- gsub(rownames(af.dm), pattern = "_", replacement = "")
af.dm.clean <- af.dm[-which(rowSums(af.dm) == 0), ]

NRAS.vaf2 <- read.csv2(file = "./data/AML/mtDNA/20220805_AML1026_2_NRAS.csv", row.names = "X")
rownames(NRAS.vaf2) <- NRAS.vaf2$bc
SF3B1.vaf2 <- read.csv2(file = "./data/AML/mtDNA/20220805_AML1026_2_SF3B1.csv", row.names = "X")
rownames(SF3B1.vaf2) <- SF3B1.vaf2$bc
SRSF2.vaf2 <- read.csv2(file = "./data/AML/mtDNA/20220805_AML1026_2_SRSF2.csv", row.names = "X")
rownames(SRSF2.vaf2) <- SRSF2.vaf2$bc

df <- as.data.frame(t(af.dm.clean[c("10685G>A", "3106C>A", "3106C>G", "3106C>T", "15615G>A", "9254A>G"), ]))
rownames(df) <- paste0("AML1026.2_", rownames(df), "-1")
df$NRAS <- NRAS.vaf2[rownames(df), "mutated"]
df$NRAS <- ifelse(df$NRAS == "mutated", 100, 0)
df$SF3B1 <- SF3B1.vaf2[rownames(df), "mutated"]
df$SF3B1 <- ifelse(df$SF3B1 == "mutated", 100, 0)
df$SRSF2 <- SRSF2.vaf2[rownames(df), "mutated"]
df$SRSF2 <- ifelse(df$SRSF2 == "mutated", 100, 0)


df <- df[intersect(rownames(df), colnames(AML1026.2)), ]
df <- cbind(df, AML1026.2@reductions$ref.umap@cell.embeddings[rownames(df), ])

p <- ggplot(df, aes(x = refUMAP_1, y = refUMAP_2, color = `10685G>A`)) +
  geom_point() +
  scale_color_gradientn(colours = c("grey90", BuenColors::jdb_palette(name = "solar_rojos"))) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )
ggsave("./AML/figures/mtDNA/UMAP/20221010_AML1026_10685G>A.png", width = 4, height = 4, dpi = 600, plot = p)

p <- ggplot() +
  geom_point(data = df[which(is.na(df$NRAS)), ], aes(x = refUMAP_1, y = refUMAP_2), color = "grey90") +
  geom_point(data = df[which(!is.na(df$NRAS)), ], aes(x = refUMAP_1, y = refUMAP_2, color = NRAS)) +
  scale_color_gradient(low = "black", high = "firebrick") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )
ggsave("./AML/figures/mtDNA/UMAP/20221010_AML1026_NRAS.png", width = 4, height = 4, dpi = 600, plot = p)

p <- ggplot() +
  geom_point(data = df[which(is.na(df$SF3B1)), ], aes(x = refUMAP_1, y = refUMAP_2), color = "grey90") +
  geom_point(data = df[which(!is.na(df$SF3B1)), ], aes(x = refUMAP_1, y = refUMAP_2, color = SF3B1)) +
  geom_point() +
  scale_color_gradient(low = "black", high = "firebrick") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )
ggsave("./AML/figures/mtDNA/UMAP/20221010_AML1026_SF3B1.png", width = 4, height = 4, dpi = 600, plot = p)


p <- ggplot() +
  geom_point(data = df[which(is.na(df$SRSF2)), ], aes(x = refUMAP_1, y = refUMAP_2), color = "grey90") +
  geom_point(data = df[which(!is.na(df$SRSF2)), ], aes(x = refUMAP_1, y = refUMAP_2, color = SRSF2)) +
  geom_point() +
  scale_color_gradient(low = "black", high = "firebrick") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )
ggsave("./AML/figures/mtDNA/UMAP/20221010_AML1026_SRSF2.png", width = 4, height = 4, dpi = 600, plot = p)

colnames(af.dm.clean) <- paste0("AML1026.2_", colnames(af.dm.clean), "-1")

cells <- unique(c(SF3B1.vaf2$bc, NRAS.vaf2$bc, SRSF2.vaf2$bc))
cells <- cells[which(cells %in% colnames(af.dm.clean))]
df <- as.data.frame(t(af.dm.clean[c("10685G>A", "3106C>A", "3106C>G", "3106C>T", "15615G>A", "9254A>G"), cells]))

df$NRAS <- NRAS.vaf2[rownames(df), "mutated"]
df$NRAS <- ifelse(df$NRAS == "mutated", 100, 0)
df$SF3B1 <- SF3B1.vaf2[rownames(df), "mutated"]
df$SF3B1 <- ifelse(df$SF3B1 == "mutated", 100, 0)
df$SRSF2 <- SRSF2.vaf2[rownames(df), "mutated"]
df$SRSF2 <- ifelse(df$SRSF2 == "mutated", 100, 0)

# df = cbind(df, AML1026.2@reductions$ref.umap@cell.embeddings[rownames(df),])
df <- t(df)
df <- df[, which(colnames(df) %in% colnames(AML1026.2)[which(AML1026.2$predicted.celltype %in% c("HSC", "LMPP", "GMP", "CD4 Memory", "CD8 Memory_2", "NK"))])]
col_fun <- circlize::colorRamp2(breaks = seq(0, 100, 100 / 8), colors = BuenColors::jdb_palette(name = "solar_rojos"))

df <- df[, sample(ncol(df))]

ha <- columnAnnotation(
  celltype = AML1026.2$predicted.celltype[colnames(df)],
  col = list(celltype = AML.combined.colors),
  simple_anno_size = unit(5, "pt"), border = T
)
svglite::svglite("./AML/figures/mtDNA/heatmaps/20220805_AML1026_mtDNA_NRAS_SF3B1_SRSF2.svg", width = 4.5, height = 2)
Heatmap(df,
  cluster_rows = F, cluster_columns = F, show_row_names = T, show_column_names = F,
  column_split = factor(AML1026.2$predicted.celltype[colnames(df)], levels = c("HSC", "LMPP", "GMP", "CD4 Memory", "CD8 Memory_2", "NK")),
  row_names_side = "left", row_names_gp = gpar(fontsize = 8), row_split = c(rep("mtDNA", 6), rep("somatic", 3)), row_title_gp = gpar(fontsize = 0),
  column_title_side = "bottom", column_title_rot = 90, column_title_gp = gpar(fontsize = 8),
  col = col_fun, border = T, use_raster = T, raster_quality = 10, na_col = "grey90", top_annotation = ha
)
dev.off()
