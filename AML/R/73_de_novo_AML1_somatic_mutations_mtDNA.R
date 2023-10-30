# explore mitochondrial and somatic mutations in AML1

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(nanoranger.R)
library(Seurat)
library(SummarizedExperiment)

FLT3_ITD <- readRDS("./data/AML/objects/20230523_FLT3-ITD.rds")

# functionality from MAESTER (https://github.com/petervangalen/MAESTER-2021)
source("./R/210215_FunctionsGeneral.R")

FLT3_ITD1_NPM1_1 <- read.csv2("./data/AML/mutations/FLT3-ITD1_NPM1.csv", sep = "\t") %>% as.data.frame()
rownames(FLT3_ITD1_NPM1_1) <- FLT3_ITD1_NPM1_1$bc
FLT3_ITD1_NPM1_2 <- read.csv2("./data/AML/mutations/FLT3-ITD1_NPM1.2.csv", sep = "\t") %>% as.data.frame()
rownames(FLT3_ITD1_NPM1_2) <- FLT3_ITD1_NPM1_2$bc
FLT3_ITD1_FLT3 <- read.csv2("./data/AML/mutations/FLT3-ITD1_FLT3-ITD.csv", sep = "\t") %>% as.data.frame()
rownames(FLT3_ITD1_FLT3) <- FLT3_ITD1_FLT3$bc

# find cells that are covered in either amplicon
cells <- unique(c(FLT3_ITD1_NPM1_1$bc, FLT3_ITD1_NPM1_2$bc))

# process maegatk output
magtk.output <- readRDS("./data/AML/mtDNA/FLT3-ITD1.maegatk/FLT3-ITD1.rds")
af.dm <- data.matrix(computeAFMutMatrix(magtk.output)) * 100
colnames(af.dm) <- paste0("FLT3-ITD1_", colnames(af.dm), "-1")
af.dm.clean <- af.dm[-which(rowSums(af.dm) == 0), ] %>% as.data.frame()
rownames(af.dm.clean) <- gsub(rownames(af.dm.clean), pattern = "_", replacement = "")
cells <- intersect(cells, colnames(af.dm))

# identify mutated and wildtype cells 
af.dm.clean <- af.dm.clean[, cells]
FLT3_ITD1_NPM1_1 <- FLT3_ITD1_NPM1_1[cells, ]
FLT3_ITD1_NPM1_2 <- FLT3_ITD1_NPM1_2[cells, ]

cells.mutated.1 <- FLT3_ITD1_NPM1_1$bc[which(FLT3_ITD1_NPM1_1$mutated == "mutated")]
cells.mutated.2 <- FLT3_ITD1_NPM1_2$bc[which(FLT3_ITD1_NPM1_2$mutated == "mutated")]

cells.non.mutated.1 <- FLT3_ITD1_NPM1_1$bc[which(FLT3_ITD1_NPM1_1$mutated == "wildtype")]
cells.non.mutated.2 <- FLT3_ITD1_NPM1_2$bc[which(FLT3_ITD1_NPM1_2$mutated == "wildtype")]

df <- data.frame(
  mutation = rownames(af.dm.clean),
  heteroplasmy.1 = rowMeans(af.dm.clean[, cells.mutated.1]),
  heteroplasmy.2 = rowMeans(af.dm.clean[, cells.mutated.2])
)
df$ratio <- df$heteroplasmy.1 / df$heteroplasmy.2

variants.1 <- c("2581A>G", "6776T>C", "13710A>T")
variants.2 <- c("3010G>A", "4733T>C", "14793A>C")

show.cells <- c(cells.mutated.1, cells.mutated.2)

# demonstrate mutual exclusivity of mtDNA mutations in both clones
ggplot() +
  geom_abline(slope = 1) +
  ggrastr::rasterize(geom_point(data = df, aes(x = heteroplasmy.1, y = heteroplasmy.2), color = "grey", size = 0.5), dpi = 600) +
  geom_point(data = df[c(variants.1, variants.2), ], aes(x = heteroplasmy.1, y = heteroplasmy.2), color = "firebrick", size = 0.5) +
  geom_label_repel(
    data = df[c(variants.1, variants.2), ], aes(x = heteroplasmy.1, y = heteroplasmy.2, label = mutation),
    label.size = 0, size = 3
  ) +
  scale_x_continuous("% heteroplasmy NPM1W288fs") +
  scale_y_continuous("% heteroplasmy NPM1W287fs") +
  theme_classic() +
  theme(
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./AML/figures/de_novo/plots/20230803_de_novo_AML1_FLT3-ITD1_mtDNA.svg", width = 2, height = 2)

# create heatmap illustrating both clones

ha <- columnAnnotation(
  NPM1_1 = FLT3_ITD1_NPM1_1[show.cells, "mutated"],
  NPM1_2 = FLT3_ITD1_NPM1_2[show.cells, "mutated"],
  FLT3_ITD = FLT3_ITD1_FLT3[show.cells, "mutated"],
  annotation_label = c("NPM1W287fs", "NPM1W288fs", "FLT3-ITD"),
  col = list(
    "NPM1_1" = c("mutated" = "firebrick", "wildtype" = "white"),
    "NPM1_2" = c("mutated" = "firebrick", "wildtype" = "white"),
    "FLT3_ITD" = c("mutated" = "firebrick", "wildtype" = "white")
  ),
  na_col = "grey90", border = T, annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10), simple_anno_size = unit(10, "pt")
)

col_fun <- circlize::colorRamp2(breaks = seq(0, 100, 100 / 8), colors = BuenColors::jdb_palette(name = "solar_extra"))
svglite::svglite("./AML/figures/de_novo/heatmaps/20230731_de_novo_AML1_FLT3-ITD1_combined_somatic_mtDNA.svg", width = 5, height = 1.5)
Heatmap(af.dm.clean[c(variants.1, variants.2), c(cells.mutated.1, cells.mutated.2)],
  show_row_names = T, show_column_names = F, cluster_rows = T, cluster_columns = F,
  show_column_dend = F, show_row_dend = F, row_names_gp = gpar(fontsize = 10),
  row_names_side = "left", border = T, na_col = "white",
  col = col_fun, use_raster = T, raster_quality = 10, top_annotation = ha,
  column_title = paste0(length(c(cells.mutated.1, cells.mutated.2)), " cells"),
  column_title_side = "bottom", column_title_gp = gpar(fontsize = 8)
)
dev.off()

### combine with numbat output

cnv.calls <- read.csv("./data/AML/numbat/FLT3_ITD1/FLT3-ITD1_cnv_calls.csv", sep = "\t") %>% as.data.frame()
cnv.calls[show.cells, ]

ha <- columnAnnotation(
  NPM1_1 = FLT3_ITD1_NPM1_1[show.cells, "mutated"],
  NPM1_2 = FLT3_ITD1_NPM1_2[show.cells, "mutated"],
  FLT3_ITD = FLT3_ITD1_FLT3[show.cells, "mutated"],
  LOH = cnv.calls[show.cells, "mutated"],
  annotation_label = c("NPM1W287fs", "NPM1W288fs", "FLT3-ITD", "loh(13)"),
  col = list(
    "NPM1_1" = c("mutated" = "firebrick", "wildtype" = "white"),
    "NPM1_2" = c("mutated" = "firebrick", "wildtype" = "white"),
    "FLT3_ITD" = c("mutated" = "firebrick", "wildtype" = "white"),
    "LOH" = c("mutated" = "firebrick", "wildtype" = "white")
  ),
  na_col = "grey90", border = T, annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10), simple_anno_size = unit(10, "pt")
)

col_fun <- circlize::colorRamp2(breaks = seq(0, 100, 100 / 8), colors = BuenColors::jdb_palette(name = "solar_extra"))

svglite::svglite("./AML/figures/de_novo/heatmaps/20230731_de_novo_AML1_FLT3-ITD1_combined_somatic_mtDNA_numbat.svg", width = 5, height = 1.5)
Heatmap(af.dm.clean[c(variants.1, variants.2), c(cells.mutated.1, cells.mutated.2)],
  show_row_names = T, show_column_names = F, cluster_rows = T, cluster_columns = F,
  show_column_dend = F, show_row_dend = F, row_names_gp = gpar(fontsize = 10),
  row_names_side = "left", border = T, na_col = "white",
  col = col_fun, use_raster = T, raster_quality = 10, top_annotation = ha,
  column_title = paste0(length(c(cells.mutated.1, cells.mutated.2)), " cells"),
  column_title_side = "bottom", column_title_gp = gpar(fontsize = 8)
)
dev.off()

### visualize individual mtDNA mutations on UMAP
af.dm.clean <- af.dm[-which(rowSums(af.dm) == 0), ] %>% as.data.frame()
rownames(af.dm.clean) <- gsub(rownames(af.dm.clean), pattern = "_", replacement = "")
df <- Embeddings(FLT3_ITD, reduction = "ref.umap") %>% as.data.frame()
cells <- intersect(colnames(af.dm.clean), rownames(df))
colnames(df) <- c("UMAP1", "UMAP2")
df <- df[cells, ]
df <- cbind(df, t(af.dm.clean[c(variants.1, variants.2), rownames(df)]))

ggplot() +
  geom_point(
    data = Embeddings(FLT3_ITD, reduction = "ref.umap") %>% as.data.frame(),
    aes(x = refUMAP_1, y = refUMAP_2), color = "grey", size = 0.5
  ) +
  geom_point(data = df, aes(x = UMAP1, y = UMAP2, color = `2581A>G`), size = 0.5) +
  scale_color_gradientn(colours = c(BuenColors::jdb_palette(name = "solar_extra"))) +
  theme_classic() +
  NoAxes() +
  theme(legend.position = "none")
ggsave("./AML/figures/de_novo/UMAP/20230801_de_novo_AML1_FLT3_ITD1_2581A>G.png", width = 4, height = 4, dpi = 600)

ggplot() +
  geom_point(
    data = Embeddings(FLT3_ITD, reduction = "ref.umap") %>% as.data.frame(),
    aes(x = refUMAP_1, y = refUMAP_2), color = "grey", size = 0.5
  ) +
  geom_point(data = df, aes(x = UMAP1, y = UMAP2, color = `6776T>C`), size = 0.5) +
  scale_color_gradientn(colours = c(BuenColors::jdb_palette(name = "solar_extra"))) +
  theme_classic() +
  NoAxes() +
  theme(legend.position = "none")
ggsave("./AML/figures/de_novo/UMAP/20230801_de_novo_AML1_FLT3_ITD1_6776T>C.png", width = 4, height = 4, dpi = 600)


ggplot() +
  geom_point(
    data = Embeddings(FLT3_ITD, reduction = "ref.umap") %>% as.data.frame(),
    aes(x = refUMAP_1, y = refUMAP_2), color = "grey", size = 0.5
  ) +
  geom_point(data = df, aes(x = UMAP1, y = UMAP2, color = `3010G>A`), size = 0.5) +
  scale_color_gradientn(colours = c(BuenColors::jdb_palette(name = "solar_extra"))) +
  theme_classic() +
  NoAxes() +
  theme(legend.position = "none")
ggsave("./AML/figures/de_novo/UMAP/20230801_de_novo_AML1_FLT3_ITD1_3010G>A.png", width = 4, height = 4, dpi = 600)

ggplot() +
  geom_point(
    data = Embeddings(FLT3_ITD, reduction = "ref.umap") %>% as.data.frame(),
    aes(x = refUMAP_1, y = refUMAP_2), color = "grey", size = 0.5
  ) +
  geom_point(data = df, aes(x = UMAP1, y = UMAP2, color = `4733T>C`), size = 0.5) +
  scale_color_gradientn(colours = c(BuenColors::jdb_palette(name = "solar_extra"))) +
  theme_classic() +
  NoAxes() +
  theme(legend.position = "none")
ggsave("./AML/figures/de_novo/UMAP/20230801_de_novo_AML1_FLT3_ITD1_4733T>C.png", width = 4, height = 4, dpi = 600)
