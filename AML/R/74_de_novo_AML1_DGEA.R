# perform differential gene expression analysis between subclones

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(nanoranger.R)
library(Seurat)
library(SummarizedExperiment)

FLT3_ITD <- readRDS("./data/AML/objects/20230523_FLT3-ITD.rds")

FLT3_ITD1_NPM1_1 <- read.csv2("./data/AML/mutations/FLT3-ITD1_NPM1.csv", sep = "\t") %>% as.data.frame()
rownames(FLT3_ITD1_NPM1_1) <- FLT3_ITD1_NPM1_1$bc
FLT3_ITD1_NPM1_2 <- read.csv2("./data/AML/mutations/FLT3-ITD1_NPM1.2.csv", sep = "\t") %>% as.data.frame()
rownames(FLT3_ITD1_NPM1_2) <- FLT3_ITD1_NPM1_2$bc

cells.mutated.1 <- FLT3_ITD1_NPM1_1$bc[which(FLT3_ITD1_NPM1_1$mutated == "mutated")]
cells.mutated.2 <- FLT3_ITD1_NPM1_2$bc[which(FLT3_ITD1_NPM1_2$mutated == "mutated")]

FLT3_ITD$clone <- "none"
FLT3_ITD$clone[cells.mutated.1] <- "clone1"
FLT3_ITD$clone[cells.mutated.2] <- "clone2"

boo <- subset(FLT3_ITD, predicted.celltype == "GMP" & clone != "none")

markers <- FindMarkers(boo, group.by = "clone", ident.1 = "clone1", ident.2 = "clone2")

markers$gene <- rownames(markers)
markers$FDR <- -log10(markers$p_val_adj)
markers$significant <- ifelse(markers$FDR > 5 & abs(markers$avg_log2FC) > 2, "yes", "no")

ggplot() +
  geom_point(data = markers, aes(x = avg_log2FC, y = FDR), color = "grey", size = 0.5) +
  geom_hline(yintercept = 5) +
  geom_vline(xintercept = c(-2, 2)) +
  geom_point(data = markers[which(markers$significant == "yes"), ], aes(x = avg_log2FC, y = FDR), color = "firebrick", size = 0.5) +
  scale_x_continuous("log2FC", limits = c(-5, 5)) +
  scale_y_continuous("-log10(FDR)") +
  theme_classic() +
  theme(
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./AML/figures/de_novo/plots/20230803_de_novo_AML1_FLT3-ITD1_volcano.svg", width = 1.8, height = 1.5)

mat <- as.matrix(GetAssayData(boo)[markers$gene[which(markers$significant == "yes")], ], slot = "scale.data")
mat <- mat[, c(colnames(boo)[which(boo$clone == "clone1")], colnames(boo)[which(boo$clone == "clone2")])]

ha <- columnAnnotation(
  clone = factor(c(
    rep("clone1", length(which(boo$clone == "clone1"))),
    rep("clone2", length(which(boo$clone == "clone2")))
  )),
  col = list("clone" = c("clone1" = "darkgreen", "clone2" = "purple")),
  simple_anno_size = unit(5, "pt"), border = T, annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 10)
)
col_fun <- circlize::colorRamp2(breaks = seq(0, 6, 6 / 8), colors = BuenColors::jdb_palette(name = "brewer_purple"))
svglite::svglite("./AML/figures/de_novo/heatmaps/20230803_de_novo_AML1_FLT3_ITD1_DGEA.svg", width = 3, height = 2.5)
Heatmap(mat,
  show_column_names = F, cluster_rows = T, cluster_columns = F, border = T,
  row_names_side = "left", row_names_gp = gpar(fontsize = 8, fontface = "italic"), show_row_dend = F,
  use_raster = T, raster_quality = 10, col = col_fun, top_annotation = ha
)
dev.off()
