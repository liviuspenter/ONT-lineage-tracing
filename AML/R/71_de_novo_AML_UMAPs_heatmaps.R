library(ComplexHeatmap)
library(ggplot2)
library(nanoranger.R)
library(Seurat)

FLT3_ITD <- readRDS("./data/AML/objects/20230523_FLT3-ITD.rds")

# de-novo AML1 (FLT3_ITD1)
FLT3_ITD1 <- read.csv2("./data/AML/mutations/FLT3-ITD1_FLT3-ITD.csv", sep = "\t")
FLT3_ITD1_NPM1_1 <- read.csv2("./data/AML/mutations/FLT3-ITD1_NPM1.csv", sep = "\t")
FLT3_ITD1_NPM1_2 <- read.csv2("./data/AML/mutations/FLT3-ITD1_NPM1.2.csv", sep = "\t")

FLT3_ITD1.mutated <- FLT3_ITD1$bc[which(FLT3_ITD1$mutated == "mutated")]
FLT3_ITD1_NPM1_1.mutated <- FLT3_ITD1_NPM1_1$bc[which(FLT3_ITD1_NPM1_1$mutated == "mutated")]
FLT3_ITD1_NPM1_2.mutated <- FLT3_ITD1_NPM1_2$bc[which(FLT3_ITD1_NPM1_2$mutated == "mutated")]

p <- DimPlot(FLT3_ITD,
             reduction = "ref.umap", cells.highlight = FLT3_ITD1_NPM1_1.mutated,
             cols.highlight = "purple", sizes.highlight = 2
) +
  NoLegend() + NoAxes()
q <- DimPlot(FLT3_ITD,
             reduction = "ref.umap", cells.highlight = FLT3_ITD1_NPM1_2.mutated,
             cols.highlight = "darkgreen", sizes.highlight = 2
) +
  NoLegend() + NoAxes()
r <- DimPlot(FLT3_ITD,
             reduction = "ref.umap", cells.highlight = FLT3_ITD1.mutated,
             cols.highlight = "purple", sizes.highlight = 2
) +
  NoLegend() + NoAxes()
ggsave("./AML/figures/de_novo/UMAP/20230526_de_novo_AML1_FLT3-ITD1_NPM1_1.png", width = 4, height = 4, dpi = 600, plot = p)
ggsave("./AML/figures/de_novo/UMAP/20230526_de_novo_AML1_FLT3-ITD1_NPM1_2.png", width = 4, height = 4, dpi = 600, plot = q)
ggsave("./AML/figures/de_novo/UMAP/20230526_de_novo_AML1_FLT3-ITD1_FLT3-ITD.png", width = 4, height = 4, dpi = 600, plot = r)

# create heatmap with both clones
boo <- merge(FLT3_ITD1_NPM1_1[, c("bc", "mutated")], FLT3_ITD1_NPM1_2[, c("bc", "mutated")], by = "bc", all = F)
boo <- merge(boo[, c("bc", "mutated.x", "mutated.y")], FLT3_ITD1[, c("bc", "mutated")], by = "bc", all = F)
colnames(boo) <- c("bc", "NPM1_1", "NPM1_2", "FLT3-ITD")
boo[is.na(boo)] <- -1
boo[boo == "mutated"] <- 1
boo[boo == "wildtype"] <- 0

col_fun <- c("1" = "red", "0" = "white", "-1" = "grey90")

svglite::svglite("./AML/figures/de_novo/heatmaps/20230526_de_novo_AML1_FLT3-ITD1.svg", width = 3, height = 1.2)
Heatmap(t(boo[order(boo$NPM1_2), c("NPM1_1", "NPM1_2", "FLT3-ITD")]),
        cluster_columns = F, col = col_fun,
        border = T, show_column_names = F, row_labels = c("NPM1W287fs", "NPM1W288fs", "FLT3-ITD"),
        row_names_side = "left", row_names_gp = gpar(fontsize = 10)
)
dev.off()

# de-novo AML3 (FLT3_ITD2)
FLT3_ITD2 <- read.csv2("./data/AML/mutations/FLT3-ITD2_FLT3-ITD.csv", sep = "\t")
FLT3_ITD2.mutated <- FLT3_ITD2$bc[which(FLT3_ITD2$mutated == "mutated")]

r <- DimPlot(FLT3_ITD,
             reduction = "ref.umap", cells.highlight = FLT3_ITD2.mutated,
             cols.highlight = "black", sizes.highlight = 2
) +
  NoLegend() + NoAxes()
ggsave("./AML/figures/de_novo/UMAP/20230526_de_novo_AML3_FLT3-ITD2_FLT3-ITD.png", width = 4, height = 4, dpi = 600, plot = r)

# de-novo AML2 (FLT3_ITD3)
FLT3_ITD3 <- read.csv2("./data/AML/mutations/FLT3-ITD3_FLT3-ITD.csv", sep = "\t")
FLT3_ITD3_NPM1 <- read.csv2("./data/AML/mutations/FLT3-ITD3_NPM1.csv", sep = "\t")
FLT3_ITD3_DNMT3A <- read.csv2("./data/AML/mutations/FLT3-ITD3_DNMT3A.csv", sep = "\t")

FLT3_ITD3.mutated <- FLT3_ITD3$bc[which(FLT3_ITD3$mutated == "mutated")]
FLT3_ITD3_NPM1.mutated <- FLT3_ITD3_NPM1$bc[which(FLT3_ITD3_NPM1$mutated == "mutated")]
FLT3_ITD3_DNMT3A.mutated <- FLT3_ITD3_DNMT3A$bc[which(FLT3_ITD3_DNMT3A$mutated == "mutated")]

p <- DimPlot(FLT3_ITD,
             reduction = "ref.umap", cells.highlight = FLT3_ITD3.mutated,
             cols.highlight = "black", sizes.highlight = 2
) +
  NoLegend() + NoAxes()
q <- DimPlot(FLT3_ITD,
             reduction = "ref.umap", cells.highlight = FLT3_ITD3_NPM1.mutated,
             cols.highlight = "black", sizes.highlight = 2
) +
  NoLegend() + NoAxes()
r <- DimPlot(FLT3_ITD,
             reduction = "ref.umap", cells.highlight = FLT3_ITD3_DNMT3A.mutated,
             cols.highlight = "black", sizes.highlight = 2
) +
  NoLegend() + NoAxes()
ggsave("./AML/figures/de_novo/UMAP/20230526_de_novo_AML2_FLT3-ITD3_FLT3-ITD.png", width = 4, height = 4, dpi = 600, plot = p)
ggsave("./AML/figures/de_novo/UMAP/20230526_de_novo_AML2_FLT3-ITD3_NPM1.png", width = 4, height = 4, dpi = 600, plot = q)
ggsave("./AML/figures/de_novo/UMAP/20230526_de_novo_AML2_FLT3-ITD3_DNMT3A.png", width = 4, height = 4, dpi = 600, plot = r)
