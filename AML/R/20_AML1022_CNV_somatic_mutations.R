# visualize somatic mutations and CNV changes in AML1022.1

library(ComplexHeatmap)
library(ggplot2)
library(nanoranger.R)
library(Seurat)

AML1022 <- readRDS("./data/AML/objects/20220729_AML1022_1.rds")

cnv.data <- read.csv2("./data/AML/numbat/AML1022.1/cnv_calls.csv", sep = "\t")
# cnv.data$bc = paste0('AML1022.1_', cnv.data$cell)
cnv.data$p_cnv <- as.numeric(cnv.data$p_cnv)
cnv.data <- cnv.data %>%
  select(cell, seg, p_cnv) %>%
  filter(seg %in% c("20b")) %>%
  tidyr::pivot_wider(names_from = "seg", values_from = "p_cnv")
colnames(cnv.data) <- c("bc", "20b")

mutations.df <- data.frame()
for (f in c("AML1022_DNMT3A.csv", "AML1022_RUNX1.csv", "AML1022_SF3B1.csv")) {
  boo <- as.data.frame(read.csv2(file = paste0("./data/AML/mutations/", f), sep = "\t"))
  boo$gene <- gsub(stringr::str_split_fixed(f, pattern = "_", n = 2)[, 2], pattern = ".csv", replacement = "")
  mutations.df <- rbind(mutations.df, boo[, c("bc", "alt", "ref", "mutated", "vaf", "gene")])
}
mutations.df$sample <- stringr::str_split_fixed(mutations.df$bc, pattern = "_", n = 2)[, 1]
mutations.df.condensed <- mutations.df %>%
  group_by(bc) %>%
  summarize(mutation = ifelse("mutated" %in% mutated, "mutated", "wildtype"))
mutations.df.condensed <- merge(mutations.df.condensed, cnv.data, by = "bc", all.x = T)
mutations.df.condensed <- merge(mutations.df.condensed,
  tidyr::pivot_wider(
    data = mutations.df[, c("bc", "mutated", "gene")],
    values_from = "mutated", names_from = "gene"
  )[, c("bc", "DNMT3A", "RUNX1", "SF3B1")],
  by = "bc", all.x = T
)
mutations.df.condensed$predicted.celltype <- AML1022$predicted.celltype[mutations.df.condensed$bc]
mutations.df.condensed$sample <- AML1022$orig.ident[mutations.df.condensed$bc]
rownames(mutations.df.condensed) <- mutations.df.condensed$bc
mutations.df.condensed <- mutations.df.condensed[intersect(mutations.df.condensed$bc, colnames(AML1022)), ]

AML.subset <- subset(AML1022, cells = mutations.df.condensed$bc)
AML.subset <- ScaleData(AML.subset, features = rownames(AML.subset))
DefaultAssay(AML.subset) <- "ADT"
AML.subset <- ScaleData(AML.subset, features = rownames(AML.subset))

marker.genes <- c("CD31", "CD33", "CD117", "HLADR", "CD14", "CD15", "CD16", "CD3", "CD4", "CD8a", "CD56")
marker.genes.RNA <- c("CD34", "HLA-DRA", "CD33", "KIT", "MPO", "CD38", "CD14", "GATA2", "ZFPM1", "FLI1", "NFE2", "PF4", "GATA1", "HBA1", "HBB", "CD3D", "CD4", "CD8A", "NCAM1", "MKI67")

mutations.df.condensed <- cbind(mutations.df.condensed, t(GetAssayData(AML.subset, assay = "ADT", slot = "scale.data")[marker.genes, ]))
mutations.df.condensed <- cbind(mutations.df.condensed, t(GetAssayData(AML.subset, assay = "RNA", slot = "scale.data")[marker.genes.RNA, ]))

mutations.df.condensed$predicted.celltype[which(mutations.df.condensed$predicted.celltype %in%
  c("CD8 Naive", "CD8 Effector_1", "CD8 Effector_2", "CD8 Memory_1", "CD8 Memory_2"))] <- "CD8"
mutations.df.condensed$predicted.celltype[which(mutations.df.condensed$predicted.celltype %in%
  c("CD4 Naive", "CD4 Memory", "Treg"))] <- "CD4"
mutations.df.condensed$predicted.celltype[which(mutations.df.condensed$predicted.celltype %in%
  c("NK", "CD56 bright NK"))] <- "NK cells"
mutations.df.condensed$predicted.celltype[which(mutations.df.condensed$predicted.celltype %in%
  c("gdT", "MAIT"))] <- "other"

mutations.df.condensed <- mutations.df.condensed[order(mutations.df.condensed$`20b`), ]

cells.1 <- mutations.df.condensed$bc[which(mutations.df.condensed$sample == "AML1022.1" &
  mutations.df.condensed$predicted.celltype %in% c("HSC", "LMPP", "GMP", "Prog_Mk", "Prog_RBC"))]
# cells.1 = c(cells.1, mutations.df.condensed$bc[which(mutations.df.condensed$sample == 'AML1022.1' &
#                                              mutations.df.condensed$predicted.celltype %in% c('CD4', 'CD8', 'NK cells'))])

ha <- columnAnnotation(
  celltype = mutations.df.condensed[cells.1, "predicted.celltype"],
  DNMT3A = mutations.df.condensed[cells.1, "DNMT3A"],
  RUNX1 = mutations.df.condensed[cells.1, "RUNX1"],
  SF3B1 = mutations.df.condensed[cells.1, "SF3B1"],
  annotation_name_gp = gpar(fontsize = 8),
  col = list(
    "celltype" = c(AML.combined.colors, "CD4" = "lightblue", "CD8" = "blue", "NK cells" = "purple"),
    "DNMT3A" = c("mutated" = "red", "wildtype" = "white"),
    "RUNX1" = c("mutated" = "red", "wildtype" = "white"),
    "SF3B1" = c("mutated" = "red", "wildtype" = "white")
  ),
  na_col = "grey90", border = T, simple_anno_size = unit(7, "pt")
)

col_fun <- circlize::colorRamp2(breaks = seq(0, 1, 1 / 8), colors = BuenColors::jdb_palette(name = "solar_rojos"))
h1 <- Heatmap(t(mutations.df.condensed[cells.1, c("20b")]),
  column_split = factor(mutations.df.condensed[cells.1, "predicted.celltype"], levels = c(names(AML.combined.colors), "CD4", "CD8", "NK cells")),
  cluster_rows = F, cluster_columns = F, show_column_names = F, border = T, top_annotation = ha, column_title_gp = gpar(fontsize = 0),
  column_title_rot = 90, na_col = "grey90", col = col_fun, row_labels = c("del(20q)"),
  row_names_side = "left", row_names_gp = gpar(fontsize = 8), use_raster = T, raster_quality = 10
)


col_fun <- circlize::colorRamp2(breaks = seq(0, 2, 2 / 8), colors = BuenColors::jdb_palette(name = "brewer_spectra"))
h2 <- Heatmap(t(mutations.df.condensed[cells.1, marker.genes]),
  column_split = factor(mutations.df.condensed[cells.1, "predicted.celltype"], levels = c(names(AML.combined.colors), "CD4", "CD8", "NK cells")),
  cluster_rows = F, cluster_columns = F, show_column_names = F, border = T, column_title_side = "bottom",
  column_title_rot = 90, na_col = "white", col = col_fun,
  row_names_side = "left", row_names_gp = gpar(fontsize = 8), use_raster = T, raster_quality = 10
)

col_fun <- circlize::colorRamp2(breaks = seq(0, 4, 4 / 8), colors = BuenColors::jdb_palette(name = "brewer_purple"))
h3 <- Heatmap(t(mutations.df.condensed[cells.1, marker.genes.RNA]),
  column_split = factor(mutations.df.condensed[cells.1, "predicted.celltype"], levels = c(names(AML.combined.colors), "CD4", "CD8", "NK cells")),
  cluster_rows = F, cluster_columns = F, show_column_names = F, border = T, column_title_side = "bottom",
  column_title_rot = 90, na_col = "white", col = col_fun,
  row_names_side = "left", row_names_gp = gpar(fontsize = 8), use_raster = T, raster_quality = 10
)

svglite::svglite("./AML/figures/AML1022/heatmaps/20221011_AML1022_1.svg", width = 6, height = 2.8)
draw(h1 %v% h3)
dev.off()
