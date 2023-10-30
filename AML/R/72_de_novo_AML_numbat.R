# process numbat output for de-novo AML1 and de-novo AML3

library(ComplexHeatmap)
library(ggplot2)
library(numbat)
library(Seurat)

FLT3_ITD <- readRDS("./data/AML/objects/20230523_FLT3-ITD.rds")

# de-novo AML1 (FLT3-ITD1) - use loh(13)
nb <- Numbat$new(out_dir = "./data/AML/numbat/FLT3_ITD1")
nb$joint_post$present <- ifelse(nb$joint_post$p_cnv > 0.8, 1, 0)
cnv.data <- as.data.frame(nb$joint_post %>% filter(seg %in% c("13a")) %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y, present))
cnv.data$p_cnv <- as.numeric(cnv.data$p_cnv)
cnv.data <- as.data.frame(cnv.data %>% select(cell, seg, p_cnv) %>% filter(seg %in% c("13a")) %>% tidyr::pivot_wider(names_from = "seg", values_from = "p_cnv"))
rownames(cnv.data) <- cnv.data$cell
colnames(cnv.data) <- c("bc", "loh13")
cnv.data$mutated <- ifelse(cnv.data$loh13 > 0.8, "mutated", "wildtype")
cnv.data$predicted.celltype <- FLT3_ITD$predicted.celltype[rownames(cnv.data)]
cnv.data$sample <- "FLT3-ITD1"
write.table(cnv.data, file = "./data/AML/numbat/FLT3_ITD1/FLT3-ITD1_cnv_calls.csv", sep = "\t", quote = F)

p <- DimPlot(FLT3_ITD,
  reduction = "ref.umap", cells.highlight = cnv.data$bc[which(cnv.data$mutated == "mutated")],
  cols.highlight = "purple", sizes.highlight = 2
) +
  NoLegend() + NoAxes()
ggsave("./AML/figures/de_novo/UMAP/20230806_de_novo_AML1_FLT3-ITD1_loh13.png", width = 4, height = 4, dpi = 600, plot = p)

# de-novo AML3 (FLT3-ITD2) - use tri(8)
nb <- Numbat$new(out_dir = "./data/AML/numbat/FLT3_ITD2/")
nb$joint_post$present <- ifelse(nb$joint_post$p_cnv > 0.8, 1, 0)
cnv.data <- as.data.frame(nb$joint_post %>% filter(seg %in% c("8a")) %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y, present))
cnv.data$p_cnv <- as.numeric(cnv.data$p_cnv)
cnv.data <- as.data.frame(cnv.data %>% select(cell, seg, p_cnv) %>% filter(seg %in% c("8a")) %>% tidyr::pivot_wider(names_from = "seg", values_from = "p_cnv"))
rownames(cnv.data) <- cnv.data$cell
colnames(cnv.data) <- c("bc", "tri8")
cnv.data$mutated <- ifelse(cnv.data$tri8 > 0.8, "mutated", "wildtype")
cnv.data$predicted.celltype <- FLT3_ITD$predicted.celltype[rownames(cnv.data)]
cnv.data$sample <- "FLT3-ITD2"
write.table(cnv.data, file = "./data/AML/numbat/FLT3_ITD2/FLT3-ITD2_cnv_calls.csv", sep = "\t", quote = F)

p <- DimPlot(FLT3_ITD,
  reduction = "ref.umap", cells.highlight = cnv.data$bc[which(cnv.data$mutated == "mutated")],
  cols.highlight = "black", sizes.highlight = 2
) +
  NoLegend() + NoAxes()
ggsave("./AML/figures/de_novo/UMAP/20230806_de_novo_AML3_FLT3-ITD2_tri8.png", width = 4, height = 4, dpi = 600, plot = p)
