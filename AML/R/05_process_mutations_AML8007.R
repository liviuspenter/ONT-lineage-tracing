# process mutations for AML1007.14

library(ComplexHeatmap)
library(dplyr)
library(nanoranger.R)
library(ggplot2)
library(Seurat)
library(SummarizedExperiment)

AML8007 <- readRDS("./data/AML/objects/20220911_AML8007.rds")
Idents(AML8007) <- "predicted.celltype"
prop.table(table(Idents(AML8007), AML8007$orig.ident), margin = 2)

boo <- as.data.frame(prop.table(table(Idents(AML8007), AML8007$orig.ident), margin = 2))
boo$Var1 <- factor(boo$Var1, levels = rev(names(AML.combined.colors)))
ggplot(boo, aes(x = Var2, y = 100 * Freq, fill = Var1)) +
  geom_col() +
  geom_hline(yintercept = 5) +
  scale_x_discrete(labels = c("Screening", "marrow CR")) +
  scale_y_continuous("%cells", limits = c(0, 100)) +
  scale_fill_manual(values = AML.combined.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./AML/figures/AML8007/plots/20220912_celltype_kinetics.svg", width = 1.1, height = 2)

# DNMT3A
DNMT3A.vaf.1 <- extract_mutation(BC.data.file = "./data/AML/pileup/AML8007_1_pileup_DNMT3A.csv.gz", REF = "C", ALT = "T", FILTER = 10)
DNMT3A.vaf.1$bc <- paste0("AML8007.1_", DNMT3A.vaf.1$bc, "-1")
DNMT3A.vaf.1$sample <- "AML8007.1"
DNMT3A.vaf.4 <- extract_mutation(BC.data.file = "./data/AML/pileup/AML8007_4_pileup_DNMT3A.csv.gz", REF = "C", ALT = "T", FILTER = 10)
DNMT3A.vaf.4$bc <- paste0("AML8007.4_", DNMT3A.vaf.4$bc, "-1")
DNMT3A.vaf.4$sample <- "AML8007.4"
DNMT3A.vaf <- as.data.frame(rbind(DNMT3A.vaf.1, DNMT3A.vaf.4))
rownames(DNMT3A.vaf) <- DNMT3A.vaf$bc
write.table(file = "./data/AML/mutations/AML8007_DNMT3A.csv", x = DNMT3A.vaf, quote = F, sep = "\t")

p <- DimPlot(subset(AML8007, orig.ident == "AML8007.1"),
  reduction = "ref.umap",
  cells.highlight = list(
    "mutated" = DNMT3A.vaf$bc[which(DNMT3A.vaf$mutated == "mutated")],
    "wildtype" = DNMT3A.vaf$bc[which(DNMT3A.vaf$mutated == "wildtype")]
  ),
  cols.highlight = c("wildtype" = "black", "mutated" = "firebrick"), sizes.highlight = 1
) +
  NoLegend() +
  NoAxes() +
  theme(plot.title = element_blank())
q <- DimPlot(subset(AML8007, orig.ident == "AML8007.4"),
  reduction = "ref.umap",
  cells.highlight = list(
    "mutated" = DNMT3A.vaf$bc[which(DNMT3A.vaf$mutated == "mutated")],
    "wildtype" = DNMT3A.vaf$bc[which(DNMT3A.vaf$mutated == "wildtype")]
  ),
  cols.highlight = c("wildtype" = "black", "mutated" = "firebrick"), sizes.highlight = 1
) +
  NoLegend() +
  NoAxes() +
  theme(plot.title = element_blank())
cowplot::plot_grid(plotlist = list(p, q))
ggsave("./AML/figures/AML8007/UMAP/AML8007_1_DNMT3A.png", width = 3, height = 3, dpi = 600, plot = p)
ggsave("./AML/figures/AML8007/UMAP/AML8007_4_DNMT3A.png", width = 3, height = 3, dpi = 600, plot = q)

# TP53_1
TP53_1.vaf.1 <- extract_mutation(BC.data.file = "./data/AML/pileup//AML8007_1_pileup_TP53_1.csv.gz", REF = "A", ALT = "T", downsample = 500000, FILTER = 10)
TP53_1.vaf.1$bc <- paste0("AML8007.1_", TP53_1.vaf.1$bc, "-1")
TP53_1.vaf.1$sample <- "AML8007.1"
TP53_1.vaf.4 <- extract_mutation(BC.data.file = "./data/AML/pileup/AML8007_4_pileup_TP53_1.csv.gz", REF = "A", ALT = "T", downsample = 500000, FILTER = 10)
TP53_1.vaf.4$bc <- paste0("AML8007.4_", TP53_1.vaf.4$bc, "-1")
TP53_1.vaf.4$sample <- "AML8007.4"
TP53_1.vaf <- as.data.frame(rbind(TP53_1.vaf.1, TP53_1.vaf.4))
rownames(TP53_1.vaf) <- TP53_1.vaf$bc
write.table(file = "./data/AML/mutations/AML8007_TP53.1.csv", x = TP53_1.vaf, quote = F, sep = "\t")

p <- DimPlot(subset(AML8007, orig.ident == "AML8007.1"),
  reduction = "ref.umap",
  cells.highlight = list(
    "mutated" = TP53_1.vaf$bc[which(TP53_1.vaf$mutated == "mutated")],
    "wildtype" = TP53_1.vaf$bc[which(TP53_1.vaf$mutated == "wildtype")]
  ),
  cols.highlight = c("wildtype" = "black", "mutated" = "firebrick"), sizes.highlight = 1
) +
  NoLegend() +
  NoAxes() +
  theme(plot.title = element_blank())
q <- DimPlot(subset(AML8007, orig.ident == "AML8007.4"),
  reduction = "ref.umap",
  cells.highlight = list(
    "mutated" = TP53_1.vaf$bc[which(TP53_1.vaf$mutated == "mutated")],
    "wildtype" = TP53_1.vaf$bc[which(TP53_1.vaf$mutated == "wildtype")]
  ),
  cols.highlight = c("wildtype" = "black", "mutated" = "firebrick"), sizes.highlight = 1
) +
  NoLegend() +
  NoAxes() +
  theme(plot.title = element_blank())
cowplot::plot_grid(plotlist = list(p, q))
ggsave("./AML/figures/AML8007/UMAP/AML8007_1_TP53_1.png", width = 3, height = 3, dpi = 600, plot = p)
ggsave("./AML/figures/AML8007/UMAP/AML8007_4_TP53_1.png", width = 3, height = 3, dpi = 600, plot = q)


# TP53_2
TP53_2.vaf.1 <- extract_mutation(BC.data.file = "./data/AML/pileup/AML8007_1_pileup_TP53_2.csv.gz", REF = "G", ALT = "A", downsample = 500000, FILTER = 10)
TP53_2.vaf.1$bc <- paste0("AML8007.1_", TP53_2.vaf.1$bc, "-1")
TP53_2.vaf.1$sample <- "AML8007.1"
TP53_2.vaf.4 <- extract_mutation(BC.data.file = "./data/AML/pileup/AML8007_4_pileup_TP53_2.csv.gz", REF = "G", ALT = "A", downsample = 500000, FILTER = 10)
TP53_2.vaf.4$bc <- paste0("AML8007.4_", TP53_2.vaf.4$bc, "-1")
TP53_2.vaf.4$sample <- "AML8007.4"
TP53_2.vaf <- as.data.frame(rbind(TP53_2.vaf.1, TP53_2.vaf.4))
rownames(TP53_2.vaf) <- TP53_2.vaf$bc
write.table(file = "./data/AML/mutations/AML8007_TP53.2.csv", x = TP53_2.vaf, quote = F, sep = "\t")

p <- DimPlot(subset(AML8007, orig.ident == "AML8007.1"),
  reduction = "ref.umap",
  cells.highlight = list(
    "mutated" = TP53_2.vaf$bc[which(TP53_2.vaf$mutated == "mutated")],
    "wildtype" = TP53_2.vaf$bc[which(TP53_2.vaf$mutated == "wildtype")]
  ),
  cols.highlight = c("wildtype" = "black", "mutated" = "firebrick"), sizes.highlight = 1
) +
  NoLegend() +
  NoAxes() +
  theme(plot.title = element_blank())
q <- DimPlot(subset(AML8007, orig.ident == "AML8007.4"),
  reduction = "ref.umap",
  cells.highlight = list(
    "mutated" = TP53_2.vaf$bc[which(TP53_2.vaf$mutated == "mutated")],
    "wildtype" = TP53_2.vaf$bc[which(TP53_2.vaf$mutated == "wildtype")]
  ),
  cols.highlight = c("wildtype" = "black", "mutated" = "firebrick"), sizes.highlight = 1
) +
  NoLegend() +
  NoAxes() +
  theme(plot.title = element_blank())
cowplot::plot_grid(plotlist = list(p, q))
ggsave("./AML/figures/AML8007/UMAP/AML8007_1_TP53_2.png", width = 3, height = 3, dpi = 600, plot = p)
ggsave("./AML/figures/AML8007/UMAP/AML8007_4_TP53_2.png", width = 3, height = 3, dpi = 600, plot = q)
