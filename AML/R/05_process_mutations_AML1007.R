# process mutations for AML1007.135

library(ComplexHeatmap)
library(dplyr)
library(nanoranger.R)
library(ggplot2)
library(Seurat)
library(SummarizedExperiment)

AML1007 <- readRDS("./data/AML/objects/20220825_AML1007.rds")
Idents(AML1007) <- "predicted.celltype"

AML1007.genotype <- read.table("./data/AML/souporcell/20210707_AML1007_chimerism.csv")
rownames(AML1007.genotype) <- AML1007.genotype$barcode
AML1007.genotype$predicted.celltype <- AML1007$predicted.celltype[AML1007.genotype$barcode]

boo <- as.data.frame(prop.table(table(Idents(AML1007), AML1007$orig.ident), margin = 2))
boo$Var1 <- factor(boo$Var1, levels = rev(names(AML.combined.colors)))
ggplot(boo[which(boo$Var2 %in% c("AML1007.1", "AML1007.3")), ], aes(x = Var2, y = 100 * Freq, fill = Var1)) +
  geom_col() +
  geom_hline(yintercept = 5) +
  scale_x_discrete(labels = c("Screening", "CR")) +
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
ggsave("./AML/figures/AML1007/plots/20220912_celltype_kinetics.svg", width = 1.1, height = 2)

# ASXL1 mutation
ASXL1.vaf.1 <- extract_indel(BC.data.file = "./data/AML/pileup/AML1007_1_pileup_ASXL1.csv.gz", REF = "A", CONSENSUS = -2, FILTER = 20)
ASXL1.vaf.1$bc <- paste0("AML1007.1_", ASXL1.vaf.1$bc, "-1")
ASXL1.vaf.1$sample <- "AML1007.1"
ASXL1.vaf.3 <- extract_indel(BC.data.file = "./data/AML/pileup/AML1007_3_pileup_ASXL1.csv.gz", REF = "A", CONSENSUS = -2, FILTER = 20)
ASXL1.vaf.3$bc <- paste0("AML1007.3_", ASXL1.vaf.3$bc, "-1")
ASXL1.vaf.3$sample <- "AML1007.3"
ASXL1.vaf <- as.data.frame(rbind(ASXL1.vaf.1, ASXL1.vaf.3))
rownames(ASXL1.vaf) <- ASXL1.vaf$bc
write.table(file = "./data/AML/mutations/AML1007_ASXL1.csv", x = ASXL1.vaf, quote = F, sep = "\t")

p <- DimPlot(subset(AML1007, orig.ident == "AML1007.1"),
  reduction = "ref.umap",
  cells.highlight = list(
    "mutated" = ASXL1.vaf$bc[which(ASXL1.vaf$mutated == "mutated")],
    "wildtype" = ASXL1.vaf$bc[which(ASXL1.vaf$mutated == "wildtype")]
  ),
  cols.highlight = c("wildtype" = "black", "mutated" = "firebrick"), sizes.highlight = 1
) +
  NoLegend() +
  NoAxes() +
  theme(plot.title = element_blank())
q <- DimPlot(subset(AML1007, orig.ident == "AML1007.3"),
  reduction = "ref.umap",
  cells.highlight = list(
    "mutated" = ASXL1.vaf$bc[which(ASXL1.vaf$mutated == "mutated")],
    "wildtype" = ASXL1.vaf$bc[which(ASXL1.vaf$mutated == "wildtype")]
  ),
  cols.highlight = c("wildtype" = "black", "mutated" = "firebrick"), sizes.highlight = 1
) +
  NoLegend() +
  NoAxes() +
  theme(plot.title = element_blank())
cowplot::plot_grid(plotlist = list(p, q))
ggsave("./AML/figures/AML1007/UMAP/AML1007_1_ASXL1.png", width = 3, height = 3, dpi = 600, plot = p)
ggsave("./AML/figures/AML1007/UMAP/AML1007_3_ASXL1.png", width = 3, height = 3, dpi = 600, plot = q)

# IDH2 mutation
IDH2.vaf.1 <- extract_mutation(BC.data.file = "./data/AML/pileup/AML1007_1_pileup_IDH2.csv.gz", REF = "C", ALT = "T", FILTER = 20)
IDH2.vaf.1$bc <- paste0("AML1007.1_", IDH2.vaf.1$bc, "-1")
IDH2.vaf.1$sample <- "AML1007.1"
IDH2.vaf.3 <- extract_mutation(BC.data.file = "./data/AML/pileup/AML1007_3_pileup_IDH2.csv.gz", REF = "C", ALT = "T", FILTER = 20)
IDH2.vaf.3$bc <- paste0("AML1007.3_", IDH2.vaf.3$bc, "-1")
IDH2.vaf.3$sample <- "AML1007.3"
IDH2.vaf <- as.data.frame(rbind(IDH2.vaf.1, IDH2.vaf.3))
rownames(IDH2.vaf) <- IDH2.vaf$bc
write.table(file = "./data/AML/mutations/AML1007_IDH2.csv", x = IDH2.vaf, quote = F, sep = "\t")

p <- DimPlot(subset(AML1007, orig.ident == "AML1007.1"),
  reduction = "ref.umap",
  cells.highlight = list(
    "mutated" = IDH2.vaf$bc[which(IDH2.vaf$mutated == "mutated")],
    "wildtype" = IDH2.vaf$bc[which(IDH2.vaf$mutated == "wildtype")]
  ),
  cols.highlight = c("wildtype" = "black", "mutated" = "firebrick"), sizes.highlight = 1
) +
  NoLegend() +
  NoAxes() +
  theme(plot.title = element_blank())
q <- DimPlot(subset(AML1007, orig.ident == "AML1007.3"),
  reduction = "ref.umap",
  cells.highlight = list(
    "mutated" = IDH2.vaf$bc[which(IDH2.vaf$mutated == "mutated")],
    "wildtype" = IDH2.vaf$bc[which(IDH2.vaf$mutated == "wildtype")]
  ),
  cols.highlight = c("wildtype" = "black", "mutated" = "firebrick"), sizes.highlight = 1
) +
  NoLegend() +
  NoAxes() +
  theme(plot.title = element_blank())
cowplot::plot_grid(plotlist = list(p, q))
ggsave("./AML/figures/AML1007/UMAP/AML1007_1_IDH2.png", width = 3, height = 3, dpi = 600, plot = p)
ggsave("./AML/figures/AML1007/UMAP/AML1007_3_IDH2.png", width = 3, height = 3, dpi = 600, plot = q)

# STAG2 mutation
STAG2.vaf.1 <- extract_indel(BC.data.file = "./data/AML/pileup/AML1007_1_pileup_STAG22.csv.gz", REF = "C", CONSENSUS = -2, FILTER = 20)
STAG2.vaf.1$bc <- paste0("AML1007.1_", STAG2.vaf.1$bc, "-1")
STAG2.vaf.1$sample <- "AML1007.1"
STAG2.vaf.3 <- extract_indel(BC.data.file = "./data/AML/pileup/AML1007_3_pileup_STAG22.csv.gz", REF = "C", CONSENSUS = -2, FILTER = 20)
STAG2.vaf.3$bc <- paste0("AML1007.3_", STAG2.vaf.3$bc, "-1")
STAG2.vaf.3$sample <- "AML1007.3"
STAG2.vaf <- as.data.frame(rbind(STAG2.vaf.1, STAG2.vaf.3))
rownames(STAG2.vaf) <- STAG2.vaf$bc
write.table(file = "./data/AML/mutations/AML1007_STAG2.csv", x = STAG2.vaf, quote = F, sep = "\t")

p <- DimPlot(subset(AML1007, orig.ident == "AML1007.1"),
  reduction = "ref.umap",
  cells.highlight = list(
    "mutated" = STAG2.vaf$bc[which(STAG2.vaf$mutated == "mutated")],
    "wildtype" = STAG2.vaf$bc[which(STAG2.vaf$mutated == "wildtype")]
  ),
  cols.highlight = c("wildtype" = "black", "mutated" = "firebrick"), sizes.highlight = 1
) +
  NoLegend() +
  NoAxes() +
  theme(plot.title = element_blank())
q <- DimPlot(subset(AML1007, orig.ident == "AML1007.3"),
  reduction = "ref.umap",
  cells.highlight = list(
    "mutated" = STAG2.vaf$bc[which(STAG2.vaf$mutated == "mutated")],
    "wildtype" = STAG2.vaf$bc[which(STAG2.vaf$mutated == "wildtype")]
  ),
  cols.highlight = c("wildtype" = "black", "mutated" = "firebrick"), sizes.highlight = 1
) +
  NoLegend() +
  NoAxes() +
  theme(plot.title = element_blank())
cowplot::plot_grid(plotlist = list(p, q))
ggsave("./AML/figures/AML1007/UMAP/AML1007_1_STAG2.png", width = 3, height = 3, dpi = 600, plot = p)
ggsave("./AML/figures/AML1007/UMAP/AML1007_3_STAG2.png", width = 3, height = 3, dpi = 600, plot = q)

# SRSF2 mutation

SRSF2.1 <- extract_mutation("./data/AML/pileup/AML1007_1_pileup_SRSF2.csv.gz",
  ALT = "A", REF = "G", FILTER = 30
)
SRSF2.1$bc <- paste0("AML1007.1_", SRSF2.1$bc, "-1")
SRSF2.1$sample <- "AML1007.1"

SRSF2.3 <- extract_mutation("./data/AML/pileup/AML1007_3_pileup_SRSF2.csv.gz",
  ALT = "A", REF = "G", FILTER = 30
)
SRSF2.3$bc <- paste0("AML1007.3_", SRSF2.3$bc, "-1")
SRSF2.3$sample <- "AML1007.3"

SRSF2.5 <- extract_mutation("./data/AML/pileup/AML1007_5_pileup_SRSF2.csv.gz",
  ALT = "A", REF = "G", FILTER = 30
)
SRSF2.5$bc <- paste0("AML1007.5_", SRSF2.5$bc, "-1")
SRSF2.5$sample <- "AML1007.5"

SRSF2.vaf <- rbind(SRSF2.1, SRSF2.3)
SRSF2.vaf <- rbind(SRSF2.vaf, SRSF2.5) %>% as.data.frame()
rownames(SRSF2.vaf) <- SRSF2$bc
write.table(SRSF2.vaf, file = "./data/AML/mutations/AML1007_SRSF2.csv", quote = F, sep = "\t")


p <- DimPlot(subset(AML1007, orig.ident == "AML1007.1"),
  reduction = "ref.umap",
  cells.highlight = list(
    "mutated" = SRSF2.vaf$bc[which(SRSF2.vaf$mutated == "mutated")],
    "wildtype" = SRSF2.vaf$bc[which(SRSF2.vaf$mutated == "wildtype")]
  ),
  cols.highlight = c("wildtype" = "black", "mutated" = "firebrick"), sizes.highlight = 1
) +
  NoLegend() +
  NoAxes() +
  theme(plot.title = element_blank())
q <- DimPlot(subset(AML1007, orig.ident == "AML1007.3"),
  reduction = "ref.umap",
  cells.highlight = list(
    "mutated" = SRSF2.vaf$bc[which(SRSF2.vaf$mutated == "mutated")],
    "wildtype" = SRSF2.vaf$bc[which(SRSF2.vaf$mutated == "wildtype")]
  ),
  cols.highlight = c("wildtype" = "black", "mutated" = "firebrick"), sizes.highlight = 1
) +
  NoLegend() +
  NoAxes() +
  theme(plot.title = element_blank())
r <- DimPlot(subset(AML1007, orig.ident == "AML1007.5"),
  reduction = "ref.umap",
  cells.highlight = list(
    "mutated" = SRSF2.vaf$bc[which(SRSF2.vaf$mutated == "mutated")],
    "wildtype" = SRSF2.vaf$bc[which(SRSF2.vaf$mutated == "wildtype")]
  ),
  cols.highlight = c("wildtype" = "black", "mutated" = "firebrick"), sizes.highlight = 1
) +
  NoLegend() +
  NoAxes() +
  theme(plot.title = element_blank())
cowplot::plot_grid(plotlist = list(p, q, r))
ggsave("./AML/figures/AML1007/UMAP/AML1007_1_SRSF2.png", width = 3, height = 3, dpi = 600, plot = p)
ggsave("./AML/figures/AML1007/UMAP/AML1007_3_SRSF2.png", width = 3, height = 3, dpi = 600, plot = q)
ggsave("./AML/figures/AML1007/UMAP/AML1007_5_SRSF2.png", width = 3, height = 3, dpi = 600, plot = r)
