# analysis of counts Illumina vs. ONT

library(dplyr)
library(ggplot2)
library(nanoranger.R)
library(Seurat)

AML1022.1 <- readRDS("./data/AML/objects/20220729_AML1022_1.rds")
AML1022.genotype <- read.table("./data/AML/souporcell/20210707_AML1022_chimerism.csv")

###
library(stringr)
DNMT3A.reads <- as.data.frame(data.table::fread("./data/AML/pileup/AML1022_1_DNMT3A.csv.gz") %>% group_by(bc) %>% tally())
DNMT3A.reads$bc <- paste0("AML1022.1_", DNMT3A.reads$bc, "-1")
rownames(DNMT3A.reads) <- DNMT3A.reads$bc
RUNX1.reads <- as.data.frame(data.table::fread("./data/AML/pileup/AML1022_1_RUNX1.csv.gz") %>% group_by(bc) %>% tally())
RUNX1.reads$bc <- paste0("AML1022.1_", RUNX1.reads$bc, "-1")
rownames(RUNX1.reads) <- RUNX1.reads$bc
SF3B1.reads <- as.data.frame(data.table::fread("./data/AML/pileup/AML1022_1_SF3B1.csv.gz") %>% group_by(bc) %>% tally())
SF3B1.reads$bc <- paste0("AML1022.1_", SF3B1.reads$bc, "-1")
rownames(SF3B1.reads) <- SF3B1.reads$bc

DNMT3A.reads$nCount_RNA <- AML1022.1$nCount_RNA[DNMT3A.reads$bc]
DNMT3A.reads$nFeature_RNA <- AML1022.1$nFeature_RNA[DNMT3A.reads$bc]
DNMT3A.reads$predicted.celltype <- AML1022.1$predicted.celltype[DNMT3A.reads$bc]
DNMT3A.reads[intersect(DNMT3A.reads$bc, colnames(AML1022.1)), "Illumina.reads"] <-
  GetAssayData(AML1022.1, slot = "counts", assay = "RNA")["DNMT3A", intersect(DNMT3A.reads$bc, colnames(AML1022.1))]
DNMT3A.reads$Illumina.reads[is.na(DNMT3A.reads$Illumina.reads)] <- 0
DNMT3A.reads$Illumina.reads <- jitter(DNMT3A.reads$Illumina.reads, amount = 0.2)
# DNMT3A.reads$Illumina.reads[which(DNMT3A.reads$Illumina.reads == 0 | is.na(DNMT3A.reads$Illumina.reads))] =
#  jitter(rep(0.1, length(which(DNMT3A.reads$Illumina.reads == 0 | is.na(DNMT3A.reads$Illumina.reads)))), amount = 0.3)

RUNX1.reads$nCount_RNA <- AML1022.1$nCount_RNA[RUNX1.reads$bc]
RUNX1.reads$nFeature_RNA <- AML1022.1$nFeature_RNA[RUNX1.reads$bc]
RUNX1.reads$predicted.celltype <- AML1022.1$predicted.celltype[RUNX1.reads$bc]
RUNX1.reads[intersect(RUNX1.reads$bc, colnames(AML1022.1)), "Illumina.reads"] <-
  GetAssayData(AML1022.1, slot = "counts", assay = "RNA")["RUNX1", intersect(RUNX1.reads$bc, colnames(AML1022.1))]
RUNX1.reads$Illumina.reads[is.na(RUNX1.reads$Illumina.reads)] <- 0
RUNX1.reads$Illumina.reads <- jitter(RUNX1.reads$Illumina.reads, amount = 0.2)


SF3B1.reads$nCount_RNA <- AML1022.1$nCount_RNA[SF3B1.reads$bc]
SF3B1.reads$nFeature_RNA <- AML1022.1$nFeature_RNA[SF3B1.reads$bc]
SF3B1.reads$predicted.celltype <- AML1022.1$predicted.celltype[SF3B1.reads$bc]
SF3B1.reads[intersect(SF3B1.reads$bc, colnames(AML1022.1)), "Illumina.reads"] <-
  GetAssayData(AML1022.1, slot = "counts", assay = "RNA")["SF3B1", intersect(SF3B1.reads$bc, colnames(AML1022.1))]
SF3B1.reads$Illumina.reads[is.na(SF3B1.reads$Illumina.reads)] <- 0
SF3B1.reads$Illumina.reads <- jitter(SF3B1.reads$Illumina.reads, amount = 0.2)

ggplot(DNMT3A.reads, aes(x = Illumina.reads, y = n)) +
  ggrastr::rasterize(geom_point(aes(color = predicted.celltype), size = 0.5), dpi = 600) +
  scale_color_manual(values = AML.combined.colors) +
  scale_x_continuous("Illumina reads", limits = c(-0.5, 25)) +
  scale_y_log10("ONT reads") +
  theme_classic() +
  ggtitle("DNMT3A") +
  theme(
    legend.position = "none",
    plot.title = element_text("Arial", size = 10, color = "black", face = "bold.italic", hjust = 0.5),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./AML/figures/AML1022/plots/20220821_AML1022_1_DNMT3A_counts.svg", width = 2, height = 1.5)

ggplot(RUNX1.reads, aes(x = Illumina.reads, y = n)) +
  ggrastr::rasterize(geom_point(aes(color = predicted.celltype), size = 0.5), dpi = 600) +
  scale_color_manual(values = AML.combined.colors) +
  scale_x_continuous("Illumina reads", limits = c(-0.5, 25)) +
  scale_y_log10("ONT reads") +
  theme_classic() +
  ggtitle("RUNX1") +
  theme(
    legend.position = "none",
    plot.title = element_text("Arial", size = 10, color = "black", face = "bold.italic", hjust = 0.5),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./AML/figures/AML1022/plots/20220821_AML1022_1_RUNX1_counts.svg", width = 2, height = 1.5)

ggplot(SF3B1.reads, aes(x = Illumina.reads, y = n)) +
  ggrastr::rasterize(geom_point(aes(color = predicted.celltype), size = 0.5), dpi = 600) +
  scale_color_manual(values = AML.combined.colors) +
  scale_x_continuous("Illumina reads", limits = c(-0.5, 25)) +
  scale_y_log10("ONT reads") +
  theme_classic() +
  ggtitle("SF3B1") +
  theme(
    legend.position = "none",
    plot.title = element_text("Arial", size = 10, color = "black", face = "bold.italic", hjust = 0.5),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./AML/figures/AML1022/plots/20220821_AML1022_1_SF3B1_counts.svg", width = 2, height = 1.5)


ggplot(RUNX1.reads, aes(x = nCount_RNA, y = n)) +
  ggrastr::rasterize(geom_point(aes(color = predicted.celltype), size = 0.5), dpi = 600) +
  scale_color_manual(values = AML.combined.colors) +
  scale_x_log10() +
  scale_y_log10("Reads RUNX1") +
  geom_smooth(color = "black") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./AML/figures/AML1022/plots/20220730_AML1022_1_RUNX1_nCount_RNA.svg", width = 2.5, height = 2)

ggplot(SF3B1.reads, aes(x = nCount_RNA, y = n)) +
  ggrastr::rasterize(geom_point(aes(color = predicted.celltype), size = 0.5), dpi = 600) +
  scale_color_manual(values = AML.combined.colors) +
  scale_x_log10() +
  scale_y_log10("Reads SF3B1") +
  geom_smooth(color = "black") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./AML/figures/AML1022/plots/20220730_AML1022_1_SF3B1_nCount_RNA.svg", width = 2.5, height = 2)


ggplot(DNMT3A.reads, aes(x = nCount_RNA, y = n)) +
  ggrastr::rasterize(geom_point(aes(color = predicted.celltype), size = 0.5), dpi = 600) +
  scale_color_manual(values = AML.combined.colors) +
  scale_x_log10() +
  scale_y_log10("Reads DNMT3A") +
  geom_smooth(color = "black") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./AML/figures/AML1022/plots/20220730_AML1022_1_DNMT3A_nCount_RNA.svg", width = 2.5, height = 2)

ggplot(RUNX1.reads, aes(x = nCount_RNA, y = n)) +
  ggrastr::rasterize(geom_point(aes(color = predicted.celltype), size = 0.5), dpi = 600) +
  scale_color_manual(values = AML.combined.colors) +
  scale_x_log10() +
  scale_y_log10("Reads RUNX1") +
  geom_smooth(color = "black") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./AML/figures/AML1022/plots/20220730_AML1022_1_RUNX1_nCount_RNA.svg", width = 2.5, height = 2)

ggplot(SF3B1.reads, aes(x = nCount_RNA, y = n)) +
  ggrastr::rasterize(geom_point(aes(color = predicted.celltype), size = 0.5), dpi = 600) +
  scale_color_manual(values = AML.combined.colors) +
  scale_x_log10() +
  scale_y_log10("Reads SF3B1") +
  geom_smooth(color = "black") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./AML/figures/AML1022/plots/20220730_AML1022_1_SF3B1_nCount_RNA.svg", width = 2.5, height = 2)
