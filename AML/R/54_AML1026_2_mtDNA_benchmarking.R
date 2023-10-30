### benchmarking versus mtscATAC-seq data

library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(SummarizedExperiment)
library(dplyr)

# functionality from MAESTER (https://github.com/petervangalen/MAESTER-2021)
source("./R/210215_FunctionsGeneral.R")

# get mtDNA variants
magtk.output <- readRDS("./data/AML/mtDNA/AML1026_2.maegatk/AML1026_2.rds")
af.dm <- data.matrix(computeAFMutMatrix(magtk.output)) * 100
rownames(af.dm) <- gsub(rownames(af.dm), pattern = "_", replacement = "")
af.dm.clean <- af.dm[-which(rowSums(af.dm) == 0), ]

# ASAP-seq data
asap.cells.by.sample <- data.table::fread("./data/AML/mtDNA/20210611_cells_by_sample.csv")
AML1026.combined.frequencies <- readRDS("./data/AML/mtDNA/20210429_AML1026_combined_mutation_frequencies.rds")
AML1026.combined.frequencies <-
  AML1026.combined.frequencies[, intersect(
    colnames(AML1026.combined.frequencies),
    asap.cells.by.sample$barcode[which(asap.cells.by.sample$sample == "AML1026_2")]
  )]

benchmarking.df <- data.frame(scRNAseq = rowMeans(af.dm.clean))
benchmarking.df$mtscATAC <- NA
benchmarking.df[rownames(AML1026.combined.frequencies), "mtscATAC"] <- rowMeans(AML1026.combined.frequencies)
benchmarking.df[is.na(benchmarking.df)] <- 10^-5

ggplot(benchmarking.df, aes(x = mtscATAC, y = scRNAseq)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10()

magtk.output.scatac <- readRDS("./data/AML/mtDNA/AML1026_12.mgatk/AML1026_12.rds")
af.dm.scatac <- data.matrix(computeAFMutMatrix(magtk.output.scatac)) * 100
rownames(af.dm.scatac) <- gsub(rownames(af.dm.scatac), pattern = "_", replacement = "")
af.dm.clean.scatac <- af.dm.scatac[-which(rowSums(af.dm.scatac) == 0), ]
af.dm.clean.scatac <- af.dm.clean.scatac[, intersect(
  colnames(af.dm.clean.scatac),
  asap.cells.by.sample$barcode[which(asap.cells.by.sample$sample == "AML1026_2")]
)]

benchmarking.df <- data.frame(scRNAseq = rowMeans(af.dm.clean))
benchmarking.df$mtscATAC <- NA
benchmarking.df[rownames(af.dm.scatac), "mtscATAC"] <- rowMeans(af.dm.scatac)
benchmarking.df[is.na(benchmarking.df)] <- 10^-5 + abs(jitter(rep(10^-5, length(benchmarking.df[is.na(benchmarking.df)])), amount = 10^-5))
benchmarking.df[benchmarking.df == 0] <- 10^-5 + abs(jitter(rep(10^-5, length(benchmarking.df[benchmarking.df == 0])), amount = 10^-5))

p <- ggplot() +
  geom_abline(slope = 1) +
  ggrastr::rasterize(geom_point(
    data = benchmarking.df[which(!rownames(benchmarking.df) %in% rownames(AML1026.combined.frequencies)), ],
    aes(x = mtscATAC, y = scRNAseq), size = 0.5, color = "grey"
  ), dpi = 600) +
  ggrastr::rasterize(geom_point(
    data = benchmarking.df[which(rownames(benchmarking.df) %in% rownames(AML1026.combined.frequencies)), ],
    aes(x = mtscATAC, y = scRNAseq), size = 0.5, color = "firebrick"
  ), dpi = 600) +
  scale_x_log10("mtscATAC-seq", limits = c(10^-5, 100), breaks = c(2 * 10^-5, 10^-4, 10^-2, 1, 100), labels = c("ND", "0.0001", "0.01", "1", "100")) +
  scale_y_log10("ONT", limits = c(10^-5, 100), breaks = c(2 * 10^-5, 10^-4, 10^-2, 1, 100), labels = c("ND", "0.0001", "0.01", "1", "100")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./AML/figures/mtDNA/plots/20220803_AML1026_scATACseq_ONT.svg", width = 2.5, height = 2.5, plot = p)
