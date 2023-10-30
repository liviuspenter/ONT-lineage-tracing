# analysis of mtDNA data: donor-recipient deconvolution

library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(SummarizedExperiment)
library(dplyr)

# functionality from MAESTER (https://github.com/petervangalen/MAESTER-2021)
source("./R/210215_FunctionsGeneral.R")
source("./R/celltypes.R")

AML.all <- readRDS("./data/AML/objects/AML.all.rds")
AML1026.2 <- readRDS("./data/AML/objects/20220711_AML1026_2.rds")

# get mtDNA variants
magtk.output <- readRDS("./data/AML/mtDNA/AML1026_2.maegatk/AML1026_2.rds")
af.dm <- data.matrix(computeAFMutMatrix(magtk.output)) * 100

rownames(af.dm) <- gsub(rownames(af.dm), pattern = "_", replacement = "")
af.dm.clean <- af.dm[-which(rowSums(af.dm) == 0), ]
# af.dm.clean = af.dm[, coverage$V1[which(coverage$V2 > 5)]]

# identify potentially informative mtDNA mutations:
# - presence in enough cells
# - variance high enough
# - mean variant allele frequency high enough

qc.df <- data.frame(
  presence = apply(af.dm.clean, 1, FUN = function(x) {
    length(which(x != 0))
  }),
  variance = rowVars(af.dm.clean),
  vmr = apply(af.dm.clean, 1, FUN = function(x) {
    litteR::iod(x)
  }),
  mean.vaf = apply(af.dm.clean, 1, FUN = function(x) {
    mean(x[which(x != 0)])
  })
)
qc.df$annotation <- ifelse((qc.df$mean.vaf > 80 & -log10(qc.df$vmr) > -1.75), "germline", "non.clonal")
qc.df$annotation[which(qc.df$annotation != "germline" & qc.df$presence > 50)] <- "clonal"

ggplot() +
  geom_point(data = qc.df, aes(x = mean.vaf, y = -log10(vmr), color = annotation)) +
  scale_color_manual(values = c("clonal" = "firebrick", "non.clonal" = "grey", "germline" = "blue")) +
  theme_classic() +
  theme(legend.position = "none")

# find germline mutations for donor-recipient deconvolution
germline.mutations <- which(qc.df$annotation %in% c("germline", "clonal") & qc.df$mean.vaf > 50)

cells.sample <- sample(ncol(af.dm.clean), size = 500)
Heatmap(af.dm.clean[germline.mutations, cells.sample],
  show_row_names = T, show_column_names = F, cluster_rows = T, cluster_columns = T,
  show_column_dend = F, show_row_dend = F, row_names_gp = gpar(fontsize = 5), row_names_side = "left"
)

recipient.variants <- c("8433T>C", "9722T>C", "15833C>T", "4011C>T")
donor.variants <- c("14905G>A", "14766C>T", "15452C>A", "15607A>G")

col_fun <- circlize::colorRamp2(breaks = seq(0, 100, 100 / 8), colors = BuenColors::jdb_palette(name = "solar_rojos"))
svglite::svglite("./AML/figures/mtDNA/heatmaps/20220803_AML1026_donor_recipient.svg", width = 4, height = 2)
Heatmap(af.dm.clean[c(donor.variants, recipient.variants), ],
  show_row_names = T, show_column_names = F, cluster_rows = T, cluster_columns = T,
  show_column_dend = F, show_row_dend = F, row_names_gp = gpar(fontsize = 5), row_names_side = "left", column_split = 2, border = T,
  col = col_fun, use_raster = T, raster_quality = 10
)
dev.off()

deconvolution.df <- data.frame(
  donor = colMeans(af.dm.clean[donor.variants, ]),
  recipient = colMeans(af.dm.clean[recipient.variants, ])
)
rownames(deconvolution.df) <- paste0("AML1026.2_", rownames(deconvolution.df), "-1")
deconvolution.df$annotation <- ifelse(deconvolution.df$donor > deconvolution.df$recipient, "donor", "recipient")
deconvolution.df$annotation[which(deconvolution.df$donor == deconvolution.df$recipient)] <- "none"
donor <- rownames(deconvolution.df)[which(deconvolution.df$annotation == "donor")]
recipient <- rownames(deconvolution.df)[which(deconvolution.df$annotation == "recipient")]

AML1026.genotype <- read.table("./data/AML/souporcell/20210707_AML1026_chimerism.csv")
rownames(AML1026.genotype) <- AML1026.genotype$barcode
AML1026.genotype$predicted.celltype <- AML1026.2$predicted.celltype[AML1026.genotype$barcode]
deconvolution.df$chimerism <- AML1026.genotype[rownames(deconvolution.df), "chimerism"]

ggplot(deconvolution.df, aes(x = donor, y = recipient)) +
  geom_hline(yintercept = c(25, 75)) +
  geom_vline(xintercept = c(25, 75)) +
  geom_point(aes(color = annotation), size = 0.5) +
  scale_color_manual(values = c("none" = "grey", "recipient" = "purple", "donor" = "orange")) +
  scale_x_continuous("% donor") +
  scale_y_continuous("% recipient") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./AML/figures/mtDNA/plots/20220803_AML1026_deconvolution_donor_recipient.svg", width = 2, height = 2)

ggplot(deconvolution.df, aes(x = donor, y = recipient)) +
  geom_abline(slope = 1) +
  geom_point(aes(color = chimerism), size = 0.5) +
  scale_color_manual(values = c("unassigned" = "grey", "recipient" = "purple", "donor" = "orange")) +
  scale_x_continuous("% donor") +
  scale_y_continuous("% recipient") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./AML/figures/mtDNA/plots/20220803_AML1026_deconvolution_donor_recipient_souporcell.svg", width = 2, height = 2)
