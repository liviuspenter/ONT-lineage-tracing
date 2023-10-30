# visualize coverage per celltype

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(nanoranger.R)

source("./R/celltypes.R")
AML1022 <- readRDS("./data/benchmarking/objects/20230812_AML1022_nanoranger_GoT.rds")

DNMT3A.GoT <- read.table(file = "./data/benchmarking/GoT_nanoranger/AML1022_pileup_DNMT3A.GoT.ONT.processed.csv")
DNMT3A.GoT$condition <- "GoT"
DNMT3A.GoT.Illumina <- read.table(file = "./data/benchmarking/GoT_nanoranger/AML1022_pileup_DNMT3A.GoT.Illumina.processed.csv")
DNMT3A.GoT.Illumina$condition <- "GoT.Illumina"
DNMT3A.normal <- read.table(file = "./data/benchmarking/GoT_nanoranger/AML1022_pileup_DNMT3A.normal.processed.csv")
DNMT3A.normal$condition <- "normal"
DNMT3A.spikein <- read.table(file = "./data/benchmarking/GoT_nanoranger/AML1022_pileup_DNMT3A.spikein.processed.csv")
DNMT3A.spikein$condition <- "spikein"
DNMT3A.data <- dplyr::bind_rows(DNMT3A.GoT, DNMT3A.GoT.Illumina, DNMT3A.normal, DNMT3A.spikein) %>%
  filter(bc %in% colnames(AML1022))
DNMT3A.data$predicted.celltype <- AML1022$predicted.celltype[DNMT3A.data$bc]

RUNX1.GoT <- read.table(file = "./data/benchmarking/GoT_nanoranger/AML1022_pileup_RUNX1.GoT.ONT.processed.csv")
RUNX1.GoT$condition <- "GoT"
RUNX1.GoT.Illumina <- read.table(file = "./data/benchmarking/GoT_nanoranger/AML1022_pileup_RUNX1.GoT.Illumina.processed.csv")
RUNX1.GoT.Illumina$condition <- "GoT.Illumina"
RUNX1.normal <- read.table(file = "./data/benchmarking/GoT_nanoranger/AML1022_pileup_RUNX1.normal.processed.csv")
RUNX1.normal$condition <- "normal"
RUNX1.spikein <- read.table(file = "./data/benchmarking/GoT_nanoranger/AML1022_pileup_RUNX1.spikein.processed.csv")
RUNX1.spikein$condition <- "spikein"
RUNX1.data <- dplyr::bind_rows(RUNX1.GoT, RUNX1.GoT.Illumina, RUNX1.normal, RUNX1.spikein) %>%
  filter(bc %in% colnames(AML1022))
RUNX1.data$predicted.celltype <- AML1022$predicted.celltype[RUNX1.data$bc]

SF3B1.GoT <- read.table(file = "./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.GoT.ONT.processed.csv")
SF3B1.GoT$condition <- "GoT"
SF3B1.GoT.Illumina <- read.table(file = "./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.GoT.Illumina.processed.csv")
SF3B1.GoT.Illumina$condition <- "GoT.Illumina"
SF3B1.normal <- read.table(file = "./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.normal.processed.csv")
SF3B1.normal$condition <- "normal"
SF3B1.spikein <- read.table(file = "./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.spikein.processed.csv")
SF3B1.spikein$condition <- "spikein"
SF3B1.data <- dplyr::bind_rows(SF3B1.GoT, SF3B1.GoT.Illumina, SF3B1.normal, SF3B1.spikein) %>%
  filter(bc %in% colnames(AML1022))
SF3B1.data$predicted.celltype <- AML1022$predicted.celltype[SF3B1.data$bc]

sample.df <- data.frame(
  barcode = colnames(AML1022),
  sample = AML1022$orig.ident,
  predicted.celltype = AML1022$predicted.celltype
) %>%
  group_by(sample, predicted.celltype) %>%
  tally() %>%
  tidyr::pivot_wider(names_from = "sample", values_from = "n") %>%
  as.data.frame()
rownames(sample.df) <- sample.df$predicted.celltype

DNMT3A.mat <- DNMT3A.data %>%
  group_by(condition, predicted.celltype) %>%
  summarize(cells = length(predicted.celltype)) %>%
  tidyr::pivot_wider(names_from = "condition", values_from = "cells") %>%
  as.data.frame()
rownames(DNMT3A.mat) <- DNMT3A.mat$predicted.celltype
DNMT3A.mat$predicted.celltype <- factor(DNMT3A.mat$predicted.celltype, levels = names(AML.combined.colors))

RUNX1.mat <- RUNX1.data %>%
  group_by(condition, predicted.celltype) %>%
  summarize(cells = length(predicted.celltype)) %>%
  tidyr::pivot_wider(names_from = "condition", values_from = "cells") %>%
  as.data.frame()
rownames(RUNX1.mat) <- RUNX1.mat$predicted.celltype
RUNX1.mat$predicted.celltype <- factor(RUNX1.mat$predicted.celltype, levels = names(AML.combined.colors))

SF3B1.mat <- SF3B1.data %>%
  group_by(condition, predicted.celltype) %>%
  summarize(cells = length(predicted.celltype)) %>%
  tidyr::pivot_wider(names_from = "condition", values_from = "cells") %>%
  as.data.frame()
rownames(SF3B1.mat) <- SF3B1.mat$predicted.celltype
SF3B1.mat$predicted.celltype <- factor(SF3B1.mat$predicted.celltype, levels = names(AML.combined.colors))

celltypes <- c(myeloid.clusters, "Prog_RBC", "Prog_Mk")

# DNMT3A
sample.df$DNMT3A.normal = DNMT3A.mat[sample.df$predicted.celltype, 'normal'] / sample.df$AML1022.1
sample.df$DNMT3A.spikein = DNMT3A.mat[sample.df$predicted.celltype, 'spikein'] / sample.df$AML1022.1G
sample.df$DNMT3A.GoT.Illumina = DNMT3A.mat[sample.df$predicted.celltype, 'GoT.Illumina'] / sample.df$AML1022.1G
sample.df$DNMT3A.GoT.ONT = DNMT3A.mat[sample.df$predicted.celltype, 'GoT'] / sample.df$AML1022.1G

max.cells = max(sample.df[,c('DNMT3A.normal', 'DNMT3A.spikein', 
                             'DNMT3A.GoT.Illumina', 'DNMT3A.GoT.ONT')], na.rm = T)
max.cells = 0.5
col_fun = circlize::colorRamp2(breaks = seq(0,max.cells, max.cells/9),  
                               colors = c('white',BuenColors::jdb_palette(name = 'samba_night')))
ha = columnAnnotation(condition = c('normal', 'GoT.Illumina', 'GoT.ONT'),
                      col = list('condition' = c('normal' = 'orange',
                                                 'GoT.Illumina' = 'black', 'GoT.ONT' = '#00aeef')),
                      simple_anno_size = unit(5, 'pt'), border=T)
svglite::svglite('./benchmarking/figures/GoT_nanoranger/heatmaps/20230912_AML1022_DNMT3A_genotyped_cells.svg', 
                 width = 3.2, height = 2.8)
ComplexHeatmap::Heatmap(sample.df[celltypes, c('DNMT3A.normal', 
                                               'DNMT3A.GoT.Illumina', 'DNMT3A.GoT.ONT')], 
                        col = col_fun, row_names_side = 'left', row_names_gp = gpar(fontsize=10),
                        column_names_gp = gpar(fontsize=10),
                        na_col = 'grey90', cluster_rows = F, cluster_columns = F, border = T,
                        top_annotation = ha, column_labels = c('normal', 'GoT Illumina', 'GoT ONT'))
dev.off()

# RUNX1
sample.df$RUNX1.normal = RUNX1.mat[sample.df$predicted.celltype, 'normal'] / sample.df$AML1022.1
sample.df$RUNX1.spikein = RUNX1.mat[sample.df$predicted.celltype, 'spikein'] / sample.df$AML1022.1G
sample.df$RUNX1.GoT.Illumina = RUNX1.mat[sample.df$predicted.celltype, 'GoT.Illumina'] / sample.df$AML1022.1G
sample.df$RUNX1.GoT.ONT = RUNX1.mat[sample.df$predicted.celltype, 'GoT'] / sample.df$AML1022.1G

max.cells = max(sample.df[,c('RUNX1.normal', 'RUNX1.spikein', 
                             'RUNX1.GoT.Illumina', 'RUNX1.GoT.ONT')], na.rm = T)
max.cells = 0.5
col_fun = circlize::colorRamp2(breaks = seq(0,max.cells, max.cells/9),  
                               colors = c('white',BuenColors::jdb_palette(name = 'samba_night')))
ha = columnAnnotation(condition = c('normal', 'GoT.Illumina', 'GoT.ONT'),
                      col = list('condition' = c('normal' = 'orange',
                                                 'GoT.Illumina' = 'black', 'GoT.ONT' = '#00aeef')),
                      simple_anno_size = unit(5, 'pt'), border=T)
svglite::svglite('./benchmarking/figures/GoT_nanoranger/heatmaps/20230912_AML1022_RUNX1_genotyped_cells.svg', 
                 width = 3.2, height = 2.8)
ComplexHeatmap::Heatmap(sample.df[celltypes, c('RUNX1.normal', 
                                               'RUNX1.GoT.Illumina', 'RUNX1.GoT.ONT')], 
                        col = col_fun, row_names_side = 'left', row_names_gp = gpar(fontsize=10),
                        column_names_gp = gpar(fontsize=10),
                        na_col = 'grey90', cluster_rows = F, cluster_columns = F, border = T,
                        top_annotation = ha, column_labels = c('normal', 'GoT Illumina', 'GoT ONT'))
dev.off()


# SF3B1
sample.df$SF3B1.normal = SF3B1.mat[sample.df$predicted.celltype, 'normal'] / sample.df$AML1022.1
sample.df$SF3B1.spikein = SF3B1.mat[sample.df$predicted.celltype, 'spikein'] / sample.df$AML1022.1G
sample.df$SF3B1.GoT.Illumina = SF3B1.mat[sample.df$predicted.celltype, 'GoT.Illumina'] / sample.df$AML1022.1G
sample.df$SF3B1.GoT.ONT = SF3B1.mat[sample.df$predicted.celltype, 'GoT'] / sample.df$AML1022.1G

max.cells = max(sample.df[,c('SF3B1.normal', 'SF3B1.spikein', 
                             'SF3B1.GoT.Illumina', 'SF3B1.GoT.ONT')], na.rm = T)
max.cells = 0.5
col_fun = circlize::colorRamp2(breaks = seq(0,max.cells, max.cells/9),  
                               colors = c('white',BuenColors::jdb_palette(name = 'samba_night')))
ha = columnAnnotation(condition = c('normal', 'GoT.Illumina', 'GoT.ONT'),
                      col = list('condition' = c('normal' = 'orange',
                                                 'GoT.Illumina' = 'black', 'GoT.ONT' = '#00aeef')),
                      simple_anno_size = unit(5, 'pt'), border=T)
svglite::svglite('./benchmarking/figures/GoT_nanoranger/heatmaps/20230912_AML1022_SF3B1_genotyped_cells.svg', 
                 width = 3.2, height = 2.8)
ComplexHeatmap::Heatmap(sample.df[celltypes, c('SF3B1.normal', 
                                               'SF3B1.GoT.Illumina', 'SF3B1.GoT.ONT')], 
                        col = col_fun, row_names_side = 'left', row_names_gp = gpar(fontsize=10),
                        column_names_gp = gpar(fontsize=10),
                        na_col = 'grey90', cluster_rows = F, cluster_columns = F, border = T,
                        top_annotation = ha, column_labels = c('normal', 'GoT Illumina', 'GoT ONT'))
dev.off()
