# number of mutated cells across the 3 de-novo AML cases according to somatic mutations and CNVs

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(nanoranger.R)
library(Seurat)
library(SummarizedExperiment)

FLT3_ITD <- readRDS("./data/AML/objects/20230523_FLT3-ITD.rds")

source("./R/celltypes.R")

# identify mutated cells
FLT3_ITD1 <- read.csv2("./data/AML/mutations/FLT3-ITD1_FLT3-ITD.csv", sep = "\t")
FLT3_ITD1_NPM1_1 <- read.csv2("./data/AML/mutations/FLT3-ITD1_NPM1.csv", sep = "\t")
FLT3_ITD1_NPM1_2 <- read.csv2("./data/AML/mutations/FLT3-ITD1_NPM1.2.csv", sep = "\t")
FLT3_ITD1.mutated <- FLT3_ITD1$bc[which(FLT3_ITD1$mutated == "mutated")]
FLT3_ITD1_NPM1_1.mutated <- FLT3_ITD1_NPM1_1$bc[which(FLT3_ITD1_NPM1_1$mutated == "mutated")]
FLT3_ITD1_NPM1_2.mutated <- FLT3_ITD1_NPM1_2$bc[which(FLT3_ITD1_NPM1_2$mutated == "mutated")]

FLT3_ITD2 <- read.csv2("./data/AML/mutations/FLT3-ITD2_FLT3-ITD.csv", sep = "\t")
FLT3_ITD2.mutated <- FLT3_ITD2$bc[which(FLT3_ITD2$mutated == "mutated")]

FLT3_ITD3 <- read.csv2("./data/AML/mutations/FLT3-ITD3_FLT3-ITD.csv", sep = "\t")
FLT3_ITD3_NPM1 <- read.csv2("./data/AML/mutations/FLT3-ITD3_NPM1.csv", sep = "\t")
FLT3_ITD3_DNMT3A <- read.csv2("./data/AML/mutations/FLT3-ITD3_DNMT3A.csv", sep = "\t")
FLT3_ITD3.mutated <- FLT3_ITD3$bc[which(FLT3_ITD3$mutated == "mutated")]
FLT3_ITD3_NPM1.mutated <- FLT3_ITD3_NPM1$bc[which(FLT3_ITD3_NPM1$mutated == "mutated")]
FLT3_ITD3_DNMT3A.mutated <- FLT3_ITD3_DNMT3A$bc[which(FLT3_ITD3_DNMT3A$mutated == "mutated")]

mutated.cells <- unique(c(
  FLT3_ITD1_NPM1_1.mutated, FLT3_ITD1_NPM1_2.mutated, FLT3_ITD2.mutated,
  FLT3_ITD3_NPM1.mutated, FLT3_ITD3_DNMT3A.mutated
))

genotyped.cells <- unique(c(
  FLT3_ITD1$bc, FLT3_ITD1_NPM1_1$bc, FLT3_ITD1_NPM1_2$bc,
  FLT3_ITD2$bc,
  FLT3_ITD3$bc, FLT3_ITD3_DNMT3A$bc, FLT3_ITD3_NPM1$bc
))

# generate statistics
df <- data.frame(
  bc = colnames(FLT3_ITD),
  predicted.celltype = as.character(FLT3_ITD$predicted.celltype)
)
rownames(df) <- df$bc
df$mutated <- ifelse(rownames(df) %in% mutated.cells, "yes", "no")
df <- df[intersect(genotyped.cells, rownames(df)), ]
df$sample <- stringr::str_split_fixed(df$bc, pattern = "_", n = 2)[, 1]

stats <- df %>%
  group_by(sample, predicted.celltype) %>%
  summarize(
    mut = length(which(mutated == "yes")),
    wt = length(which(mutated == "no"))
  )
stats$cells <- stats$mut + stats$wt
stats$mut.freq <- stats$mut / (stats$mut + stats$wt)
stats$predicted.celltype <- factor(stats$predicted.celltype, levels = names(nanoranger.R::AML.combined.colors))

# plot
ggplot(
  stats[which(stats$predicted.celltype %in% c(myeloid.clusters, "Prog_Mk", "Prog_RBC")), ],
  aes(x = predicted.celltype, y = 100 * mut.freq)
) +
  stat_summary(geom = "crossbar", width = 0.5, fun = median, fun.min = median, fun.max = median) +
  geom_point(aes(fill = predicted.celltype), shape = 21) +
  scale_fill_manual(values = AML.combined.colors) +
  scale_y_continuous("% mutated") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank()
  )
ggsave("./AML/figures/de_novo/plots/20230806_de_novo_AML_celltypes_mutation_rate.svg", width = 2, height = 2)

# find cells with CNV calls
cnv.data.1 <- read.csv(file = "./data/AML/numbat/FLT3_ITD1/FLT3-ITD1_cnv_calls.csv", sep = "\t")
cnv.data.2 <- read.csv(file = "./data/AML/numbat/FLT3_ITD2/FLT3-ITD2_cnv_calls.csv", sep = "\t")
cnv.data <- rbind(
  cnv.data.1[, c("bc", "mutated", "predicted.celltype")],
  cnv.data.2[, c("bc", "mutated", "predicted.celltype")]
)

# statistics
df <- data.frame(
  bc = colnames(FLT3_ITD),
  predicted.celltype = as.character(FLT3_ITD$predicted.celltype)
)
rownames(df) <- df$bc
df$mutated <- ifelse(rownames(df) %in% cnv.data$bc[which(cnv.data$mutated == "mutated")], "yes", "no")
df <- df[intersect(cnv.data$bc, rownames(df)), ]
df$sample <- stringr::str_split_fixed(df$bc, pattern = "_", n = 2)[, 1]

stats <- df %>%
  group_by(sample, predicted.celltype) %>%
  summarize(
    mut = length(which(mutated == "yes")),
    wt = length(which(mutated == "no"))
  )
stats$cells <- stats$mut + stats$wt
stats$mut.freq <- stats$mut / (stats$mut + stats$wt)
stats$predicted.celltype <- factor(stats$predicted.celltype, levels = names(nanoranger.R::AML.combined.colors))

ggplot(
  stats[which(stats$predicted.celltype %in% c(myeloid.clusters, "Prog_Mk", "Prog_RBC")), ],
  aes(x = predicted.celltype, y = 100 * mut.freq)
) +
  stat_summary(geom = "crossbar", width = 0.5, fun = median, fun.min = median, fun.max = median) +
  geom_point(aes(fill = predicted.celltype), shape = 21) +
  scale_fill_manual(values = AML.combined.colors) +
  scale_y_continuous("% CNV", limits = c(0, 100)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank()
  )
ggsave("./AML/figures/de_novo/plots/20230806_de_novo_AML_celltypes_CNV_rate.svg", width = 2, height = 2)
