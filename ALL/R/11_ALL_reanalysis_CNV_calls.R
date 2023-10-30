# reanalysis of data from Witkowski et al., Cancer Cell 2020 and Caron et al., Scientific Reports 2020
#
# perform CNV calls from numbat output

library(dplyr)
library(ggplot2)
library(numbat)
library(Seurat)
library(numbat)

ALL.validation <- readRDS(file = "./data/ALL/objects/20221027_ALL_validation.rds")

df <- as.data.frame(ALL.validation@reductions$ref.umap@cell.embeddings)
df$mutated <- NA
df$sample <- ALL.validation$orig.ident

# ETV001 - use +16, +17 and +21
nb <- Numbat$new(out_dir = "./data/ALL/reanalysis/ETV001.numbat/ETV001_internal_ref/")
nb$joint_post$present <- ifelse(nb$joint_post$p_cnv > 0.8, 1, 0)
cnv.data <- as.data.frame(nb$joint_post %>% filter(seg %in% c("16a", "17b", "21a")) %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y, present))
cnv.data$p_cnv <- as.numeric(cnv.data$p_cnv)
cnv.data <- as.data.frame(cnv.data %>% select(cell, seg, p_cnv) %>% filter(seg %in% c("16a", "17b", "21a")) %>% tidyr::pivot_wider(names_from = "seg", values_from = "p_cnv"))
rownames(cnv.data) <- cnv.data$cell
colnames(cnv.data) <- c("bc", "amp16", "amp17", "amp21")
cnv.data$mutated <- ifelse(rowMeans(cnv.data[, c("amp16", "amp17", "amp21")]) > 0.8, "mutated", "wildtype")
cnv.data$predicted.celltype <- ALL.validation$predicted.celltype[rownames(cnv.data)]
cnv.data$sample <- "ETV001"
write.table(cnv.data, file = "./data/ALL/numbat/ETV001_cnv_calls.csv", sep = "\t", quote = F)

cnv.data <- read.csv2("./data/ALL/numbat/ETV001_cnv_calls.csv", sep = "\t")
df[rownames(cnv.data), "mutated"] <- cnv.data$mutated

ggplot() +
  geom_point(data = df, aes(x = refUMAP_1, y = refUMAP_2), color = "grey90") +
  geom_point(data = df[which(df$sample == "ETV001" & df$mutated == "wildtype"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "black") +
  geom_point(data = df[which(df$sample == "ETV001" & df$mutated == "mutated"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "red") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
ggsave("./ALL/figures/UMAP/20221028_ETV001_mutation.png", width = 4, height = 4, dpi = 600)

cnv.data %>%
  group_by(predicted.celltype) %>%
  summarize(
    mut.freq = length(which(mutated == "mutated")) / length(predicted.celltype),
    cells = length(predicted.celltype)
  )


# ETV002 - use del(6q) and del(11p)
nb <- Numbat$new(out_dir = "./data/ALL/reanalysis/ETV002.numbat/ETV002_internal_ref/")
nb$joint_post$present <- ifelse(nb$joint_post$p_cnv > 0.8, 1, 0)
cnv.data <- as.data.frame(nb$joint_post %>% filter(seg %in% c("6b", "11a")) %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y, present))
cnv.data$p_cnv <- as.numeric(cnv.data$p_cnv)
cnv.data <- as.data.frame(cnv.data %>% select(cell, seg, p_cnv) %>% filter(seg %in% c("6b", "11a")) %>% tidyr::pivot_wider(names_from = "seg", values_from = "p_cnv"))
rownames(cnv.data) <- cnv.data$cell
colnames(cnv.data) <- c("bc", "del6", "del11")
cnv.data$mutated <- ifelse(rowMeans(cnv.data[, c("del6", "del11")]) > 0.8, "mutated", "wildtype")
cnv.data$predicted.celltype <- ALL.validation$predicted.celltype[rownames(cnv.data)]
cnv.data$sample <- "ETV002"
write.table(cnv.data, file = "./data/ALL/numbat/ETV002_cnv_calls.csv", sep = "\t", quote = F)

cnv.data <- read.csv2("./data/ALL/numbat/ETV002_cnv_calls.csv", sep = "\t")
df[rownames(cnv.data), "mutated"] <- cnv.data$mutated

ggplot() +
  geom_point(data = df, aes(x = refUMAP_1, y = refUMAP_2), color = "grey90") +
  geom_point(data = df[which(df$sample == "ETV002" & df$mutated == "wildtype"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "black") +
  geom_point(data = df[which(df$sample == "ETV002" & df$mutated == "mutated"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "red") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
ggsave("./ALL/figures/UMAP/20221028_ETV002_mutation.png", width = 4, height = 4, dpi = 600)

cnv.data %>%
  group_by(predicted.celltype) %>%
  summarize(
    mut.freq = length(which(mutated == "mutated")) / length(predicted.celltype),
    cells = length(predicted.celltype)
  )

# ETV003 - use loh(6q), amp(10p), del(10q), del(15q)
nb <- Numbat$new(out_dir = "./data/ALL/reanalysis/ETV003.numbat/ETV003_internal_ref/")
nb$joint_post$present <- ifelse(nb$joint_post$p_cnv > 0.8, 1, 0)
cnv.data <- as.data.frame(nb$joint_post %>% filter(seg %in% c("6b", "10a", "10c", "15b")) %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y, present))
cnv.data$p_cnv <- as.numeric(cnv.data$p_cnv)
cnv.data <- as.data.frame(cnv.data %>% select(cell, seg, p_cnv) %>% filter(seg %in% c("6b", "10a", "10c", "15b")) %>% tidyr::pivot_wider(names_from = "seg", values_from = "p_cnv"))
rownames(cnv.data) <- cnv.data$cell
colnames(cnv.data) <- c("bc", "del6", "amp10", "del10", "del15")
cnv.data$mutated <- ifelse(rowMeans(cnv.data[, c("del6", "amp10", "del10", "del15")]) > 0.8, "mutated", "wildtype")
cnv.data$predicted.celltype <- ALL.validation$predicted.celltype[rownames(cnv.data)]
cnv.data$sample <- "ETV003"
write.table(cnv.data, file = "./data/ALL/numbat/ETV003_cnv_calls.csv", sep = "\t", quote = F)

cnv.data <- read.csv2("./data/ALL/numbat/ETV003_cnv_calls.csv", sep = "\t")
df[rownames(cnv.data), "mutated"] <- cnv.data$mutated

ggplot() +
  geom_point(data = df, aes(x = refUMAP_1, y = refUMAP_2), color = "grey90") +
  geom_point(data = df[which(df$sample == "ETV003" & df$mutated == "wildtype"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "black") +
  geom_point(data = df[which(df$sample == "ETV003" & df$mutated == "mutated"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "red") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
ggsave("./ALL/figures/UMAP/20221028_ETV003_mutation.png", width = 4, height = 4, dpi = 600)

cnv.data %>%
  group_by(predicted.celltype) %>%
  summarize(
    mut.freq = length(which(mutated == "mutated")) / length(predicted.celltype),
    cells = length(predicted.celltype)
  )


# PH001 - use loh(9p), loh(15p)
nb <- Numbat$new(out_dir = "./data/ALL/reanalysis/PH001.numbat/PH001_internal_ref//")
nb$joint_post$present <- ifelse(nb$joint_post$p_cnv > 0.8, 1, 0)
cnv.data <- as.data.frame(nb$joint_post %>% filter(seg %in% c("9a", "15b")) %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y, present))
cnv.data$p_cnv <- as.numeric(cnv.data$p_cnv)
cnv.data <- as.data.frame(cnv.data %>% select(cell, seg, p_cnv) %>% filter(seg %in% c("9a", "15b")) %>% tidyr::pivot_wider(names_from = "seg", values_from = "p_cnv"))
rownames(cnv.data) <- cnv.data$cell
colnames(cnv.data) <- c("bc", "del9", "loh15")
cnv.data$mutated <- ifelse(rowMeans(cnv.data[, c("del9", "loh15")]) > 0.8, "mutated", "wildtype")
# cnv.data$mutated = ifelse(cnv.data[,c('loh15')] > 0.8, 'mutated', 'wildtype')
cnv.data$predicted.celltype <- ALL.validation$predicted.celltype[rownames(cnv.data)]
cnv.data$sample <- "PH001"
write.table(cnv.data, file = "./data/ALL/numbat/PH001_cnv_calls.csv", sep = "\t", quote = F)

cnv.data <- read.csv2("./data/ALL/numbat/PH001_cnv_calls.csv", sep = "\t")
df[rownames(cnv.data), "mutated"] <- cnv.data$mutated

ggplot() +
  geom_point(data = df, aes(x = refUMAP_1, y = refUMAP_2), color = "grey90") +
  geom_point(data = df[which(df$sample == "PH001" & df$mutated == "wildtype"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "black") +
  geom_point(data = df[which(df$sample == "PH001" & df$mutated == "mutated"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "red") +
  scale_color_manual(values = c("mutated" = "red", "wildtype" = "black")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
ggsave("./ALL/figures/UMAP/20221028_PH001_mutation.png", width = 4, height = 4, dpi = 600)


# PH002 - use +2 +5 +6 +14 +18 +21
nb <- Numbat$new(out_dir = "./data/ALL/reanalysis/PH002.numbat/PH002/")
nb$joint_post$present <- ifelse(nb$joint_post$p_cnv > 0.8, 1, 0)
cnv.data <- as.data.frame(nb$joint_post %>% filter(seg %in% c("2b", "5a", "6a", "14a", "18a", "21a")) %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y, present))
cnv.data$p_cnv <- as.numeric(cnv.data$p_cnv)
cnv.data <- as.data.frame(cnv.data %>% select(cell, seg, p_cnv) %>% filter(seg %in% c("2b", "5a", "6a", "14a", "18a", "21a")) %>% tidyr::pivot_wider(names_from = "seg", values_from = "p_cnv"))
rownames(cnv.data) <- cnv.data$cell
colnames(cnv.data) <- c("bc", "amp2", "amp5", "amp6", "amp14", "amp18", "amp21")
cnv.data$mutated <- ifelse(rowMeans(cnv.data[, c("amp2", "amp5", "amp6", "amp14", "amp14", "amp21")]) > 0.8, "mutated", "wildtype")
cnv.data$predicted.celltype <- ALL.validation$predicted.celltype[rownames(cnv.data)]
cnv.data$sample <- "PH002"
write.table(cnv.data, file = "./data/ALL/numbat/PH002_cnv_calls.csv", sep = "\t", quote = F)

cnv.data <- read.csv2("./data/ALL/numbat/PH002_cnv_calls.csv", sep = "\t")
df[rownames(cnv.data), "mutated"] <- cnv.data$mutated

ggplot() +
  geom_point(data = df, aes(x = refUMAP_1, y = refUMAP_2), color = "grey90") +
  geom_point(data = df[which(df$sample == "PH002" & df$mutated == "wildtype"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "black") +
  geom_point(data = df[which(df$sample == "PH002" & df$mutated == "mutated"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "red") +
  scale_color_manual(values = c("mutated" = "red", "wildtype" = "black")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
ggsave("./ALL/figures/UMAP/20221028_PH002_mutation.png", width = 4, height = 4, dpi = 600)

# ETV6-RUNX1_1 - use +10 +16
nb <- Numbat$new(out_dir = "./data/ALL/reanalysis/ETV6-RUNX1_1.numbat/ETV6-RUNX1_1/")
nb$joint_post$present <- ifelse(nb$joint_post$p_cnv > 0.8, 1, 0)
cnv.data <- as.data.frame(nb$joint_post %>% filter(seg %in% c("10a", "16a")) %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y, present))
cnv.data$p_cnv <- as.numeric(cnv.data$p_cnv)
cnv.data <- as.data.frame(cnv.data %>% select(cell, seg, p_cnv) %>% filter(seg %in% c("10a", "16a")) %>% tidyr::pivot_wider(names_from = "seg", values_from = "p_cnv"))
rownames(cnv.data) <- cnv.data$cell
colnames(cnv.data) <- c("bc", "amp10", "amp16")
cnv.data$mutated <- ifelse(rowMeans(cnv.data[, c("amp10", "amp16")]) > 0.8, "mutated", "wildtype")
cnv.data$predicted.celltype <- ALL.validation$predicted.celltype[rownames(cnv.data)]
cnv.data$sample <- "ETV6-RUNX1_1"
write.table(cnv.data, file = "./data/ALL/numbat/ETV6-RUNX1_1_cnv_calls.csv", sep = "\t", quote = F)

cnv.data <- read.csv2("./data/ALL/numbat/ETV6-RUNX1_1_cnv_calls.csv", sep = "\t")
df[rownames(cnv.data), "mutated"] <- cnv.data$mutated

ggplot() +
  geom_point(data = df, aes(x = refUMAP_1, y = refUMAP_2), color = "grey90") +
  geom_point(data = df[which(df$sample == "ETV6-RUNX1_1" & df$mutated == "wildtype"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "black") +
  geom_point(data = df[which(df$sample == "ETV6-RUNX1_1" & df$mutated == "mutated"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "red") +
  scale_color_manual(values = c("mutated" = "red", "wildtype" = "black")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
ggsave("./ALL/figures/UMAP/20221028_ETV-RUNX1_1_mutation.png", width = 4, height = 4, dpi = 600)


# ETV6-RUNX1_2 - use del(5q), del(11q), del(12p)
nb <- Numbat$new(out_dir = "./data/ALL/reanalysis/ETV6-RUNX1_2.numbat/ETV6-RUNX1_2/")
nb$joint_post$present <- ifelse(nb$joint_post$p_cnv > 0.8, 1, 0)
cnv.data <- as.data.frame(nb$joint_post %>% filter(seg %in% c("5c", "11b", "12a")) %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y, present))
cnv.data$p_cnv <- as.numeric(cnv.data$p_cnv)
cnv.data <- as.data.frame(cnv.data %>% select(cell, seg, p_cnv) %>% filter(seg %in% c("5c", "11b", "12a")) %>% tidyr::pivot_wider(names_from = "seg", values_from = "p_cnv"))
rownames(cnv.data) <- cnv.data$cell
colnames(cnv.data) <- c("bc", "del5", "del11", "del12")
cnv.data$mutated <- ifelse(rowMeans(cnv.data[, c("del11", "del12")]) > 0.8, "mutated", "wildtype")
cnv.data$predicted.celltype <- ALL.validation$predicted.celltype[rownames(cnv.data)]
cnv.data$sample <- "ETV6-RUNX1_2"
write.table(cnv.data, file = "./data/ALL/numbat/ETV6-RUNX1_2_cnv_calls.csv", sep = "\t", quote = F)

cnv.data <- read.csv2("./data/ALL/numbat/ETV6-RUNX1_2_cnv_calls.csv", sep = "\t")
df[rownames(cnv.data), "mutated"] <- cnv.data$mutated

ggplot() +
  geom_point(data = df, aes(x = refUMAP_1, y = refUMAP_2), color = "grey90") +
  geom_point(data = df[which(df$sample == "ETV6-RUNX1_2" & df$mutated == "wildtype"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "black") +
  geom_point(data = df[which(df$sample == "ETV6-RUNX1_2" & df$mutated == "mutated"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "red") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
ggsave("./ALL/figures/UMAP/20221028_ETV-RUNX1_2_mutation.png", width = 4, height = 4, dpi = 600)

# ETV6-RUNX1_3 - use del(8q) del(12p)
#
# this CVN calling did not work
nb <- Numbat$new(out_dir = "./data/ALL/reanalysis/ETV6-RUNX1_3.numbat/ETV6-RUNX1_3/")
nb$joint_post$present <- ifelse(nb$joint_post$p_cnv > 0.8, 1, 0)
cnv.data <- as.data.frame(nb$joint_post %>% filter(seg %in% c("8a", "12a")) %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y, present))
cnv.data$p_cnv <- as.numeric(cnv.data$p_cnv)
cnv.data <- as.data.frame(cnv.data %>% select(cell, seg, p_cnv) %>% filter(seg %in% c("8a", "12a")) %>% tidyr::pivot_wider(names_from = "seg", values_from = "p_cnv"))
rownames(cnv.data) <- cnv.data$cell
colnames(cnv.data) <- c("bc", "del8", "del12")
cnv.data$mutated <- ifelse(rowMeans(cnv.data[, c("del8", "del12")]) > 0.8, "mutated", "wildtype")
cnv.data$predicted.celltype <- ALL.validation$predicted.celltype[rownames(cnv.data)]
cnv.data$sample <- "ETV6-RUNX1_3"
write.table(cnv.data, file = "./data/ALL/numbat/ETV6-RUNX1_3_cnv_calls.csv", sep = "\t", quote = F)

cnv.data <- read.csv2("./data/ALL/numbat/ETV6-RUNX1_3_cnv_calls.csv", sep = "\t")
df[rownames(cnv.data), "mutated"] <- cnv.data$mutated

ggplot() +
  geom_point(data = df, aes(x = refUMAP_1, y = refUMAP_2), color = "grey90") +
  geom_point(data = df[which(df$sample == "ETV6-RUNX1_3" & df$mutated == "wildtype"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "black") +
  geom_point(data = df[which(df$sample == "ETV6-RUNX1_3" & df$mutated == "mutated"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "red") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
ggsave("./ALL/figures/UMAP/20221028_ETV-RUNX1_3_mutation.png", width = 4, height = 4, dpi = 600)

# ETV6-RUNX1_4 - use del(9p), del(12p), amp(17q), del(18q)
nb <- Numbat$new(out_dir = "./data/ALL/reanalysis/ETV6-RUNX1_4.numbat/ETV6-RUNX1_4/")
nb$joint_post$present <- ifelse(nb$joint_post$p_cnv > 0.8, 1, 0)
cnv.data <- as.data.frame(nb$joint_post %>% filter(seg %in% c("9a", "12a", "17b", "18b")) %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y, present))
cnv.data$p_cnv <- as.numeric(cnv.data$p_cnv)
cnv.data <- as.data.frame(cnv.data %>% select(cell, seg, p_cnv) %>% filter(seg %in% c("9a", "12a", "17b", "18b")) %>% tidyr::pivot_wider(names_from = "seg", values_from = "p_cnv"))
rownames(cnv.data) <- cnv.data$cell
colnames(cnv.data) <- c("bc", "del9", "del12", "amp17", "del18")
cnv.data$mutated <- ifelse(rowMeans(cnv.data[, c("del9", "del12", "amp17", "del18")]) > 0.8, "mutated", "wildtype")
cnv.data$predicted.celltype <- ALL.validation$predicted.celltype[rownames(cnv.data)]
cnv.data$sample <- "ETV6-RUNX1_4"
write.table(cnv.data, file = "./data/ALL/numbat/ETV6-RUNX1_4_cnv_calls.csv", sep = "\t", quote = F)

cnv.data <- read.csv2("./data/ALL/numbat/ETV6-RUNX1_4_cnv_calls.csv", sep = "\t")
df[rownames(cnv.data), "mutated"] <- cnv.data$mutated

ggplot() +
  geom_point(data = df, aes(x = refUMAP_1, y = refUMAP_2), color = "grey90") +
  geom_point(data = df[which(df$sample == "ETV6-RUNX1_4" & df$mutated == "wildtype"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "black") +
  geom_point(data = df[which(df$sample == "ETV6-RUNX1_4" & df$mutated == "mutated"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "red") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
ggsave("./ALL/figures/UMAP/20221028_ETV-RUNX1_4_mutation.png", width = 4, height = 4, dpi = 600)


# HHD1 - use +4, +8, +9, +14, +17, +18, +21
nb <- Numbat$new(out_dir = "./data/ALL/reanalysis/HHD_1.numbat/HHD_1/")
nb$joint_post$present <- ifelse(nb$joint_post$p_cnv > 0.8, 1, 0)
cnv.data <- as.data.frame(nb$joint_post %>% filter(seg %in% c("4a", "8a", "9b", "14a", "17a", "18a", "21a")) %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y, present))
cnv.data$p_cnv <- as.numeric(cnv.data$p_cnv)
cnv.data <- as.data.frame(cnv.data %>% select(cell, seg, p_cnv) %>% filter(seg %in% c("4a", "8a", "9b", "14a", "17a", "18a", "21a")) %>% tidyr::pivot_wider(names_from = "seg", values_from = "p_cnv"))
rownames(cnv.data) <- cnv.data$cell
colnames(cnv.data) <- c("bc", "amp4", "amp8", "amp9", "amp14", "amp17", "amp18", "amp21")
cnv.data$mutated <- ifelse(rowMeans(cnv.data[, c("amp4", "amp8", "amp9", "amp14", "amp17", "amp18", "amp21")]) > 0.8, "mutated", "wildtype")
cnv.data$predicted.celltype <- ALL.validation$predicted.celltype[rownames(cnv.data)]
cnv.data$sample <- "HHD_1"
write.table(cnv.data, file = "./data/ALL/numbat/HHD_1_cnv_calls.csv", sep = "\t", quote = F)

cnv.data <- read.csv2("./data/ALL/numbat/HHD_1_cnv_calls.csv", sep = "\t")
df[rownames(cnv.data), "mutated"] <- cnv.data$mutated

ggplot() +
  geom_point(data = df, aes(x = refUMAP_1, y = refUMAP_2), color = "grey90") +
  geom_point(data = df[which(df$sample == "HHD_1" & df$mutated == "wildtype"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "black") +
  geom_point(data = df[which(df$sample == "HHD_1" & df$mutated == "mutated"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "red") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
ggsave("./ALL/figures/UMAP/20221028_HHD_1_mutation.png", width = 4, height = 4, dpi = 600)


# HHD_2 - use +4, +6, +14, +18
nb <- Numbat$new(out_dir = "./data/ALL/reanalysis/HHD_2.numbat/HHD_2/")
nb$joint_post$present <- ifelse(nb$joint_post$p_cnv > 0.8, 1, 0)
cnv.data <- as.data.frame(nb$joint_post %>% filter(seg %in% c("4a", "6a", "14a", "18a")) %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y, present))
cnv.data$p_cnv <- as.numeric(cnv.data$p_cnv)
cnv.data <- as.data.frame(cnv.data %>% select(cell, seg, p_cnv) %>% filter(seg %in% c("4a", "6a", "14a", "18a")) %>% tidyr::pivot_wider(names_from = "seg", values_from = "p_cnv"))
rownames(cnv.data) <- cnv.data$cell
colnames(cnv.data) <- c("bc", "amp4", "amp8", "amp14", "amp18")
cnv.data$mutated <- ifelse(rowMeans(cnv.data[, c("amp4", "amp8", "amp14", "amp18")]) > 0.8, "mutated", "wildtype")
cnv.data$predicted.celltype <- ALL.validation$predicted.celltype[rownames(cnv.data)]
cnv.data$sample <- "HHD_2"
write.table(cnv.data, file = "./data/ALL/numbat/HHD_2_cnv_calls.csv", sep = "\t", quote = F)

cnv.data <- read.csv2("./data/ALL/numbat/HHD_2_cnv_calls.csv", sep = "\t")
df[rownames(cnv.data), "mutated"] <- cnv.data$mutated

ggplot() +
  geom_point(data = df, aes(x = refUMAP_1, y = refUMAP_2), color = "grey90") +
  geom_point(data = df[which(df$sample == "HHD_2" & df$mutated == "wildtype"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "black") +
  geom_point(data = df[which(df$sample == "HHD_2" & df$mutated == "mutated"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "red") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
ggsave("./ALL/figures/UMAP/20221028_HHD_2_mutation.png", width = 4, height = 4, dpi = 600)



# PRE-T_1 - use amp(9p), del(11q)
nb <- Numbat$new(out_dir = "./data/ALL/reanalysis/PRE-T_1.numbat/PRE-T_1/")
nb$joint_post$present <- ifelse(nb$joint_post$p_cnv > 0.8, 1, 0)
cnv.data <- as.data.frame(nb$joint_post %>% filter(seg %in% c("9a", "11a")) %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y, present))
cnv.data$p_cnv <- as.numeric(cnv.data$p_cnv)
cnv.data <- as.data.frame(cnv.data %>% select(cell, seg, p_cnv) %>% filter(seg %in% c("9a", "11a")) %>% tidyr::pivot_wider(names_from = "seg", values_from = "p_cnv"))
rownames(cnv.data) <- cnv.data$cell
colnames(cnv.data) <- c("bc", "amp9", "del11")
cnv.data$mutated <- ifelse(rowMeans(cnv.data[, c("amp9", "del11")]) > 0.8, "mutated", "wildtype")
cnv.data$predicted.celltype <- ALL.validation$predicted.celltype[rownames(cnv.data)]
cnv.data$sample <- "PRE-T_1"
write.table(cnv.data, file = "./data/ALL/numbat/PRE-T_1_cnv_calls.csv", sep = "\t", quote = F)

cnv.data <- read.csv2("./data/ALL/numbat/PRE-T_1_cnv_calls.csv", sep = "\t")
df[rownames(cnv.data), "mutated"] <- cnv.data$mutated

ggplot() +
  geom_point(data = df, aes(x = refUMAP_1, y = refUMAP_2), color = "grey90") +
  geom_point(data = df[which(df$sample == "PRE-T_1" & df$mutated == "wildtype"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "black") +
  geom_point(data = df[which(df$sample == "PRE-T_1" & df$mutated == "mutated"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "red") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
ggsave("./ALL/figures/UMAP/20221028_PRE-T_1_mutation.png", width = 4, height = 4, dpi = 600)


# PRE-T_2 - use del(6q), loh(9p)
nb <- Numbat$new(out_dir = "./data/ALL/reanalysis/PRE-T_2.numbat/PRE-T_2/")
nb$joint_post$present <- ifelse(nb$joint_post$p_cnv > 0.8, 1, 0)
cnv.data <- as.data.frame(nb$joint_post %>% filter(seg %in% c("6b", "9a")) %>% select(cell, CHROM, seg, cnv_state, p_cnv, p_cnv_x, p_cnv_y, present))
cnv.data$p_cnv <- as.numeric(cnv.data$p_cnv)
cnv.data <- as.data.frame(cnv.data %>% select(cell, seg, p_cnv) %>% filter(seg %in% c("6b", "9a")) %>% tidyr::pivot_wider(names_from = "seg", values_from = "p_cnv"))
rownames(cnv.data) <- cnv.data$cell
colnames(cnv.data) <- c("bc", "del9", "loh9")
cnv.data$mutated <- ifelse(rowMeans(cnv.data[, c("del9", "loh9")]) > 0.8, "mutated", "wildtype")
cnv.data$predicted.celltype <- ALL.validation$predicted.celltype[rownames(cnv.data)]
cnv.data$sample <- "PRE-T_2"
write.table(cnv.data, file = "./data/ALL/numbat/PRE-T_2_cnv_calls.csv", sep = "\t", quote = F)

cnv.data <- read.csv2("./data/ALL/numbat/PRE-T_2_cnv_calls.csv", sep = "\t")
df[rownames(cnv.data), "mutated"] <- cnv.data$mutated

ggplot() +
  geom_point(data = df, aes(x = refUMAP_1, y = refUMAP_2), color = "grey90") +
  geom_point(data = df[which(df$sample == "PRE-T_2" & df$mutated == "wildtype"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "black") +
  geom_point(data = df[which(df$sample == "PRE-T_2" & df$mutated == "mutated"), ], aes(x = refUMAP_1, y = refUMAP_2), color = "red") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )
ggsave("./ALL/figures/UMAP/20221028_PRE-T_2_mutation.png", width = 4, height = 4, dpi = 600)
