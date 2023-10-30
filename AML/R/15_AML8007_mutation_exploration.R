# visualize mutations for AML1007.14

library(ComplexHeatmap)
library(dplyr)
library(nanoranger.R)
library(ggplot2)
library(Seurat)
library(SummarizedExperiment)

source("./R/celltypes.R")

AML8007 <- readRDS("./data/AML/objects/20220911_AML8007.rds")
Idents(AML8007) <- "predicted.celltype"

DNMT3A.vaf <- read.table("./data/AML/mutations/AML8007_DNMT3A.csv", sep = "\t") %>% as.data.frame()
TP53_1.vaf <- read.table("./data/AML/mutations/AML8007_TP53.1.csv", sep = "\t") %>% as.data.frame()
TP53_2.vaf <- read.table("./data/AML/mutations/AML8007_TP53.2.csv", sep = "\t") %>% as.data.frame()

####
mutations.8007 <- merge(DNMT3A.vaf[, c("bc", "mutated")], TP53_1.vaf[, c("bc", "mutated")], by = "bc", all = T)
colnames(mutations.8007) <- c("bc", "DNMT3A", "TP53.1")
mutations.8007 <- merge(mutations.8007, TP53_2.vaf[, c("bc", "mutated")], by = "bc", all = T)
colnames(mutations.8007)[4] <- "TP53.2"
rownames(mutations.8007) <- mutations.8007$bc
mutations.8007 <- mutations.8007[intersect(colnames(AML8007), mutations.8007$bc), ]
mutations.8007$predicted.celltype <- AML8007$predicted.celltype[mutations.8007$bc]
mutations.8007 <- mutations.8007[which(mutations.8007$predicted.celltype %in% names(which(table(mutations.8007$predicted.celltype) > 5))), ]

celltype <- factor(mutations.8007[which(grepl("AML8007.1", rownames(mutations.8007))), "predicted.celltype", ],
  levels = names(AML.combined.colors)[which(names(AML.combined.colors) %in%
    mutations.8007[which(grepl("AML8007.1", rownames(mutations.8007))), "predicted.celltype", ])]
)
ha <- columnAnnotation(celltype = celltype, col = list("celltype" = AML.combined.colors), border = T, simple_anno_size = unit(5, "pt"))
svglite::svglite("./AML/figures/AML8007/heatmaps/20220912_AML8007_1.svg", width = 5, height = 1.5)
ComplexHeatmap::Heatmap(
  t(mutations.8007[
    which(grepl("AML8007.1", rownames(mutations.8007))),
    c("DNMT3A", "TP53.1", "TP53.2")
  ]),
  show_column_names = F, col = c("mutated" = "firebrick", "wildtype" = "white"), na_col = "white",
  column_split = celltype, column_title_side = "bottom", column_title_rot = 90, column_title_gp = gpar(fontsize = 8),
  row_names_side = "left", row_names_gp = gpar(fontsize = 8),
  row_labels = c("DNMT3AV296M", "TP53C176S", "TP53R282W"),
  top_annotation = ha, border = T, raster_quality = 10, use_raster = T
)
dev.off()


celltype <- factor(mutations.8007[which(grepl("AML8007.4", rownames(mutations.8007))), "predicted.celltype", ],
  levels = names(AML.combined.colors)[which(names(AML.combined.colors) %in%
    mutations.8007[which(grepl("AML8007.4", rownames(mutations.8007))), "predicted.celltype", ])]
)
ha <- columnAnnotation(celltype = celltype, col = list("celltype" = AML.combined.colors), border = T, simple_anno_size = unit(5, "pt"))
svglite::svglite("./AML/figures/AML8007/heatmaps/20220912_AML8007_4.svg", width = 5, height = 1.5)
ComplexHeatmap::Heatmap(
  t(mutations.8007[
    which(grepl("AML8007.4", rownames(mutations.8007))),
    c("DNMT3A", "TP53.1", "TP53.2")
  ]),
  show_column_names = F, col = c("mutated" = "firebrick", "wildtype" = "white"), na_col = "white",
  column_split = celltype, column_title_side = "bottom", column_title_rot = 90, column_title_gp = gpar(fontsize = 8),
  row_names_side = "left", row_names_gp = gpar(fontsize = 8),
  row_labels = c("DNMT3AV296M", "TP53C176S", "TP53R282W"),
  top_annotation = ha, border = T, raster_quality = 10, use_raster = T
)
dev.off()

### statistics
mutations.8007$sample <- stringr::str_split_fixed(mutations.8007$bc, pattern = "_", n = 2)[, 1]
statistics.8007 <- mutations.8007 %>%
  group_by(sample, predicted.celltype) %>%
  summarize(
    DNMT3A.mut = length(which(DNMT3A == "mutated")),
    DNMT3A.wt = length(which(DNMT3A == "wildtype")),
    TP53.1.mut = length(which(TP53.1 == "mutated")),
    TP53.1.wt = length(which(TP53.1 == "wildtype")),
    TP53.2.mut = length(which(TP53.2 == "mutated")),
    TP53.2.wt = length(which(TP53.2 == "wildtype"))
  )

statistics.8007$DNMT3A.mut.freq <- statistics.8007$DNMT3A.mut / (statistics.8007$DNMT3A.mut + statistics.8007$DNMT3A.wt)
statistics.8007$TP53.1.mut.freq <- statistics.8007$TP53.1.mut / (statistics.8007$TP53.1.mut + statistics.8007$TP53.1.wt)
statistics.8007$TP53.2.mut.freq <- statistics.8007$TP53.2.mut / (statistics.8007$TP53.2.mut + statistics.8007$TP53.2.wt)

# DNMT3A
ggplot(statistics.8007[which(statistics.8007$predicted.celltype %in% c(myeloid.clusters, "Prog_Mk", "Prog_RBC")), ], aes(x = sample, y = 100 * DNMT3A.mut.freq, color = predicted.celltype)) +
  geom_point() +
  geom_line(aes(group = predicted.celltype)) +
  scale_x_discrete(labels = c("Screening", "marrow CR")) +
  scale_y_continuous("%DNMT3AV296M", limits = c(0, 100)) +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./AML/figures/AML8007/plots/20220912_AML8007_DNMT3A_myeloid_frequency.svg", width = 1.1, height = 2)

ggplot(statistics.8007[which(statistics.8007$predicted.celltype %in% TNK.clusters), ], aes(x = sample, y = 100 * DNMT3A.mut.freq, color = predicted.celltype)) +
  geom_point() +
  geom_line(aes(group = predicted.celltype)) +
  scale_x_discrete(labels = c("Screening", "marrow CR")) +
  scale_y_continuous("%DNMT3AV296M", limits = c(0, 100)) +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./AML/figures/AML8007/plots/20220912_AML8007_DNMT3A_TNK_frequency.svg", width = 1.1, height = 2)

# TP53.1
ggplot(statistics.8007[which(statistics.8007$predicted.celltype %in% c(myeloid.clusters, "Prog_Mk", "Prog_RBC")), ], aes(x = sample, y = 100 * TP53.1.mut.freq, color = predicted.celltype)) +
  geom_point() +
  geom_line(aes(group = predicted.celltype)) +
  scale_x_discrete(labels = c("Screening", "marrow CR")) +
  scale_y_continuous("%TP53C176S", limits = c(0, 100)) +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./AML/figures/AML8007/plots/20220912_AML8007_TP53.1_myeloid_frequency.svg", width = 1.1, height = 2)

ggplot(statistics.8007[which(statistics.8007$predicted.celltype %in% TNK.clusters), ], aes(x = sample, y = 100 * TP53.1.mut.freq, color = predicted.celltype)) +
  geom_point() +
  geom_line(aes(group = predicted.celltype)) +
  scale_x_discrete(labels = c("Screening", "marrow CR")) +
  scale_y_continuous("%TP53C176S", limits = c(0, 100)) +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./AML/figures/AML8007/plots/20220912_AML8007_TP53.1_TNK_frequency.svg", width = 1.1, height = 2)


# TP53.2
ggplot(statistics.8007[which(statistics.8007$predicted.celltype %in% c(myeloid.clusters, "Prog_Mk", "Prog_RBC")), ], aes(x = sample, y = 100 * TP53.2.mut.freq, color = predicted.celltype)) +
  geom_point() +
  geom_line(aes(group = predicted.celltype)) +
  scale_x_discrete(labels = c("Screening", "marrow CR")) +
  scale_y_continuous("%TP53R282W", limits = c(0, 100)) +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./AML/figures/AML8007/plots/20220912_AML8007_TP53.2_myeloid_frequency.svg", width = 1.1, height = 2)

ggplot(statistics.8007[which(statistics.8007$predicted.celltype %in% TNK.clusters), ], aes(x = sample, y = 100 * TP53.2.mut.freq, color = predicted.celltype)) +
  geom_point() +
  geom_line(aes(group = predicted.celltype)) +
  scale_x_discrete(labels = c("Screening", "marrow CR")) +
  scale_y_continuous("%TP53R282W", limits = c(0, 100)) +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./AML/figures/AML8007/plots/20220912_AML8007_TP53.2_TNK_frequency.svg", width = 1.1, height = 2)
