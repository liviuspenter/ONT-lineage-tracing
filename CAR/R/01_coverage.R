# plot coverage of CD28 and CD247 amplicon

library(ggplot2)
library(Seurat)

technology.colors <- c("Both" = "#283350", "ONT" = "#f93800", "Illumina" = "#ffb500")

so <- readRDS(file = "./data/20220819_CAR_so.rds")

CAR.coverage <- data.table::fread("./data/CAR/CAR_coverage_CD28.csv")
CAR.coverage$transcript <- ifelse(CAR.coverage$V4 < 1310, "CAR", "CD28")
ggplot(CAR.coverage, aes(x = V4, y = 100 * V5 / max(V5), fill = transcript)) +
  ggrastr::rasterize(geom_col(), dpi = 600) +
  scale_x_continuous("position") +
  scale_y_continuous("% coverage") +
  scale_fill_manual(values = c("CD28" = "#00aeef", "CAR" = "firebrick")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./CAR/figures/plots/20220819_coverage_CD28.svg", width = 2, height = 2)

CAR.coverage <- data.table::fread("./data/CAR/CAR_coverage_CD247.csv")
CAR.coverage$transcript <- ifelse(CAR.coverage$V4 > 1132, "CD247", "CAR")
CAR.coverage$transcript[which(CAR.coverage$V4 %in% seq(810, 1131))] <- "CD28"
ggplot(CAR.coverage, aes(x = V4, y = 100 * V5 / max(V5), fill = transcript)) +
  ggrastr::rasterize(geom_col(), dpi = 600) +
  scale_x_continuous("position") +
  scale_y_continuous("% coverage") +
  scale_fill_manual(values = c("CD28" = "darkgreen", "CAR" = "firebrick", "CD247" = "black")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./CAR/figures/plots/20220822_coverage_CD247.svg", width = 3, height = 2)
