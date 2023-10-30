# plot BCRABL1 coverage for different libraries

library(dplyr)
library(ggplot2)

BCR.exon.start <- c(145, 1876, 2058, 2163, 2349, 2457, 2518, 2571, 2712, 2834, 3003, 3123, 3199, 3304)
BCR.exon.end <- c(1874, 2057, 2162, 2348, 2456, 2517, 2570, 2711, 2833, 3002, 3122, 3198, 3303, 3378)
annotation.df <- data.frame(start = BCR.exon.start, end = BCR.exon.end, name = paste("exon", seq(1, 14)))

ALL1.coverage <- data.table::fread("./data/ALL/coverage/ALL1_BCRABL1_coverage.csv")
ALL1.coverage$sample <- "ALL1"
ALL2.coverage <- data.table::fread("./data/ALL/coverage/ALL2_BCRABL1_coverage.csv")
ALL2.coverage$sample <- "ALL2"
K562.coverage <- data.table::fread("./data/ALL/coverage/K562_BCRABL1_coverage.csv")
K562.coverage$sample <- "K562"
CML_E3.coverage <- data.table::fread("./data/ALL/coverage/CML_E3_BCRABL1_coverage.csv")
CML_E3.coverage$sample <- "CML_E3"

BCRABL1.coverage <- dplyr::bind_rows(ALL1.coverage, ALL2.coverage, K562.coverage, CML_E3.coverage)
colnames(BCRABL1.coverage) <- c("gene", "start", "end", "pos", "coverage", "sample")

BCRABL1.coverage <- BCRABL1.coverage %>%
  group_by(sample, gene) %>%
  mutate(coverage = coverage / max(coverage))

cov <- ggplot() +
  geom_rect(data = annotation.df, aes(xmin = start, xmax = end, ymin = 0, ymax = 100, fill = name)) +
  scale_fill_manual(values = BuenColors::jdb_palette(name = "corona", n = 14)) +
  scale_x_continuous(limits = c(0, 3400)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )
p <- ggplot(BCRABL1.coverage[which(BCRABL1.coverage$sample == "ALL1" & grepl("BCR", BCRABL1.coverage$gene)), ], aes(x = pos, y = 100 * as.numeric(coverage))) +
  ggrastr::rasterize(geom_col(fill = "black"), dpi = 600) +
  scale_x_continuous(limits = c(0, 3400)) +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )
q <- ggplot(BCRABL1.coverage[which(BCRABL1.coverage$sample == "ALL2" & grepl("BCR", BCRABL1.coverage$gene)), ], aes(x = pos, y = 100 * as.numeric(coverage))) +
  ggrastr::rasterize(geom_col(fill = "firebrick"), dpi = 600) +
  scale_x_continuous(limits = c(0, 3400)) +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )
r <- ggplot(BCRABL1.coverage[which(BCRABL1.coverage$sample == "K562" & grepl("BCR", BCRABL1.coverage$gene)), ], aes(x = pos, y = 100 * as.numeric(coverage))) +
  ggrastr::rasterize(geom_col(fill = "blue"), dpi = 600) +
  scale_x_continuous(limits = c(0, 3400)) +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  )
s <- ggplot(BCRABL1.coverage[which(BCRABL1.coverage$sample == "CML_E3" & grepl("BCR", BCRABL1.coverage$gene)), ], aes(x = pos, y = 100 * as.numeric(coverage))) +
  ggrastr::rasterize(geom_col(fill = "grey"), dpi = 600) +
  scale_x_continuous("position on BCR mRNA", limits = c(0, 3400)) +
  scale_y_continuous(breaks = c(0, 50, 100)) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_blank(),
    axis.text.x = element_blank()
  )

svglite::svglite("./ALL/figures/coverage/20221007_coverage_BCRABL1.svg", width = 2, height = 3)
cowplot::plot_grid(plotlist = list(cov, p, q, r, s), ncol = 1)
dev.off()
