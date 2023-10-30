# QC plots

library(dplyr)
library(ggplot2)
library(nanoranger.R)

AML1022 <- readRDS("./data/benchmarking/objects/20230812_AML1022_nanoranger_GoT.rds")

DNMT3A.GoT.ONT <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022_pileup_DNMT3A.GoT.ONT.csv.gz")
DNMT3A.GoT.Illumina <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022_pileup_DNMT3A.GoT.Illumina.csv.gz")
DNMT3A.spikein <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022_pileup_DNMT3A.spikein.csv.gz")

DNMT3A.GoT.ONT <- DNMT3A.GoT.ONT %>%
  group_by(bc) %>%
  tally()
DNMT3A.GoT.Illumina <- DNMT3A.GoT.Illumina %>%
  mutate(bc = gsub(bc, pattern = "-1", replacement = "")) %>%
  group_by(bc) %>%
  tally()
DNMT3A.spikein <- DNMT3A.spikein %>%
  group_by(bc) %>%
  tally()

df <- merge(DNMT3A.GoT.ONT, DNMT3A.GoT.Illumina, by = "bc", all = T)
df <- merge(df, DNMT3A.spikein, by = "bc", all = T)

colnames(df) <- c("bc", "GoT.ONT", "GoT.Illumina", "spikein")
df$bc <- paste0("AML1022.1G_", df$bc, "-1")

df$high.qc <- ifelse(df$bc %in% colnames(AML1022), "yes", "no")
df[is.na(df)] <- 0
df[df == 0] <- jitter(rep(0.1, length(which(df == 0))), amount = 0.02)

ggplot(df, aes(x = GoT.ONT, y = spikein, color = high.qc)) +
  ggrastr::rasterize(geom_point(size = 0.5), dpi = 600) +
  scale_x_log10("reads GoT ONT") +
  scale_y_log10("reads nanoranger") +
  scale_color_manual(values = c("yes" = "firebrick", "no" = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230916_DNMT3A_reads_GoT_nanoranger.svg", width = 2.5, height = 2.5)

DNMT3A.spikein.processed <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_DNMT3A.spikein.processed.csv") %>%
  as.data.frame()
DNMT3A.GoT.ONT.processed <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_DNMT3A.GoT.ONT.processed.csv") %>%
  as.data.frame()

DNMT3A.GoT.ONT <- DNMT3A.GoT.ONT[which(paste0("AML1022.1G_", DNMT3A.GoT.ONT$bc, "-1")
%in% DNMT3A.GoT.ONT.processed$bc), ]
DNMT3A.GoT.ONT$rank <- rank(-DNMT3A.GoT.ONT$n, ties.method = "random")
DNMT3A.GoT.ONT$nanoranger <- ifelse(paste0("AML1022.1G_", DNMT3A.GoT.ONT$bc, "-1")
%in% DNMT3A.spikein.processed$bc, "yes", "no")

ggplot() +
  geom_line(data = DNMT3A.GoT.ONT, aes(x = rank, y = n), color = "grey90") +
  ggrastr::rasterize(geom_point(data = DNMT3A.GoT.ONT, aes(x = rank, y = n, color = nanoranger), size = 0.5), dpi = 600) +
  scale_x_log10("rank GoT ONT") +
  scale_y_log10("reads GoT ONT") +
  scale_color_manual(values = c("yes" = "orange", "no" = "#00aeef")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230916_DNMT3A_kneeplot_GOT_ONT.svg", width = 2.5, height = 2.5)

ggplot(DNMT3A.GoT.ONT, aes(x = nanoranger, y = n, color = nanoranger)) +
  ggrastr::rasterize(geom_jitter(size = 0.5), dpi = 600) +
  scale_x_discrete(labels = c("detected", "not detected"), limits = c("yes", "no")) +
  scale_y_log10("reads GoT ONT") +
  scale_color_manual(values = c("yes" = "orange", "no" = "#00aeef")) +
  geom_violin(color = "black") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230916_DNMT3A_violin_plot_GOT_ONT.svg", width = 2.5, height = 2.5)

### RUNX1
RUNX1.GoT.ONT <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022_pileup_RUNX1.GoT.ONT.csv.gz")
RUNX1.GoT.Illumina <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022_pileup_RUNX1.GoT.Illumina.csv.gz")
RUNX1.spikein <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022_pileup_RUNX1.spikein.csv.gz")

RUNX1.GoT.ONT <- RUNX1.GoT.ONT %>%
  group_by(bc) %>%
  tally()
RUNX1.GoT.Illumina <- RUNX1.GoT.Illumina %>%
  mutate(bc = gsub(bc, pattern = "-1", replacement = "")) %>%
  group_by(bc) %>%
  tally()
RUNX1.spikein <- RUNX1.spikein %>%
  group_by(bc) %>%
  tally()

df <- merge(RUNX1.GoT.ONT, RUNX1.GoT.Illumina, by = "bc", all = T)
df <- merge(df, RUNX1.spikein, by = "bc", all = T)

colnames(df) <- c("bc", "GoT.ONT", "GoT.Illumina", "spikein")
df$bc <- paste0("AML1022.1G_", df$bc, "-1")

df$high.qc <- ifelse(df$bc %in% colnames(AML1022), "yes", "no")
df[is.na(df)] <- 0
df[df == 0] <- jitter(rep(0.1, length(which(df == 0))), amount = 0.02)

ggplot(df, aes(x = GoT.ONT, y = spikein, color = high.qc)) +
  ggrastr::rasterize(geom_point(size = 0.5), dpi = 600) +
  scale_x_log10("reads GoT ONT", limits = c(0.08, 15000), breaks = c(1, 100, 10000), labels = c(1, 100, 10000)) +
  scale_y_log10("reads nanoranger") +
  scale_color_manual(values = c("yes" = "firebrick", "no" = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230916_RUNX1_reads_GoT_nanoranger.svg", width = 2.5, height = 2.5)


RUNX1.spikein.processed <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_RUNX1.spikein.processed.csv") %>%
  as.data.frame()
RUNX1.GoT.ONT.processed <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_RUNX1.GoT.ONT.processed.csv") %>%
  as.data.frame()


RUNX1.GoT.ONT <- RUNX1.GoT.ONT[which(paste0("AML1022.1G_", RUNX1.GoT.ONT$bc, "-1")
%in% RUNX1.GoT.ONT.processed$bc), ]
RUNX1.GoT.ONT$rank <- rank(-RUNX1.GoT.ONT$n, ties.method = "random")
RUNX1.GoT.ONT$nanoranger <- ifelse(paste0("AML1022.1G_", RUNX1.GoT.ONT$bc, "-1")
%in% RUNX1.spikein.processed$bc, "yes", "no")

ggplot() +
  geom_line(data = RUNX1.GoT.ONT, aes(x = rank, y = n), color = "grey90") +
  ggrastr::rasterize(geom_point(data = RUNX1.GoT.ONT, aes(x = rank, y = n, color = nanoranger), size = 0.5), dpi = 600) +
  scale_x_log10("rank GoT ONT") +
  scale_y_log10("reads GoT ONT") +
  scale_color_manual(values = c("yes" = "orange", "no" = "#00aeef")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230916_RUNX1_kneeplot_GOT_ONT.svg", width = 2.5, height = 2.5)

ggplot(RUNX1.GoT.ONT, aes(x = nanoranger, y = n, color = nanoranger)) +
  ggrastr::rasterize(geom_jitter(size = 0.5), dpi = 600) +
  scale_x_discrete(labels = c("detected", "not detected"), limits = c("yes", "no")) +
  scale_y_log10("reads GoT ONT") +
  scale_color_manual(values = c("yes" = "orange", "no" = "#00aeef")) +
  geom_violin(color = "black") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230916_RUNX1_violin_plot_GOT_ONT.svg", width = 2.5, height = 2.5)


#### SF3B1
SF3B1.GoT.ONT <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.GoT.ONT.csv.gz")
SF3B1.GoT.Illumina <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.GoT.Illumina.csv.gz")
SF3B1.spikein <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.spikein.csv.gz")

SF3B1.GoT.ONT <- SF3B1.GoT.ONT %>%
  group_by(bc) %>%
  tally()
SF3B1.GoT.Illumina <- SF3B1.GoT.Illumina %>%
  mutate(bc = gsub(bc, pattern = "-1", replacement = "")) %>%
  group_by(bc) %>%
  tally()
SF3B1.spikein <- SF3B1.spikein %>%
  group_by(bc) %>%
  tally()

df <- merge(SF3B1.GoT.ONT, SF3B1.GoT.Illumina, by = "bc", all = T)
df <- merge(df, SF3B1.spikein, by = "bc", all = T)

colnames(df) <- c("bc", "GoT.ONT", "GoT.Illumina", "spikein")
df$bc <- paste0("AML1022.1G_", df$bc, "-1")

df$high.qc <- ifelse(df$bc %in% colnames(AML1022), "yes", "no")
df[is.na(df)] <- 0
df[df == 0] <- jitter(rep(0.1, length(which(df == 0))), amount = 0.02)

ggplot(df, aes(x = GoT.ONT, y = spikein, color = high.qc)) +
  ggrastr::rasterize(geom_point(size = 0.5), dpi = 600) +
  scale_x_log10("reads GoT ONT", limits = c(0.08, 2050), breaks = c(1, 10, 100, 1000), labels = c(1, 10, 100, 1000)) +
  scale_y_log10("reads nanoranger") +
  scale_color_manual(values = c("yes" = "firebrick", "no" = "grey")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230916_SF3B1_reads_GoT_nanoranger.svg", width = 2.5, height = 2.5)


SF3B1.spikein.processed <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.spikein.processed.csv") %>%
  as.data.frame()
SF3B1.GoT.ONT.processed <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.GoT.ONT.processed.csv") %>%
  as.data.frame()

SF3B1.GoT.ONT <- SF3B1.GoT.ONT[which(paste0("AML1022.1G_", SF3B1.GoT.ONT$bc, "-1")
%in% SF3B1.GoT.ONT.processed$bc), ]
SF3B1.GoT.ONT$rank <- rank(-SF3B1.GoT.ONT$n, ties.method = "random")
SF3B1.GoT.ONT$nanoranger <- ifelse(paste0("AML1022.1G_", SF3B1.GoT.ONT$bc, "-1")
%in% SF3B1.spikein.processed$bc, "yes", "no")

ggplot() +
  geom_line(data = SF3B1.GoT.ONT, aes(x = rank, y = n), color = "grey90") +
  ggrastr::rasterize(geom_point(data = SF3B1.GoT.ONT, aes(x = rank, y = n, color = nanoranger), size = 0.5), dpi = 600) +
  scale_x_log10("rank GoT ONT") +
  scale_y_log10("reads GoT ONT") +
  scale_color_manual(values = c("yes" = "orange", "no" = "#00aeef")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230916_SF3B1_kneeplot_GOT_ONT.svg", width = 2.5, height = 2.5)

ggplot(SF3B1.GoT.ONT, aes(x = nanoranger, y = n, color = nanoranger)) +
  ggrastr::rasterize(geom_jitter(size = 0.5), dpi = 600) +
  scale_x_discrete(labels = c("detected", "not detected"), limits = c("yes", "no")) +
  scale_y_log10("reads GoT ONT") +
  scale_color_manual(values = c("yes" = "orange", "no" = "#00aeef")) +
  geom_violin(color = "black") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230916_SF3B1_violin_plot_GOT_ONT.svg", width = 2.5, height = 2.5)
