# compare minimal read length of cell barcodes found with Illumina or ONT or both

library(dplyr)
library(ggplot2)
library(nanoranger.R)

AML1022 <- readRDS("./data/benchmarking/objects/20230812_AML1022_nanoranger_GoT.rds")

DNMT3A.GoT.ONT <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022_pileup_DNMT3A.GoT.ONT.csv.gz")
DNMT3A.GoT.ONT <- DNMT3A.GoT.ONT %>%
  group_by(bc) %>%
  reframe(
    n = length(bc),
    read.len = mean(length),
    min.read.len = min(length),
    max.read.len = max(length)
  )
DNMT3A.GoT.ONT$bc <- paste0("AML1022.1G_", DNMT3A.GoT.ONT$bc, "-1")
DNMT3A.GoT.ONT <- DNMT3A.GoT.ONT[which(DNMT3A.GoT.ONT$bc %in% colnames(AML1022)), ]

DNMT3A.GoT.Illumina <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022_pileup_DNMT3A.GoT.Illumina.csv.gz")
DNMT3A.GoT.Illumina <- DNMT3A.GoT.Illumina %>%
  group_by(bc) %>%
  tally() %>%
  as.data.frame()
DNMT3A.GoT.Illumina$bc <- paste0("AML1022.1G_", DNMT3A.GoT.Illumina$bc)
DNMT3A.GoT.Illumina <- DNMT3A.GoT.Illumina[which(DNMT3A.GoT.Illumina$bc %in% colnames(AML1022)), ]

DNMT3A.GoT.ONT$overlap <- ifelse(DNMT3A.GoT.ONT$bc %in% DNMT3A.GoT.Illumina$bc, "yes", "no")

ggplot(DNMT3A.GoT.ONT, aes(x = min.read.len, y = n, color = overlap)) +
  geom_point(size = 0.5) +
  scale_x_continuous("minimal read length", limits = c(0, max(DNMT3A.GoT.ONT$min.read.len))) +
  scale_y_log10("reads GoT ONT") +
  scale_color_manual(values = c("no" = "#00aeef", "yes" = "black")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230815_DNMT3A_minimal_read_length_ONT_Illumina.svg", width = 1.8, height = 1.8)


RUNX1.GoT.ONT <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022_pileup_RUNX1.GoT.ONT.csv.gz")
RUNX1.GoT.ONT <- RUNX1.GoT.ONT %>%
  group_by(bc) %>%
  reframe(
    n = length(bc),
    read.len = mean(length),
    min.read.len = min(length),
    max.read.len = max(length)
  )
RUNX1.GoT.ONT$bc <- paste0("AML1022.1G_", RUNX1.GoT.ONT$bc, "-1")
RUNX1.GoT.ONT <- RUNX1.GoT.ONT[which(RUNX1.GoT.ONT$bc %in% colnames(AML1022)), ]

RUNX1.GoT.Illumina <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022_pileup_RUNX1.GoT.Illumina.csv.gz")
RUNX1.GoT.Illumina <- RUNX1.GoT.Illumina %>%
  group_by(bc) %>%
  tally() %>%
  as.data.frame()
RUNX1.GoT.Illumina$bc <- paste0("AML1022.1G_", RUNX1.GoT.Illumina$bc)
RUNX1.GoT.Illumina <- RUNX1.GoT.Illumina[which(RUNX1.GoT.Illumina$bc %in% colnames(AML1022)), ]

RUNX1.GoT.ONT$overlap <- ifelse(RUNX1.GoT.ONT$bc %in% RUNX1.GoT.Illumina$bc, "yes", "no")

ggplot(RUNX1.GoT.ONT, aes(x = min.read.len, y = n, color = overlap)) +
  geom_point(size = 0.5) +
  scale_x_continuous("minimal read length", limits = c(0, max(RUNX1.GoT.ONT$min.read.len))) +
  scale_y_log10("reads GoT ONT") +
  scale_color_manual(values = c("no" = "#00aeef", "yes" = "black")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230815_RUNX1_minimal_read_length_ONT_Illumina.svg", width = 1.8, height = 1.8)


SF3B1.GoT.ONT <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.GoT.ONT.csv.gz")
SF3B1.GoT.ONT <- SF3B1.GoT.ONT %>%
  group_by(bc) %>%
  reframe(
    n = length(bc),
    read.len = mean(length),
    min.read.len = min(length),
    max.read.len = max(length)
  )
SF3B1.GoT.ONT$bc <- paste0("AML1022.1G_", SF3B1.GoT.ONT$bc, "-1")
SF3B1.GoT.ONT <- SF3B1.GoT.ONT[which(SF3B1.GoT.ONT$bc %in% colnames(AML1022)), ]

SF3B1.GoT.Illumina <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.GoT.Illumina.csv.gz")
SF3B1.GoT.Illumina <- SF3B1.GoT.Illumina %>%
  group_by(bc) %>%
  tally() %>%
  as.data.frame()
SF3B1.GoT.Illumina$bc <- paste0("AML1022.1G_", SF3B1.GoT.Illumina$bc)
SF3B1.GoT.Illumina <- SF3B1.GoT.Illumina[which(SF3B1.GoT.Illumina$bc %in% colnames(AML1022)), ]

SF3B1.GoT.ONT$overlap <- ifelse(SF3B1.GoT.ONT$bc %in% SF3B1.GoT.Illumina$bc, "yes", "no")

ggplot(SF3B1.GoT.ONT, aes(x = min.read.len, y = n, color = overlap)) +
  geom_point(size = 0.5) +
  scale_x_continuous("minimal read length", limits = c(0, max(SF3B1.GoT.ONT$min.read.len))) +
  scale_y_log10("reads GoT ONT") +
  scale_color_manual(values = c("no" = "#00aeef", "yes" = "black")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230815_SF3B1_minimal_read_length_ONT_Illumina.svg", width = 1.8, height = 1.8)
