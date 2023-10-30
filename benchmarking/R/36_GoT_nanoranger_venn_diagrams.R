# visualize overlap of cell barcodes between conditions

library(dplyr)
library(ggplot2)
library(Seurat)
library(ggVennDiagram)

AML1022 <- readRDS("./data/benchmarking/objects/20230812_AML1022_nanoranger_GoT.rds")

DNMT3A.normal <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_DNMT3A.normal.processed.csv") %>%
  as.data.frame()
DNMT3A.spikein <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_DNMT3A.spikein.processed.csv") %>%
  as.data.frame()
DNMT3A.GoT.Illumina <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_DNMT3A.GoT.Illumina.processed.csv") %>%
  as.data.frame()
DNMT3A.GoT.ONT <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_DNMT3A.GoT.ONT.processed.csv") %>%
  as.data.frame()

DNMT3A.normal$condition <- "normal"
DNMT3A.spikein$condition <- "spikein"
DNMT3A.GoT.ONT$condition <- "GoT.ONT"
DNMT3A.GoT.Illumina$condition <- "GoT.Illumina"

DNMT3A.data <- dplyr::bind_rows(DNMT3A.GoT.ONT, DNMT3A.GoT.Illumina, DNMT3A.normal, DNMT3A.spikein) %>%
  filter(bc %in% colnames(AML1022)) %>%
  filter(!is.na(vaf))

ggVennDiagram(
  x = list(
    DNMT3A.data$bc[which(DNMT3A.data$condition == "GoT.Illumina")],
    DNMT3A.data$bc[which(DNMT3A.data$condition == "spikein")],
    DNMT3A.data$bc[which(DNMT3A.data$condition == "GoT.ONT")]
  ),
  category.names = c("GoT.Illumina", "nanoranger", "GoT.ONT"),
  set_color = c("red", "black", "grey"), label_size = 2, set_size = 2, label = "count",
) +
  scale_fill_distiller(palette = "PRGn") +
  scale_color_manual(values = c("spikein" = "red", "GoT.Illumina" = "black", "GoT.ONT" = "grey"))
ggsave("./benchmarking/figures/GoT_nanoranger/venndiagrams/20230823_DNMT3A_venn.svg", width = 2, height = 2)



RUNX1.normal <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_RUNX1.normal.processed.csv") %>%
  as.data.frame()
RUNX1.spikein <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_RUNX1.spikein.processed.csv") %>%
  as.data.frame()
RUNX1.GoT.Illumina <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_RUNX1.GoT.Illumina.processed.csv") %>%
  as.data.frame()
RUNX1.GoT.ONT <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_RUNX1.GoT.ONT.processed.csv") %>%
  as.data.frame()

RUNX1.normal$condition <- "normal"
RUNX1.spikein$condition <- "spikein"
RUNX1.GoT.ONT$condition <- "GoT.ONT"
RUNX1.GoT.Illumina$condition <- "GoT.Illumina"

RUNX1.data <- dplyr::bind_rows(RUNX1.GoT.ONT, RUNX1.GoT.Illumina, RUNX1.normal, RUNX1.spikein) %>%
  filter(bc %in% colnames(AML1022)) %>%
  filter(!is.na(vaf))

ggVennDiagram(
  x = list(
    RUNX1.data$bc[which(RUNX1.data$condition == "GoT.Illumina")],
    RUNX1.data$bc[which(RUNX1.data$condition == "spikein")],
    RUNX1.data$bc[which(RUNX1.data$condition == "GoT.ONT")]
  ),
  category.names = c("GoT.Illumina", "nanoranger", "GoT.ONT"),
  set_color = c("red", "black", "grey"), label_size = 2, set_size = 2, label = "count",
) +
  scale_fill_distiller(palette = "PRGn") +
  scale_color_manual(values = c("spikein" = "red", "GoT.Illumina" = "black", "GoT.ONT" = "grey"))
ggsave("./benchmarking/figures/GoT_nanoranger/venndiagrams/20230823_RUNX1_venn.svg", width = 2, height = 2)



SF3B1.normal <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.normal.processed.csv") %>%
  as.data.frame()
SF3B1.spikein <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.spikein.processed.csv") %>%
  as.data.frame()
SF3B1.GoT.Illumina <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.GoT.Illumina.processed.csv") %>%
  as.data.frame()
SF3B1.GoT.ONT <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.GoT.ONT.processed.csv") %>%
  as.data.frame()

SF3B1.normal$condition <- "normal"
SF3B1.spikein$condition <- "spikein"
SF3B1.GoT.ONT$condition <- "GoT.ONT"
SF3B1.GoT.Illumina$condition <- "GoT.Illumina"

SF3B1.data <- dplyr::bind_rows(SF3B1.GoT.ONT, SF3B1.GoT.Illumina, SF3B1.normal, SF3B1.spikein) %>%
  filter(bc %in% colnames(AML1022)) %>%
  filter(!is.na(vaf))

ggVennDiagram(
  x = list(
    SF3B1.data$bc[which(SF3B1.data$condition == "GoT.Illumina")],
    SF3B1.data$bc[which(SF3B1.data$condition == "spikein")],
    SF3B1.data$bc[which(SF3B1.data$condition == "GoT.ONT")]
  ),
  category.names = c("GoT.Illumina", "nanoranger", "GoT.ONT"),
  set_color = c("red", "black", "grey"), label_size = 2, set_size = 2, label = "count",
) +
  scale_fill_distiller(palette = "PRGn") +
  scale_color_manual(values = c("spikein" = "red", "GoT.Illumina" = "black", "GoT.ONT" = "grey"))
ggsave("./benchmarking/figures/GoT_nanoranger/venndiagrams/20230823_SF3B1_venn.svg", width = 2, height = 2)
