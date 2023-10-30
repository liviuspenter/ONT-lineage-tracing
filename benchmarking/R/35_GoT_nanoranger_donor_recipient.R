# compare mutations calls with donor/recipient annotation

library(dplyr)
library(ggplot2)
library(Seurat)

AML1022 <- readRDS("./data/benchmarking/objects/20230812_AML1022_nanoranger_GoT.rds")

AML1022.1.souporcell <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022.1/souporcell/clusters.tsv") %>%
  as.data.frame()
AML1022.1.souporcell$barcode <- paste0("AML1022.1_", AML1022.1.souporcell$barcode)
AML1022.1.spikein.souporcell <- data.table::fread("./data/benchmarking/GoT_nanoranger/AML1022.1_GoT_with_amplicon/souporcell/clusters.tsv") %>%
  as.data.frame()
AML1022.1.spikein.souporcell$barcode <- paste0("AML1022.1G_", AML1022.1.spikein.souporcell$barcode)

# AML1022.1 - donor
DimPlot(AML1022,
  reduction = "ref.umap",
  cells.highlight = AML1022.1.souporcell$barcode[which(AML1022.1.souporcell$assignment == "0")]
)
# AML1022.1 - recipient
DimPlot(AML1022,
  reduction = "ref.umap",
  cells.highlight = AML1022.1.souporcell$barcode[which(AML1022.1.souporcell$assignment == "1")]
)

# AML1022.1G - recipient
DimPlot(AML1022,
  reduction = "ref.umap",
  cells.highlight = AML1022.1.spikein.souporcell$barcode[which(AML1022.1.spikein.souporcell$assignment == "0")]
)
# AML1022.1G - donor
DimPlot(AML1022,
  reduction = "ref.umap",
  cells.highlight = AML1022.1.spikein.souporcell$barcode[which(AML1022.1.spikein.souporcell$assignment == "1")]
)

AML1022.1.souporcell$chimerism <- "none"
AML1022.1.souporcell$chimerism[which(AML1022.1.souporcell$assignment == "0")] <- "donor"
AML1022.1.souporcell$chimerism[which(AML1022.1.souporcell$assignment == "1")] <- "recipient"
AML1022.1.spikein.souporcell$chimerism <- "none"
AML1022.1.spikein.souporcell$chimerism[which(AML1022.1.spikein.souporcell$assignment == "0")] <- "recipient"
AML1022.1.spikein.souporcell$chimerism[which(AML1022.1.spikein.souporcell$assignment == "1")] <- "donor"
rownames(AML1022.1.souporcell) <- AML1022.1.souporcell$barcode
rownames(AML1022.1.spikein.souporcell) <- AML1022.1.spikein.souporcell$barcode
souporcell.data <- rbind(AML1022.1.souporcell, AML1022.1.spikein.souporcell)

write.table(
  stringr::str_split_fixed(AML1022.1.souporcell$barcode[which(AML1022.1.souporcell$chimerism == "recipient")],
    pattern = "_", n = 2
  )[, 2],
  file = "./data/benchmarking/GoT_nanoranger/AML1022.1/barcodes.recipient.tsv", quote = F, sep = "\t",
  row.names = F, col.names = F
)

write.table(
  stringr::str_split_fixed(AML1022.1.spikein.souporcell$barcode[which(AML1022.1.spikein.souporcell$chimerism == "recipient")],
    pattern = "_", n = 2
  )[, 2],
  file = "./data/benchmarking/GoT_nanoranger/AML1022.1_GoT_with_amplicon/barcodes.recipient.tsv", quote = F, sep = "\t",
  row.names = F, col.names = F
)

write.table(AML1022.1.souporcell,
  file = "./data/benchmarking/GoT_nanoranger/AML1022.1/souporcell/AML1022.1.souporcell.processed.csv",
  sep = "\t", quote = F
)
write.table(AML1022.1.spikein.souporcell,
  file = "./data/benchmarking/GoT_nanoranger/AML1022.1_GoT_with_amplicon/souporcell/AML1022.1G.souporcell.processed.csv",
  sep = "\t", quote = F
)


### correlate with SF3B1 genotyping data
SF3B1.normal <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.normal.processed.csv") %>%
  as.data.frame()
SF3B1.spikein <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.spikein.processed.csv") %>%
  as.data.frame()
SF3B1.spikein.additive <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.spikein.additive.processed.csv") %>%
  as.data.frame()
SF3B1.GoT.Illumina <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.GoT.Illumina.processed.csv") %>%
  as.data.frame()
SF3B1.GoT.ONT <- read.table("./data/benchmarking/GoT_nanoranger/AML1022_pileup_SF3B1.GoT.ONT.processed.csv") %>%
  as.data.frame()

SF3B1.normal$condition <- "normal"
SF3B1.spikein$condition <- "spikein"
SF3B1.spikein.additive$condition <- "spikein.additive"
SF3B1.GoT.ONT$condition <- "GoT.ONT"
SF3B1.GoT.Illumina$condition <- "GoT.Illumina"

SF3B1.data <- dplyr::bind_rows(SF3B1.GoT.ONT, SF3B1.GoT.Illumina, SF3B1.normal, SF3B1.spikein, SF3B1.spikein.additive) %>%
  filter(bc %in% colnames(AML1022)) %>%
  filter(!is.na(vaf))

SF3B1.data$predicted.celltype <- AML1022$predicted.celltype[SF3B1.data$bc]
SF3B1.data$chimerism <- souporcell.data[SF3B1.data$bc, "chimerism"]


ggplot(
  SF3B1.data[which(SF3B1.data$condition == "normal" &
    SF3B1.data$chimerism %in% c("recipient", "donor")), ],
  aes(x = chimerism, y = 100 * vaf)
) +
  scale_y_continuous("% apparent VAF") +
  geom_jitter(size = 0.5, color = "orange") +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230823_AML1022_SF3B1_donor_recipient_normal.svg", width = 1.2, height = 2)

ggplot(
  SF3B1.data[which(SF3B1.data$condition == "spikein" &
    SF3B1.data$chimerism %in% c("recipient", "donor")), ],
  aes(x = chimerism, y = 100 * vaf)
) +
  scale_y_continuous("% apparent VAF") +
  geom_jitter(size = 0.5, color = "red") +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230823_AML1022_SF3B1_donor_recipient_spikein.svg", width = 1.2, height = 2)

ggplot(
  SF3B1.data[which(SF3B1.data$condition == "spikein.additive" &
    SF3B1.data$chimerism %in% c("recipient", "donor")), ],
  aes(x = chimerism, y = 100 * vaf)
) +
  scale_y_continuous("% apparent VAF") +
  geom_jitter(size = 0.5, color = "firebrick") +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230823_AML1022_SF3B1_donor_recipient_spikein.additive.svg", width = 1.2, height = 2)

ggplot(
  SF3B1.data[which(SF3B1.data$condition == "GoT.Illumina" &
    SF3B1.data$chimerism %in% c("recipient", "donor")), ],
  aes(x = chimerism, y = 100 * vaf)
) +
  scale_y_continuous("% apparent VAF") +
  geom_jitter(size = 0.5, color = "black") +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230823_AML1022_SF3B1_donor_recipient_GoT.Illumina.svg", width = 1.2, height = 2)

ggplot(
  SF3B1.data[which(SF3B1.data$condition == "GoT.ONT" &
    SF3B1.data$chimerism %in% c("recipient", "donor")), ],
  aes(x = chimerism, y = 100 * vaf)
) +
  scale_y_continuous("% apparent VAF") +
  geom_jitter(size = 0.5, color = "#00aeef") +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
ggsave("./benchmarking/figures/GoT_nanoranger/plots/20230823_AML1022_SF3B1_donor_recipient_GoT.ONT.svg", width = 1.2, height = 2)

View(SF3B1.data %>% filter(chimerism %in% c("recipient", "donor")) %>% group_by(condition, chimerism, mutated) %>% tally())
