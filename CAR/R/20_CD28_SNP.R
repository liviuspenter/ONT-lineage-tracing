# compare CD28 SNP versus CAR/wildtype CD28

library(ggplot2)
library(nanoranger.R)
library(Seurat)

### compare raw reads at level of cell barcode and UMI
CAR.reads <- data.table::fread("./data/CAR/97_6_read_length.csv.gz")
SNP.reads <- data.table::fread("./data/CAR/97_6_pileup_Yescarta.csv.gz")

CAR.reads.reduced <- CAR.reads %>%
  group_by(bc, umi) %>%
  summarize(ref.start = mean(ref_start)) %>%
  as.data.frame()
CAR.reads.reduced$bc.umi <- paste0(CAR.reads.reduced$bc, ".", CAR.reads.reduced$umi)

SNP.reads.reduced <- SNP.reads %>%
  group_by(bc, umi) %>%
  summarize(read = names(which.max(table(base)))) %>%
  as.data.frame()
SNP.reads.reduced$bc.umi <- paste0(SNP.reads.reduced$bc, ".", SNP.reads.reduced$umi)

df <- merge(CAR.reads.reduced[, c("bc.umi", "ref.start")],
  SNP.reads.reduced[, c("bc.umi", "read")],
  by = "bc.umi", all = F
)

ggplot() +
  geom_histogram(data = df[which(df$read == "G"), ], aes(x = 1572 - ref.start), fill = "firebrick") +
  geom_histogram(data = df[which(df$read == "T"), ], aes(x = 1572 - ref.start), fill = "#00aeef") +
  scale_x_continuous("length of transcript", limits = c(0, 1600), breaks = c(0, 500, 1000, 1500)) +
  scale_y_continuous("transcripts") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./CAR/figures/plots/20230807_CAR_transcript_length_vs_SNP.svg", width = 2, height = 2)

### compare at the level of cells

# CAR versus CD28 wildtype
CAR.ONT <- extract_fusion_gene("./data/CAR/97_6_CART_CD28.csv.gz", WILDTYPE = "CD28", FUSION = "CARTmod")
CAR.ONT$bc <- paste0(CAR.ONT$bc, "-1")

# CD28 SNP
CD28.SNP <- extract_mutation("./data/CAR/97_6_pileup_Yescarta.csv.gz", ALT = "T", REF = "G")
CD28.SNP$bc <- paste0(CD28.SNP$bc, "-1")
write.table(CD28.SNP, file = "./data/CAR/97_6_pileup_Yescarta.processed.csv", sep = "\t", quote = F)

df <- merge(CAR.ONT[, c("bc", "CARTmod", "CD28", "mutated")],
  CD28.SNP[, c("bc", "alt", "ref", "mutated")],
  by = "bc", all = F
)

colnames(df) <- c("bc", "CARTmod", "CD28", "CAR", "alt", "ref", "CD28.SNP")

# CAR not detected and CD28 germline SNP: 948
length(which(df$CAR == "wildtype" & df$CD28.SNP == "mutated"))
# CAR detected and CD28 CAR SNP: 2,699
length(which(df$CAR == "mutated" & df$CD28.SNP == "wildtype"))
# CAR not detected but CD28 CAR SNP: 16
length(which(df$CAR == "wildtype" & df$CD28.SNP == "wildtype"))
# CAR detected but CD28 germline SNP: 379
length(which(df$CAR == "mutated" & df$CD28.SNP == "mutated"))

# wildtype CD28 versus CD28 germline SNP
ggplot(df, aes(x = CD28, y = alt)) +
  geom_abline(slope = 1) +
  ggpointdensity::geom_pointdensity(size = 0.5) +
  scale_color_viridis_c() +
  scale_x_continuous("wildtype CD28 transcripts", limits = c(0, 25)) +
  scale_y_continuous("transcripts with germline CD28 SNP", limits = c(0, 25)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./CAR/figures/plots/20230807_CD28_wildtype_vs_CD28_germline_SNP.svg", width = 2, height = 3)

# CAR versus CD28 CAR SNP
ggplot(df, aes(x = CARTmod, y = ref)) +
  geom_abline(slope = 1) +
  ggpointdensity::geom_pointdensity(size = 0.5) +
  scale_color_viridis_c() +
  scale_x_continuous("CAR transcripts", limits = c(0, 100)) +
  scale_y_continuous("transcripts with CAR CD28 SNP", limits = c(0, 100)) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./CAR/figures/plots/20230807_CAR_transcript_vs_CD28_CAR_SNP.svg", width = 2, height = 3)
