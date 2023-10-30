# process mutations for AML1010.1

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(nanoranger.R)
library(Seurat)
library(SummarizedExperiment)

AML1010.1 <- readRDS("./data/AML/objects/20220711_AML1010_1.rds")
AML1010.genotype <- read.table("./data/AML/souporcell/20210707_AML1010_chimerism.csv")
rownames(AML1010.genotype) <- AML1010.genotype$barcode

NRAS.vaf <- extract_mutation(BC.data.file = "./data/AML/pileup/AML1010_1_pileup_NRAS.csv.gz", ALT = "T", REF = "C")
NRAS.vaf$bc <- paste0("AML1010.1_", NRAS.vaf$bc, "-1")
NRAS.vaf <- NRAS.vaf[which(NRAS.vaf$bc %in% colnames(AML1010.1)), ]
NRAS.vaf$chimerism <- AML1010.genotype[NRAS.vaf$bc, "chimerism"]
NRAS.vaf$predicted.celltype <- AML1010.1$predicted.celltype[NRAS.vaf$bc]
write.table(NRAS.vaf, file = "./data/AML/mutations/AML1010_NRAS.csv", sep = "\t", quote = F)

NRAS.mutated.cells <- NRAS.vaf$bc[which(NRAS.vaf$mutated == "mutated")]
NRAS.wildtype.cells <- NRAS.vaf$bc[which(NRAS.vaf$mutated == "wildtype")]

NRAS.statistics <- NRAS.vaf %>%
  group_by(predicted.celltype) %>%
  summarize(
    NRAS.mut = length(which(mutated == "mutated")),
    NRAS.wt = length(which(mutated == "wildtype"))
  )
NRAS.statistics$freq <- NRAS.statistics$NRAS.mut / (NRAS.statistics$NRAS.mut + NRAS.statistics$NRAS.wt)
NRAS.statistics$predicted.celltype <- factor(NRAS.statistics$predicted.celltype, levels = names(AML.combined.colors))

ggplot(
  NRAS.statistics[which(NRAS.statistics$predicted.celltype %in% c("HSC", "LMPP", "GMP", "CD14 Mono", "CD16 Mono", "cDC2", "Prog_DC")), ],
  aes(x = predicted.celltype, y = 100 * freq)
) +
  geom_col(aes(fill = predicted.celltype), color = "black") +
  scale_y_continuous("% NRASG12D") +
  scale_fill_manual(values = AML.combined.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank()
  )
ggsave("./AML/figures/AML1010/plots/20220803_AML1010_1_NRAS_vaf.svg", width = 1.5, height = 2)
