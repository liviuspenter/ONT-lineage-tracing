library(ComplexHeatmap)
library(ggplot2)
library(nanoranger.R)
library(Seurat)

FLT3_ITD <- readRDS("./data/AML/objects/20230523_FLT3-ITD.rds")

### de-novo AML1 (FLT3_ITD1)
FLT3_ITD1 <- extract_indel("./data/AML/pileup/FLT3_ITD1_pileup_FLT3.csv.gz", REF = "T", CONSENSUS = 24)
FLT3_ITD1$bc <- paste0("FLT3-ITD1_", FLT3_ITD1$bc, "-1")
FLT3_ITD1$sample <- "FLT3_ITD1"
write.table(FLT3_ITD1, "./data/AML/mutations/FLT3-ITD1_FLT3-ITD.csv", quote = F, sep = "\t")
FLT3_ITD1.mutated <- FLT3_ITD1$bc[which(FLT3_ITD1$mutated == "mutated")]

FLT3_ITD1_NPM1_1 <- extract_indel("./data/AML/pileup/FLT3_ITD1_pileup_NPM1_1.csv.gz", REF = "C", CONSENSUS = 4)
FLT3_ITD1_NPM1_1$bc <- paste0("FLT3-ITD1_", FLT3_ITD1_NPM1_1$bc, "-1")
FLT3_ITD1_NPM1_1$sample <- "FLT3_ITD1"
write.table(FLT3_ITD1_NPM1_1, "./data/AML/mutations/FLT3-ITD1_NPM1.csv", quote = F, sep = "\t")
FLT3_ITD1_NPM1_1.mutated <- FLT3_ITD1_NPM1_1$bc[which(FLT3_ITD1_NPM1_1$mutated == "mutated")]

FLT3_ITD1_NPM1_2 <- extract_indel("./data/AML/pileup/FLT3_ITD1_pileup_NPM1_2.csv.gz", REF = "C", CONSENSUS = 4)
FLT3_ITD1_NPM1_2$bc <- paste0("FLT3-ITD1_", FLT3_ITD1_NPM1_2$bc, "-1")
FLT3_ITD1_NPM1_2$sample <- "FLT3_ITD1"
write.table(FLT3_ITD1_NPM1_2, "./data/AML/mutations/FLT3-ITD1_NPM1.2.csv", quote = F, sep = "\t")

FLT3_ITD1_NPM1_2$predicted.celltype <- FLT3_ITD$predicted.celltype[FLT3_ITD1_NPM1_2$bc]
FLT3_ITD1_NPM1_2.mutated <- FLT3_ITD1_NPM1_2$bc[which(FLT3_ITD1_NPM1_2$mutated == "mutated")]

### de-novo AML3 (FLT3_ITD2)
FLT3_ITD2 <- extract_indel("./data/AML/pileup/FLT3_ITD2_pileup_FLT3.csv.gz", REF = "T", CONSENSUS = 21)
FLT3_ITD2$bc <- paste0("FLT3-ITD2_", FLT3_ITD2$bc, "-1")
FLT3_ITD2$sample <- "FLT3_ITD2"
write.table(FLT3_ITD2, "./data/AML/mutations/FLT3-ITD2_FLT3-ITD.csv", quote = F, sep = "\t")

FLT3_ITD2.mutated <- FLT3_ITD2$bc[which(FLT3_ITD2$mutated == "mutated")]

### de-novo AML2 (FLT3_ITD3)
FLT3_ITD3 <- extract_indel("./data/AML/pileup/FLT3_ITD3_pileup_FLT3.csv.gz", REF = "T", CONSENSUS = 24)
FLT3_ITD3$bc <- paste0("FLT3-ITD3_", FLT3_ITD3$bc, "-1")
FLT3_ITD3$sample <- "FLT3_ITD3"
write.table(FLT3_ITD3, "./data/AML/mutations/FLT3-ITD3_FLT3-ITD.csv", quote = F, sep = "\t")

FLT3_ITD3.mutated <- FLT3_ITD3$bc[which(FLT3_ITD3$mutated == "mutated")]

FLT3_ITD3_NPM1_1 <- extract_indel("./data/AML/pileup/FLT3_ITD3_pileup_NPM1_1.csv.gz", REF = "C", CONSENSUS = 4)
FLT3_ITD3_NPM1_1$bc <- paste0("FLT3-ITD3_", FLT3_ITD3_NPM1_1$bc, "-1")
FLT3_ITD3_NPM1_1$sample <- "FLT3_ITD3"
write.table(FLT3_ITD3_NPM1_1, "./data/AML/mutations/FLT3-ITD3_NPM1.csv", quote = F, sep = "\t")

FLT3_ITD3_NPM1_1.mutated <- FLT3_ITD3_NPM1_1$bc[which(FLT3_ITD3_NPM1_1$mutated == "mutated")]

FLT3_ITD3_DNMT3A <- extract_mutation("./data/AML/pileup/FLT3_ITD3_pileup_DNMT3A.csv.gz", REF = "C", ALT = "T")
FLT3_ITD3_DNMT3A$bc <- paste0("FLT3-ITD3_", FLT3_ITD3_DNMT3A$bc, "-1")
FLT3_ITD3_DNMT3A$sample <- "FLT3_ITD3"
write.table(FLT3_ITD3_DNMT3A, "./data/AML/mutations/FLT3-ITD3_DNMT3A.csv", quote = F, sep = "\t")

FLT3_ITD3_DNMT3A.mutated <- FLT3_ITD3_DNMT3A$bc[which(FLT3_ITD3_DNMT3A$mutated == "mutated")]
