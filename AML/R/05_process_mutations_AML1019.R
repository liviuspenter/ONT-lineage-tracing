# process mutations for AML1019.1

library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(nanoranger.R)
library(Seurat)
library(SummarizedExperiment)
library(dplyr)

ASXL1 = as.data.frame(extract_mutation(BC.data.file = './data/AML/pileup/AML1019_pileup_ASXL1.csv.gz', REF = 'C', ALT = 'T', FILTER = 10))
ASXL1$bc = paste0('AML1019.1_',ASXL1$bc,'-1')
rownames(ASXL1) = ASXL1$bc
write.table(file = './data/AML/mutations/AML1019_ASXL1.csv', x = ASXL1, quote = F, sep = '\t')

RUNX1 = as.data.frame(extract_mutation(BC.data.file = './data/AML/pileup/AML1019_pileup_RUNX1.csv.gz', REF = 'G', ALT = 'A', FILTER = 10))
RUNX1$bc = paste0('AML1019.1_',RUNX1$bc,'-1')
rownames(RUNX1) = RUNX1$bc
write.table(file = './data/AML/mutations/AML1019_RUNX1.csv', x = RUNX1, quote = F, sep = '\t')


