# process mutations for AML1026.2

library(ComplexHeatmap)
library(ggplot2)
library(nanoranger.R)
library(Seurat)
library(SummarizedExperiment)

AML1026.2 = readRDS('./data/AML/objects/20220711_AML1026_2.rds')

# NRAS
NRAS.vaf2 = extract_mutation(BC.data.file = './data/AML/pileup/AML1026_2_NRAS.csv.gz', ALT = 'G', REF = 'C', FILTER = 10)
NRAS.vaf2$bc = paste0('AML1026.2_', NRAS.vaf2$bc, '-1')
write.table(file = './data/AML/mutations/AML1026_NRAS.csv', x = NRAS.vaf2, quote = F, sep = '\t')

NRAS.vaf2 = NRAS.vaf2[which(NRAS.vaf2$bc %in% colnames(AML1026.2)),]
NRAS.mutated.cells.2 = NRAS.vaf2$bc[which(NRAS.vaf2$mutated == 'mutated')]
NRAS.wildtype.cells.2 = NRAS.vaf2$bc[which(NRAS.vaf2$mutated == 'wildtype')]
NRAS.vaf2$predicted.celltype = AML1026.2$predicted.celltype[NRAS.vaf2$bc]

# SF3B1
SF3B1.vaf2 = extract_mutation(BC.data.file = './data/AML/pileup/AML1026_2_SF3B1.csv.gz', ALT = 'T', REF = 'C', FILTER = 10)
SF3B1.vaf2$bc = paste0('AML1026.2_', SF3B1.vaf2$bc, '-1')
write.table(file = './data/AML/mutations/AML1026_SF3B1.csv', x = SF3B1.vaf2, quote = F, sep = '\t')

SF3B1.vaf2 = SF3B1.vaf2[which(SF3B1.vaf2$bc %in% colnames(AML1026.2)),]
SF3B1.mutated.cells.2 = SF3B1.vaf2$bc[which(SF3B1.vaf2$mutated == 'mutated')]
SF3B1.wildtype.cells.2 = SF3B1.vaf2$bc[which(SF3B1.vaf2$mutated == 'wildtype')]
SF3B1.vaf2$predicted.celltype = AML1026.2$predicted.celltype[SF3B1.vaf2$bc]

### SRSF2
SRSF2.vaf2 = extract_mutation('./data/Mehdi/AML1026_2_mito_ASXL1_NRAS_SF3B1_SRSF2/AML1026_2_SRSF2.csv.gz', 
                         ALT = 'T', REF = 'G')
SRSF2.vaf2$bc = paste0('AML1026.2_',SRSF2.vaf2$bc, '-1')
write.table(SRSF2.vaf2, file = './data/AML/mutations/AML1026_SRSF2.csv', quote = F, sep = '\t')

SRSF2.vaf2 = read.csv2('./data/AML/mutations/AML1026_SRSF2.csv', dec = ',', row.names = 1) %>% as.data.frame()

SRSF2.vaf2 = SRSF2.vaf2[which(SRSF2.vaf2$bc %in% colnames(AML1026.2)),]
SRSF2.mutated.cells.2 = SRSF2.vaf2$bc[which(SRSF2.vaf2$mutated == 'mutated')]
SRSF2.wildtype.cells.2 = SRSF2.vaf2$bc[which(SRSF2.vaf2$mutated == 'wildtype')]
SRSF2.vaf2$predicted.celltype = AML1026.2$predicted.celltype[SRSF2.vaf2$bc]

# write output for mtDNA analysis
write.csv2(SF3B1.vaf2, file = './data/AML/mtDNA/20220805_AML1026_2_SF3B1.csv', quote = F)
write.csv2(NRAS.vaf2, file = './data/AML/mtDNA/20220805_AML1026_2_NRAS.csv', quote = F)
write.csv2(SRSF2.vaf2, file = './data/AML/mtDNA/20220805_AML1026_2_SRSF2.csv', quote = F)

