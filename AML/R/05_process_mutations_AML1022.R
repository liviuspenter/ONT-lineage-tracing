# process mutations for AML1022.1

library(ComplexHeatmap)
library(ggplot2)
library(nanoranger.R)
library(Seurat)
library(SummarizedExperiment)

AML1022.1 = readRDS('./data/AML/objects/20220729_AML1022_1.rds')
AML1022.genotype = read.table('./data/AML/souporcell/20210707_AML1022_chimerism.csv')
rownames(AML1022.genotype) = AML1022.genotype$barcode

DNMT3A.vaf = extract_mutation(BC.data.file = './data/AML/pileup/AML1022_1_DNMT3A.csv.gz', ALT = 'T', REF = 'C', FILTER = 10)
DNMT3A.vaf$bc = paste0('AML1022.1_', DNMT3A.vaf$bc, '-1')
write.table(file = './data/AML/mutations/AML1022_DNMT3A.csv', x = DNMT3A.vaf, quote = F, sep = '\t')

DNMT3A.vaf = DNMT3A.vaf[which(DNMT3A.vaf$bc %in% colnames(AML1022.1)),]
DNMT3A.vaf$chimerism = AML1022.genotype[DNMT3A.vaf$bc, 'chimerism']

RUNX1.vaf = extract_mutation(BC.data.file = './data/AML/pileup/AML1022_1_RUNX1.csv.gz', ALT = 'C', REF = 'A', FILTER = 10)
RUNX1.vaf$bc = paste0('AML1022.1_', RUNX1.vaf$bc, '-1')
write.table(file = './data/AML/mutations/AML1022_RUNX1.csv', x = RUNX1.vaf, quote = F, sep = '\t')
RUNX1.vaf = RUNX1.vaf[which(RUNX1.vaf$bc %in% colnames(AML1022.1)),]
RUNX1.vaf$chimerism = AML1022.genotype[RUNX1.vaf$bc, 'chimerism']

SF3B1.vaf = extract_mutation(BC.data.file = './data/AML/pileup/AML1022_1_SF3B1.csv.gz', ALT = 'C', REF = 'T', FILTER = 10)
SF3B1.vaf$bc = paste0('AML1022.1_', SF3B1.vaf$bc, '-1')
write.table(file = './data/AML/mutations/AML1022_SF3B1.csv', x = SF3B1.vaf, quote = F, sep = '\t')
SF3B1.vaf = SF3B1.vaf[which(SF3B1.vaf$bc %in% colnames(AML1022.1)),]
SF3B1.vaf$chimerism = AML1022.genotype[SF3B1.vaf$bc, 'chimerism']

