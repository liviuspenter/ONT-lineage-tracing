# process mutations for AML3003.14

library(ComplexHeatmap)
library(ggplot2)
library(nanoranger.R)
library(Seurat)
library(SummarizedExperiment)

AML3003.14 = readRDS('./data/AML/objects/202208_AML3003_14.rds')

# ASXL1 mutation
ASXL1.vaf.1 = extract_indel(BC.data.file = './data/AML/pileup/AML3003_1_pileup_ASXL1.csv.gz', REF = 'T', CONSENSUS = 1)
ASXL1.vaf.1$bc = paste0('AML3003.1_', ASXL1.vaf.1$bc, '-1')
ASXL1.vaf.1$sample = 'AML3003.1'
ASXL1.vaf.4 = extract_indel(BC.data.file = './data/AML/pileup/AML3003_4_pileup_ASXL1.csv.gz', REF = 'T', CONSENSUS = 1)
ASXL1.vaf.4$bc = paste0('AML3003.4_', ASXL1.vaf.4$bc, '-1')
ASXL1.vaf.4$sample = 'AML3003.4'
ASXL1.vaf = rbind(ASXL1.vaf.1, ASXL1.vaf.4)
write.table(file = './data/AML/mutations/AML3003_ASXL1.csv', x = ASXL1.vaf, quote = F, sep = '\t')

# STAG2 mutation
STAG2.vaf.1 = extract_mutation(BC.data.file = './data/AML/pileup/AML3003_1_pileup_STAG2.csv.gz', REF = 'C', ALT = 'T', FILTER = 10)
STAG2.vaf.1$bc = paste0('AML3003.1_', STAG2.vaf.1$bc, '-1')
STAG2.vaf.1$sample = 'AML3003.1'
STAG2.vaf.4 = extract_mutation(BC.data.file = './data/AML/pileup/AML3003_4_pileup_STAG2.csv.gz', REF = 'C', ALT = 'T', FILTER = 10)
STAG2.vaf.4$bc = paste0('AML3003.4_', STAG2.vaf.4$bc, '-1')
STAG2.vaf.4$sample = 'AML3003.4'
STAG2.vaf = rbind(STAG2.vaf.1, STAG2.vaf.4)
write.table(file = './data/AML/mutations/AML3003_STAG2.csv', x = STAG2.vaf, quote = F, sep = '\t')

p=DimPlot(subset(AML3003.14, orig.ident == 'AML3003.1'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = STAG2.vaf$bc[which(STAG2.vaf$mutated == 'mutated')], 
                                 'wildtype' = STAG2.vaf$bc[which(STAG2.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
q=DimPlot(subset(AML3003.14, orig.ident == 'AML3003.4'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = STAG2.vaf$bc[which(STAG2.vaf$mutated == 'mutated')], 
                                 'wildtype' = STAG2.vaf$bc[which(STAG2.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
ggsave('./AML/figures/AML3003/UMAP/AML3003_1_STAG2.png', width = 3, height = 3, dpi = 600, plot = p)
ggsave('./AML/figures/AML3003/UMAP/AML3003_4_STAG2.png', width = 3, height = 3, dpi = 600, plot = q)



# this mutation doesn't exist but was described in the patient's history
TET2.vaf.1 = extract_indel(BC.data.file = './data/AML/pileup//AML3003_1_pileup_TET2.csv.gz', REF = 'G', CONSENSUS = 1)
TET2.vaf.1$bc = paste0('AML3003.1_', TET2.vaf.1$bc, '-1')
TET2.vaf.1$sample = 'AML3003.1'
TET2.vaf.4 = extract_indel(BC.data.file = './data/AML/pileup//AML3003_4_pileup_TET2.csv.gz', REF = 'G', CONSENSUS = 1)
TET2.vaf.4$bc = paste0('AML3003.4_', TET2.vaf.4$bc, '-1')
TET2.vaf.4$sample = 'AML3003.4'

DimPlot(AML3003.14, reduction = 'ref.umap', cells.highlight = TET2.vaf.1$bc[which(TET2.vaf.1$mutated == 'mutated')])
DimPlot(AML3003.14, reduction = 'ref.umap', cells.highlight = TET2.vaf.4$bc[which(TET2.vaf.4$mutated == 'wildtype')])

# TET2.2 mutation
TET2.2.vaf.1 = extract_length_diff(BC.data.file = './data/AML/pileup/AML3003_1_pileup_TET2.2_length_difference.csv.gz', POS = 105272892)
TET2.2.vaf.1$bc = paste0('AML3003.1_', TET2.2.vaf.1$bc, '-1')
TET2.2.vaf.1$sample = 'AML3003.1'
TET2.2.vaf.4 = extract_length_diff(BC.data.file = './data/AML/pileup/AML3003_4_pileup_TET2.2_length_difference.csv.gz', POS = 105272892)
TET2.2.vaf.4$bc = paste0('AML3003.4_', TET2.2.vaf.4$bc, '-1')
TET2.2.vaf.4$sample = 'AML3003.4'
TET2.2.vaf = rbind(TET2.2.vaf.1, TET2.2.vaf.4)
write.table(file = './data/AML/mutations/AML3003_TET2.2.csv', x = TET2.2.vaf, quote = F, sep = '\t')

p=DimPlot(subset(AML3003.14, orig.ident == 'AML3003.1'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = TET2.2.vaf$bc[which(TET2.2.vaf$mutated == 'mutated')], 
                                 'wildtype' = TET2.2.vaf$bc[which(TET2.2.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
q=DimPlot(subset(AML3003.14, orig.ident == 'AML3003.4'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = TET2.2.vaf$bc[which(TET2.2.vaf$mutated == 'mutated')], 
                                 'wildtype' = TET2.2.vaf$bc[which(TET2.2.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
ggsave('./AML/figures/AML3003/UMAP/AML3003_1_TET2.2.png', width = 3, height = 3, dpi = 600, plot = p)
ggsave('./AML/figures/AML3003/UMAP/AML3003_4_TET2.2.png', width = 3, height = 3, dpi = 600, plot = q)


# TET2.3 mutation
TET2.3.vaf.1 = extract_indel(BC.data.file = './data/AML/pileup/AML3003_1_pileup_TET2.3.csv.gz', REF = 'G', FILTER = 10)
TET2.3.vaf.1$bc = paste0('AML3003.1_', TET2.3.vaf.1$bc, '-1')
TET2.3.vaf.1$sample = 'AML3003.1'
TET2.3.vaf.4 = extract_indel(BC.data.file = './data/AML/pileup/AML3003_4_pileup_TET2.3.csv.gz', REF = 'G', FILTER = 10)
TET2.3.vaf.4$bc = paste0('AML3003.4_', TET2.3.vaf.4$bc, '-1')
TET2.3.vaf.4$sample = 'AML3003.4'
TET2.3.vaf = rbind(TET2.3.vaf.1, TET2.3.vaf.4)
write.table(file = './data/AML/mutations/AML3003_TET2.3.csv', x = TET2.3.vaf, quote = F, sep = '\t')

p=DimPlot(subset(AML3003.14, orig.ident == 'AML3003.1'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = TET2.3.vaf$bc[which(TET2.3.vaf$mutated == 'mutated')], 
                                 'wildtype' = TET2.3.vaf$bc[which(TET2.3.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
q=DimPlot(subset(AML3003.14, orig.ident == 'AML3003.4'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = TET2.3.vaf$bc[which(TET2.3.vaf$mutated == 'mutated')], 
                                 'wildtype' = TET2.3.vaf$bc[which(TET2.3.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
ggsave('./AML/figures/AML3003/UMAP/AML3003_1_TET2.3.png', width = 3, height = 3, dpi = 600, plot = p)
ggsave('./AML/figures/AML3003/UMAP/AML3003_4_TET2.3.png', width = 3, height = 3, dpi = 600, plot = q)


# TET2.4 mutation
TET2.4.vaf.1 = extract_mutation(BC.data.file = './data/AML/pileup/AML3003_1_pileup_TET2.4.csv.gz', REF = 'A', ALT = 'G', FILTER = 10)
TET2.4.vaf.1$bc = paste0('AML3003.1_', TET2.4.vaf.1$bc, '-1')
TET2.4.vaf.1$sample = 'AML3003.1'
TET2.4.vaf.4 = extract_mutation(BC.data.file = './data/AML/pileup/AML3003_4_pileup_TET2.4.csv.gz', REF = 'A', ALT = 'G', FILTER = 10)
TET2.4.vaf.4$bc = paste0('AML3003.4_', TET2.4.vaf.4$bc, '-1')
TET2.4.vaf.4$sample = 'AML3003.4'
TET2.4.vaf = rbind(TET2.4.vaf.1, TET2.4.vaf.4)
write.table(file = './data/AML/mutations/AML3003_TET2.4.csv', x = TET2.4.vaf, quote = F, sep = '\t')

p=DimPlot(subset(AML3003.14, orig.ident == 'AML3003.1'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = TET2.4.vaf$bc[which(TET2.4.vaf$mutated == 'mutated')], 
                                 'wildtype' = TET2.4.vaf$bc[which(TET2.4.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 2) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
q=DimPlot(subset(AML3003.14, orig.ident == 'AML3003.4'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = TET2.4.vaf$bc[which(TET2.4.vaf$mutated == 'mutated')], 
                                 'wildtype' = TET2.4.vaf$bc[which(TET2.4.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 2) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
ggsave('./AML/figures/AML3003/UMAP/AML3003_1_TET2.4.png', width = 3, height = 3, dpi = 600, plot = p)
ggsave('./AML/figures/AML3003/UMAP/AML3003_4_TET2.4.png', width = 3, height = 3, dpi = 600, plot = q)


# U2AF1 mutation not found
U2AF1.vaf.1 = extract_mutation(BC.data.file = './data/AML/pileup/AML3003_1_pileup_U2AF1.csv.gz', REF = 'G', ALT = 'A')
U2AF1.vaf.1$bc = paste0('AML3003.1_', U2AF1.vaf.1$bc, '-1')
U2AF1.vaf.1$sample = 'AML3003.1'
U2AF1.vaf.4 = extract_mutation(BC.data.file = './data/AML/pileup/AML3003_4_pileup_U2AF1.csv.gz', REF = 'G', ALT = 'A')
U2AF1.vaf.4$bc = paste0('AML3003.4_', U2AF1.vaf.4$bc, '-1')
U2AF1.vaf.4$sample = 'AML3003.4'
U2AF1.vaf = rbind(U2AF1.vaf.1, U2AF1.vaf.4)
#write.table(file = './data/AML/mutations/AML3003_U2AF1.csv', x = U2AF1.vaf, quote = F, sep = '\t')

DimPlot(AML3003.14, reduction = 'ref.umap', cells.highlight = U2AF1.vaf.1$bc[which(U2AF1.vaf.1$mutated == 'wildtype')])
DimPlot(AML3003.14, reduction = 'ref.umap', cells.highlight = U2AF1.vaf.4$bc[which(U2AF1.vaf.4$mutated == 'wildtype')])
