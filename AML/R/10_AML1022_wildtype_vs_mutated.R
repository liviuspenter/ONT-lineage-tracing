# analysis of mutations in AML1022.1 versus donor-recipient genotype 

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(nanoranger.R)
library(Seurat)
library(SummarizedExperiment)

AML1022.1 = readRDS('./data/AML/objects/20220729_AML1022_1.rds')
AML1022.genotype = read.table('./data/AML/souporcell/20210707_AML1022_chimerism.csv')
rownames(AML1022.genotype) = AML1022.genotype$barcode

DNMT3A.vaf = read.csv2('./data/AML/mutations/AML1022_DNMT3A.csv', sep = '\t')
RUNX1.vaf = read.csv2('./data/AML/mutations/AML1022_RUNX1.csv', sep = '\t')
SF3B1.vaf = read.csv2('./data/AML/mutations/AML1022_SF3B1.csv', sep = '\t')

DNMT3A.vaf = DNMT3A.vaf[which(DNMT3A.vaf$bc %in% colnames(AML1022.1)),]
DNMT3A.vaf$chimerism = AML1022.genotype[DNMT3A.vaf$bc, 'chimerism']

RUNX1.vaf = RUNX1.vaf[which(RUNX1.vaf$bc %in% colnames(AML1022.1)),]
RUNX1.vaf$chimerism = AML1022.genotype[RUNX1.vaf$bc, 'chimerism']

SF3B1.vaf = SF3B1.vaf[which(SF3B1.vaf$bc %in% colnames(AML1022.1)),]
SF3B1.vaf$chimerism = AML1022.genotype[SF3B1.vaf$bc, 'chimerism']


length(which((DNMT3A.vaf$alt + DNMT3A.vaf$ref) == 1)) / nrow(DNMT3A.vaf)
length(which((RUNX1.vaf$alt + RUNX1.vaf$ref) == 1)) / nrow(RUNX1.vaf)
length(which((SF3B1.vaf$alt + SF3B1.vaf$ref) == 1)) / nrow(SF3B1.vaf)

DNMT3A.mutated.cells = DNMT3A.vaf$bc[which(DNMT3A.vaf$mutated == 'mutated')]
DNMT3A.wildtype.cells = DNMT3A.vaf$bc[which(DNMT3A.vaf$mutated == 'wildtype')]
DNMT3A.vaf$predicted.celltype = AML1022.1$predicted.celltype[DNMT3A.vaf$bc]

RUNX1.mutated.cells = RUNX1.vaf$bc[which(RUNX1.vaf$mutated == 'mutated')]
RUNX1.wildtype.cells = RUNX1.vaf$bc[which(RUNX1.vaf$mutated == 'wildtype')]
RUNX1.vaf$predicted.celltype = AML1022.1$predicted.celltype[RUNX1.vaf$bc]

SF3B1.mutated.cells = SF3B1.vaf$bc[which(SF3B1.vaf$mutated == 'mutated')]
SF3B1.wildtype.cells = SF3B1.vaf$bc[which(SF3B1.vaf$mutated == 'wildtype')]
SF3B1.vaf$predicted.celltype = AML1022.1$predicted.celltype[SF3B1.vaf$bc]

p=DimPlot(AML1022.1, reduction = 'ref.umap', cells.highlight = list('mutated' = DNMT3A.mutated.cells, 'wildtype' = DNMT3A.wildtype.cells),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'orange'), sizes.highlight = 0.5, pt.size = 0.5) +
  NoAxes() + NoLegend()
ggsave('./AML/figures/AML1022/UMAP/20220730_AML1022_1_DNMT3A.png', width = 2, height = 2, plot = p)

p=DimPlot(AML1022.1, reduction = 'ref.umap', cells.highlight = list('mutated' = DNMT3A.mutated.cells, 'wildtype' = RUNX1.wildtype.cells),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 0.5, pt.size = 0.5) +
  NoAxes() + NoLegend()
ggsave('./AML/figures/AML1022/UMAP/20220730_AML1022_1_RUNX1.png', width = 2, height = 2, plot = p)

p=DimPlot(AML1022.1, reduction = 'ref.umap', cells.highlight = list('mutated' = DNMT3A.mutated.cells, 'wildtype' = SF3B1.wildtype.cells),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'blue'), sizes.highlight = 0.5, pt.size = 0.5) +
  NoAxes() + NoLegend()
ggsave('./AML/figures/AML1022/UMAP/20220730_AML1022_1_SF3B1.png', width = 2, height = 2, plot = p)

boo=DNMT3A.vaf %>% group_by(predicted.celltype) %>% summarize(mutated = length(which(mutated == 'mutated')),
                                                              recipient = length(which(chimerism == 'recipient')),
                                                              donor = length(which(chimerism == 'donor')))
ggplot(boo, aes(x=100*mutated/(recipient + donor), y=100-100*recipient/(recipient + donor), fill=predicted.celltype)) + 
  geom_point(aes(size=donor+recipient), color='black', shape=21, alpha=0.8) +
  scale_x_continuous('% DNMT3AR882H', limits = c(0,100)) +
  scale_y_continuous('% donor', limits = c(0,100)) +
  scale_fill_manual(values = AML.combined.colors, guide = 'none') +
  scale_size_continuous('# cells') + 
  theme_classic() +
  theme(legend.position = 'bottom',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./AML/figures/AML1022/plots/20220729_AML1022_1_DNMT3A_VAF_chimerism.svg', width = 2, height = 2)

boo=RUNX1.vaf %>% group_by(predicted.celltype) %>% summarize(mutated = length(which(mutated == 'mutated')),
                                                             recipient = length(which(chimerism == 'recipient')),
                                                             donor = length(which(chimerism == 'donor')))
ggplot(boo, aes(x=100*mutated/(recipient + donor), y=100-100*recipient/(recipient + donor), fill=predicted.celltype)) + 
  geom_point(aes(size=donor+recipient), color='black', shape=21) +
  scale_x_continuous('% RUNX1I177S', limits = c(0,100)) +
  scale_y_continuous('% donor', limits = c(0,100)) +
  scale_fill_manual(values = AML.combined.colors, guide = 'none') +
  scale_size_continuous('# cells') + 
  theme_classic() +
  theme(legend.position = 'bottom',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./AML/figures/AML1022/plots/20220729_AML1022_1_RUNX1_VAF_chimerism.svg', width = 2, height = 2)

boo=SF3B1.vaf %>% group_by(predicted.celltype) %>% summarize(mutated = length(which(mutated == 'mutated')),
                                                             recipient = length(which(chimerism == 'recipient')),
                                                             donor = length(which(chimerism == 'donor')))
ggplot(boo, aes(x=100*mutated/(recipient + donor), y=100-100*recipient/(recipient + donor), fill=predicted.celltype)) + 
  geom_point(aes(size=donor+recipient), color='black', shape=21) +
  scale_x_continuous('% SF3B1K700E', limits = c(0,100)) +
  scale_y_continuous('% donor', limits = c(0,100)) +
  scale_fill_manual(values = AML.combined.colors, guide = 'none') +
  scale_size_continuous('# cells') + 
  theme_classic() +
  theme(legend.position = 'bottom',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./AML/figures/AML1022/plots/20220729_AML1022_1_SF3B1_VAF_chimerism.svg', width = 2, height = 2)

ggplot() + 
  geom_histogram(data=SF3B1.vaf, aes(x=100*as.numeric(vaf)), fill='blue', alpha=1) +
  geom_histogram(data=RUNX1.vaf, aes(x=100*as.numeric(vaf)), fill='firebrick', alpha=1) +
  geom_histogram(data=DNMT3A.vaf, aes(x=100*as.numeric(vaf)), fill='orange', alpha=1) +
  scale_x_continuous('% VAF') + 
  theme_classic() +
  theme(legend.position = 'bottom',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./AML/figures/AML1022/plots/20220730_AML1022_1_distribution_vaf.svg', width = 2, height = 1.5)

ggplot() + 
  geom_histogram(data=SF3B1.vaf, aes(x=alt+ref), fill='blue') +
  geom_histogram(data=RUNX1.vaf, aes(x=alt+ref), fill='firebrick') +
  geom_histogram(data=DNMT3A.vaf, aes(x=alt+ref), fill='orange') +
  scale_x_continuous('UMI per cell barcode', breaks = c(0,1,5,10)) + 
  scale_y_continuous('cell barcodes') +
  theme_classic() +
  theme(legend.position = 'bottom',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./AML/figures/AML1022/plots/20220730_AML1022_1_distribution_cbc.svg', width = 2, height = 1.5)
