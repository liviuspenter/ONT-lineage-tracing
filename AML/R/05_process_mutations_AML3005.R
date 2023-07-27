# process mutations for AML3005.13

library(ComplexHeatmap)
library(ggplot2)
library(nanoranger.R)
library(Seurat)
library(SummarizedExperiment)

AML3005.13 = readRDS('./data/AML/objects/20220718_AML3005_13.rds')

cluster.colors = c('Progenitor' = 'firebrick', 'Erythropoiesis' = 'grey', 'B lymphopoiesis' = 'darkgreen', 
                   'NK cell' = 'black', 'CD8 T cell' = 'darkblue', 'CD4 T cell' = 'lightblue')

p=DimPlot(AML3005.13, reduction = 'umap', group.by = 'manual.cluster', cols = cluster.colors) +
  NoLegend() + 
  NoAxes() +
  theme(plot.title = element_blank())
ggsave('./AML/figures/AML3005/UMAP/20220719_AML3005_13_UMAP_celltypes.png', width = 4, height = 4, dpi = 600, plot = p)

Idents(AML3005.13) = 'manual.cluster'
AML3005.13$manual.cluster = factor(AML3005.13$manual.cluster, 
                                   levels = c('Progenitor', 'Erythropoiesis', 'B lymphopoiesis', 'NK cell', 'CD8 T cell', 'CD4 T cell'))
AML.df = as.data.frame(prop.table(table(Idents(AML3005.13), AML3005.13$orig.ident), margin = 2))
AML.df$Var1 = factor(AML.df$Var1, levels = rev(names(cluster.colors)))
ggplot(AML.df, aes(x=Var2, y=100*Freq, fill=Var1)) + 
  geom_col() + 
  scale_fill_manual(values = cluster.colors) +
  scale_x_discrete(labels = c('Screening', 'CR')) + 
  scale_y_continuous('% cells') +
  geom_hline(yintercept = 5) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML3005/plots/20220719_AML3005_13_clusters_dynamics.svg', width = 1.3, height = 2)

# DNMT3A mutation
DNMT3A.vaf = extract_indel(BC.data.file = './data/AML/pileup/AML3005_1_pileup_DNMT3A.csv.gz', REF = 'G', FILTER = 10)
DNMT3A.vaf$bc = paste0('AML3005.1_', DNMT3A.vaf$bc, '-1')
DNMT3A.vaf$sample = 'AML3005.1'

DNMT3A.vaf2 = extract_indel(BC.data.file = './data/AML/pileup/AML3005_3_pileup_DNMT3A.csv.gz', REF = 'G', FILTER = 10)
DNMT3A.vaf2$bc = paste0('AML3005.3_', DNMT3A.vaf2$bc, '-1')
DNMT3A.vaf2$sample = 'AML3005.3'

DNMT3A.combined = rbind(DNMT3A.vaf, DNMT3A.vaf2)
write.table(DNMT3A.combined, file = './data/AML/mutations/AML3005_DNMT3A.csv', sep = '\t', quote = F)

p=DimPlot(AML3005.13, reduction = 'umap', 
          cells.highlight = list('mutated' = DNMT3A.combined$bc[which(DNMT3A.combined$mutated == 'mutated')], 
                                 'wildtype' = DNMT3A.combined$bc[which(DNMT3A.combined$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 2) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
ggsave('./AML/figures/AML3005/UMAP/20220719_AML3005_13_UMAP_DNMT3A.png', width = 4, height = 4, dpi = 600, plot = p)

DNMT3A.combined$predicted.celltype = AML3005.13$predicted.celltype[DNMT3A.combined$bc]
DNMT3A.combined$manual.cluster = AML3005.13$manual.cluster[DNMT3A.combined$bc]

DNMT3A.condensed = DNMT3A.combined %>% filter(bc %in% colnames(AML3005.13)) %>%
  group_by(manual.cluster) %>% 
  summarize(mutated.1 = length(which(mutated == 'mutated' & sample == 'AML3005.1')),
            wildtype.1 = length(which(mutated == 'wildtype' & sample == 'AML3005.1')),
            mutated.3 = length(which(mutated == 'mutated' & sample == 'AML3005.3')),
            wildtype.3 = length(which(mutated == 'wildtype' & sample == 'AML3005.3')))
DNMT3A.condensed$mutated.freq.1 = DNMT3A.condensed$mutated.1 / (DNMT3A.condensed$mutated.1 + DNMT3A.condensed$wildtype.1)
DNMT3A.condensed$mutated.freq.3 = DNMT3A.condensed$mutated.3 / (DNMT3A.condensed$mutated.3 + DNMT3A.condensed$wildtype.3)

ggplot(data=reshape2::melt(DNMT3A.condensed), aes(x=variable, y=100*value, color=manual.cluster)) +
  geom_point(size=0.5) +
  geom_line(aes(group=manual.cluster)) + 
  scale_x_discrete(limits = c('mutated.freq.1', 'mutated.freq.3'), labels = c('Screening', 'CR')) +
  scale_y_continuous(limits = c(0,100), '% cells') +
  scale_color_manual(values = cluster.colors) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML3005/plots/20220719_AML3005_13_DNMT3A_dynamics.svg', width = 1.1, height = 2)


# U2AF1 mutation
U2AF1.vaf = extract_mutation(BC.data.file = './data/AML/pileup/AML3005_1_pileup_U2AF1.csv.gz', ALT = 'T', REF = 'G', FILTER = 10)
U2AF1.vaf$bc = paste0('AML3005.1_', U2AF1.vaf$bc, '-1')
U2AF1.vaf$sample = 'AML3005.1'

U2AF1.vaf2 = extract_mutation(BC.data.file = './data/AML/pileup/AML3005_3_pileup_U2AF1.csv.gz', ALT = 'T', REF = 'G', FILTER = 10)
U2AF1.vaf2$bc = paste0('AML3005.3_', U2AF1.vaf2$bc, '-1')
U2AF1.vaf2$sample = 'AML3005.3'

U2AF1.combined = rbind(U2AF1.vaf, U2AF1.vaf2)
write.table(U2AF1.combined, file = './data/AML/mutations/AML3005_U2AF1.csv', sep = '\t', quote = F)        

p=DimPlot(AML3005.13, reduction = 'umap', 
        cells.highlight = list('mutated' = U2AF1.combined$bc[which(U2AF1.combined$mutated == 'mutated')], 
                               'wildtype' = U2AF1.combined$bc[which(U2AF1.combined$mutated == 'wildtype')]),
        cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 2) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
ggsave('./AML/figures/AML3005/UMAP/20220719_AML3005_13_UMAP_U2AF1.png', width = 4, height = 4, dpi = 600, plot = p)

U2AF1.combined$predicted.celltype = AML3005.13$predicted.celltype[U2AF1.combined$bc]
U2AF1.combined$manual.cluster = AML3005.13$manual.cluster[U2AF1.combined$bc]

U2AF1.condensed = U2AF1.combined %>% filter(bc %in% colnames(AML3005.13)) %>% 
  group_by(manual.cluster) %>% 
  summarize(mutated.1 = length(which(mutated == 'mutated' & sample == 'AML3005.1')),
            wildtype.1 = length(which(mutated == 'wildtype' & sample == 'AML3005.1')),
            mutated.3 = length(which(mutated == 'mutated' & sample == 'AML3005.3')),
            wildtype.3 = length(which(mutated == 'wildtype' & sample == 'AML3005.3')))
U2AF1.condensed$mutated.freq.1 = U2AF1.condensed$mutated.1 / (U2AF1.condensed$mutated.1 + U2AF1.condensed$wildtype.1)
U2AF1.condensed$mutated.freq.3 = U2AF1.condensed$mutated.3 / (U2AF1.condensed$mutated.3 + U2AF1.condensed$wildtype.3)

ggplot(data=reshape2::melt(U2AF1.condensed), aes(x=variable, y=100*value, color=manual.cluster)) +
  geom_point(size=0.5) +
  geom_line(aes(group=manual.cluster)) + 
  scale_x_discrete(limits = c('mutated.freq.1', 'mutated.freq.3'), labels = c('Screening', 'CR')) +
  scale_y_continuous(limits = c(0,100), '% cells') +
  scale_color_manual(values = cluster.colors) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML3005/plots/20220719_AML3005_13_U2AF1_dynamics.svg', width = 1.1, height = 2)


# TP53 mutation
TP53.vaf = extract_mutation(BC.data.file = './data/AML/pileup/AML3005_1/AML3005_1_pileup_TP53.csv.gz', ALT = 'T', REF = 'C', FILTER = 10)
TP53.vaf$bc = paste0('AML3005.1_', TP53.vaf$bc, '-1')
TP53.vaf$sample = 'AML3005.1'

TP53.vaf2 = extract_mutation(BC.data.file = './data/AML/pileup/AML3005_3/AML3005_3_pileup_TP53.csv.gz', ALT = 'T', REF = 'C', FILTER = 10)
TP53.vaf2$bc = paste0('AML3005.3_', TP53.vaf2$bc, '-1')
TP53.vaf2$sample = 'AML3005.3'

TP53.combined = rbind(TP53.vaf, TP53.vaf2)
write.table(TP53.combined, file = './data/AML/mutations/AML3005_TP53.csv', sep = '\t', quote = F)

p=DimPlot(AML3005.13, reduction = 'umap', 
          cells.highlight = list('mutated' = TP53.combined$bc[which(TP53.combined$mutated == 'mutated')], 
                                 'wildtype' = TP53.combined$bc[which(TP53.combined$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 2) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
ggsave('./AML/figures/AML3005/UMAP/20220719_AML3005_13_UMAP_TP53.png', width = 4, height = 4, dpi = 600, plot = p)

TP53.combined$predicted.celltype = AML3005.13$predicted.celltype[TP53.combined$bc]
TP53.combined$manual.cluster = AML3005.13$manual.cluster[TP53.combined$bc]

TP53.condensed = TP53.combined %>% filter(bc %in% colnames(AML3005.13)) %>%
  group_by(manual.cluster) %>% 
  summarize(mutated.1 = length(which(mutated == 'mutated' & sample == 'AML3005.1')),
            wildtype.1 = length(which(mutated == 'wildtype' & sample == 'AML3005.1')),
            mutated.3 = length(which(mutated == 'mutated' & sample == 'AML3005.3')),
            wildtype.3 = length(which(mutated == 'wildtype' & sample == 'AML3005.3')))
TP53.condensed$mutated.freq.1 = TP53.condensed$mutated.1 / (TP53.condensed$mutated.1 + TP53.condensed$wildtype.1)
TP53.condensed$mutated.freq.3 = TP53.condensed$mutated.3 / (TP53.condensed$mutated.3 + TP53.condensed$wildtype.3)


ggplot(data=reshape2::melt(TP53.condensed), aes(x=variable, y=100*value, color=manual.cluster)) +
  geom_point(size=2) +
  geom_line(aes(group=manual.cluster)) + 
  scale_x_discrete(limits = c('mutated.freq.1', 'mutated.freq.3'), labels = c('Screening', 'CR')) +
  scale_y_continuous(limits = c(0,100), '% cells') +
  scale_color_manual(values = cluster.colors) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML3005/plots/20220719_AML3005_13_TP53_dynamics.svg', width = 1.1, height = 2)


