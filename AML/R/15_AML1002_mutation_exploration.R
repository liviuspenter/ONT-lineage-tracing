# visualize mutations for AML1002.13

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(nanoranger.R)
library(Seurat)
library(SummarizedExperiment)

source('./R/celltypes.R')

AML1002 = readRDS('./data/AML/objects/20220911_AML1002.rds')
Idents(AML1002) = 'predicted.celltype'

IDH2.vaf = read.table(file = './data/AML/mutations/AML1002_IDH2.csv', sep = '\t')
TET2.vaf = read.table(file = './data/AML/mutations/AML1002_TET2.csv', sep = '\t')
TP53_1.vaf = read.table(file = './data/AML/mutations/AML1002_TP53.1.csv', sep = '\t')
TP53_2.vaf = read.table(file = './data/AML/mutations/AML1002_TP53.2.csv', sep = '\t')
U2AF1.vaf = read.table(file = './data/AML/mutations/AML1002_U2AF1.csv', sep = '\t')

####
mutations.1002 = merge(IDH2.vaf[,c('bc', 'mutated')], TET2.vaf[,c('bc', 'mutated')], by='bc', all=T)
colnames(mutations.1002) = c('bc', 'IDH2', 'TET2')
mutations.1002 = merge(mutations.1002, TP53_1.vaf[,c('bc', 'mutated')], by='bc', all=T)
colnames(mutations.1002)[4] = 'TP53.1'
mutations.1002 = merge(mutations.1002, TP53_2.vaf[,c('bc', 'mutated')], by='bc', all=T)
colnames(mutations.1002)[5] = 'TP53.2'
mutations.1002 = merge(mutations.1002, U2AF1.vaf[,c('bc', 'mutated')], by='bc', all=T)
colnames(mutations.1002)[6] = 'U2AF1'
rownames(mutations.1002) = mutations.1002$bc
mutations.1002 = mutations.1002[intersect(colnames(AML1002), mutations.1002$bc),]
mutations.1002$predicted.celltype = AML1002$predicted.celltype[mutations.1002$bc]
mutations.1002 = mutations.1002[which(mutations.1002$predicted.celltype %in% names(which(table(mutations.1002$predicted.celltype) > 5))),]


celltype = factor(mutations.1002[which(grepl('AML1002.1', rownames(mutations.1002))),'predicted.celltype',],
                  levels = names(AML.combined.colors)[which(names(AML.combined.colors) %in% 
                                                              mutations.1002[which(grepl('AML1002.1', rownames(mutations.1002))),'predicted.celltype',])])
ha = columnAnnotation(celltype = celltype, col = list('celltype' = AML.combined.colors), border=T, simple_anno_size = unit(5, 'pt'))
svglite::svglite('./AML/figures/AML1002/heatmaps/20220912_AML1002_1.svg', width = 5, height = 1.5)
ComplexHeatmap::Heatmap(t(mutations.1002[which(grepl('AML1002.1', rownames(mutations.1002))),
                                         c('IDH2', 'TET2', 'TP53.1', 'TP53.2', 'U2AF1')]), 
                        show_column_names = F, col = c('mutated' = 'firebrick', 'wildtype' = 'white'), na_col = 'white',
                        column_split = celltype, column_title_side = 'bottom', column_title_rot = 90, column_title_gp = gpar(fontsize=8),
                        row_names_side = 'left', row_names_gp = gpar(fontsize=8), 
                        row_labels = c('IDH2R140Q', 'TET2I1873T', 'TP53H179R', 'TP53P278S', 'U2AF1S34Y'), 
                        top_annotation = ha, border = T, raster_quality = 10, use_raster = T)
dev.off()

celltype = factor(mutations.1002[which(grepl('AML1002.3', rownames(mutations.1002))),'predicted.celltype',],
                  levels = names(AML.combined.colors)[which(names(AML.combined.colors) %in% 
                                                              mutations.1002[which(grepl('AML1002.3', rownames(mutations.1002))),'predicted.celltype',])])
ha = columnAnnotation(celltype = celltype, col = list('celltype' = AML.combined.colors), border=T, simple_anno_size = unit(5, 'pt'))
svglite::svglite('./AML/figures/AML1002/heatmaps/20220912_AML1002_3.svg', width = 5, height = 1.5)
ComplexHeatmap::Heatmap(t(mutations.1002[which(grepl('AML1002.3', rownames(mutations.1002))),
                                         c('IDH2', 'TET2', 'TP53.1', 'TP53.2', 'U2AF1')]), 
                        show_column_names = F, col = c('mutated' = 'firebrick', 'wildtype' = 'white'), na_col = 'white',
                        column_split = celltype, column_title_side = 'bottom', column_title_rot = 90, column_title_gp = gpar(fontsize=8),
                        row_names_side = 'left', row_names_gp = gpar(fontsize=8), 
                        row_labels = c('IDH2R140Q', 'TET2I1873T', 'TP53H179R', 'TP53P278S', 'U2AF1S34Y'), 
                        top_annotation = ha, border = T, raster_quality = 10, use_raster = T)
dev.off()



### statistics
mutations.1002$sample = stringr::str_split_fixed(mutations.1002$bc, pattern = '_', n=2)[,1]
statistics.1002 = mutations.1002 %>% 
  group_by(sample, predicted.celltype) %>% 
  summarize(IDH2.mut = length(which(IDH2 == 'mutated')),
            IDH2.wt = length(which(IDH2 == 'wildtype')),
            TET2.mut = length(which(TET2 == 'mutated')),
            TET2.wt = length(which(TET2 == 'wildtype')),
            TP53.1.mut = length(which(TP53.1 == 'mutated')),
            TP53.1.wt = length(which(TP53.1 == 'wildtype')),
            TP53.2.mut = length(which(TP53.2 == 'mutated')),
            TP53.2.wt = length(which(TP53.2 == 'wildtype')),
            U2AF1.mut = length(which(U2AF1 == 'mutated')),
            U2AF1.wt = length(which(U2AF1 == 'wildtype')))

statistics.1002$IDH2.mut.freq = statistics.1002$IDH2.mut / (statistics.1002$IDH2.mut + statistics.1002$IDH2.wt)
statistics.1002$TET2.mut.freq = statistics.1002$TET2.mut / (statistics.1002$TET2.mut + statistics.1002$TET2.wt)
statistics.1002$TP53.1.mut.freq = statistics.1002$TP53.1.mut / (statistics.1002$TP53.1.mut + statistics.1002$TP53.1.wt)
statistics.1002$TP53.2.mut.freq = statistics.1002$TP53.2.mut / (statistics.1002$TP53.2.mut + statistics.1002$TP53.2.wt)
statistics.1002$U2AF1.mut.freq = statistics.1002$U2AF1.mut / (statistics.1002$U2AF1.mut + statistics.1002$U2AF1.wt)

ggplot(statistics.1002[which(statistics.1002$predicted.celltype %in% myeloid.clusters),], aes(x=sample, y=100*TET2.mut.freq, color=predicted.celltype)) + 
  geom_point() + 
  geom_line(aes(group=predicted.celltype)) + 
  scale_x_discrete(labels = c('Screening', 'SD')) + 
  scale_y_continuous('%TET2S34Y',limits = c(0,100)) +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML1002/plots/20220912_AML1002_TET2_myeloid_frequency.svg', width = 1.1, height = 2)

ggplot(statistics.1002[which(statistics.1002$predicted.celltype %in% TNK.clusters),], aes(x=sample, y=100*TET2.mut.freq, color=predicted.celltype)) + 
  geom_point() + 
  geom_line(aes(group=predicted.celltype)) + 
  scale_x_discrete(labels = c('Screening', 'SD')) + 
  scale_y_continuous('%TET2S34Y',limits = c(0,100)) +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML1002/plots/20220912_AML1002_TET2_TNK_frequency.svg', width = 1.1, height = 2)

# TP53.1
ggplot(statistics.1002[which(statistics.1002$predicted.celltype %in% myeloid.clusters),], aes(x=sample, y=100*TP53.1.mut.freq, color=predicted.celltype)) + 
  geom_point() + 
  geom_line(aes(group=predicted.celltype)) + 
  scale_x_discrete(labels = c('Screening', 'SD')) + 
  scale_y_continuous('%TP53H179R',limits = c(0,100)) +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML1002/plots/20220912_AML1002_TP53.1_myeloid_frequency.svg', width = 1.1, height = 2)

ggplot(statistics.1002[which(statistics.1002$predicted.celltype %in% TNK.clusters),], aes(x=sample, y=100*TP53.1.mut.freq, color=predicted.celltype)) + 
  geom_point() + 
  geom_line(aes(group=predicted.celltype)) + 
  scale_x_discrete(labels = c('Screening', 'SD')) + 
  scale_y_continuous('%TP53H179R',limits = c(0,100)) +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML1002/plots/20220912_AML1002_TP53.1_TNK_frequency.svg', width = 1.1, height = 2)

# TP53.2
ggplot(statistics.1002[which(statistics.1002$predicted.celltype %in% myeloid.clusters),], aes(x=sample, y=100*TP53.2.mut.freq, color=predicted.celltype)) + 
  geom_point() + 
  geom_line(aes(group=predicted.celltype)) + 
  scale_x_discrete(labels = c('Screening', 'SD')) + 
  scale_y_continuous('%TP53P278S',limits = c(0,100)) +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML1002/plots/20220912_AML1002_TP53.2_myeloid_frequency.svg', width = 1.1, height = 2)

ggplot(statistics.1002[which(statistics.1002$predicted.celltype %in% TNK.clusters),], aes(x=sample, y=100*TP53.2.mut.freq, color=predicted.celltype)) + 
  geom_point() + 
  geom_line(aes(group=predicted.celltype)) + 
  scale_x_discrete(labels = c('Screening', 'SD')) + 
  scale_y_continuous('%TP53P278S',limits = c(0,100)) +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML1002/plots/20220912_AML1002_TP53.2_TNK_frequency.svg', width = 1.1, height = 2)

# U2AF1
ggplot(statistics.1002[which(statistics.1002$predicted.celltype %in% myeloid.clusters),], aes(x=sample, y=100*U2AF1.mut.freq, color=predicted.celltype)) + 
  geom_point() + 
  geom_line(aes(group=predicted.celltype)) + 
  scale_x_discrete(labels = c('Screening', 'SD')) + 
  scale_y_continuous('%U2AF1S34Y',limits = c(0,100)) +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML1002/plots/20220912_AML1002_U2AF1_myeloid_frequency.svg', width = 1.1, height = 2)

ggplot(statistics.1002[which(statistics.1002$predicted.celltype %in% TNK.clusters),], aes(x=sample, y=100*U2AF1.mut.freq, color=predicted.celltype)) + 
  geom_point() + 
  geom_line(aes(group=predicted.celltype)) + 
  scale_x_discrete(labels = c('Screening', 'SD')) + 
  scale_y_continuous('%U2AF1S34Y',limits = c(0,100)) +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML1002/plots/20220912_AML1002_U2AF1_TNK_frequency.svg', width = 1.1, height = 2)