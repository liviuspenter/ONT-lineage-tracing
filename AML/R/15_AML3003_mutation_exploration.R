# visualize mutations for AML3003.14

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(nanoranger.R)
library(Seurat)
library(SummarizedExperiment)

source('./R/celltypes.R')

AML3003.14 = readRDS('./data/AML/objects/202208_AML3003_14.rds')

ASXL1.vaf = read.table(file = './data/AML/mutations/AML3003_ASXL1.csv', sep = '\t')
STAG2.vaf = read.table(file = './data/AML/mutations/AML3003_STAG2.csv', sep = '\t')
TET2.2.vaf = read.table(file = './data/AML/mutations/AML3003_TET2.2.csv', sep = '\t')
TET2.3.vaf = read.table(file = './data/AML/mutations/AML3003_TET2.3.csv', sep = '\t')
TET2.4.vaf = read.table(file = './data/AML/mutations/AML3003_TET2.4.csv', sep = '\t')

### combine ASXL1, STAG2, TET2.2, TET2.3 and TET2.4 mutations
mutations.3003 = merge(ASXL1.vaf[,c('bc', 'mutated')], STAG2.vaf[,c('bc', 'mutated')], by='bc', all=T)
colnames(mutations.3003) = c('bc', 'ASXL1', 'STAG2')
mutations.3003 = merge(mutations.3003, TET2.2.vaf[,c('bc', 'mutated')], by='bc', all=T)
colnames(mutations.3003)[4] = 'TET2.2'
mutations.3003 = merge(mutations.3003, TET2.3.vaf[,c('bc', 'mutated')], by='bc', all=T)
colnames(mutations.3003)[5] = 'TET2.3'
mutations.3003 = merge(mutations.3003, TET2.4.vaf[,c('bc', 'mutated')], by='bc', all=T)
colnames(mutations.3003)[6] = 'TET2.4'
rownames(mutations.3003) = mutations.3003$bc
mutations.3003 = mutations.3003[intersect(colnames(AML3003.14), mutations.3003$bc),]
mutations.3003$predicted.celltype = AML3003.14$predicted.celltype[mutations.3003$bc]
mutations.3003 = mutations.3003[which(mutations.3003$predicted.celltype %in% names(which(table(mutations.3003$predicted.celltype) > 5))),]

celltype = factor(mutations.3003[which(grepl('AML3003.1', rownames(mutations.3003))),'predicted.celltype',],
                  levels = names(AML.combined.colors)[which(names(AML.combined.colors) %in% 
                                                              mutations.3003[which(grepl('AML3003.1', rownames(mutations.3003))),'predicted.celltype',])])
ha = columnAnnotation(celltype = celltype, col = list('celltype' = AML.combined.colors), border=T, simple_anno_size = unit(5, 'pt'))
svglite::svglite('./AML/figures/AML3003/heatmaps/20220816_AML3003_1.svg', width = 5, height = 1.5)
ComplexHeatmap::Heatmap(t(mutations.3003[which(grepl('AML3003.1', rownames(mutations.3003))),
                                         c('ASXL1', 'STAG2', 'TET2.2', 'TET2.3', 'TET2.4')]), 
                        show_column_names = F, col = c('mutated' = 'firebrick', 'wildtype' = 'white'), na_col = 'white',
                        column_split = celltype, column_title_side = 'bottom', column_title_rot = 90, column_title_gp = gpar(fontsize=8),
                        row_names_side = 'left', row_names_gp = gpar(fontsize=8), 
                        row_labels = c('ASXL1Y591*', 'STAG2Q275*', 'TET2N1504Kfs', 'TET2N488Mfs', 'TET2I1762V'), 
                        top_annotation = ha, border = T)
dev.off()

celltype = factor(mutations.3003[which(grepl('AML3003.4', rownames(mutations.3003))),'predicted.celltype',],
                  levels = names(AML.combined.colors)[which(names(AML.combined.colors) %in% 
                                                              mutations.3003[which(grepl('AML3003.4', rownames(mutations.3003))),'predicted.celltype',])])
ha = columnAnnotation(celltype = celltype, col = list('celltype' = AML.combined.colors), border=T, simple_anno_size = unit(5, 'pt'))
svglite::svglite('./AML/figures/AML3003/heatmaps/20220816_AML3003_4.svg', width = 5, height = 1.5)
ComplexHeatmap::Heatmap(t(mutations.3003[which(grepl('AML3003.4', rownames(mutations.3003))),
                                         c('ASXL1', 'STAG2', 'TET2.2', 'TET2.3', 'TET2.4')]), 
                        show_column_names = F, col = c('mutated' = 'firebrick', 'wildtype' = 'white'), na_col = 'white',
                        column_split = celltype, column_title_side = 'bottom', column_title_rot = 90, column_title_gp = gpar(fontsize=8),
                        row_names_side = 'left', row_names_gp = gpar(fontsize=8), 
                        row_labels = c('ASXL1Y591*', 'STAG2Q275*', 'TET2N1504Kfs', 'TET2N488Mfs', 'TET2I1762V'), 
                        top_annotation = ha, border = T)
dev.off()

### statistics
mutations.3003$sample = stringr::str_split_fixed(mutations.3003$bc, pattern = '_', n=2)[,1]
statistics.3003 = mutations.3003 %>% group_by(sample, predicted.celltype) %>% summarize(ASXL1.mut = length(which(ASXL1 == 'mutated')),
                                                                                        ASXL1.wt = length(which(ASXL1 == 'wildtype')),
                                                                                        STAG2.mut = length(which(STAG2 == 'mutated')),
                                                                                        STAG2.wt = length(which(STAG2 == 'wildtype')),
                                                                                        TET2.2.mut = length(which(TET2.2 == 'mutated')),
                                                                                        TET2.2.wt = length(which(TET2.2 == 'wildtype')),
                                                                                        TET2.3.mut = length(which(TET2.3 == 'mutated')),
                                                                                        TET2.3.wt = length(which(TET2.3 == 'wildtype')),
                                                                                        TET2.4.mut = length(which(TET2.4 == 'mutated')),
                                                                                        TET2.4.wt = length(which(TET2.4 == 'wildtype')))

statistics.3003$ASXL1.mut.freq = statistics.3003$ASXL1.mut / (statistics.3003$ASXL1.mut + statistics.3003$ASXL1.wt)
statistics.3003$STAG2.mut.freq = statistics.3003$STAG2.mut / (statistics.3003$STAG2.mut + statistics.3003$STAG2.wt)
statistics.3003$TET2.2.mut.freq = statistics.3003$TET2.2.mut / (statistics.3003$TET2.2.mut + statistics.3003$TET2.2.wt)
statistics.3003$TET2.3.mut.freq = statistics.3003$TET2.3.mut / (statistics.3003$TET2.3.mut + statistics.3003$TET2.3.wt)
statistics.3003$TET2.4.mut.freq = statistics.3003$TET2.4.mut / (statistics.3003$TET2.4.mut + statistics.3003$TET2.4.wt)

ggplot(statistics.3003[which(statistics.3003$predicted.celltype %in% myeloid.clusters),], aes(x=sample, y=100*STAG2.mut.freq, color=predicted.celltype)) + 
  geom_point() + 
  geom_line(aes(group=predicted.celltype)) + 
  scale_x_discrete(labels = c('Screening', 'CR')) + 
  scale_y_continuous('%STAG2Q275*',limits = c(0,100)) +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./figures/plots/3003/20220816_AML3003_1_STAG2_frequency.svg', width = 1.1, height = 2)

ggplot(statistics.3003[which(statistics.3003$predicted.celltype %in% myeloid.clusters),], aes(x=sample, y=100*TET2.3.mut.freq, color=predicted.celltype)) + 
  geom_point() + 
  geom_line(aes(group=predicted.celltype)) + 
  scale_x_discrete(labels = c('Screening', 'CR')) + 
  scale_y_continuous('%TET2N488Mfs',limits = c(0,100)) +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./figures/plots/3003/20220816_AML3003_1_TET2.3_frequency.svg', width = 1.1, height = 2)

statistics.3003 = mutations.3003 %>% 
  group_by(predicted.celltype) %>% 
  summarize(ASXL1.mut = length(which(ASXL1 == 'mutated')),
            ASXL1.wt = length(which(ASXL1 == 'wildtype')),
            STAG2.mut = length(which(STAG2 == 'mutated')),
            STAG2.wt = length(which(STAG2 == 'wildtype')),
            TET2.2.mut = length(which(TET2.2 == 'mutated')),
            TET2.2.wt = length(which(TET2.2 == 'wildtype')),
            TET2.3.mut = length(which(TET2.3 == 'mutated')),
            TET2.3.wt = length(which(TET2.3 == 'wildtype')),
            TET2.4.mut = length(which(TET2.4 == 'mutated')),
            TET2.4.wt = length(which(TET2.4 == 'wildtype')))

statistics.3003$ASXL1.mut.freq = statistics.3003$ASXL1.mut / (statistics.3003$ASXL1.mut + statistics.3003$ASXL1.wt)
statistics.3003$STAG2.mut.freq = statistics.3003$STAG2.mut / (statistics.3003$STAG2.mut + statistics.3003$STAG2.wt)
statistics.3003$TET2.2.mut.freq = statistics.3003$TET2.2.mut / (statistics.3003$TET2.2.mut + statistics.3003$TET2.2.wt)
statistics.3003$TET2.3.mut.freq = statistics.3003$TET2.3.mut / (statistics.3003$TET2.3.mut + statistics.3003$TET2.3.wt)
statistics.3003$TET2.4.mut.freq = statistics.3003$TET2.4.mut / (statistics.3003$TET2.4.mut + statistics.3003$TET2.4.wt)
statistics.3003$predicted.celltype = factor(statistics.3003$predicted.celltype, levels = names(AML.combined.colors))

ggplot(statistics.3003[which(statistics.3003$predicted.celltype %in% TNK.clusters),], aes(x=predicted.celltype, y=100*STAG2.mut.freq, fill=predicted.celltype)) + 
  geom_col() + 
  scale_y_continuous('%STAG2Q275',limits = c(0,100)) +
  scale_fill_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML3003/plots/20220816_AML3003_STAG2_frequency_immune_cells.svg', width = 1.5, height = 2)

ggplot(statistics.3003[which(statistics.3003$predicted.celltype %in% TNK.clusters),], aes(x=predicted.celltype, y=100*TET2.3.mut.freq, fill=predicted.celltype)) + 
  geom_col() + 
  scale_y_continuous('%TET2N488Mfs',limits = c(0,100)) +
  scale_fill_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML3003/plots/20220816_AML3003_TET2.3_frequency_immune_cells.svg', width = 1.5, height = 2)

### STAG2 mutated versus non-mutated
STAG2.vaf$predicted.celltype = AML3003.14$predicted.celltype[STAG2.vaf$bc]
STAG2.vaf = as.data.frame(STAG2.vaf)
rownames(STAG2.vaf) = STAG2.vaf$bc
boo = subset(AML3003.14, cells = which(colnames(AML3003.14) %in% STAG2.vaf$bc[which(STAG2.vaf$predicted.celltype == 'LMPP' & !is.na(STAG2.vaf$mutated))]))
boo$mutated = STAG2.vaf[colnames(boo), 'mutated']
STAG2.signature = FindMarkers(boo, group.by = 'mutated', ident.1 = 'wildtype', ident.2 = 'mutated')
STAG2.signature$gene = rownames(STAG2.signature)
highlight.genes = STAG2.signature$gene[which(-log10(STAG2.signature$p_val) > 2.5 & abs(STAG2.signature$avg_log2FC) > 0.5)]
STAG2.signature$highlight = ifelse(STAG2.signature$gene %in% highlight.genes, 'yes', 'no')
p=ggplot(STAG2.signature, aes(x=avg_log2FC, y=-log10(p_val))) + 
  ggrastr::rasterize(geom_point(aes(color=highlight), size=0.5), dpi=600) + 
  geom_hline(yintercept = 2.5) + 
  geom_vline(xintercept = c(-0.5,0.5)) +
  scale_x_continuous('Log2FC') + 
  scale_y_continuous('-log10(p value)') + 
  scale_color_manual(values = c('no' = 'grey', 'yes' = 'black')) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'))
ggsave('./AML/figures/AML3003/volcano/20220816_AML3003_LMPP_STAG2_volcano.svg', width = 2, height = 2, plot = p)

p=DoHeatmap(boo, group.by = 'mutated', features = highlight.genes, disp.max = 2, 
            group.colors = c('mutated' = 'firebrick', 'wildtype' = 'black'), label = F) + 
  scale_fill_gradientn(colours = BuenColors::jdb_palette(name = 'brewer_celsius'), na.value = 'black') + 
  theme(axis.text = element_text('Arial', size=8, color='black'))
ggsave('./AML/figures/AML3003/heatmaps/20220816_AML3003_LMPP_STAG2_mutated.svg', width = 3, height = 3.5, plot = p)
