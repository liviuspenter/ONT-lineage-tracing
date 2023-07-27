# combined analysis of somatic mutations across all cases

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(Seurat)

source('./R/celltypes.R')

AML.combined = readRDS('./data/AML/objects/AML.all.rds')

files = list.files('./data/AML/mutations/')

mutations.df = data.frame()
for (f in files) {
  boo = as.data.frame(read.csv2(file = paste0('./data/AML/mutations/',f), sep = '\t'))
  mutations.df = rbind(mutations.df, boo[,c('bc', 'alt', 'ref', 'mutated', 'vaf')])
}
mutations.df$sample = stringr::str_split_fixed(mutations.df$bc, pattern = '_', n=2)[,1]

mutations.df = mutations.df[which(mutations.df$bc %in% colnames(AML.combined)),]
mutations.df.condensed = mutations.df %>% group_by(bc) %>% summarize(mutation = ifelse('mutated' %in% mutated, 'mutated', 'wildtype'))
mutations.df.condensed$patient = stringr::str_split_fixed(mutations.df.condensed$bc, pattern = '\\.', n=2)[,1]
mutations.df.condensed$sample = stringr::str_split_fixed(mutations.df.condensed$bc, pattern = '_', n=2)[,1]

AML.subset = subset(AML.combined, cells = colnames(AML.combined)[which(colnames(AML.combined) %in% unique(mutations.df$bc))])

p=DimPlot(AML.subset, sizes.highlight = 0.5, 
          reduction = 'ref.umap', 
          cells.highlight = list('mutated' = mutations.df.condensed$bc[which(mutations.df.condensed$mutation == 'mutated')],
                                 'wildtype' = mutations.df.condensed$bc[which(mutations.df.condensed$mutation == 'wildtype')])) +
  scale_color_manual(values = c('wildtype' = 'black', 'mutated' = 'red')) +
  NoLegend() + 
  NoAxes()
ggsave('./AML/figures/combined/UMAP/20221005_AML_mutations.png', width = 4, height = 4, dpi = 600)

mutations.df.condensed$predicted.celltype = AML.combined$predicted.celltype[mutations.df.condensed$bc]
boo = mutations.df.condensed %>% group_by(predicted.celltype, patient) %>%
  summarize(mut.freq = length(which(mutation == 'mutated')) / length(predicted.celltype))
boo = as.data.frame(tidyr::pivot_wider(data = boo, names_from = 'predicted.celltype', values_from = 'mut.freq'))
rownames(boo) = boo$patient
boo = boo[,-1]
boo = t(boo)
boo = boo[names(AML.combined.colors),]
boo = boo[c(myeloid.clusters, 'Prog_Mk', 'Prog_RBC', TNK.clusters),]
col_fun = circlize::colorRamp2(breaks = seq(0,1,1/8), colors = BuenColors::jdb_palette(name = 'solar_rojos'))
cell.number = as.data.frame(mutations.df.condensed %>% group_by(predicted.celltype) %>% summarize(n = n()))
rownames(cell.number) = cell.number$predicted.celltype
n = cell.number[c(myeloid.clusters, 'Prog_Mk', 'Prog_RBC', TNK.clusters), 'n']
names(n) = c(myeloid.clusters, 'Prog_Mk', 'Prog_RBC', TNK.clusters)

ha = rowAnnotation(barplot = anno_barplot(x = n), annotation_label = '# cells', annotation_name_gp = gpar(fontsize=8))
ha2 = rowAnnotation(celltype = c(myeloid.clusters, 'Prog_Mk', 'Prog_RBC', TNK.clusters), 
                    col = list('celltype' = AML.combined.colors[c(myeloid.clusters, 'Prog_Mk', 'Prog_RBC', TNK.clusters)]), border=T,
                    annotation_name_gp = gpar(fontsize=8), simple_anno_size = unit(5, 'pt'))

svglite::svglite('./AML/figures/combined/heatmaps/20221004_mutations_celltype.svg', width = 4.7, height = 3.3)
Heatmap(boo, na_col = 'white', cluster_rows = F, cluster_columns = T, show_column_dend = F, col = col_fun, row_names_side = 'left', border=T,
        row_split = factor(c(rep('myeloid', 8), rep('Prog_RBC_MK', 2), rep('TNK', 12)), levels = c('myeloid', 'Prog_RBC_MK', 'TNK')),
        right_annotation = ha, left_annotation = ha2, row_names_gp = gpar(fontsize=8), row_title_gp = gpar(fontsize=0), column_names_gp = gpar(fontsize=8)) 
dev.off()

# save for table
boo = mutations.df.condensed %>% group_by(predicted.celltype, patient) %>%
  summarize(mut.cells = length(which(mutation == 'mutated')),
            cells = length(predicted.celltype))
write.csv2(boo, file = './data/AML/mutations/20230211_combined_mutations.csv', quote = F)

### percent mutated versus number of cells
boo = mutations.df.condensed %>% group_by(predicted.celltype, patient) %>%
  summarize(cells = length(predicted.celltype), 
            mut.freq = length(which(mutation == 'mutated')) / length(predicted.celltype))

ggplot(boo[which(boo$cells > 4),], aes(x=predicted.celltype, y=100*mut.freq, color=predicted.celltype)) +
  geom_jitter(size=0.5) + 
  stat_summary(geom = 'crossbar', width = 0.5, size=0.5, fun.min = median, fun.max = median, fun = median, color='black') + 
  scale_x_discrete(limits = c(myeloid.clusters, 'Prog_Mk', 'Prog_RBC', TNK.clusters)) + 
  scale_y_continuous('% cells') + 
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/combined/plots/20221025_mutation_frequencies.svg', width = 2.8, height = 2)
