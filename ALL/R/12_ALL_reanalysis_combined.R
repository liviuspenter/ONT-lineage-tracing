# reanalysis of data from Witkowski et al., Cancer Cell 2020 and Caron et al., Scientific Reports 2020
# 
# perform CNV calls from numbat output

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(numbat)
library(Seurat)
library(numbat)

### generate combined UMAP and heatmap
ALL.validation = readRDS(file='./data/ALL/objects/20221027_ALL_validation.rds')

myeloid.clusters = c('HSC', 'LMPP', 'GMP', 'CD14 Mono', 'CD16 Mono', 'pDC', 'cDC2', 'Prog_DC')

TNK.clusters = c('CD8 Naive','CD8 Effector_1','CD8 Effector_2','CD8 Memory_1','CD8 Memory_2','CD56 bright NK','NK','gdT','MAIT','CD4 Naive','CD4 Memory','Treg')

df = as.data.frame(ALL.validation@reductions$ref.umap@cell.embeddings)
df$mutated = NA
df$sample = ALL.validation$orig.ident

use.samples = c('ETV001', 'ETV002', 'ETV003', 'PH001', 'PH002', 'ETV6-RUNX1_1', 'ETV6-RUNX1_2', 'ETV6-RUNX1_4', 'HHD_1', 'HHD_2', 'PRE-T_1', 'PRE-T_2')
cnv.data.combined = data.frame()
for (s in use.samples) {
  cnv.data = read.csv2(paste0('./data/ALL/numbat/', s, '_cnv_calls.csv'), sep = '\t')
  cnv.data.combined = rbind(cnv.data.combined, cnv.data[,c('bc', 'mutated', 'predicted.celltype', 'sample')])
}
rownames(cnv.data.combined) = cnv.data.combined$bc

df$mutated = cnv.data.combined[rownames(df), 'mutated']
df = df[intersect(rownames(df), rownames(cnv.data.combined)),]

ggplot() + 
  geom_point(data=df[which(!df$sample %in% c('PRE-T_1', 'PRE-T_2')),], aes(x=refUMAP_1, y=refUMAP_2, color=mutated), size=0.5) +
  scale_color_manual(values = c('mutated' = 'red', 'wildtype' = 'black')) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
ggsave('./ALL/figures/UMAP/20221029_B_ALL_validation.png', width = 4, height = 4, dpi = 600)

ggplot() + 
  geom_point(data=df[which(df$sample %in% c('PRE-T_1', 'PRE-T_2')),], aes(x=refUMAP_1, y=refUMAP_2, color=mutated), size=0.5) +
  scale_color_manual(values = c('mutated' = 'red', 'wildtype' = 'black')) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())
ggsave('./ALL/figures/UMAP/20221029_T_ALL_validation.png', width = 4, height = 4, dpi = 600)

boo = as.data.frame(cnv.data.combined %>% group_by(sample, predicted.celltype) %>% 
                      summarize(mut.freq = length(which(mutated == 'mutated')) / length(predicted.celltype)) %>%
                      tidyr::pivot_wider(names_from = 'predicted.celltype', values_from = 'mut.freq'))
rownames(boo) = boo$sample
boo = boo[,-1]
boo = boo[,c(myeloid.clusters, 'Prog_Mk', 'Prog_RBC', TNK.clusters, 'Prog_B 1', 'Prog_B 2', 'Naive B','Memory B', 'Plasmablast')]
boo = boo[c('ETV001', 'ETV002', 'ETV003', 'PH001', 'PH002', 'ETV6-RUNX1_1', 'ETV6-RUNX1_2', 'ETV6-RUNX1_4', 'HHD_1', 'HHD_2'),]

cell.number = as.data.frame(cnv.data.combined %>% group_by(predicted.celltype) %>% summarize(n = n()))
rownames(cell.number) = cell.number$predicted.celltype
n = cell.number[c(myeloid.clusters, 'Prog_Mk', 'Prog_RBC', TNK.clusters, 'Prog_B 1', 'Prog_B 2', 'Naive B','Memory B', 'Plasmablast'), 'n']

ha = rowAnnotation(barplot = anno_barplot(x = n), annotation_label = '# cells', annotation_name_gp = gpar(fontsize=8))
ha2 = rowAnnotation(celltype = c(myeloid.clusters, 'Prog_Mk', 'Prog_RBC', TNK.clusters, 'Prog_B 1', 'Prog_B 2', 'Naive B','Memory B', 'Plasmablast'), 
                    col = list('celltype' = nanoranger.R::AML.combined.colors[c(myeloid.clusters, 'Prog_Mk', 'Prog_RBC', TNK.clusters, 'Prog_B 1', 'Prog_B 2', 'Naive B','Memory B', 'Plasmablast')]), border=T,
                    annotation_name_gp = gpar(fontsize=8), simple_anno_size = unit(5, 'pt'))

svglite::svglite('./ALL/figures/heatmaps/20221029_ALL_CNV_celltype.svg', width = 5, height = 4)
col_fun = circlize::colorRamp2(breaks = seq(0,1,1/8), colors = c(BuenColors::jdb_palette(name = 'solar_rojos')))
Heatmap(t(boo), cluster_rows = F, cluster_columns = F, col = col_fun, na_col = 'white', border=T, row_names_side = 'left',
        right_annotation = ha, left_annotation = ha2, row_names_gp = gpar(fontsize=8), row_title_gp = gpar(fontsize=0), column_names_gp = gpar(fontsize=8),
        row_split = factor(c(rep('myeloid', length(myeloid.clusters)),
                             rep('MK_RBC', 2),
                             rep('TNK', length(TNK.clusters)),
                             rep('B', 5)), c('B', 'myeloid', 'MK_RBC', 'TNK')))
dev.off()
