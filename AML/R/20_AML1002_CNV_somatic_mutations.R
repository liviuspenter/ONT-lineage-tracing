# visualize somatic mutations and CNV changes in AML1002.13

library(ComplexHeatmap)
library(dplyr)
library(nanoranger.R)
library(ggplot2)
library(Seurat)

AML1002 = readRDS('./data/AML/objects/20220930_AML1002_all.rds')

cnv.data = read.csv2('./data/AML/numbat/AML1002.1/cnv_calls.csv', sep = '\t')
cnv.data$bc = paste0('AML1002.1_', cnv.data$cell)
cnv.data$p_cnv = as.numeric(cnv.data$p_cnv)
cnv.data = cnv.data %>% select(bc, seg, p_cnv) %>% filter(seg %in% c('7b')) %>% tidyr::pivot_wider(names_from = 'seg', values_from = 'p_cnv')

cnv.data.2 = read.csv2('./data/AML/numbat/AML1002.3.all/cnv_calls.csv', sep = '\t')
cnv.data.2$bc = paste0('AML1002.3_', cnv.data.2$cell)
cnv.data.2$p_cnv = as.numeric(cnv.data.2$p_cnv)
cnv.data.2 = cnv.data.2 %>% select(bc, seg, p_cnv) %>% filter(seg %in% c('7b')) %>% tidyr::pivot_wider(names_from = 'seg', values_from = 'p_cnv')

cnv.data = rbind(cnv.data, cnv.data.2)

mutations.df = data.frame()
for (f in c('AML1002_IDH2.csv', 'AML1002_TP53.1.csv', 'AML1002_TP53.2.csv', 'AML1002_TET2.csv', 'AML1002_U2AF1.csv')) {
  boo = as.data.frame(read.csv2(file = paste0('./data/AML/mutations/',f), sep = '\t'))
  boo$gene = gsub(stringr::str_split_fixed(f, pattern = '_', n=2)[,2], pattern = '.csv', replacement = '')
  mutations.df = rbind(mutations.df, boo[,c('bc', 'alt', 'ref', 'mutated', 'vaf', 'gene')])
}
mutations.df$sample = stringr::str_split_fixed(mutations.df$bc, pattern = '_', n=2)[,1]
mutations.df.condensed = mutations.df %>% group_by(bc) %>% summarize(mutation = ifelse('mutated' %in% mutated, 'mutated', 'wildtype'))
mutations.df.condensed = merge(mutations.df.condensed, cnv.data, by = 'bc', all.x=T)
mutations.df.condensed = merge(mutations.df.condensed, 
                               tidyr::pivot_wider(data=mutations.df[,c('bc', 'mutated', 'gene')], 
                                                  values_from = 'mutated', names_from = 'gene')[,c('bc','IDH2', 'TP53.1', 'TP53.2', 'TET2', 'U2AF1')], 
      by = 'bc', all.x = T)
mutations.df.condensed$predicted.celltype = AML1002$predicted.celltype[mutations.df.condensed$bc]
mutations.df.condensed$sample = AML1002$orig.ident[mutations.df.condensed$bc]
rownames(mutations.df.condensed) = mutations.df.condensed$bc
mutations.df.condensed = mutations.df.condensed[intersect(mutations.df.condensed$bc, colnames(AML1002)),]

AML.subset = subset(AML1002, cells = mutations.df.condensed$bc)
AML.subset = ScaleData(AML.subset, features = rownames(AML.subset))

marker.genes.RNA = c('CD34', 'HLA-DRA','CD33', 'KIT', 'MPO', 'CD38', 'CD14', 'GATA2','ZFPM1', 'FLI1', 'NFE2','PF4', 'GATA1','HBA1', 'HBB','CD3D', 'CD4', 'CD8A', 'NCAM1', 'MKI67')

mutations.df.condensed = cbind(mutations.df.condensed, t(GetAssayData(AML.subset, assay = 'RNA', slot = 'scale.data')[marker.genes.RNA,]))

mutations.df.condensed$predicted.celltype[which(mutations.df.condensed$predicted.celltype %in% 
                                                  c('CD8 Naive', 'CD8 Effector_1', 'CD8 Effector_2', 'CD8 Memory_1', 'CD8 Memory_2'))] = 'CD8'
mutations.df.condensed$predicted.celltype[which(mutations.df.condensed$predicted.celltype %in% 
                                                  c('CD4 Naive', 'CD4 Memory', 'Treg'))] = 'CD4'
mutations.df.condensed$predicted.celltype[which(mutations.df.condensed$predicted.celltype %in% 
                                                  c('NK', 'CD56 bright NK'))] = 'NK cells'
mutations.df.condensed$predicted.celltype[which(mutations.df.condensed$predicted.celltype %in% 
                                                  c('gdT', 'MAIT'))] = 'other'

cells.1 = mutations.df.condensed$bc[which(mutations.df.condensed$sample == 'AML1002.1' & 
                                            mutations.df.condensed$predicted.celltype %in% c('HSC', 'LMPP', 'GMP', 'Prog_Mk','Prog_RBC'))]
cells.1 = c(cells.1, mutations.df.condensed$bc[which(mutations.df.condensed$sample == 'AML1002.1' & 
                                              mutations.df.condensed$predicted.celltype %in% c('CD4', 'CD8', 'NK cells'))])
cells.3 = mutations.df.condensed$bc[which(mutations.df.condensed$sample == 'AML1002.3' & 
                                            mutations.df.condensed$predicted.celltype %in% c('HSC', 'LMPP', 'GMP', 'Prog_Mk','Prog_RBC'))]
cells.3 = c(cells.3, mutations.df.condensed$bc[which(mutations.df.condensed$sample == 'AML1002.3' & 
                                              mutations.df.condensed$predicted.celltype %in% c('CD4', 'CD8', 'NK cells'))])

cells.3 = cells.3[which(cells.3 %in% cnv.data$bc)]

mutations.df.condensed$TP53.mutated = ifelse(mutations.df.condensed$TP53.1 == 'mutated' | mutations.df.condensed$TP53.2 == 'mutated', 'mutated', 'wildtype')
ggplot(mutations.df.condensed[which(mutations.df.condensed$predicted.celltype %in% c('NK cells')),], aes(x=`7b`, y=NCAM1)) + geom_point(aes(color=TP53.mutated))
ggplot(mutations.df.condensed[which(mutations.df.condensed$predicted.celltype %in% c('CD4')),], aes(x=`7b`, y=CD4)) + geom_point(aes(color=TP53.mutated))

ha = columnAnnotation(celltype = mutations.df.condensed[cells.1, 'predicted.celltype'],
                      IDH2 = mutations.df.condensed[cells.1, 'IDH2'],
                      TP53.1 = mutations.df.condensed[cells.1, 'TP53.1'],
                      TP53.2 = mutations.df.condensed[cells.1, 'TP53.2'],
                      TET2 = mutations.df.condensed[cells.1, 'TET2'],
                      U2AF1 = mutations.df.condensed[cells.1, 'U2AF1'],
                      annotation_name_gp= gpar(fontsize = 8),
                      col = list('celltype' = c(AML.combined.colors, 'CD4' = 'lightblue', 'CD8' = 'blue', 'NK cells' = 'purple'),
                                 'IDH2' = c('mutated' = 'red', 'wildtype' = 'white'),
                                 'TP53.1' = c('mutated' = 'red', 'wildtype' = 'white'),
                                 'TP53.2' = c('mutated' = 'red', 'wildtype' = 'white'),
                                 'TET2' = c('mutated' = 'red', 'wildtype' = 'white'),
                                 'U2AF1' = c('mutated' = 'red', 'wildtype' = 'white')), 
                      na_col = 'grey90', border=T, simple_anno_size = unit(7, 'pt'))

col_fun = circlize::colorRamp2(breaks = seq(0,1,1/8), colors = BuenColors::jdb_palette(name = 'solar_rojos'))
h1=Heatmap(t(mutations.df.condensed[cells.1,c('7b')]), 
        column_split = factor(mutations.df.condensed[cells.1, 'predicted.celltype'], levels = c(names(AML.combined.colors), 'CD4', 'CD8', 'NK cells')), 
        cluster_rows = F, cluster_columns = F, show_column_names = F, border = T, top_annotation = ha, column_title_gp = gpar(fontsize=0),
        column_title_rot = 90, na_col = 'grey90', col = col_fun, row_labels = c('del(7q)'), 
        row_names_side = 'left', row_names_gp = gpar(fontsize=8), use_raster = T, raster_quality = 10)

col_fun = circlize::colorRamp2(breaks = seq(0,1,1/8), colors = BuenColors::jdb_palette(name = 'brewer_purple'))
h2=Heatmap(t(mutations.df.condensed[cells.1,marker.genes.RNA]), 
           column_split = factor(mutations.df.condensed[cells.1, 'predicted.celltype'], levels = c(names(AML.combined.colors), 'CD4', 'CD8', 'NK cells')), 
           cluster_rows = F, cluster_columns = F, show_column_names = F, border = T, column_title_side = 'bottom', 
           column_title_rot = 90, na_col = 'white', col = col_fun, 
           row_names_side = 'left', row_names_gp = gpar(fontsize=8), use_raster = T, raster_quality = 10)

svglite::svglite('./AML/figures/AML1002/heatmaps/20220930_AML1002_1.svg', width = 6, height = 3.8)
draw(h1 %v% h2) 
dev.off()

###

ha = columnAnnotation(celltype = mutations.df.condensed[cells.3, 'predicted.celltype'],
                      IDH2 = mutations.df.condensed[cells.3, 'IDH2'],
                      TP53.1 = mutations.df.condensed[cells.3, 'TP53.1'],
                      TP53.2 = mutations.df.condensed[cells.3, 'TP53.2'],
                      TET2 = mutations.df.condensed[cells.3, 'TET2'],
                      U2AF1 = mutations.df.condensed[cells.3, 'U2AF1'],
                      annotation_name_gp= gpar(fontsize = 8),
                      col = list('celltype' = c(AML.combined.colors, 'CD4' = 'lightblue', 'CD8' = 'blue', 'NK cells' = 'purple'),
                                 'IDH2' = c('mutated' = 'red', 'wildtype' = 'white'),
                                 'TP53.1' = c('mutated' = 'red', 'wildtype' = 'white'),
                                 'TP53.2' = c('mutated' = 'red', 'wildtype' = 'white'),
                                 'TET2' = c('mutated' = 'red', 'wildtype' = 'white'),
                                 'U2AF1' = c('mutated' = 'red', 'wildtype' = 'white')), 
                      na_col = 'grey90', border=T, simple_anno_size = unit(7, 'pt'))

col_fun = circlize::colorRamp2(breaks = seq(0,1,1/8), colors = BuenColors::jdb_palette(name = 'solar_rojos'))
h1=Heatmap(t(mutations.df.condensed[cells.3,c('7b')]), 
           column_split = factor(mutations.df.condensed[cells.3, 'predicted.celltype'], levels = c(names(AML.combined.colors), 'CD4', 'CD8', 'NK cells')), 
           cluster_rows = F, cluster_columns = F, show_column_names = F, border = T, top_annotation = ha, column_title_gp = gpar(fontsize=0),
           column_title_rot = 90, na_col = 'grey90', col = col_fun, row_labels = c('del(7q)'), 
           row_names_side = 'left', row_names_gp = gpar(fontsize=8), use_raster = T, raster_quality = 10)

col_fun = circlize::colorRamp2(breaks = seq(0,1,1/8), colors = BuenColors::jdb_palette(name = 'brewer_purple'))
h2=Heatmap(t(mutations.df.condensed[cells.3,marker.genes.RNA]), 
           column_split = factor(mutations.df.condensed[cells.3, 'predicted.celltype'], levels = c(names(AML.combined.colors), 'CD4', 'CD8', 'NK cells')), 
           cluster_rows = F, cluster_columns = F, show_column_names = F, border = T, column_title_side = 'bottom', 
           column_title_rot = 90, na_col = 'white', col = col_fun, 
           row_names_side = 'left', row_names_gp = gpar(fontsize=8), use_raster = T, raster_quality = 10)

svglite::svglite('./AML/figures/AML1002/heatmaps/20220930_AML1002_3.svg', width = 6, height = 3.8)
draw(h1 %v% h2) 
dev.off()

###

ha = columnAnnotation(celltype = mutations.df.condensed[c(cells.1, cells.3), 'predicted.celltype'],
                      IDH2 = mutations.df.condensed[c(cells.1, cells.3), 'IDH2'],
                      TP53.1 = mutations.df.condensed[c(cells.1, cells.3), 'TP53.1'],
                      TP53.2 = mutations.df.condensed[c(cells.1, cells.3), 'TP53.2'],
                      TET2 = mutations.df.condensed[c(cells.1, cells.3), 'TET2'],
                      U2AF1 = mutations.df.condensed[c(cells.1, cells.3), 'U2AF1'],
                      annotation_name_gp= gpar(fontsize = 8),
                      col = list('celltype' = c(AML.combined.colors, 'CD4' = 'lightblue', 'CD8' = 'blue', 'NK cells' = 'purple'),
                                 'IDH2' = c('mutated' = 'red', 'wildtype' = 'white'),
                                 'TP53.1' = c('mutated' = 'red', 'wildtype' = 'white'),
                                 'TP53.2' = c('mutated' = 'red', 'wildtype' = 'white'),
                                 'TET2' = c('mutated' = 'red', 'wildtype' = 'white'),
                                 'U2AF1' = c('mutated' = 'red', 'wildtype' = 'white')), 
                      na_col = 'grey90', border=T, simple_anno_size = unit(7, 'pt'))

col_fun = circlize::colorRamp2(breaks = seq(0,1,1/8), colors = BuenColors::jdb_palette(name = 'solar_rojos'))
h1=Heatmap(t(mutations.df.condensed[c(cells.1, cells.3),c('7b')]), 
           column_split = factor(mutations.df.condensed[c(cells.1, cells.3), 'predicted.celltype'], levels = c(names(AML.combined.colors), 'CD4', 'CD8', 'NK cells')), 
           cluster_rows = F, cluster_columns = F, show_column_names = F, border = T, top_annotation = ha, column_title_gp = gpar(fontsize=0),
           column_title_rot = 90, na_col = 'grey90', col = col_fun, row_labels = c('del(7q)'), 
           row_names_side = 'left', row_names_gp = gpar(fontsize=8), use_raster = T, raster_quality = 10)

col_fun = circlize::colorRamp2(breaks = seq(0,4,4/8), colors = BuenColors::jdb_palette(name = 'brewer_purple'))
h2=Heatmap(t(mutations.df.condensed[c(cells.1, cells.3),marker.genes.RNA]), 
           column_split = factor(mutations.df.condensed[c(cells.1, cells.3), 'predicted.celltype'], levels = c(names(AML.combined.colors), 'CD4', 'CD8', 'NK cells')), 
           cluster_rows = F, cluster_columns = F, show_column_names = F, border = T, column_title_side = 'bottom', 
           column_title_rot = 90, na_col = 'white', col = col_fun, 
           row_names_side = 'left', row_names_gp = gpar(fontsize=8), use_raster = T, raster_quality = 10)

svglite::svglite('./AML/figures/AML1002/heatmaps/20220930_AML1002_13.svg', width = 8, height = 3)
draw(h1 %v% h2, use_raster = T, raster_quality = 10) 
dev.off()

###
boo = 
  mutations.df.condensed %>% group_by(predicted.celltype) %>% 
  summarize(IDH2.freq = length(which(IDH2 == 'mutated')) / length(which(!is.na(IDH2))),
            IDH2.count = length(which(IDH2 == 'mutated')),
            TET2.freq = length(which(TET2 == 'mutated')) / length(which(!is.na(TET2))),
            TET2.count = length(which(IDH2 == 'mutated')),
            U2AF1.freq = length(which(U2AF1 == 'mutated')) / length(which(!is.na(U2AF1))),
            U2AF1.count = length(which(U2AF1 == 'mutated')),
            TP53.1.freq = length(which(TP53.1 == 'mutated')) / length(which(!is.na(TP53.1))),
            TP53.1.count = length(which(TP53.1 == 'mutated')),
            TP53.2.freq = length(which(TP53.2 == 'mutated')) / length(which(!is.na(TP53.2))),
            TP53.2.count = length(which(TP53.2 == 'mutated'))) %>%
  tidyr::pivot_longer(names_to = 'mutation', cols = c(IDH2.freq, TET2.freq, U2AF1.freq, TP53.1.freq, TP53.2.freq))
boo[is.na(boo)] = 0
ggplot(boo[which(boo$predicted.celltype %in% c('HSC', 'LMPP', 'GMP', 'Prog_Mk','Prog_RBC', 'CD4', 'CD8', 'NK cells')),], 
       aes(y=mutation, x=100*value, fill=factor(predicted.celltype, levels = rev(c('HSC', 'LMPP', 'GMP', 'Prog_Mk','Prog_RBC', 'CD4', 'CD8', 'NK cells'))))) +
  geom_col(position = 'dodge') +
  scale_y_discrete(limits = rev(c('IDH2.freq', 'TP53.1.freq', 'TP53.2.freq', 'TET2.freq', 'U2AF1.freq')),
                   labels = rev(c('IDH2R140Q', 'TP53H179R', 'TP53P278S', 'TET2I1873T', 'U2AF1S34Y'))) + 
  scale_x_continuous('% mutated') + 
  scale_fill_manual(values = c(AML.combined.colors, 'CD4' = 'lightblue', 'CD8' = 'blue', 'NK cells' = 'purple')) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'),
        #axis.text.x = element_text('Arial', size=10, angle=90, hjust=1, vjust=0.5),
        axis.title.y = element_blank())
ggsave('./AML/figures/AML1002/plots/20221005_AML1002_celltype_mutation_freq.svg', width = 2.5, height = 2.5)
