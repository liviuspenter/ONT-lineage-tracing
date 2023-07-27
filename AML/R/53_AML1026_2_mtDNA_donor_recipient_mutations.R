# analysis of mtDNA data: donor-recipient deconvolution

library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(Seurat)
library(SummarizedExperiment)
library(dplyr)

# functionality from MAESTER (https://github.com/petervangalen/MAESTER-2021)
source('./R/210215_FunctionsGeneral.R')
source('./R/celltypes.R')

AML.all = readRDS('./data/AML/objects/AML.all.rds')
AML1026.2 = readRDS('./data/AML/objects/20220711_AML1026_2.rds')

# get mtDNA variants
magtk.output = readRDS('./data/AML/mtDNA/AML1026_2.maegatk/AML1026_2.rds')
af.dm = data.matrix(computeAFMutMatrix(magtk.output))*100

rownames(af.dm) = gsub(rownames(af.dm), pattern = '_', replacement = '')
af.dm.clean = af.dm[-which(rowSums(af.dm) == 0),]


colnames(af.dm.clean) = paste0('AML1026.2_', colnames(af.dm.clean), '-1')
TNK.cells = colnames(AML1026.2)[which(AML1026.2$predicted.celltype %in% TNK.clusters)]
TNK.cells = TNK.cells[which(TNK.cells %in% colnames(af.dm.clean))]

mtDNA.variants = names(sort(rowMeans(af.dm.clean[, TNK.cells]), decreasing = T))[1:50]
col_fun = circlize::colorRamp2(breaks = seq(0,10, 10/8), colors = BuenColors::jdb_palette(name = 'solar_rojos'))
Heatmap(af.dm.clean[mtDNA.variants, TNK.cells], show_row_dend = F, show_column_dend = F, show_row_names = T, row_names_gp = gpar(fontsize=6),
        show_column_names = F, col = col_fun)

#AML1026.2 = AML1026.2[intersect(colnames(AML1026.2), colnames(af.dm.clean))]
AML1026.2 = subset(AML1026.2, cells = intersect(colnames(AML1026.2), colnames(af.dm.clean)))
AML1026.2$boo = af.dm.clean['9719C>T', intersect(colnames(AML1026.2), colnames(af.dm.clean))]

mtDNA.df = t(af.dm.clean)
mtDNA.df = mtDNA.df[intersect(colnames(AML1026.2), rownames(mtDNA.df)),]
cells = intersect(colnames(AML1026.2), rownames(mtDNA.df))
cells = cells[which(cells %in% rownames(AML1026.genotype)[which(AML1026.genotype$chimerism == 'recipient')])]
mtDNA.df = cbind(mtDNA.df[cells,], AML1026.2@reductions$ref.umap@cell.embeddings[cells,])
mtDNA.df = as.data.frame(mtDNA.df)
mtDNA.df$predicted.celltype = AML1026.2$predicted.celltype[cells]

# find mutations unique to cell types
mutations = rownames(af.dm.clean)
mutations.df = data.frame()
for (celltype in unique(mtDNA.df$predicted.celltype)) {
  message(celltype)
  cells.celltype = rownames(mtDNA.df)[which(mtDNA.df$predicted.celltype == celltype)]
  boo = mtDNA.df[cells.celltype, mutations]
  mutations.df = rbind(mutations.df, colMeans(boo))
}
colnames(mutations.df) = mutations
rownames(mutations.df) = unique(mtDNA.df$predicted.celltype)
mutations.df = as.data.frame(t(mutations.df))
mutations.df$highlight = 'none'
mutations.df$highlight[which(rownames(mutations.df) %in% c('10685G>A', '3106C>A',   '3106C>G', '3106C>T', '15615G>A','9254A>G'))] = 'highlight'
mutations.df$label = rownames(mutations.df)

# HSC: 10685G>A, 3106C>A, 9663G>C
ggplot() +
  scale_x_log10('% HSC',limits = c(1,100)) + 
  scale_y_log10('% recipient CD4+ T cells',limits = c(1,100)) + 
  scale_color_manual(values = c('none' = 'grey', 'highlight' = 'firebrick')) + 
  geom_abline(slope = 1) + 
  geom_point(data=mutations.df, aes(x=`CD4 Memory`, y=HSC, color=highlight), size=0.5) +
  geom_label_repel(data=mutations.df[which(mutations.df$highlight == 'highlight'),], aes(x=`CD4 Memory`, y=HSC, label=label), size=2, label.size = 0) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./AML/figures/mtDNA/plots/20220804_AML1026_HSC_CD4_Tcells_mtDNA_mutations.svg', width = 2.5, height = 2.5)




HSC.cells = colnames(AML1026.2)[which(AML1026.2$predicted.celltype == 'HSC' & 
                                        colnames(AML1026.2) %in% AML1026.genotype$barcode[which(AML1026.genotype$chimerism == 'recipient')])]
HSC.cells = HSC.cells[which(HSC.cells %in% colnames(af.dm.clean))]
LMPP.cells = colnames(AML1026.2)[which(AML1026.2$predicted.celltype == 'LMPP' & 
                                         colnames(AML1026.2) %in% AML1026.genotype$barcode[which(AML1026.genotype$chimerism == 'recipient')])]
LMPP.cells = LMPP.cells[which(LMPP.cells %in% colnames(af.dm.clean))]
GMP.cells = colnames(AML1026.2)[which(AML1026.2$predicted.celltype == 'GMP' & 
                                        colnames(AML1026.2) %in% AML1026.genotype$barcode[which(AML1026.genotype$chimerism == 'recipient')])]
GMP.cells = GMP.cells[which(GMP.cells %in% colnames(af.dm.clean))]
CD4.recipient = colnames(AML1026.2)[which(AML1026.2$predicted.celltype == 'CD4 Memory' & 
                                            colnames(AML1026.2) %in% AML1026.genotype$barcode[which(AML1026.genotype$chimerism == 'recipient')])]
CD4.recipient = CD4.recipient[which(CD4.recipient %in% colnames(af.dm.clean))]
NK.recipient = colnames(AML1026.2)[which(AML1026.2$predicted.celltype == 'NK' & 
                                           colnames(AML1026.2) %in% AML1026.genotype$barcode[which(AML1026.genotype$chimerism == 'recipient')])]
NK.recipient = NK.recipient[which(NK.recipient %in% colnames(af.dm.clean))]
CD4.donor = colnames(AML1026.2)[which(AML1026.2$predicted.celltype == 'CD4 Memory' & 
                                        colnames(AML1026.2) %in% AML1026.genotype$barcode[which(AML1026.genotype$chimerism == 'donor')])]
CD4.donor = CD4.donor[which(CD4.donor %in% colnames(af.dm.clean))]
NK.donor = colnames(AML1026.2)[which(AML1026.2$predicted.celltype == 'NK' & 
                                       colnames(AML1026.2) %in% AML1026.genotype$barcode[which(AML1026.genotype$chimerism == 'donor')])]
NK.donor = NK.donor[which(NK.donor %in% colnames(af.dm.clean))]


ha = columnAnnotation(celltype = c(rep('HSC', length(HSC.cells)),
                                   rep('LMPP', length(LMPP.cells)),
                                   rep('GMP', length(GMP.cells)),
                                   rep('CD4.recipient', length(CD4.recipient)),
                                   rep('NK.recipient', length(NK.recipient)),
                                   rep('CD4.donor', length(CD4.donor)),
                                   rep('NK.donor', length(NK.donor))),
                      col = list(celltype = c('HSC' = as.character(AML.combined.colors['HSC']),
                                              'LMPP' = as.character(AML.combined.colors['LMPP']),
                                              'GMP' = as.character(AML.combined.colors['GMP']),
                                              'CD4.recipient' = as.character(AML.combined.colors['CD4 Memory']),
                                              'NK.recipient' = as.character(AML.combined.colors['NK']),
                                              'CD4.donor' = as.character(AML.combined.colors['CD4 Memory']),
                                              'NK.donor' = as.character(AML.combined.colors['NK']))),
                      simple_anno_size = unit(5, 'pt'), border=T)
col_fun = circlize::colorRamp2(breaks = seq(0,100,100/8), colors = BuenColors::jdb_palette(name = 'solar_rojos'))
svglite::svglite('./AML/figures/mtDNA/heatmaps/20220804_AML1026_mtDNA.svg', width = 5.5, height = 2.5)
Heatmap(af.dm.clean[c('10685G>A', '3106C>A',   '3106C>G', '3106C>T', '15615G>A','9254A>G', donor.variants, recipient.variants),
                    c(HSC.cells, GMP.cells, LMPP.cells, CD4.recipient, NK.recipient, CD4.donor, NK.donor)], top_annotation = ha, 
        column_split = factor(c(rep('HSC', length(HSC.cells)),
                                rep('LMPP', length(LMPP.cells)),
                                rep('GMP', length(GMP.cells)),
                                rep('CD4.recipient', length(CD4.recipient)),
                                rep('NK.recipient', length(NK.recipient)),
                                rep('CD4.donor', length(CD4.donor)),
                                rep('NK.donor', length(NK.donor))), 
                              levels = c('HSC', 'LMPP', 'GMP', 'CD4.recipient', 'NK.recipient', 'CD4.donor', 'NK.donor')),
        show_column_names = F, cluster_columns = F, cluster_rows = F, row_names_side = 'left', row_names_gp = gpar(fontsize=8),
        column_title_side = 'bottom', column_title_rot = 90, column_title_gp = gpar(fontsize=8), 
        col = col_fun, border=T, use_raster = T, raster_quality = 10)
dev.off()

