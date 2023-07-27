# reanalysis of beat AML data 

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(Seurat)

### van Galen signature in AML scRNA-seq dataset
AML.all = readRDS('./data/AML/objects/AML.all.rds')
signature.van.Galen = readxl::read_excel('./data/AML/signatures/20211013_AML_signature_van_galen.xlsx')

AML.all = AddModuleScore(AML.all, features = list('vG.HSC' = signature.van.Galen$HSC), ctrl = 5, name = 'vG.HSC')
AML.all = AddModuleScore(AML.all, features = list('vG.progenitor' = signature.van.Galen$progenitor), ctrl = 5, name = 'vG.progenitor')
AML.all = AddModuleScore(AML.all, features = list('vG.GMP' = signature.van.Galen$GMP), ctrl = 5, name = 'vG.GMP')
AML.all = AddModuleScore(AML.all, features = list('vG.promono' = signature.van.Galen$promono), ctrl = 5, name = 'vG.promono')
AML.all = AddModuleScore(AML.all, features = list('vG.mono' = signature.van.Galen$mono), ctrl = 5, name = 'vG.mono')
AML.all = AddModuleScore(AML.all, features = list('vG.cDC' = signature.van.Galen$cDC), ctrl = 5, name = 'vG.cDC')

AML.all = AddModuleScore(AML.all, features = list('Erythroid' = c('GATA1', 'CA1', 'HBA1', 'HBB', 'ALAD')), ctrl = 5, name = 'Erythroid')
AML.all = AddModuleScore(AML.all, features = list('Megakaryocytic' = c('GATA2', 'ZFPM1', 'FLI1', 'NFE2', 'PF4')), ctrl = 5, name = 'Megakaryocytic')

p=FeaturePlot(AML.all, reduction = 'ref.umap', features = c('vG.HSC1'), min.cutoff = 'q05', max.cutoff = 'q95') + 
  NoAxes() + NoLegend() + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))
ggsave('./AML/figures/beataml/UMAP/20221017_AML_all_HSC.png', width = 4, height = 4, dpi = 600, plot = p)

p=FeaturePlot(AML.all, reduction = 'ref.umap', features = c('vG.progenitor1'), min.cutoff = 'q05', max.cutoff = 'q95') + 
  NoAxes() + NoLegend() + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))
ggsave('./AML/figures/beataml/UMAP/20221017_AML_all_progenitor.png', width = 4, height = 4, dpi = 600, plot = p)

p=FeaturePlot(AML.all, reduction = 'ref.umap', features = c('vG.GMP1'), min.cutoff = 'q05', max.cutoff = 'q95') + 
  NoAxes() + NoLegend() + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))
ggsave('./AML/figures/beataml/UMAP/20221017_AML_all_GMP.png', width = 4, height = 4, dpi = 600, plot = p)

p=FeaturePlot(AML.all, reduction = 'ref.umap', features = c('vG.promono1'), min.cutoff = 'q05', max.cutoff = 'q95') + 
  NoAxes() + NoLegend() + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))
ggsave('./AML/figures/beataml/UMAP/20221017_AML_all_Promono.png', width = 4, height = 4, dpi = 600, plot = p)

p=FeaturePlot(AML.all, reduction = 'ref.umap', features = c('vG.mono1'), min.cutoff = 'q05', max.cutoff = 'q95') + 
  NoAxes() + NoLegend() + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))
ggsave('./AML/figures/beataml/UMAP/20221017_AML_all_Mono.png', width = 4, height = 4, dpi = 600, plot = p)

p=FeaturePlot(AML.all, reduction = 'ref.umap', features = c('vG.cDC1'), min.cutoff = 'q05', max.cutoff = 'q95') + 
  NoAxes() + NoLegend() + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))
ggsave('./AML/figures/beataml/UMAP/20221017_AML_all_cDC.png', width = 4, height = 4, dpi = 600, plot = p)

p=FeaturePlot(AML.all, reduction = 'ref.umap', features = c('Erythroid1'), min.cutoff = 'q05', max.cutoff = 'q95') + 
  NoAxes() + NoLegend() + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))
ggsave('./AML/figures/beataml/UMAP/20221017_AML_all_Erythroid.png', width = 4, height = 4, dpi = 600, plot = p)

p=FeaturePlot(AML.all, reduction = 'ref.umap', features = c('Megakaryocytic1'), min.cutoff = 'q05', max.cutoff = 'q95') + 
  NoAxes() + NoLegend() + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))
ggsave('./AML/figures/beataml/UMAP/20221017_AML_all_Megakaryocytic.png', width = 4, height = 4, dpi = 600, plot = p)

### load beataml dataset
beataml.RNA = as.data.frame(data.table::fread('./data/AML/beataml/beataml_waves1to4_norm_exp_dbgap.txt'))
beataml.meta = as.data.frame(readxl::read_excel('./data/AML/beataml/beataml_wv1to4_clinical.xlsx'))
beataml.meta = beataml.meta[which(beataml.meta$dbgap_rnaseq_sample %in% intersect(beataml.meta$dbgap_rnaseq_sample, colnames(beataml.RNA))),]
rownames(beataml.meta) = beataml.meta$dbgap_rnaseq_sample

# exclude remission samples
beataml.meta = beataml.meta[-which(beataml.meta$diseaseStageAtSpecimenCollection == 'Remission'),]

rownames(beataml.RNA) = beataml.RNA$display_label
beataml.RNA = beataml.RNA[,-c(1,2,3,4)]

# use only samples with metadata
beataml.RNA = beataml.RNA[,rownames(beataml.meta)]
beataml.RNA.t = as.data.frame(t(beataml.RNA))[rownames(beataml.meta),]

combined.data = cbind(beataml.meta, beataml.RNA.t[,c('HBA1', 'HBB', 'CD34', 'CD14', 'MECOM')])

# calculate van Galen score for beataml
combined.data$vG.HSC = rowMeans(beataml.RNA.t[,signature.van.Galen$HSC[which(signature.van.Galen$HSC %in% colnames(beataml.RNA.t))]])
combined.data$vG.progenitor = rowMeans(beataml.RNA.t[,signature.van.Galen$progenitor[which(signature.van.Galen$progenitor %in% colnames(beataml.RNA.t))]])
combined.data$vG.GMP = rowMeans(beataml.RNA.t[,signature.van.Galen$GMP[which(signature.van.Galen$GMP %in% colnames(beataml.RNA.t))]])
combined.data$vG.promono = rowMeans(beataml.RNA.t[,signature.van.Galen$promono[which(signature.van.Galen$promono %in% colnames(beataml.RNA.t))]])
combined.data$vG.mono = rowMeans(beataml.RNA.t[,signature.van.Galen$mono[which(signature.van.Galen$mono %in% colnames(beataml.RNA.t))]])
combined.data$vG.cDC = rowMeans(beataml.RNA.t[,signature.van.Galen$cDC[which(signature.van.Galen$cDC %in% colnames(beataml.RNA.t))]])
combined.data$Erythroid = rowMeans(beataml.RNA.t[,c('GATA1', 'CA1', 'HBA1', 'HBB', 'ALAD')])
combined.data$Megakaryocytic = rowMeans(beataml.RNA.t[,c('GATA2', 'ZFPM1', 'FLI1', 'NFE2', 'PF4')])

ggplot(data=combined.data[which(!combined.data$specificDxAtAcquisition %in% c('unknown','Unknown') ),], 
       aes(x=reorder(specificDxAtAcquisition, vG.HSC, FUN = median), y=vG.HSC)) + geom_boxplot() + geom_jitter(size=0.5) + 
  scale_y_continuous('HSC score') + 
  theme_classic() +
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.title.x = element_blank())
ggsave('./AML/figures/beataml/plots/20221016_HSC_score.svg', width = 3.5, height = 5)

ggplot(data=combined.data[which(!combined.data$specificDxAtAcquisition %in% c('unknown','Unknown') ),], 
       aes(x=reorder(specificDxAtAcquisition, vG.progenitor, FUN = median), y=vG.progenitor)) + geom_boxplot() + geom_jitter(size=0.5) + 
  scale_y_continuous('Progenitor score') + 
  theme_classic() +
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.title.x = element_blank())
ggsave('./AML/figures/beataml/plots/20221016_Progenitor_score.svg', width = 3.5, height = 5)

ggplot(data=combined.data[which(!combined.data$specificDxAtAcquisition %in% c('unknown','Unknown') ),], 
       aes(x=reorder(specificDxAtAcquisition, vG.GMP, FUN = median), y=vG.GMP)) + geom_boxplot() + geom_jitter(size=0.5) + 
  scale_y_continuous('GMP score') + 
  theme_classic() +
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.title.x = element_blank())
ggsave('./AML/figures/beataml/plots/20221016_GMP_score.svg', width = 3.5, height = 5)

ggplot(data=combined.data[which(!combined.data$specificDxAtAcquisition %in% c('unknown','Unknown') ),], 
       aes(x=reorder(specificDxAtAcquisition, vG.promono, FUN = median), y=vG.promono)) + geom_boxplot() + geom_jitter(size=0.5) + 
  scale_y_continuous('Promono score') + 
  theme_classic() +
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.title.x = element_blank())
ggsave('./AML/figures/beataml/plots/20221016_Promono_score.svg', width = 3.5, height = 5)

ggplot(data=combined.data[which(!combined.data$specificDxAtAcquisition %in% c('unknown','Unknown') ),], 
       aes(x=reorder(specificDxAtAcquisition, vG.mono, FUN = median), y=vG.mono)) + geom_boxplot() + geom_jitter(size=0.5) + 
  scale_y_continuous('Mono score') + 
  theme_classic() +
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.title.x = element_blank())
ggsave('./AML/figures/beataml/plots/20221016_Mono_score.svg', width = 3.5, height = 5)

ggplot(data=combined.data[which(!combined.data$specificDxAtAcquisition %in% c('unknown','Unknown') ),], 
       aes(x=reorder(specificDxAtAcquisition, vG.cDC, FUN = median), y=vG.cDC)) + geom_boxplot() + geom_jitter(size=0.5) + 
  scale_y_continuous('cDC score') + 
  theme_classic() +
  theme(axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.title.x = element_blank())
ggsave('./AML/figures/beataml/plots/20221016_cDC_score.svg', width = 3.5, height = 5)

ggplot(data=combined.data[which(!combined.data$specificDxAtAcquisition %in% c('unknown','Unknown') ),], 
       aes(x=reorder(specificDxAtAcquisition, Erythroid, FUN = median), y=Erythroid)) + geom_boxplot(outlier.size = 0) + geom_jitter(size=0.5, aes(color=Erythroid)) + 
  scale_y_continuous('Erythroid score') + 
  scale_color_gradientn(colors = BuenColors::jdb_palette(name = 'brewer_red')) + 
  theme_classic() +
  theme(legend.position = 'none', 
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.title.x = element_blank())
ggsave('./AML/figures/beataml/plots/20221016_Erythroid_score.svg', width = 3.5, height = 5)

ggplot(data=combined.data[which(!combined.data$specificDxAtAcquisition %in% c('unknown','Unknown') ),], 
       aes(x=reorder(specificDxAtAcquisition, Megakaryocytic, FUN = median), y=Megakaryocytic)) + geom_boxplot(outlier.size = 0) + geom_jitter(size=0.5, aes(color=Megakaryocytic)) + 
  scale_y_continuous('Megakaryocytic score') + 
  scale_color_gradientn(colors = BuenColors::jdb_palette(name = 'brewer_purple')) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
        axis.title.x = element_blank())
ggsave('./AML/figures/beataml/plots/20221016_Megakaryocytic_score.svg', width = 3.5, height = 5)


mat = as.data.frame(combined.data %>% filter(!specificDxAtAcquisition %in% c('unknown', 'Unknown')) %>%
                      group_by(specificDxAtAcquisition) %>% 
                      summarize(vG.HSC = mean(vG.HSC),
                                vG.progenitor = mean(vG.progenitor),
                                vG.GMP = mean(vG.GMP),
                                vG.promono = mean(vG.promono),
                                vG.mono = mean(vG.mono),
                                vG.cDC = mean(vG.cDC),
                                Erythroid = mean(Erythroid),
                                Megakaryocytic = mean(Megakaryocytic)))
rownames(mat) = mat$specificDxAtAcquisition
mat = mat[,-1]

col_fun = circlize::colorRamp2(breaks = seq(-4, 4, 8/8), colors = BuenColors::jdb_palette(name = 'brewer_yes'))
svglite::svglite('./AML/figures/beataml/heatmaps/20221017_beataml_expression_profiles.svg', width = 4.5, height = 4)
Heatmap(scale(mat), cluster_columns = F, cluster_rows = T, show_row_dend = F, row_names_side = 'left', border=T,
        row_names_gp = gpar(fontsize=8), column_names_gp = gpar(fontsize=8), col = col_fun,
        column_split = factor(c(rep('vG', 6), rep('Penter', 2)), levels = c('vG', 'Penter')), column_title_gp = gpar(fontsize=0),
        column_labels = c('HSC', 'Progenitor', 'GMP', 'Promono', 'Mono', 'cDC', 'Erythroid', 'Megakaryocytic'))
dev.off()

# analyses that didn't make it into the paper

ggplot(data=combined.data, aes(x=as.numeric(hemoglobin), y=as.numeric(HBB))) + geom_point(aes(color=priorMDS), size=0.5)  +
  scale_x_continuous('hemoglobin g/dl') + 
  scale_y_continuous('HBB expression') + 
  scale_color_manual(values = c('y' = 'black', 'n' = 'grey')) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./AML/figures/beataml/plots/20221017_beataml_HBB_hemoglobin.svg', width = 1.5, height = 1.5)

ggplot(data=combined.data[!is.na(combined.data$`%.Blasts.in.PB`),], aes(x=as.numeric(`%.Blasts.in.PB`), y=as.numeric(HBB))) + geom_point(aes(color=priorMDS), size=0.5)  +
  scale_x_continuous('% blasts in PB') + 
  scale_y_continuous('HBB expression') + 
  scale_color_manual(values = c('y' = 'black', 'n' = 'grey')) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./AML/figures/beataml/plots/20221017_beataml_HBB_blasts_PB.svg', width = 1.5, height = 1.5)

ggplot(data=combined.data[!is.na(combined.data$`%.Blasts.in.BM`),], aes(x=as.numeric(`%.Blasts.in.BM`), y=as.numeric(HBB))) + geom_point(aes(color=priorMDS), size=0.5)  +
  scale_x_continuous('% blasts in BM') + 
  scale_y_continuous('HBB expression') + 
  scale_color_manual(values = c('y' = 'black', 'n' = 'grey')) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./AML/figures/beataml/plots/20221017_beataml_HBB_blasts_BM.svg', width = 1.5, height = 1.5)

ggplot(data=combined.data[!is.na(combined.data$wbcCount),], aes(x=as.numeric(`%.Blasts.in.BM`), y=as.numeric(HBB))) + geom_point(aes(color=priorMDS), size=0.5)  +
  scale_x_continuous('WBC /nl') + 
  scale_y_continuous('HBB expression') + 
  scale_color_manual(values = c('y' = 'black', 'n' = 'grey')) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./AML/figures/beataml/plots/20221017_beataml_HBB_WBC.svg', width = 1.5, height = 1.5)

so = CreateSeuratObject(counts = beataml.RNA)
so = FindVariableFeatures(so)
so = NormalizeData(so)
so = ScaleData(so)
so = RunPCA(so)
so = RunUMAP(so, dims = 1:10)
so = FindNeighbors(so)
so = FindClusters(so, resolution = 1)

so = AddMetaData(so, metadata = beataml.meta)

so = AddModuleScore(so, features = list('vG.HSC' = c('NPTX2', 'H1F0', 'EMP1', 'CALCRL', 'TPSD1', 'TPT1', 'CRHBP', 'TSC22D1', 'DST', 
                                                   'NRIP1', 'ABCB1', 'ZBTB20', 'TPSB2', 'KMT2A', 'MEF2C', 'ST3GAL1', 'TMEM25', 'C20orf203', 'GNG11', 'HOPX')), ctrl = 5, name = 'vG.HSC')
so = AddModuleScore(so, features = list('vG.progenitor' = intersect(signature.van.Galen$progenitor, rownames(so))), ctrl = 5, name = 'vG.progenitor')
so = AddModuleScore(so, features = list('vG.GMP' = c('PRTN3', 'MPO', 'CALR', 'CLEC5A', 'ELANE', 'TRH', 'CEBPE', 'NUCB2', 'CSF3R', 'CD38', 'IGFBP2',
                                                     'PRRT4', 'SNHG5', 'FABP5', 'CLEC11A', 'SERPINB1', 'AZU1', 'HNRNPDL', 'HSPB1', 'C12orf57', 'FGFR1')), ctrl = 5, name = 'vG.GMP')
so = AddModuleScore(so, features = list('vG.promono' = c('DEFB1', 'RNASE2', 'MS4A3', 'SERPINB10', 'SESN3', 'ZFR', 'MRPL33', 'CTSG', 'SLC44A1', 'SLPI', 'FUT4',
                                                         'SRGN', 'CD70', 'PLD3', 'RETN', 'TP53INP2', 'HSPA5', 'RNASE3', 'CCL23', 'EMB', 'ATP8B4', 'CLU',
                                                         'FAM107B', 'KBTBD11', 'CSTA', 'ANKRD28', 'PIWIL4', 'RNVU1-6')), ctrl = 5, name = 'vG.promono')
so = AddModuleScore(so, features = list('vG.mono' = c('FCN1', 'S100A12', 'MAFB', 'S100A9', 'PLBD1', 'SERPINA1', 'BCL2A1', 'THBS1', 'PSAP', 'S100A8', 'FPR1', 'C5AR1',
                                                   'CD14', 'NAMPT', 'VNN2', 'CTSS', 'DUSP1', 'CEBPB', 'CR1', 'NFKBIA', 'SLC11A1', 'LILRB3', 'BCL6', 
                                                   'CYP1B1', 'TNFAIP2', 'AQP9', 'TLR4', 'APOBEC3A')), ctrl = 5, nbin = 5, name = 'vG.mono')
so = AddModuleScore(so, features = list('vG.cDC' = c('HLA-DRB5', 'CST3', 'SAMHD1', 'NAPSB', 'FCER1A', 'HLA-DRB1', 'PKIB', 'HLA-DRA', 'HLA-DRB6', 'CPVL',
                                                     'HLA-DPB1', 'HLA-DQA1', 'HLA-DPA1', 'CLEC4A', 'TMSB10', 'CAP1', 'HLA-DQB1', 'CRIP1', 'CLEC10A', 'GPX1',
                                                     'ITGB7', 'HLA-DQB2', 'DBI', 'FTH1P3', 'ACTB', 'HLA-DQA2', 'S100B', 'ALDH2')), ctrl = 5, name = 'vG.cDC')
so = AddModuleScore(so, features = list('Erythroid' = c('GATA1', 'HBA1', 'HBB', 'CA1', 'ALAD')), ctrl = 5, name = 'Erythroid')
so = AddModuleScore(so, features = list('Megakaryocytic' = c('GATA2', 'ZFPM1', 'FLI1', 'NFE2', 'PF4')), ctrl = 5, name = 'Megakaryocytic')

saveRDS('./data/beataml/20221016_beataml.rds', object = so)

p=FeaturePlot(so, c('vG.HSC1'), min.cutoff = 'q02', max.cutoff = 'q98') + NoAxes() + NoLegend() + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))
ggsave('./AML/figures/beataml/UMAP/20221016_beataml_UMAP_HSC.png', width = 4, height = 4, dpi = 600, plot = p)

p=FeaturePlot(so, c('vG.progenitor1'), min.cutoff = 'q02', max.cutoff = 'q98') + NoAxes() + NoLegend() + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))
ggsave('./AML/figures/beataml/UMAP/20221016_beataml_UMAP_Progenitor.png', width = 4, height = 4, dpi = 600, plot = p)

p=FeaturePlot(so, c('vG.GMP1'), min.cutoff = 'q02', max.cutoff = 'q98') + NoAxes() + NoLegend() + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))
ggsave('./AML/figures/beataml/UMAP/20221016_beataml_UMAP_GMP.png', width = 4, height = 4, dpi = 600, plot = p)

p=FeaturePlot(so, c('vG.promono1'), min.cutoff = 'q02', max.cutoff = 'q98') + NoAxes() + NoLegend() + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))
ggsave('./AML/figures/beataml/UMAP/20221016_beataml_UMAP_Promono.png', width = 4, height = 4, dpi = 600, plot = p)

p=FeaturePlot(so, c('vG.mono1'), min.cutoff = 'q02', max.cutoff = 'q98') + NoAxes() + NoLegend() + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))
ggsave('./AML/figures/beataml/UMAP/20221016_beataml_UMAP_Mono.png', width = 4, height = 4, dpi = 600, plot = p)

p=FeaturePlot(so, c('vG.cDC1'), min.cutoff = 'q02', max.cutoff = 'q98') + NoAxes() + NoLegend() + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))
ggsave('./AML/figures/beataml/UMAP/20221016_beataml_UMAP_cDC.png', width = 4, height = 4, dpi = 600, plot = p)

p=FeaturePlot(so, c('Erythroid1'), min.cutoff = 'q02', max.cutoff = 'q98') + NoAxes() + NoLegend() + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))
ggsave('./AML/figures/beataml/UMAP/20221016_beataml_UMAP_Erythroid.png', width = 4, height = 4, dpi = 600, plot = p)

p=FeaturePlot(so, c('Megakaryocytic1'), min.cutoff = 'q02', max.cutoff = 'q98') + NoAxes() + NoLegend() + theme(plot.title = element_blank()) + 
  scale_color_gradientn(colours = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos')))
ggsave('./AML/figures/beataml/UMAP/20221016_beataml_UMAP_Megakaryocytic.png', width = 4, height = 4, dpi = 600, plot = p)


plot.list = list()
i = 1
for (dx in unique(so$specificDxAtAcquisition)) {
  if (length(which(so$specificDxAtAcquisition == dx)) < 4) {
    next()
  }
  p = DimPlot(so, cells.highlight = colnames(so)[which(so$specificDxAtAcquisition == dx)]) + 
    NoLegend() + NoAxes() + ggtitle(dx) + theme(plot.title = element_text('Arial', size=8, color='black', hjust=0.5))
  plot.list[[i]] = p
  i = i+1
  ggsave(paste0('./AML/figures/beataml/UMAP/20221016_UMAP_beataml_', dx, '.svg'), width = 4, height = 4, dpi = 600, plot = ggrastr::rasterize(p, dpi = 600))
}
cowplot::plot_grid(plotlist = plot.list)

marker.genes.RNA = c('CD34', 'HLA-DRA','CD33', 'KIT', 'MPO', 'CD38', 'CD14', 'GATA2','ZFPM1', 'FLI1', 'NFE2','PF4', 'GATA1','HBA1', 'HBB','CA1', 'ALAD','CD3D', 'CD4', 'CD8A', 'NCAM1', 'MKI67') 

col_fun = circlize::colorRamp2(breaks = seq(-5, 15, 20/8), colors = BuenColors::jdb_palette(name = 'brewer_yes'))
col_fun2 = circlize::colorRamp2(breaks = seq(0, 20, 20/8), colors = BuenColors::jdb_palette(name = 'brewer_red'))

ha = HeatmapAnnotation(
  hemoglobin = anno_simple(as.numeric(beataml.meta$hemoglobin), col = col_fun2, simple_anno_size = unit(5, 'pt'), border=T),
  annotation_name_side = "left", annotation_name_gp = gpar(fontsize=8))

lgd_hemoglobin = Legend(title = "hemoglobin", col_fun = col_fun2, at = c(0, 5, 10, 15, 20), 
                    labels = c("0", "5", "10", "15", "20"))

ht = Heatmap(beataml.RNA[marker.genes.RNA, beataml.meta$dbgap_rnaseq_sample], cluster_rows = F, show_row_dend = F, show_column_dend = F, show_row_names = T, show_column_names = F,
        top_annotation = ha, col = col_fun, row_names_side = 'left', row_names_gp = gpar(fontsize=8),
        border=T, raster_quality = 10, use_raster = T)

svglite::svglite('./AML/figures/beataml/heatmaps/20221017_beataml_gene_expression.svg', width = 4, height = 3)
draw(ht, annotation_legend_list = lgd_hemoglobin)
dev.off()

ht = Heatmap(beataml.RNA[marker.genes.RNA, beataml.meta$dbgap_rnaseq_sample], cluster_rows = F, show_row_dend = F, 
             show_column_dend = F, show_row_names = T, show_column_names = F, column_km = 3,
             top_annotation = ha, col = col_fun, row_names_side = 'left', row_names_gp = gpar(fontsize=8),
             border=T, raster_quality = 10, use_raster = T)

svglite::svglite('./AML/figures/beataml/heatmaps/20221024_beataml_gene_expression_clustered.svg', width = 4, height = 3)
draw(ht, annotation_legend_list = lgd_hemoglobin)
dev.off()


ha = columnAnnotation(specificDxAtInclusion = beataml.meta$specificDxAtInclusion)
Heatmap(beataml.RNA[c(signature.van.Galen$HSC, 
                      signature.van.Galen$progenitor,
                      signature.van.Galen$GMP,
                      signature.van.Galen$promono,
                      signature.van.Galen$mono,
                      signature.van.Galen$cDC), beataml.meta$dbgap_rnaseq_sample], cluster_rows = F, show_row_dend = F, show_column_dend = F, 
        show_row_names = T, show_column_names = F, top_annotation = ha, column_split = beataml.meta$specificDxAtAcquisition)

DimPlot(so, group.by = 'isDenovo', cols = c('FALSE' = 'grey', 'unknown' = 'grey', 'TRUE' = 'black'))
DimPlot(so, group.by = 'priorMDS', cols = c('n' = 'grey', 'y' = 'black'))
DimPlot(so, group.by = 'isTransformed', cols = c('FALSE' = 'grey', 'unknown' = 'grey', 'TRUE' = 'black'))

plot.list = list()
i = 1
for (dx in unique(so$specificDxAtAcquisition)) {
  if (length(which(so$specificDxAtAcquisition == dx)) < 4) {
    next()
  }
  p = DimPlot(so, cells.highlight = colnames(so)[which(so$specificDxAtAcquisition == dx)]) + 
    NoLegend() + NoAxes() + ggtitle(dx) + theme(plot.title = element_text('Arial', size=8, color='black', hjust=0.5))
  plot.list[[i]] = p
  i = i+1
}
cowplot::plot_grid(plotlist = plot.list)

svglite::svglite('./AML/figures/beataml/solar_rojos.svg', width = 2, height = 0.5)
BuenColors::jdb_palette(name = 'solar_rojos', type = 'continuous')
dev.off()

###
# correlation with inhibitor data
###

beataml.inhibitor.data = as.data.frame(data.table::fread('./data/beataml/beataml_probit_curve_fits_v4_dbgap.txt') )
beataml.inhibitor.data = beataml.inhibitor.data[which(beataml.inhibitor.data$dbgap_rnaseq_sample %in% 
                                                        intersect(beataml.inhibitor.data$dbgap_rnaseq_sample, combined.data$dbgap_rnaseq_sample)),]
boo = combined.data[intersect(beataml.inhibitor.data$dbgap_rnaseq_sample, combined.data$dbgap_rnaseq_sample),]

df = data.frame()
for (inhibitor in unique(beataml.inhibitor.data$inhibitor)) {
  message(inhibitor)
  inhibitor.subset = beataml.inhibitor.data[which(beataml.inhibitor.data$inhibitor == inhibitor),]
  if (nrow(inhibitor.subset) < 50) {
    next()
  }
  boo = combined.data[inhibitor.subset$dbgap_rnaseq_sample,]
  df = rbind(df, data.frame(inhibitor = inhibitor, 
             HSC = as.numeric(cor.test(inhibitor.subset$auc, boo$vG.HSC)$estimate),
             Progenitor = as.numeric(cor.test(inhibitor.subset$auc, boo$vG.progenitor)$estimate),
             GMP = as.numeric(cor.test(inhibitor.subset$auc, boo$vG.GMP)$estimate),
             Promono = as.numeric(cor.test(inhibitor.subset$auc, boo$vG.promono)$estimate),
             Mono = as.numeric(cor.test(inhibitor.subset$auc, boo$vG.mono)$estimate),
             cDC = as.numeric(cor.test(inhibitor.subset$auc, boo$vG.cDC)$estimate),
             Erythroid = as.numeric(cor.test(inhibitor.subset$auc, boo$Erythroid)$estimate),
             Megakaryocytic = as.numeric(cor.test(inhibitor.subset$auc, boo$Megakaryocytic)$estimate)))
}
rownames(df) = df$inhibitor
df = df[,-1]
col_fun = circlize::colorRamp2(breaks = seq(-0.4,0.4,0.8/8), colors = rev(BuenColors::jdb_palette(name = 'brewer_celsius')))
Heatmap(t(df), col = col_fun)

###
# survival analysis
###

combined.data$vG.HSC.rel = combined.data$vG.HSC / max(combined.data$vG.HSC)
combined.data$vG.progenitor.rel = combined.data$vG.progenitor / max(combined.data$vG.progenitor)
combined.data$vG.GMP.rel = combined.data$vG.GMP / max(combined.data$vG.GMP)
combined.data$vG.promono.rel = combined.data$vG.promono / max(combined.data$vG.promono)
combined.data$vG.mono.rel = combined.data$vG.mono / max(combined.data$vG.mono)
combined.data$vG.cDC.rel = combined.data$vG.cDC / max(combined.data$vG.cDC)
combined.data$Erythroid.rel = combined.data$Erythroid / max(combined.data$Erythroid)
combined.data$Megakaryocytic.rel = combined.data$Megakaryocytic / max(combined.data$Megakaryocytic)

sample.identity = combined.data[,c('vG.HSC.rel', 'vG.progenitor.rel', 'vG.GMP.rel', 'vG.promono.rel', 
                                   'vG.mono.rel', 'vG.cDC.rel', 'Erythroid.rel', 'Megakaryocytic.rel')] %>% mutate(identity = names(.)[max.col(.)]) 


boo = beataml.meta %>% group_by(dbgap_rnaseq_sample) %>% summarize(dbgap_rnaseq_sample, overallSurvival = unique(overallSurvival), vitalStatus)
boo$status = NA
boo$status[which(boo$vitalStatus == 'Alive')] = 0
boo$status[which(boo$vitalStatus == 'Dead')] = 1
boo$identity = sample.identity[boo$dbgap_rnaseq_sample, 'identity']

library(ggfortify)
library(survival)

km_fit <- survfit(Surv(overallSurvival, status) ~ identity, data=boo)
autoplot(km_fit)