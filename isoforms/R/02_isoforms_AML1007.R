# visualize differentially expressed PTPRC exons in AML1007

library(dplyr)
library(nanoranger.R)
library(ggplot2)
library(Seurat)

AML1007 = readRDS('./data/AML/objects/20220825_AML1007.rds')

knee_plot(data.table::fread('./data/isoforms/AML1007.1/AML1007_1_exons.gz') %>% filter(gene == 'PTPRC'))
knee_plot(data.table::fread('./data/isoforms/AML1007.3/AML1007_3_exons.gz') %>% filter(gene == 'PTPRC'))
knee_plot(data.table::fread('./data/isoforms/AML1007.5/AML1007_5_exons.gz') %>% filter(gene == 'PTPRC'))

PTPRC.1 = extract_isoforms('./data/isoforms/AML1007.1/AML1007_1_exons.gz', GENE = 'PTPRC', EXON = 'exon4', filter = 10)
PTPRC.1$sample = 'AML1007.1'
PTPRC.3 = extract_isoforms('./data/isoforms/AML1007.3/AML1007_3_exons.gz', GENE = 'PTPRC', EXON = 'exon4', filter = 10)
PTPRC.3$sample = 'AML1007.3'
PTPRC.5 = extract_isoforms('./data/isoforms/AML1007.5/AML1007_5_exons.gz', GENE = 'PTPRC', EXON = 'exon4', filter = 10)
PTPRC.5$sample = 'AML1007.5'

PTPRC = rbind(PTPRC.1, PTPRC.3)
PTPRC = rbind(PTPRC, PTPRC.5)
PTPRC$bc = paste0(PTPRC$sample, '_', PTPRC$bc, '-1')

PTPRC = PTPRC[which(PTPRC$bc %in% colnames(AML1007)),]
PTPRC$predicted.celltype = AML1007$predicted.celltype[PTPRC$bc]
PTPRC$CD45RA = GetAssayData(AML1007, assay = 'ADT')['CD45RA', PTPRC$bc]
PTPRC$CD45RO = GetAssayData(AML1007, assay = 'ADT')['CD45RO', PTPRC$bc]
PTPRC$UMAP1 = AML1007@reductions$umap@cell.embeddings[PTPRC$bc, 'UMAP_1']
PTPRC$UMAP2 = AML1007@reductions$umap@cell.embeddings[PTPRC$bc, 'UMAP_2']

ggplot() + 
  geom_point(data=PTPRC[order(PTPRC$CD45RA),], aes(x=UMAP1, y=UMAP2, color=CD45RA), size=0.5) + 
  scale_color_gradientn(colors = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos'))) + 
  theme_classic() + 
  theme(legend.position = 'none') +
  NoAxes()
ggsave('./isoforms/figures/AML1007/UMAP/20230102_AML1007_CD45RA.png', width = 4, height = 4, dpi = 600)

ggplot() + 
  geom_point(data=PTPRC[order(PTPRC$detected),], aes(x=UMAP1, y=UMAP2, color=detected/(detected+not.detected)), size=0.5) + 
  scale_color_gradientn(colors = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos'))) + 
  theme_classic() + 
  theme(legend.position = 'none') + 
  NoAxes()
ggsave('./isoforms/figures/AML1007/UMAP/20230102_AML1007_PTPRC_exon4.png', width = 4, height = 4, dpi = 600)

PTPRC$detected.rel = PTPRC$detected/(PTPRC$detected+PTPRC$not.detected)

ggplot(PTPRC, aes(x=CD45RA, y=detected.rel)) + geom_point(size=0.5) +
  geom_smooth(method = 'lm') +
  scale_x_continuous(limits = c(0,3)) +
  theme_classic()

cor.test(PTPRC$detected.rel, PTPRC$CD45RA)

PTPRC.1007 = PTPRC

PTPRC.statistics.1007 = PTPRC.1007 %>% group_by(sample, predicted.celltype) %>% 
  summarize(cells.detected = length(which(detected > 0)),
            cells.not.detected = length(which(not.detected > 0)),
            cells = length(predicted.celltype),
            CD45RA = mean(CD45RA))
PTPRC.statistics.1007$detected.rel = PTPRC.statistics.1007$cells.detected / (PTPRC.statistics.1007$cells.detected + PTPRC.statistics.1007$cells.not.detected)
PTPRC.statistics.1007$predicted.celltype.sample = paste0(PTPRC.statistics.1007$predicted.celltype, '.', PTPRC.statistics.1007$sample)

ggplot(PTPRC.statistics.1007[which(PTPRC.statistics.1007$predicted.celltype %in% c('CD4 Naive', 'CD4 Memory', 'CD8 Naive', 'CD8 Memory_2')),], 
       aes(x=predicted.celltype.sample, y=detected.rel, fill=predicted.celltype)) + 
  geom_col(color='black') +
  scale_x_discrete(limits = paste0(rep(c('CD4 Naive', 'CD4 Memory', 'CD8 Naive', 'CD8 Memory_2'), each=3), '.',c('AML1007.1', 'AML1007.3', 'AML1007.5'))) +
  scale_y_continuous(' PTPRC exon 4') + 
  scale_fill_manual(values = AML.combined.colors) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
ggsave('./isoforms/figures/AML1007/plots/20230103_PTPRC_exon4_kinetics.svg', width = 2, height = 2)


ggplot(PTPRC.statistics.1007[which(PTPRC.statistics.1007$predicted.celltype %in% c('CD4 Naive', 'CD4 Memory', 'CD8 Naive', 'CD8 Memory_2')),], 
       aes(x=predicted.celltype.sample, y=CD45RA, fill=predicted.celltype)) + 
  geom_col(color='black') +
  scale_x_discrete(limits = paste0(rep(c('CD4 Naive', 'CD4 Memory', 'CD8 Naive', 'CD8 Memory_2'), each=3), '.',c('AML1007.1', 'AML1007.3', 'AML1007.5'))) +
  scale_y_continuous('CD45RA (CITE-seq)') + 
  scale_fill_manual(values = AML.combined.colors) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
ggsave('./isoforms/figures/AML1007/plots/20230103_CD45RA_kinetics.svg', width = 2, height = 2)


p=DimPlot(AML1007, cells = colnames(AML1007)[which(AML1007$orig.ident == 'AML1007.1')], 
        reduction = 'ref.umap', cells.highlight = list('detected' = PTPRC.1007$bc[which(PTPRC.1007$detected > 0)],
                                                   'not.detected' = PTPRC.1007$bc[which(PTPRC.1007$not.detected > 0)]),
        cols.highlight = c('detected' = 'firebrick', 'not.detected' = 'blue')) + NoLegend()

q=DimPlot(AML1007, cells = colnames(AML1007)[which(AML1007$orig.ident == 'AML1007.3')], 
        reduction = 'ref.umap', cells.highlight = list('detected' = PTPRC.1007$bc[which(PTPRC.1007$detected > 0)],
                                                   'not.detected' = PTPRC.1007$bc[which(PTPRC.1007$not.detected > 0)]),
        cols.highlight = c('detected' = 'firebrick', 'not.detected' = 'blue'))

r=DimPlot(AML1007, cells = colnames(AML1007)[which(AML1007$orig.ident == 'AML1007.5')], 
        reduction = 'ref.umap', cells.highlight = list('detected' = PTPRC.1007$bc[which(PTPRC.1007$detected > 0)],
                                                   'not.detected' = PTPRC.1007$bc[which(PTPRC.1007$not.detected > 0)]),
        cols.highlight = c('detected' = 'firebrick', 'not.detected' = 'blue'))
