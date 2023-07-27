# PTPRC exon 4 expression in TILs versus matched PBMC of melanoma patient

library(dplyr)
library(ggplot2)
library(nanoranger.R)
library(Seurat)

TIL3 = readRDS('./data/TCR/objects/20220820_TCR3.rds')
PBMC8 = readRDS('./data/TCR/objects/20220820_TCR8.rds')

knee_plot(data.table::fread('./data/isoforms/TIL3/TIL3_exons.gz') %>% filter(gene == 'PTPRC'))
knee_plot(data.table::fread('./data/isoforms/PBMC8/PBMC8_exons.gz') %>% filter(gene == 'PTPRC'))

PTPRC.TIL3 = extract_isoforms('./data/isoforms/TIL3/TIL3_exons.gz', GENE = 'PTPRC', EXON = 'exon4', filter = 5)
PTPRC.TIL3 = PTPRC.TIL3[which(PTPRC.TIL3$bc %in% colnames(TIL3)),]

PTPRC.PBMC8 = extract_isoforms('./data/isoforms/PBMC8/PBMC8_exons.gz', GENE = 'PTPRC', EXON = 'exon4', filter = 5)
PTPRC.PBMC8 = PTPRC.PBMC8[which(PTPRC.PBMC8$bc %in% colnames(PBMC8)),]

PTPRC.TIL3$CD45RA = GetAssayData(TIL3, assay = 'ADT')['CD45RA', PTPRC.TIL3$bc]
PTPRC.TIL3$CD45RO = GetAssayData(TIL3, assay = 'ADT')['CD45RO', PTPRC.TIL3$bc]
PTPRC.TIL3$UMAP1 = TIL3@reductions$umap@cell.embeddings[PTPRC.TIL3$bc, 'UMAP_1']
PTPRC.TIL3$UMAP2 = TIL3@reductions$umap@cell.embeddings[PTPRC.TIL3$bc, 'UMAP_2']

PTPRC.PBMC8$CD45RA = GetAssayData(PBMC8, assay = 'ADT')['CD45RA', PTPRC.PBMC8$bc]
PTPRC.PBMC8$CD45RO = GetAssayData(PBMC8, assay = 'ADT')['CD45RO', PTPRC.PBMC8$bc]
PTPRC.PBMC8$UMAP1 = PBMC8@reductions$umap@cell.embeddings[PTPRC.PBMC8$bc, 'UMAP_1']
PTPRC.PBMC8$UMAP2 = PBMC8@reductions$umap@cell.embeddings[PTPRC.PBMC8$bc, 'UMAP_2']

ggplot() +
  geom_point(data=PTPRC.TIL3[order(PTPRC.TIL3$CD45RA),], aes(x=UMAP1, y=UMAP2, color=CD45RA), size=0.5) + 
  scale_color_gradientn(colors = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos'))) + 
  theme_classic() + 
  theme(legend.position = 'none') +
  NoAxes()
ggsave('./isoforms/figures/melanoma/UMAP/20230102_TIL3_CD45RA.png', width = 4, height = 4, dpi = 600)

ggplot() +
  geom_point(data=PTPRC.TIL3[order(PTPRC.TIL3$CD45RA),], aes(x=UMAP1, y=UMAP2, color=CD45RA), size=0.5) + 
  scale_color_gradientn(colors = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos'))) 
ggsave('./isoforms/figures/melanoma/UMAP/20230102_TIL3_CD45RA.svg', width = 4, height = 4, dpi = 600)

ggplot() + 
  geom_point(data=PTPRC.TIL3[order(PTPRC.TIL3$detected/(PTPRC.TIL3$detected+PTPRC.TIL3$not.detected)),], 
             aes(x=UMAP1, y=UMAP2, color=detected/(detected+not.detected)), size=0.5) + 
  scale_color_gradientn(colors = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos'))) + 
  theme_classic() + 
  theme(legend.position = 'none') + 
  NoAxes()
ggsave('./isoforms/figures/melanoma/UMAP/20230102_TIL3_PTPRC_exon4.png', width = 4, height = 4, dpi = 600)


ggplot() + 
  geom_point(data=PTPRC.TIL3[order(PTPRC.TIL3$detected/(PTPRC.TIL3$detected+PTPRC.TIL3$not.detected)),], 
             aes(x=UMAP1, y=UMAP2, color=detected/(detected+not.detected)), size=0.5) + 
  scale_color_gradientn(colors = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos'))) 
ggsave('./isoforms/figures/melanoma/UMAP/20230102_TIL3_PTPRC_exon4.svg', width = 4, height = 4, dpi = 600)

ggplot() +
  geom_point(data=PTPRC.PBMC8[order(PTPRC.PBMC8$CD45RA),], aes(x=UMAP1, y=UMAP2, color=CD45RA), size=0.5) + 
  scale_color_gradientn(colors = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos'))) + 
  theme_classic() + 
  theme(legend.position = 'none') +
  NoAxes()
ggsave('./isoforms/figures/melanoma/UMAP/20230102_PBMC8_CD45RA.png', width = 4, height = 4, dpi = 600)

ggplot() + 
  geom_point(data=PTPRC.PBMC8[order(PTPRC.PBMC8$detected/(PTPRC.PBMC8$detected + PTPRC.PBMC8$not.detected)),], 
             aes(x=UMAP1, y=UMAP2, color=detected/(detected+not.detected)), size=0.5) + 
  scale_color_gradientn(colors = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos'))) + 
  theme_classic() + 
  theme(legend.position = 'none') + 
  NoAxes()
ggsave('./isoforms/figures/melanoma/UMAP/20230102_PBMC8_PTPRC_exon4.png', width = 4, height = 4, dpi = 600)
