# visualize differentially expressed PTPRC exons in AML3005

library(dplyr)
library(ggplot2)
library(nanoranger.R)
library(Seurat)

AML3005 = readRDS('./data/AML/objects/20220718_AML3005_13.rds')

# read isoform data
PTPRC.1 = extract_isoforms('./data/isoforms/AML3005.1/AML3005_1_exons.gz', GENE = 'PTPRC', EXON = 'exon4', filter = 10)
PTPRC.1$sample = 'AML3005.1'
PTPRC.3 = extract_isoforms('./data/isoforms/AML3005.3/AML3005_3_exons.gz', GENE = 'PTPRC', EXON = 'exon4', filter = 5)
PTPRC.3$sample = 'AML3005.3'
PTPRC = rbind(PTPRC.1, PTPRC.3)
PTPRC$bc = paste0(PTPRC$sample, '_', PTPRC$bc, '-1')
PTPRC = PTPRC[which(PTPRC$bc %in% colnames(AML3005)),]
PTPRC$predicted.celltype = AML3005$predicted.celltype[PTPRC$bc]
PTPRC$UMAP1 = AML3005@reductions$umap@cell.embeddings[PTPRC$bc, 'UMAP_1']
PTPRC$UMAP2 = AML3005@reductions$umap@cell.embeddings[PTPRC$bc, 'UMAP_2']
PTPRC$detected.rel = PTPRC$detected/(PTPRC$detected+PTPRC$not.detected)

PTPRC.3005 = PTPRC
PTPRC.statistics.3005 = PTPRC.3005 %>% group_by(sample, predicted.celltype) %>% 
  summarize(cells.detected = length(which(detected > 0)),
            cells.not.detected = length(which(not.detected > 0)),
            cells = length(predicted.celltype))
PTPRC.statistics.3005$detected.rel = PTPRC.statistics.3005$cells.detected / (PTPRC.statistics.3005$cells.detected + PTPRC.statistics.3005$cells.not.detected)
PTPRC.statistics.3005$predicted.celltype.sample = paste0(PTPRC.statistics.3005$predicted.celltype, '.', PTPRC.statistics.3005$sample)

# expression of PTPRC exon 4 before and after ipilimumab
ggplot(PTPRC.statistics.3005[which(PTPRC.statistics.3005$predicted.celltype %in% c('CD4 Naive', 'CD4 Memory', 'CD8 Naive', 'CD8 Memory_2')),], 
       aes(x=predicted.celltype.sample, y=100*detected.rel, fill=predicted.celltype)) + 
  geom_col(color='black') +
  scale_x_discrete(limits = paste0(rep(c('CD4 Naive', 'CD4 Memory', 'CD8 Naive', 'CD8 Memory_2'), each=2), '.',c('AML3005.1', 'AML3005.3'))) +
  scale_y_continuous('% PTPRC exon 4') + 
  scale_fill_manual(values = AML.combined.colors) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
ggsave('./isoforms/figures/AML3005/plots/20230103_PTPRC_exon4_kinetics.svg', width = 2, height = 2)

PTPRC.1 = PTPRC[which(PTPRC$sample == 'AML3005.1'),]
ggplot() + 
  geom_point(data=PTPRC.1[order(PTPRC.1$detected),], aes(x=UMAP1, y=UMAP2, color=detected/(detected+not.detected)), size=0.5) + 
  scale_color_gradientn(colors = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos'))) + 
  theme_classic() + 
  theme(legend.position = 'none') + 
  NoAxes()
ggsave('./isoforms/figures/AML3005/UMAP/20230102_AML3005_1_PTPRC_exon4.png', width = 4, height = 4, dpi = 600)

PTPRC.2 = PTPRC[which(PTPRC$sample == 'AML3005.3'),]
ggplot() + 
  geom_point(data=PTPRC.2[order(PTPRC.2$detected),], aes(x=UMAP1, y=UMAP2, color=detected/(detected+not.detected)), size=0.5) + 
  scale_color_gradientn(colors = c('grey90', BuenColors::jdb_palette(name = 'solar_rojos'))) + 
  theme_classic() + 
  theme(legend.position = 'none') + 
  NoAxes()
ggsave('./isoforms/figures/AML3005/UMAP/20230102_AML3005_3_PTPRC_exon4.png', width = 4, height = 4, dpi = 600)


# compare targeted versus whole-transcriptome coverage of PTPRC and CTLA-4
PTPRC.1 = data.table::fread('./data/isoforms/AML3005.1/AML3005_1_WT_exons.gz')
PTPRC.WT = data.table::fread('./data/isoforms/AML3005.1/AML3005_1_exons.gz')
PTPRC = rbind(PTPRC.1, PTPRC.WT)
PTPRC$gene.exon = paste0(PTPRC$gene, '.', PTPRC$exon)

PTPRC = PTPRC %>% group_by(sample, gene.exon, bc) %>% summarize(n = n())
PTPRC = PTPRC %>% arrange(desc(n)) %>% group_by(sample, gene.exon) %>% mutate(rank = row_number())

ggplot(PTPRC[which(PTPRC$gene.exon == 'PTPRC.exon4'),], aes(x=rank, y=n)) + 
  geom_line(aes(group=sample)) + 
  scale_x_log10('cell rank') +
  scale_y_log10('reads PTPRC') +
  theme_classic() +
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./isoforms/figures/AML3005/plots/20230104_AML3003_1_PTPRC_comparison.svg', width = 2, height = 1.5)

ggplot(PTPRC[which(PTPRC$gene.exon == 'CTLA4.exon2'),], aes(x=rank, y=n)) + 
  geom_line(aes(group=sample)) + 
  scale_x_log10('cell rank') +
  scale_y_log10('reads CTLA-4') +
  theme_classic() +
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))  
ggsave('./isoforms/figures/AML3005/plots/20230104_AML3003_1_CTLA4_comparison.svg', width = 2, height = 1.5)