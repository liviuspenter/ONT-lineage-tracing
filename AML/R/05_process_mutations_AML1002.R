# process mutations for AML1002.13

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)
library(nanoranger.R)
library(Seurat)
library(SummarizedExperiment)

AML1002 = readRDS('./data/AML/objects/20220911_AML1002.rds')
Idents(AML1002) = 'predicted.celltype'

boo=as.data.frame(prop.table(table(Idents(AML1002), AML1002$orig.ident), margin = 2))
boo$Var1 = factor(boo$Var1, levels = rev(names(AML.combined.colors)))
ggplot(boo, aes(x=Var2, y=100*Freq, fill=Var1)) + 
  geom_col() +
  geom_hline(yintercept = 5) + 
  scale_x_discrete(labels = c('Screening', 'SD')) + 
  scale_y_continuous('%cells', limits = c(0,100)) +
  scale_fill_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=8, color='black'),
        axis.title = element_text('Arial', size=8, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML1002/plots/20220912_celltype_kinetics.svg', width = 1.1, height = 2)

# IDH1 mutations - very few mutated cells 
IDH1_1.vaf.1 = extract_mutation(BC.data.file = './data/AML/pileup/AML1002_1_pileup_IDH1_1.csv.gz', REF = 'C', ALT = 'T', FILTER = 10)
IDH1_1.vaf.1$bc = paste0('AML1002.1_', IDH1_1.vaf.1$bc, '-1')
IDH1_1.vaf.1$sample = 'AML1002.1'
IDH1_2.vaf.1 = extract_mutation(BC.data.file = './data/AML/pileup/AML1002_1_pileup_IDH1_2.csv.gz', REF = 'G', ALT = 'A', FILTER = 10)
IDH1_2.vaf.1$bc = paste0('AML1002.1_', IDH1_2.vaf.1$bc, '-1')
IDH1_2.vaf.1$sample = 'AML1002.1'

IDH1_1.vaf.3 = extract_mutation(BC.data.file = './data/AML/pileup/AML1002_3_pileup_IDH1_1.csv.gz', REF = 'C', ALT = 'T', FILTER = 10)
IDH1_1.vaf.3$bc = paste0('AML1002.3_', IDH1_1.vaf.3$bc, '-1')
IDH1_1.vaf.3$sample = 'AML1002.3'
IDH1_2.vaf.3 = extract_mutation(BC.data.file = './data/AML/pileup/AML1002_3_pileup_IDH1_2.csv.gz', REF = 'G', ALT = 'A', FILTER = 10)
IDH1_2.vaf.3$bc = paste0('AML1002.3_', IDH1_2.vaf.3$bc, '-1')
IDH1_2.vaf.3$sample = 'AML1002.3'

IDH1.vaf = as.data.frame(rbind(IDH1_2.vaf.1, IDH1_2.vaf.3))
#IDH1.vaf$predicted.celltype = AML1002$predicted.celltype[IDH1.vaf$bc]
rownames(IDH1.vaf) = IDH1.vaf$bc
write.table(file = './data/AML/mutations/AML1002_IDH1.csv', x = IDH1.vaf, quote = F, sep = '\t')


p=DimPlot(subset(AML1002, orig.ident == 'AML1002.1'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = IDH1_1.vaf.1$bc[which(IDH1_1.vaf.1$mutated == 'mutated')], 
                                 'wildtype' = IDH1_1.vaf.1$bc[which(IDH1_1.vaf.1$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
q=DimPlot(subset(AML1002, orig.ident == 'AML1002.3'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = IDH1_1.vaf.3$bc[which(IDH1_1.vaf.3$mutated == 'mutated')], 
                                 'wildtype' = IDH1_1.vaf.3$bc[which(IDH1_1.vaf.3$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
cowplot::plot_grid(plotlist = list(p,q))
ggsave('./AML/figures/AML1002/UMAP/AML1002_1_IDH1_1.png', width = 3, height = 3, dpi = 600, plot = p)
ggsave('./AML/figures/AML1002/UMAP/AML1002_3_IDH1_1.png', width = 3, height = 3, dpi = 600, plot = q)


p=DimPlot(subset(AML1002, orig.ident == 'AML1002.1'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = IDH1_2.vaf.1$bc[which(IDH1_2.vaf.1$mutated == 'mutated')], 
                                 'wildtype' = IDH1_2.vaf.1$bc[which(IDH1_2.vaf.1$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
q=DimPlot(subset(AML1002, orig.ident == 'AML1002.3'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = IDH1_2.vaf.3$bc[which(IDH1_2.vaf.3$mutated == 'mutated')], 
                                 'wildtype' = IDH1_2.vaf.3$bc[which(IDH1_2.vaf.3$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
cowplot::plot_grid(plotlist = list(p,q))
ggsave('./AML/figures/AML1002/UMAP/AML1002_1_IDH1_2.png', width = 3, height = 3, dpi = 600, plot = p)
ggsave('./AML/figures/AML1002/UMAP/AML1002_3_IDH1_2.png', width = 3, height = 3, dpi = 600, plot = q)


# IDH2 mutation
IDH2.vaf.1 = extract_mutation(BC.data.file = './data/AML/pileup/AML1002_1_pileup_IDH2.csv.gz', REF = 'C', ALT = 'T', FILTER = 10)
IDH2.vaf.1$bc = paste0('AML1002.1_', IDH2.vaf.1$bc, '-1')
IDH2.vaf.1$sample = 'AML1002.1'
IDH2.vaf.3 = extract_mutation(BC.data.file = './data/AML/pileup/AML1002_3_pileup_IDH2.csv.gz', REF = 'C', ALT = 'T', FILTER = 10)
IDH2.vaf.3$bc = paste0('AML1002.3_', IDH2.vaf.3$bc, '-1')
IDH2.vaf.3$sample = 'AML1002.3'
IDH2.vaf = as.data.frame(rbind(IDH2.vaf.1, IDH2.vaf.3))
rownames(IDH2.vaf) = IDH2.vaf$bc
write.table(file = './data/AML/mutations/AML1002_IDH2.csv', x = IDH2.vaf, quote = F, sep = '\t')

p=DimPlot(subset(AML1002, orig.ident == 'AML1002.1'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = IDH2.vaf$bc[which(IDH2.vaf$mutated == 'mutated')], 
                                 'wildtype' = IDH2.vaf$bc[which(IDH2.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
q=DimPlot(subset(AML1002, orig.ident == 'AML1002.3'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = IDH2.vaf$bc[which(IDH2.vaf$mutated == 'mutated')], 
                                 'wildtype' = IDH2.vaf$bc[which(IDH2.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
cowplot::plot_grid(plotlist = list(p,q))
ggsave('./AML/figures/AML1002/UMAP/AML1002_1_IDH2.png', width = 3, height = 3, dpi = 600, plot = p)
ggsave('./AML/figures/AML1002/UMAP/AML1002_3_IDH2.png', width = 3, height = 3, dpi = 600, plot = q)

### TET2
TET2.vaf.1 = extract_mutation(BC.data.file = './data/AML/pileup/AML1002_1_pileup_TET2.csv.gz', REF = 'T', ALT = 'C', FILTER = 10)
TET2.vaf.1$bc = paste0('AML1002.1_', TET2.vaf.1$bc, '-1')
TET2.vaf.1$sample = 'AML1002.1'
TET2.vaf.3 = extract_mutation(BC.data.file = './data/AML/pileup/AML1002_3_pileup_TET2.csv.gz', REF = 'T', ALT = 'C', FILTER = 10)
TET2.vaf.3$bc = paste0('AML1002.3_', TET2.vaf.3$bc, '-1')
TET2.vaf.3$sample = 'AML1002.3'
TET2.vaf = as.data.frame(rbind(TET2.vaf.1, TET2.vaf.3))
rownames(TET2.vaf) = TET2.vaf$bc
write.table(file = './data/AML/mutations/AML1002_TET2.csv', x = TET2.vaf, quote = F, sep = '\t')
TET2.vaf$predicted.celltype = AML1002$predicted.celltype[TET2.vaf$bc]

p=DimPlot(subset(AML1002, orig.ident == 'AML1002.1'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = TET2.vaf$bc[which(TET2.vaf$mutated == 'mutated')], 
                                 'wildtype' = TET2.vaf$bc[which(TET2.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
q=DimPlot(subset(AML1002, orig.ident == 'AML1002.3'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = TET2.vaf$bc[which(TET2.vaf$mutated == 'mutated')], 
                                 'wildtype' = TET2.vaf$bc[which(TET2.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
cowplot::plot_grid(plotlist = list(p,q))
ggsave('./AML/figures/AML1002/UMAP/AML1002_1_TET2.png', width = 3, height = 3, dpi = 600, plot = p)
ggsave('./AML/figures/AML1002/UMAP/AML1002_3_TET2.png', width = 3, height = 3, dpi = 600, plot = q)


### TP53_1
TP53_1.vaf.1 = extract_mutation(BC.data.file = './data/AML/pileup/AML1002_1_pileup_TP53_1.csv.gz', REF = 'T', ALT = 'C', FILTER = 10)
TP53_1.vaf.1$bc = paste0('AML1002.1_', TP53_1.vaf.1$bc, '-1')
TP53_1.vaf.1$sample = 'AML1002.1'
TP53_1.vaf.3 = extract_mutation(BC.data.file = './data/AML/pileup/AML1002_3_pileup_TP53_1.csv.gz', REF = 'T', ALT = 'C', FILTER = 10)
TP53_1.vaf.3$bc = paste0('AML1002.3_', TP53_1.vaf.3$bc, '-1')
TP53_1.vaf.3$sample = 'AML1002.3'
TP53_1.vaf = as.data.frame(rbind(TP53_1.vaf.1, TP53_1.vaf.3))
rownames(TP53_1.vaf) = TP53_1.vaf$bc
write.table(file = './data/AML/mutations/AML1002_TP53.1.csv', x = TP53_1.vaf, quote = F, sep = '\t')

p=DimPlot(subset(AML1002, orig.ident == 'AML1002.1'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = TP53_1.vaf$bc[which(TP53_1.vaf$mutated == 'mutated')], 
                                 'wildtype' = TP53_1.vaf$bc[which(TP53_1.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
q=DimPlot(subset(AML1002, orig.ident == 'AML1002.3'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = TP53_1.vaf$bc[which(TP53_1.vaf$mutated == 'mutated')], 
                                 'wildtype' = TP53_1.vaf$bc[which(TP53_1.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
cowplot::plot_grid(plotlist = list(p,q))
ggsave('./AML/figures/AML1002/UMAP/AML1002_1_TP53_1.png', width = 3, height = 3, dpi = 600, plot = p)
ggsave('./AML/figures/AML1002/UMAP/AML1002_3_TP53_1.png', width = 3, height = 3, dpi = 600, plot = q)


### TP53_2
TP53_2.vaf.1 = extract_mutation(BC.data.file = './data/AML/pileup/AML1002_1_pileup_TP53_2.csv.gz', REF = 'G', ALT = 'A', FILTER = 10)
TP53_2.vaf.1$bc = paste0('AML1002.1_', TP53_2.vaf.1$bc, '-1')
TP53_2.vaf.1$sample = 'AML1002.1'
TP53_2.vaf.3 = extract_mutation(BC.data.file = './data/AML/pileup/AML1002_3_pileup_TP53_2.csv.gz', REF = 'G', ALT = 'A', FILTER = 10)
TP53_2.vaf.3$bc = paste0('AML1002.3_', TP53_2.vaf.3$bc, '-1')
TP53_2.vaf.3$sample = 'AML1002.3'
TP53_2.vaf = as.data.frame(rbind(TP53_2.vaf.1, TP53_2.vaf.3))
rownames(TP53_2.vaf) = TP53_2.vaf$bc
write.table(file = './data/AML/mutations/AML1002_TP53.2.csv', x = TP53_2.vaf, quote = F, sep = '\t')

p=DimPlot(subset(AML1002, orig.ident == 'AML1002.1'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = TP53_2.vaf$bc[which(TP53_2.vaf$mutated == 'mutated')], 
                                 'wildtype' = TP53_2.vaf$bc[which(TP53_2.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
q=DimPlot(subset(AML1002, orig.ident == 'AML1002.3'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = TP53_2.vaf$bc[which(TP53_2.vaf$mutated == 'mutated')], 
                                 'wildtype' = TP53_2.vaf$bc[which(TP53_2.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
cowplot::plot_grid(plotlist = list(p,q))
ggsave('./AML/figures/AML1002/UMAP/AML1002_1_TP53_2.png', width = 3, height = 3, dpi = 600, plot = p)
ggsave('./AML/figures/AML1002/UMAP/AML1002_3_TP53_2.png', width = 3, height = 3, dpi = 600, plot = q)


### U2AF1
U2AF1.vaf.1 = extract_mutation(BC.data.file = './data/AML/pileup/AML1002_1_pileup_U2AF1.csv.gz', REF = 'G', ALT = 'T', FILTER = 10, downsample = 500000)
U2AF1.vaf.1$bc = paste0('AML1002.1_', U2AF1.vaf.1$bc, '-1')
U2AF1.vaf.1$sample = 'AML1002.1'
U2AF1.vaf.3 = extract_mutation(BC.data.file = './data/AML/pileup/AML1002_3_pileup_U2AF1.csv.gz', REF = 'G', ALT = 'T', FILTER = 10, downsample = 500000)
U2AF1.vaf.3$bc = paste0('AML1002.3_', U2AF1.vaf.3$bc, '-1')
U2AF1.vaf.3$sample = 'AML1002.3'
U2AF1.vaf = as.data.frame(rbind(U2AF1.vaf.1, U2AF1.vaf.3))
rownames(U2AF1.vaf) = U2AF1.vaf$bc
write.table(file = './data/AML/mutations/AML1002_U2AF1.csv', x = U2AF1.vaf, quote = F, sep = '\t')

p=DimPlot(subset(AML1002, orig.ident == 'AML1002.1'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = U2AF1.vaf$bc[which(U2AF1.vaf$mutated == 'mutated')], 
                                 'wildtype' = U2AF1.vaf$bc[which(U2AF1.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
q=DimPlot(subset(AML1002, orig.ident == 'AML1002.3'), reduction = 'ref.umap', 
          cells.highlight = list('mutated' = U2AF1.vaf$bc[which(U2AF1.vaf$mutated == 'mutated')], 
                                 'wildtype' = U2AF1.vaf$bc[which(U2AF1.vaf$mutated == 'wildtype')]),
          cols.highlight = c('wildtype' = 'black', 'mutated' = 'firebrick'), sizes.highlight = 1) +
  NoLegend() + 
  NoAxes() + 
  theme(plot.title = element_blank())
cowplot::plot_grid(plotlist = list(p,q))
ggsave('./AML/figures/AML1002/UMAP/AML1002_1_U2AF1.png', width = 3, height = 3, dpi = 600, plot = p)
ggsave('./AML/figures/AML1002/UMAP/AML1002_3_U2AF1.png', width = 3, height = 3, dpi = 600, plot = q)


