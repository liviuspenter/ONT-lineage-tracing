# analysis of mtDNA data

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

# analysis of mtDNA coverage
coverage = as.data.frame(data.table::fread('./data/AML/mtDNA//AML1026_2.maegatk/AML1026_2.depthTable.txt'))
coverage$V1 = paste0('AML1026.2_', coverage$V1, '-1')
rownames(coverage) = coverage$V1
coverage$predicted.celltype = 'none'
coverage[intersect(coverage$V1, colnames(AML1026.2)), 'predicted.celltype'] = AML1026.2$predicted.celltype[intersect(coverage$V1, colnames(AML1026.2))]
coverage$predicted.celltype = factor(coverage$predicted.celltype, 
                                     levels = (coverage %>% group_by(predicted.celltype) %>% 
                                                 summarize(coverage = mean(V2)) %>% 
                                                 arrange(desc(coverage)))$predicted.celltype)

mt.genes = c(rownames(AML1026.2)[grepl('^MT-', rownames(AML1026.2))], 'MTRNR2L1', 'MTRNR2L4', 'MTRNR2L6', 'MTRNR2L8', 'MTRNR2L12')
counts = colSums(GetAssayData(AML1026.2, slot = 'counts')[mt.genes, intersect(coverage$V1, colnames(AML1026.2))])

coverage$UMIs = NA
coverage[names(counts), 'UMIs'] = counts

ggplot(coverage[which(coverage$predicted.celltype != 'none'),], aes(x=predicted.celltype, y=V2, fill=predicted.celltype)) + 
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = AML.combined.colors) +
  scale_y_continuous('mtDNA coverage', breaks = c(0,50,100,150), limits = c(0,100)) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank())
ggsave('./AML/figures/mtDNA/20220711_1026_mtDNA_coverage_by_cell.svg', width = 2.5, height = 2.5)

p=ggplot(coverage, aes(x=UMIs, y=V2)) + 
  ggrastr::rasterize(geom_point(aes(color=predicted.celltype), size=0.5), dpi=600) + 
  scale_color_manual(values = AML.combined.colors) +
  scale_x_log10('mtDNA UMIs Illumina') + 
  scale_y_log10('mtDNA coverage ONT', breaks = c(1,10,50,100,200), limits = c(1,250)) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./AML/figures/mtDNA/20220711_1026_mtDNA_coverage_vs_UMIs.svg', width = 2, height = 2, plot = p)
