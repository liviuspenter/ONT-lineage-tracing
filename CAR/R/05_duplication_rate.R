# plot duplication rate of CD28 transcripts

library(dplyr)
library(ggplot2)
library(nanoranger.R)
library(Seurat)

technology.colors = c('Both' = '#283350', 'ONT' = '#f93800', 'Illumina' = '#ffb500')

so = readRDS(file='./data/CAR/objects/20220819_CAR_so.rds')

### read length versus duplication 
BC.data = data.table::fread('./data/CAR/97_6_CART_CD28.csv.gz')
BC.list = BC.data %>% group_by(bc, umi) %>% summarize(n = n()) %>% filter (n > 4)

# run starcode and identify UMI clusters
UMIs.collapsed = data.frame()
for (bc in unique(BC.list$bc)) {
  message(bc)
  filehandle.in = tempfile()
  filehandle.out = tempfile()
  write.table(x=BC.list[which(BC.list$bc == bc), c('umi', 'n')], file = filehandle.in, sep = '\t', row.names = F, quote = F, col.names = F)
  
  system(paste0(STARCODE, ' -d 3 -i ', filehandle.in, ' -o ', filehandle.out, ' --print-clusters'))
  starcode.output = as.data.frame(read.csv2(filehandle.out, sep = '\t', header = F))
  colnames(starcode.output) = c('umi', 'n', 'umi.non.collapsed')
  starcode.output$bc = bc
  
  UMIs.collapsed = rbind(UMIs.collapsed, starcode.output)
  
  unlink(filehandle.in)
  unlink(filehandle.out)
}

BC.data = data.table::fread('./data/CAR/97_6_CART_CD28.csv.gz')
#BC.data$bc = paste0(BC.data$bc, '-1')
BC.data.condensed = BC.data %>% group_by(gene, bc, umi) %>% summarize(n = length(gene),
                                                                      qlen.median = median(qlen))
BC.data.condensed$bc.umi = paste0(BC.data.condensed$bc, '.', BC.data.condensed$umi)
UMIs.collapsed$bc.umi = paste0(UMIs.collapsed$bc, '.', UMIs.collapsed$umi)
rownames(UMIs.collapsed) = UMIs.collapsed$bc.umi

BC.data.condensed = BC.data.condensed[which(BC.data.condensed$bc.umi %in% UMIs.collapsed$bc.umi),]
BC.data.condensed$count = UMIs.collapsed[BC.data.condensed$bc.umi, 'n']

p=ggplot() + 
  geom_point(data=BC.data.condensed[which(BC.data.condensed$gene == 'CD28'),], aes(x=qlen.median, y=count), color='darkgreen') + 
  scale_x_continuous('transcript length',limits = c(0,700)) +
  scale_y_log10('reads per transcript') + 
  theme_classic()

q=ggplot() + 
  #geom_point(data=BC.data.condensed[which(BC.data.condensed$gene == 'CARTmod'),], aes(x=qlen.median, y=count), color='firebrick', size=0.5) + 
  ggrastr::rasterize(ggpointdensity::geom_pointdensity(data=BC.data.condensed[which(BC.data.condensed$gene == 'CARTmod'),], 
                                                       aes(x=qlen.median, y=count), size=0.5), dpi=600) +
  viridis::scale_color_viridis() +
  scale_x_continuous('transcript length',limits = c(0,2000)) +
  scale_y_log10('reads per transcript') + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./CAR/figures/plots/20220819_CAR_transcript_length_duplication_rate.svg', width = 2.5, height = 2, plot = q)
cowplot::plot_grid(plotlist = list(p,q)) 

q=ggplot() + 
  #geom_point(data=BC.data.condensed[which(BC.data.condensed$gene == 'CARTmod'),], aes(x=qlen.median, y=count), color='firebrick', size=0.5) + 
  ggrastr::rasterize(ggpointdensity::geom_pointdensity(data=BC.data.condensed[which(BC.data.condensed$gene == 'CARTmod'),], 
                                                       aes(x=-qlen.median, y=count), size=0.5), dpi=600) +
  viridis::scale_color_viridis() +
  scale_x_continuous('transcript end from primer',limits = c(-2000,0), breaks = c(-2000, -1500, -1000,-500,0), labels = c('2000', '1500', '1000', '500', '0')) +
  scale_y_log10('reads per transcript') + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./CAR/figures/plots/20220821_CAR_transcript_length_duplication_rate.svg', width = 2.5, height = 2, plot = q)

q=ggplot() + 
  #geom_point(data=BC.data.condensed[which(BC.data.condensed$gene == 'CARTmod'),], aes(x=qlen.median, y=count), color='firebrick', size=0.5) + 
  ggrastr::rasterize(ggpointdensity::geom_pointdensity(data=BC.data.condensed[which(BC.data.condensed$gene == 'CARTmod'),], 
                                                       aes(x=-qlen.median, y=count), size=0.5), dpi=600) +
  viridis::scale_color_viridis() +
  scale_x_continuous('transcript end from primer',limits = c(-2000,0), breaks = c(-2000, -1500, -1000,-500,0), labels = c('2000', '1500', '1000', '500', '0')) +
  scale_y_log10('reads per transcript') + 
  theme_classic() +
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./CAR/figures/plots/20220821_CAR_transcript_length_duplication_rate_with_legend.svg', width = 2.5, height = 2, plot = q)

#
so = readRDS('./data/CAR/objects/20220819_CAR_so.rds')
CAR.ONT = nanoranger.R::extract_fusion_gene('./data/CAR/97_6_CART_CD28.csv.gz', WILDTYPE = 'CD28', FUSION = 'CARTmod')
CAR.ONT$bc = paste0(CAR.ONT$bc, '-1')
CAR.ONT = CAR.ONT[which(CAR.ONT$bc %in% colnames(so)),]
rownames(CAR.ONT) = CAR.ONT$bc

# look at UMIs
ggplot(CAR.ONT[which(CAR.ONT$CARTmod != 0),], aes(x=CARTmod)) + geom_histogram(binwidth = 1) +
  scale_x_continuous('CAR UMIs per cell', limits = c(0,100)) + 
  scale_y_continuous('cells') + 
  theme_classic() +
  theme(legend.position = 'none', 
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./CAR/figures/plots/20220826_UMI_histogram.svg', width = 2, height = 2)
so.subset = readRDS('./data/CAR/objects/20220819_CAR_so_subset.rds')
so.subset$CAR.UMI = CAR.ONT[colnames(so.subset), 'CARTmod']
p=FeaturePlot(so.subset, 'CAR.UMI', min.cutoff = 'q05', max.cutoff = 'q95', cols = c('lightblue', 'firebrick')) + 
  NoAxes() + NoLegend() + theme(plot.title = element_blank())
ggsave('./CAR/figures/UMAP/20220826_CAR_UMIs.png', width = 4, height = 4, dpi = 600, plot = p)

p=FeaturePlot(so.subset, 'CAR.UMI', min.cutoff = 'q05', max.cutoff = 'q95', cols = c('lightblue', 'firebrick')) + 
  NoAxes() + theme(plot.title = element_blank())
ggsave('./CAR/figures/UMAP/20220826_CAR_UMIs.svg', width = 4, height = 4, dpi = 600, plot = p)

### compare short read data: CAR detected versus CD28 only

# read CAR ONT data
CAR.data.ONT = as.data.frame(data.table::fread('./data/CAR/97_6_CART_CD28.csv.gz') %>% 
                               filter(gene == 'CARTmod') %>% 
                               group_by(bc) %>% 
                               summarize(CARTmod = length(gene)))
CAR.data.ONT$bc = paste0(CAR.data.ONT$bc, '-1')
rownames(CAR.data.ONT) = CAR.data.ONT$bc
CAR.data.ONT = CAR.data.ONT[intersect(rownames(CAR.data.ONT), colnames(so)),]

# read Illumina data
CAR.data.Illumina = as.data.frame(data.table::fread('./data/CAR/97_6_CART_illumina.csv'))
colnames(CAR.data.Illumina)[1] = 'bc'
rownames(CAR.data.Illumina) = CAR.data.Illumina$bc
CAR.data.Illumina = CAR.data.Illumina[intersect(CAR.data.Illumina$bc, colnames(so)),]

# join data
CAR.data.combined = merge(CAR.data.ONT[,c('bc', 'CARTmod')], CAR.data.Illumina[,c('bc', 'CART_count')], by = 'bc', all=T)
CAR.data.combined[is.na(CAR.data.combined)] = 0
CAR.data.combined$CARTmod.plot = CAR.data.combined$CARTmod
CAR.data.combined$CART_count.plot = CAR.data.combined$CART_count
CAR.data.combined$CARTmod.plot[which(CAR.data.combined$CARTmod == 0)] = jitter(rep(0.1, length(which(CAR.data.combined$CARTmod == 0))), amount = 0.01)
CAR.data.combined$CART_count.plot[which(CAR.data.combined$CART_count == 0)] = jitter(rep(0.1, length(which(CAR.data.combined$CART_count == 0))), amount = 0.01)
CAR.data.combined = CAR.data.combined %>% arrange(-CARTmod) %>% mutate(ONT.rank = row_number())
CAR.data.combined = CAR.data.combined %>% arrange(-CART_count) %>% mutate(Illumina.rank = row_number())

CAR.cells = CAR.data.combined$bc[which(CAR.data.combined$CART_count != 0)]
boo = GetAssayData(so, slot = 'counts')['CD28',]
CD28.cells = names(which(boo != 0))
CD28.cells = CD28.cells[-which(CD28.cells %in% CAR.cells)]

so.subset2 = subset(so, cells = which(colnames(so) %in% c(CAR.cells, CD28.cells)))
so.subset2 <- subset(so.subset2, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 20)
so.subset2 = NormalizeData(so.subset2)
so.subset2 = FindVariableFeatures(so.subset2)
so.subset2 = ScaleData(so.subset2)
so.subset2 = RunPCA(so.subset2)
so.subset2 = FindNeighbors(so.subset2)
so.subset2 = FindClusters(so.subset2, resolution = 0.3)
so.subset2 = RunUMAP(so.subset2, dims = 1:20)
so.subset2$CAR = ifelse(colnames(so.subset2) %in% CAR.cells, 'detected', 'not.detected')
saveRDS(file='./data/20221029_CAR_so_subset_illumina.rds', so.subset2)
so.subset2 = readRDS('./data/CAR/objects/20221029_CAR_so_subset_illumina.rds')
markers = FindMarkers(so.subset2, group.by = 'CAR', ident.1 = 'detected', ident.2 = 'not.detected')
markers$gene = rownames(markers)
markers$highlight = ifelse(-log10(markers$p_val_adj) > 4 & abs(markers$avg_log2FC) > 0.5, 'yes', 'no')

ggplot(markers, aes(x=avg_log2FC, y=-log10(p_val_adj), color=highlight)) + 
  geom_point() + 
  ggrepel::geom_label_repel(data=markers[which(markers$highlight == 'yes'),], aes(x=avg_log2FC, y=-log10(p_val_adj), color=highlight, label=gene), 
                            size = 3, label.size = 0) + 
  geom_hline(yintercept = 4) + 
  geom_vline(xintercept = c(0.5, -0.5)) + 
  scale_x_continuous('Log2FC',limits = c(-1.2, 1.2)) + 
  scale_y_continuous('-log10(FDR)') + 
  scale_color_manual(values = c('yes' = 'black', 'no' = 'grey')) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./CAR/figures/plots/20221029_volcano_CAR_CD28_Illumina.svg', width = 2, height = 2)
