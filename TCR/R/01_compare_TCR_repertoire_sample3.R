library(dplyr)
library(ggplot2)
library(Seurat)

# create Seurat object
seurat.data = Read10X(data.dir = paste0('./data/TCR/Pool80_1-40/filtered_feature_bc_matrix/'))
rownames(x = seurat.data[['Antibody Capture']]) = gsub(pattern = '*_CITEseq', replacement = '', rownames(seurat.data[['Antibody Capture']])) 
so <- CreateSeuratObject(counts = seurat.data$`Gene Expression`, project = 'TCR3', min.cells = 3, min.features = 200)
so[["ADT"]] <- CreateAssayObject(seurat.data$`Antibody Capture`[,colnames(x = so)])
so = NormalizeData(so, assay = 'ADT', normalization.method = 'CLR')
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- subset(so, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 20)
so = NormalizeData(so)
so = FindVariableFeatures(so)
so = ScaleData(so)
so = RunPCA(so)
so = RunUMAP(so, dims = 1:30)
so = FindNeighbors(so)
so = FindClusters(so, resolution = 0.1)
so = RenameCells(so, new.names = gsub(colnames(so), pattern = '-1', replacement = '')) 
saveRDS('./data/TCR/objects/20220820_TCR3.rds', object = so)
so = readRDS('./data/TCR/objects/20220820_TCR3.rds')

# read CD8+ TIL data from Oliveira et al., Nature 2021 - needs to be obtained directly from Giacomo Oliveira
Giacomo.so = readRDS('/Users/shaka87/dfci/10026/10x/data/other_datasets/tils.CD8.R0.6.harmonized.20200422.rds')
Giacomo.barcodes = colnames(Giacomo.so)[which(Giacomo.so$sample == 'TIL_pre_CD45pos_CD3pos-1')]
Giacomo.barcodes = stringr::str_split_fixed(Giacomo.barcodes, pattern = '_', n=2)[,2]
Giacomo.barcodes = stringr::str_split_fixed(Giacomo.barcodes, pattern = '-', n=2)[,1]

# get nanoranger TCR data
ONT.data = as.data.frame(data.table::fread('./data/TCR/20220610_TCR3_clone_assignments.csv', header = T))
rownames(ONT.data) = ONT.data$V1
colnames(ONT.data)[1] = 'barcode'

# get Illumina TCR data
Illumina.data = data.table::fread('./data/TCR/Pool80_1-47/all_contig_annotations.csv')
Illumina.data$barcode = gsub(Illumina.data$barcode, pattern = '-1', replacement = '')
#rownames(Illumina.data) = Illumina.data$barcode

# combine data
Illumina.reads = Illumina.data %>% group_by(barcode, chain) %>% summarize(reads = sum(reads))
ONT.reads = rbind(ONT.data %>% group_by(barcode) %>% summarize(chain = 'TRB', reads = sum(as.numeric(strsplit(TRB_counts, split = '_')[[1]]))),
                  ONT.data %>% group_by(barcode) %>% summarize(chain = 'TRA', reads = sum(as.numeric(strsplit(TRA_counts, split = '_')[[1]]))))

length(intersect(Illumina.reads$barcode, ONT.reads$barcode))
length(setdiff(Illumina.reads$barcode, ONT.reads$barcode))

combined.reads = rbind(merge(Illumina.reads[which(Illumina.reads$chain == 'TRA'),], 
                             ONT.reads[which(ONT.reads$chain == 'TRA'),], 
                             by='barcode', all=T),
                       merge(Illumina.reads[which(Illumina.reads$chain == 'TRB'),], 
                             ONT.reads[which(ONT.reads$chain == 'TRB'),], 
                             by='barcode', all=T))
colnames(combined.reads) = c('barcode', 'Illumina.chain', 'Illumina.reads', 'ONT.chain', 'ONT.reads')
combined.reads$Illumina.reads[which(is.na(combined.reads$Illumina.reads))] = 0
combined.reads$ONT.reads[which(is.na(combined.reads$ONT.reads))] = 0
combined.reads = combined.reads[-which(combined.reads$Illumina.reads == 0 & combined.reads$ONT.reads == 0),]

combined.reads$Illumina.reads.log10 = combined.reads$Illumina.reads
combined.reads$Illumina.reads.log10[which(combined.reads$Illumina.reads == 0)] = 
  0.5+abs(jitter(combined.reads$Illumina.reads.log10[which(combined.reads$Illumina.reads == 0)], factor = 5))
combined.reads$ONT.reads.log10 = combined.reads$ONT.reads
combined.reads$ONT.reads.log10[which(combined.reads$ONT.reads == 0)] = 
  0.5+abs(jitter(combined.reads$ONT.reads.log10[which(combined.reads$ONT.reads == 0)], factor = 5))
combined.reads$Illumina.chain[which(is.na(combined.reads$Illumina.chain))] = combined.reads$ONT.chain[which(is.na(combined.reads$Illumina.chain))]
combined.reads$ONT.chain[which(is.na(combined.reads$ONT.chain))] = combined.reads$Illumina.chain[which(is.na(combined.reads$ONT.chain))]

p=ggplot(combined.reads[sample(nrow(combined.reads), nrow(combined.reads), replace = F),], 
         aes(x=Illumina.reads.log10, y=ONT.reads.log10, color=Illumina.chain)) + 
  ggrastr::rasterize(geom_point(size=0.5, alpha=0.3, shape=1), dpi=600) + 
  scale_x_log10('Illumina reads',breaks = c(0.55, 1, 10, 100, 1000, 10000), labels = c('ND', '1', '10', '100', '1000', '10000'), limits = c(0.5, 25000)) + 
  scale_y_log10('ONT reads',breaks = c(0.55, 1, 10, 100, 1000, 10000), labels = c('ND', '1', '10', '100', '1000', '10000'), limits = c(0.5, 25000)) + 
  scale_color_manual(values = c('TRA' = 'green', 'TRB' = 'darkgreen')) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=8, color='black'))
ggsave('./TCR/figures/plots/20220728_Illumina_ONT_reads_sample3.svg', width = 2, height = 2, plot = p)


combined.reads = combined.reads %>% group_by(ONT.chain) %>% arrange(-ONT.reads) %>% mutate(rank.ONT = row_number())
combined.reads = combined.reads %>% group_by(Illumina.chain) %>% arrange(-Illumina.reads) %>% mutate(rank.Illumina = row_number())

# compare library structures
ggplot() + 
  geom_line(data=combined.reads[which(combined.reads$ONT.chain == 'TRA' & combined.reads$ONT.reads != 0),], 
            aes(x=rank.ONT, y=ONT.reads), color='blue') + 
  geom_line(data=combined.reads[which(combined.reads$Illumina.chain == 'TRA' & combined.reads$Illumina.reads != 0),], 
            aes(x=rank.Illumina, y=Illumina.reads), color='blue', linetype = 'dashed') + 
  scale_x_log10('cell rank') + 
  scale_y_log10('reads',limits = c(1,30000)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./TCR/figures/plots/20220820_kneeplot_TRA.svg', width = 2.5, height = 2)

ggplot() + 
  geom_line(data=combined.reads[which(combined.reads$ONT.chain == 'TRB' & combined.reads$ONT.reads != 0),], 
            aes(x=rank.ONT, y=ONT.reads), color='blue') + 
  geom_line(data=combined.reads[which(combined.reads$Illumina.chain == 'TRB' & combined.reads$Illumina.reads != 0),], 
            aes(x=rank.Illumina, y=Illumina.reads), color='orange', linetype = 'dashed') + 
  scale_x_log10('cell rank') + 
  scale_y_log10('reads',limits = c(1,30000)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./TCR/figures/plots/20220820_kneeplot_TRB.svg', width = 2.5, height = 2)

### compare CDR3 sequences between Illumina and nanoranger data

ONT.data$CDR3b = lapply(ONT.data$TRBs, FUN = function(x) {stringr::str_split(x, pattern = '_')[[1]]})
ONT.data$CDR3a = lapply(ONT.data$TRAs, FUN = function(x) {stringr::str_split(x, pattern = '_')[[1]]})

CDR3b.count.ONT = as.data.frame(table(unlist(ONT.data$CDR3b)))
CDR3b.count.ONT$percentage = CDR3b.count.ONT$Freq / sum(CDR3b.count.ONT$Freq)
CDR3b.count.Illumina = as.data.frame(table(Illumina.data$cdr3[which(grepl('TRB',Illumina.data$v_gene))]))
CDR3b.count.Illumina$percentage = CDR3b.count.Illumina$Freq / sum(CDR3b.count.Illumina$Freq)
CDR3b.combined = merge(CDR3b.count.Illumina, CDR3b.count.ONT, by = 'Var1', all=T)

CDR3b.combined$Illumina.freq = CDR3b.combined$percentage.x
CDR3b.combined$Illumina.freq[which(is.na(CDR3b.combined$Illumina.freq))] = 
  abs(jitter(rep(0.00001, length(which(is.na(CDR3b.combined$Illumina.freq)))), factor = 5))

CDR3b.combined$ONT.freq = CDR3b.combined$percentage.y
CDR3b.combined$ONT.freq[which(is.na(CDR3b.combined$ONT.freq))] = 
  abs(jitter(rep(0.00001, length(which(is.na(CDR3b.combined$ONT.freq)))), factor = 5))

ggplot(CDR3b.combined[which(CDR3b.combined$Var1 != 'None'),], aes(x=Illumina.freq, y=ONT.freq)) + 
  geom_point(size=0.5) +
  scale_x_log10('% Illumina',breaks = c(0.00001, 0.001, 0.01, 0.1, 1, 10), labels = c('ND', '0.001','0.01', '0.1', '1', '10'), limits = c(0.000005, 10)) + 
  scale_y_log10('% ONT',breaks = c(0.00001, 0.001, 0.01, 0.1, 1, 10), labels = c('ND', '0.001','0.01', '0.1', '1', '10'), limits = c(0.000005, 10)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./TCR/figures/plots/20220728_Illumina_ONT_CDR3b_sample3.svg', width = 2.5, height = 2.5)

CDR3b.combined$Freq.x[which(is.na(CDR3b.combined$Freq.x))] = 0
CDR3b.combined$Freq.y[which(is.na(CDR3b.combined$Freq.y))] = 0

# outlier:
# CASSHRGDYEQYF not found in Illumina vs. CASSHRGGDYEQYF
# CISVPAELEEETQYF not found in ONT

Illumina.cdr3a = Illumina.data %>% group_by(barcode) %>% filter(chain == 'TRA' & cdr3 != 'None') %>% summarize(cdr3 = paste(sort(cdr3), collapse = '_'))
Illumina.cdr3a$number = unlist(lapply(Illumina.cdr3a$cdr3, FUN = function(x) {length(stringr::str_split(x, pattern = '_')[[1]])}))
ONT.cdr3a = ONT.data %>% group_by(barcode) %>% filter(TRAs != 'None') %>% summarize(cdr3 = paste(sort(stringr::str_split(TRAs, pattern = '_')[[1]]), collapse = '_'))
ONT.cdr3a$number = unlist(lapply(ONT.cdr3a$cdr3, FUN = function(x) {length(stringr::str_split(x, pattern = '_')[[1]])}))
combined.cdr3a = merge(Illumina.cdr3a, ONT.cdr3a, by = 'barcode', all=T)

CDR3a.count.ONT = as.data.frame(table(unlist(ONT.data$CDR3a)))
CDR3a.count.ONT$percentage = CDR3a.count.ONT$Freq / sum(CDR3a.count.ONT$Freq)
CDR3a.count.Illumina = as.data.frame(table(Illumina.data$cdr3[which(grepl('TRA',Illumina.data$v_gene))]))
CDR3a.count.Illumina$percentage = CDR3a.count.Illumina$Freq / sum(CDR3a.count.Illumina$Freq)
CDR3a.combined = merge(CDR3a.count.Illumina, CDR3a.count.ONT, by = 'Var1', all=T)

CDR3a.combined$Illumina.freq = CDR3a.combined$percentage.x
CDR3a.combined$Illumina.freq[which(is.na(CDR3a.combined$Illumina.freq))] = 
  abs(jitter(rep(0.00001, length(which(is.na(CDR3a.combined$Illumina.freq)))), factor = 5))

CDR3a.combined$ONT.freq = CDR3a.combined$percentage.y
CDR3a.combined$ONT.freq[which(is.na(CDR3a.combined$ONT.freq))] = 
  abs(jitter(rep(0.00001, length(which(is.na(CDR3a.combined$ONT.freq)))), factor = 5))

ggplot(CDR3a.combined[which(CDR3a.combined$Var1 != 'None'),], aes(x=Illumina.freq, y=ONT.freq)) + 
  geom_point(size=0.5) +
  scale_x_log10('% Illumina',breaks = c(0.00001, 0.001, 0.01, 0.1, 1, 10), labels = c('ND', '0.001','0.01', '0.1', '1', '10'), limits = c(0.000005, 10)) + 
  scale_y_log10('% ONT',breaks = c(0.00001, 0.001, 0.01, 0.1, 1, 10), labels = c('ND', '0.001','0.01', '0.1', '1', '10'), limits = c(0.000005, 10)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./TCR/figures/plots/20220728_Illumina_ONT_CDR3a_sample3.svg', width = 2.5, height = 2.5)

CDR3a.combined$Freq.x[which(is.na(CDR3a.combined$Freq.x))] = 0
CDR3a.combined$Freq.y[which(is.na(CDR3a.combined$Freq.y))] = 0

# not found in Illumina: CALGERGGGSQGNLIF
# not found in Illumina: LVCEGTMGFQGAQKLVF
# not foundin ONT: CALSESGFGNEKLTF

length(which(combined.cdr3a$cdr3.x == combined.cdr3a$cdr3.y))

Illumina.cdr3b = Illumina.data %>% group_by(barcode) %>% filter(chain == 'TRB' & cdr3 != 'None') %>% summarize(cdr3 = paste(sort(cdr3), collapse = '_'))
Illumina.cdr3b$number = unlist(lapply(Illumina.cdr3b$cdr3, FUN = function(x) {length(stringr::str_split(x, pattern = '_')[[1]])}))
ONT.cdr3b = ONT.data %>% group_by(barcode) %>% filter(TRBs != 'None') %>% summarize(cdr3 = paste(sort(stringr::str_split(TRBs, pattern = '_')[[1]]), collapse = '_'))
ONT.cdr3b$number = unlist(lapply(ONT.cdr3b$cdr3, FUN = function(x) {length(stringr::str_split(x, pattern = '_')[[1]])}))

# compare number of cells with CDR3a/b across nanoranger and Illumina

Illumina.cdr3b = Illumina.data %>% filter(chain == 'TRB' & cdr3 != 'None') %>% group_by(barcode)
Illumina.cdr3b = as.data.frame(sort(table(Illumina.cdr3b$cdr3)))
Illumina.cdr3b$Percentage = Illumina.cdr3b$Freq / length(unique(Illumina.data$barcode))

Illumina.cdr3b$ONT.count = sapply(Illumina.cdr3b$Var1, FUN = function(x){length(which(grepl(x,ONT.data$TRBs)))})
Illumina.cdr3b$ONT.freq = Illumina.cdr3b$ONT.count / length(unique(ONT.data$barcode))
Illumina.cdr3b$ONT.count[which(Illumina.cdr3b$ONT.count == 0)] = 0.1
ggplot(Illumina.cdr3b, aes(x=Freq, y=ONT.count)) +
  ggrastr::rasterize(geom_point(size=0.5, color = 'blue'), dpi=600) +
  scale_x_log10('# Illumina',breaks = c(1, 10, 100, 1000), labels = c('1', '10', '100', '1000'), limits = c(1, 2000)) + 
  scale_y_log10('# ONT',breaks = c(0.1, 1, 10, 100, 1000), labels = c('ND', '1', '10', '100', '1000'), limits = c(0.1, 2000)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./TCR/figures/plots/20220728_Illumina_ONT_CDR3b_count_sample3.svg', width = 2.5, height = 2.5)

Illumina.cdr3a = Illumina.data %>% filter(chain == 'TRA' & cdr3 != 'None') %>% group_by(barcode)
Illumina.cdr3a = as.data.frame(sort(table(Illumina.cdr3a$cdr3)))
Illumina.cdr3a$Percentage = Illumina.cdr3a$Freq / length(unique(Illumina.data$barcode))

Illumina.cdr3a$ONT.count = sapply(Illumina.cdr3a$Var1, FUN = function(x){length(which(grepl(x,ONT.data$TRAs)))})
Illumina.cdr3a$ONT.freq = Illumina.cdr3a$ONT.count / length(unique(ONT.data$barcode))
Illumina.cdr3a$ONT.count[which(Illumina.cdr3a$ONT.count == 0)] = 0.1
ggplot(Illumina.cdr3a, aes(x=Freq, y=ONT.count)) +
  ggrastr::rasterize(geom_point(size=0.5, color='orange'), dpi=600) +
  scale_x_log10('# Illumina',breaks = c(1, 10, 100, 1000), labels = c('1', '10', '100', '1000'), limits = c(1, 2000)) + 
  scale_y_log10('# ONT',breaks = c(0.1, 1, 10, 100, 1000), labels = c('ND', '1', '10', '100', '1000'), limits = c(0.1, 2000)) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'))
ggsave('./TCR/figures/plots/20220728_Illumina_ONT_cdr3a_count_sample3.svg', width = 2.5, height = 2.5)

boo = rbind(Illumina.cdr3a, Illumina.cdr3b)
boo$chain = c(rep('TRA', nrow(Illumina.cdr3a)), rep('TRB', nrow(Illumina.cdr3b)))
boo = boo[sample(nrow(boo), nrow(boo)),]

# combined visualization for CDR3a and CDR3b
ggplot() +
  ggrastr::rasterize(geom_point(data=boo, aes(x=Freq, y=ONT.count, color=chain), size=0.5), dpi=600) +
  scale_x_log10('Illumina cells',breaks = c(1, 10, 100, 1000), labels = c('1', '10', '100', '1000'), limits = c(1, 2000)) + 
  scale_y_log10('ONT cells',breaks = c(0.1, 1, 10, 100, 1000), labels = c('ND', '1', '10', '100', '1000'), limits = c(0.1, 2000)) +
  scale_color_manual(values = c('TRA' = 'green', 'TRB' = 'darkgreen')) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=8, color='black'))
ggsave('./TCR/figures/plots/20220820_Illumina_ONT_cdr3ab_count_sample3.svg', width = 2, height = 2)

### plot statistics
ggplot(data.frame(n = c(length(which(colnames(so) %in% unique(ONT.data$barcode))),
                 length(which(!colnames(so) %in% unique(ONT.data$barcode))),
                 length(which(colnames(so) %in% unique(Illumina.data$barcode))),
                 length(which(!colnames(so) %in% unique(Illumina.data$barcode)))),
           technology = c('ONT', 'ONT', 'Illumina', 'Illumina'),
           status = c('detected', 'not.detected', 'detected', 'not.detected'))) +
  geom_col(aes(x=technology, y=n, fill=status), position = 'dodge') +
  scale_y_continuous('cells') +
  scale_fill_manual(values = c('not.detected' = 'grey', 'detected' = 'black')) + 
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text('Arial', size=10, color='black'),
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank())
ggsave('./TCR/figures/plots/20220820_TCR3_cells_detected.svg', width = 1.3, height = 2)

# UMAPs
p=DimPlot(so, group.by = 'detected', cols = c('yes' = 'black', 'no' = 'grey')) + NoAxes() + NoLegend() + theme(plot.title = element_blank())
ggsave('./TCR/figures/UMAP/20220820_TCR3_detected.png', width = 3, height = 3, dpi = 600)
p=FeaturePlot(so, features = 'CD3E', cols = c('grey', 'darkgreen'), min.cutoff = 'q05', max.cutoff = 'q95', order = T) + NoAxes() + NoLegend() + theme(plot.title = element_blank())
ggsave('./TCR/figures/UMAP/20220820_TCR3_CD3E.png', width = 3, height = 3, dpi = 600)
p=FeaturePlot(so, features = 'CD14', cols = c('grey', 'darkgreen'), min.cutoff = 'q05', max.cutoff = 'q95', order = T) + NoAxes() + NoLegend() + theme(plot.title = element_blank())
ggsave('./TCR/figures/UMAP/20220820_TCR3_CD14.png', width = 3, height = 3, dpi = 600)

so$manual.cluster = 'none'
so$manual.cluster[which(so$seurat_clusters %in% c(0,1,3))] = 'T cell'
so$manual.cluster[which(so$seurat_clusters %in% c(2))] = 'Mono'

p=DimPlot(so, group.by = 'manual.cluster', cols = c('T cell' = 'green', 'Mono' = 'blue')) + NoAxes() + NoLegend() + theme(plot.title = element_blank())
ggsave('./TCR/figures/UMAP/20220820_TCR3_manual_annotation.png', width = 3, height = 3, dpi = 600)
