# analysis of BCR-ABL1 in primary ALL samples

library(ggplot2)
library(Seurat)

# bone marrow reference
bm = readRDS(file = './data/reference/bm.reference.rds')
bm[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(object = bm[["spca.annoy.neighbors"]], 
                                               file = "./data/reference/reftmp.idx")
Misc(bm[["wnn.umap"]], slot="model")$num_precomputed_nns <- 1

# create Seurat object
seurat.data = Read10X(data.dir = paste0('./data/ALL/Pool115-2_18/filtered_feature_bc_matrix/'))
so <- CreateSeuratObject(counts = seurat.data, project = 'ALL.1', min.cells = 3, min.features = 200)
so[['percent.mt']] <- PercentageFeatureSet(so, pattern = '^MT-')
so = RenameCells(so, add.cell.id = 'ALL.1')

seurat.data = Read10X(data.dir = paste0('./data/ALL/Pool115-2_19/filtered_feature_bc_matrix/'))
so2 <- CreateSeuratObject(counts = seurat.data, project = 'ALL.2', min.cells = 3, min.features = 200)
so2[['percent.mt']] <- PercentageFeatureSet(so2, pattern = '^MT-')
so2 = RenameCells(so2, add.cell.id = 'ALL.2')

seurat.data = Read10X(data.dir = paste0('./data/ALL/Pool116_1/filtered_feature_bc_matrix/'))
so3 <- CreateSeuratObject(counts = seurat.data, project = 'ALL.3', min.cells = 3, min.features = 200)
so3[['percent.mt']] <- PercentageFeatureSet(so3, pattern = '^MT-')
so3 = RenameCells(so3, add.cell.id = 'ALL.3')

seurat.data = Read10X(data.dir = paste0('./data/ALL/Pool116_2/filtered_feature_bc_matrix/'))
so4 <- CreateSeuratObject(counts = seurat.data, project = 'ALL.4', min.cells = 3, min.features = 200)
so4[['percent.mt']] <- PercentageFeatureSet(so4, pattern = '^MT-')
so4 = RenameCells(so4, add.cell.id = 'ALL.4')

so = merge(so, list(so2, so3, so4))
so = NormalizeData(so)
so = FindVariableFeatures(so)
so = ScaleData(so)
so = RunPCA(so)
so = RunUMAP(so, dims = 1:30)
so = FindNeighbors(so)
so = FindClusters(so, resolution = 1)

anchors <- FindTransferAnchors(reference = bm, query = so, k.filter = NA, reference.reduction = "spca", 
                               reference.neighbors = "spca.annoy.neighbors", dims = 1:50)
ALL <- MapQuery(anchorset = anchors, so, reference = bm, 
                      refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
                      reference.reduction = "spca", reduction.model = "wnn.umap")

ALL$manual.cluster = 'ALL'
ALL$manual.cluster[which(ALL$seurat_clusters %in% c(2,18))] = 'NK'
ALL$manual.cluster[which(ALL$seurat_clusters %in% c(5,15))] = 'CD8'
ALL$manual.cluster[which(ALL$seurat_clusters %in% c(0,10))] = 'CD4'
ALL$manual.cluster[which(ALL$seurat_clusters %in% c(20,12,7,27,22,11,23))] = 'Ery'
ALL$manual.cluster[which(ALL$seurat_clusters %in% c(14,21))] = 'B cell'
ALL$manual.cluster[which(ALL$seurat_clusters %in% c(16))] = 'Progenitor'
ALL$manual.cluster[which(ALL$seurat_clusters %in% c(9))] = 'Mono'
saveRDS('./data/ALL/objects/20220913_ALL.rds', object = ALL)

ALL = readRDS('./data/ALL/objects/20220913_ALL.rds')

markers = FindMarkers(ALL, group.by = 'manual.cluster', ident.1 = 'B cell', ident.2 = 'ALL', logfc.threshold = 0.5)
markers$gene = rownames(markers)
markers$FDR = -log10(markers$p_val_adj)
markers$FDR[which(markers$p_val_adj == 0)] = max(markers$FDR[which(markers$FDR != Inf)])
markers$significant = ifelse(abs(markers$avg_log2FC) > 2.5 & markers$FDR > 250, 'yes', 'no')
highlight.genes = c('SOCS2', 'DNTT', 'CD34', 'H2AFY', 'MS4A1', 'FCER2')
ggplot() + 
  ggrastr::rasterize(geom_point(data=markers[which(!markers$gene %in% highlight.genes),], aes(x=avg_log2FC, y=FDR), color = 'grey', size=0.5), dpi=600) + 
  ggrastr::rasterize(geom_point(data=markers[which(markers$gene %in% highlight.genes),], aes(x=avg_log2FC, y=FDR), color = 'magenta', size=0.5), dpi=600) + 
  ggrepel::geom_label_repel(data=markers[which(markers$gene %in% highlight.genes),], aes(x=avg_log2FC, y=FDR, label=gene), size=3, label.size = 0) + 
  scale_x_continuous('log2(FC)') + 
  scale_y_continuous('-log10(FDR)',limits = c(0,350)) + 
  theme_classic() + 
  theme(axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./ALL/figures/plots/20220913_volcano_ALL_Bcells.svg', width = 3, height = 2)


# read BCR-ABL1 data
ALL1.genotype = nanoranger.R::extract_fusion_gene('./data/ALL/ALL1/ALL1_BCRABL1.csv.gz', WILDTYPE = 'ABL1', FUSION = 'BCR', filter = 10)
ALL1.genotype$bc = paste0('ALL.1_', ALL1.genotype$bc, '-1')
rownames(ALL1.genotype) = ALL1.genotype$bc

ALL2.genotype = nanoranger.R::extract_fusion_gene('./data/ALL/ALL2/ALL2_BCRABL1.csv.gz', WILDTYPE = 'ABL1', FUSION = 'BCR', filter = 10)
ALL2.genotype$bc = paste0('ALL.2_', ALL2.genotype$bc, '-1')
rownames(ALL2.genotype) = ALL2.genotype$bc

ALL.cells = c(rownames(ALL1.genotype)[which(ALL1.genotype$mutated == 'mutated')], 
              rownames(ALL2.genotype)[which(ALL2.genotype$mutated == 'mutated')])

df = as.data.frame(ALL@reductions$umap@cell.embeddings)
df$manual.cluster = ALL$manual.cluster[rownames(df)]
df$sample = ALL$orig.ident[rownames(df)]
ggplot(df, aes(x=UMAP_1, y=UMAP_2)) + 
  geom_point(aes(color=manual.cluster), size=0.5) + 
  #geom_point(data=df[ALL.cells,], color='black', size=2) + 
  scale_color_manual(values = c('Ery' = 'grey', 'Progenitor' = 'yellow', 'ALL' = 'darkorange', 'B cell' = 'firebrick', 'CD4' = 'lightblue', 'CD8' = 'darkblue', 
                                'NK' = 'purple', 'Mono' = 'darkgreen')) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
ggsave('./ALL/figures/UMAP/20220913_ALL.png', dpi = 600, width = 4, height = 4)

ggplot(df, aes(x=UMAP_1, y=UMAP_2)) + 
  geom_point(aes(color=sample), size=0.5) + 
  scale_color_manual(values = c('ALL.1' = 'black', 'ALL.2' = 'firebrick', 'ALL.3' = 'lightblue', 'ALL.4' = 'lightgreen')) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
ggsave('./ALL/figures/UMAP/20221026_ALL_samples.png', dpi = 600, width = 4, height = 4)

ggplot(df, aes(x=UMAP_1, y=UMAP_2)) + 
  geom_point(color='grey', size=0.5) + 
  geom_point(data=df[rownames(ALL1.genotype)[which(ALL1.genotype$mutated == 'mutated')],], color='black', size=2) + 
  geom_point(data=df[rownames(ALL2.genotype)[which(ALL2.genotype$mutated == 'mutated')],], color='firebrick', size=2) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
ggsave('./ALL/figures/UMAP/20220913_ALL_mutation.png', dpi = 600, width = 4, height = 4)

DimPlot(ALL, reduction = 'ref.umap', cells.highlight = list('ALL.1' = rownames(ALL1.genotype)[which(ALL1.genotype$mutated == 'mutated')],
                                                            'ALL.2' = rownames(ALL2.genotype)[which(ALL2.genotype$mutated == 'mutated')]),
        cols.highlight = c('ALL.1' = 'black', 'ALL.2' = 'firebrick'), sizes.highlight = 1) + 
  NoLegend() + 
  NoAxes()
ggsave('./ALL/figures/UMAP/20221026_ALL_mutation_ref_UMAP.png', width = 4, height = 4, dpi = 600)

p = FeaturePlot(ALL, 'DNTT', cols = c('grey', BuenColors::jdb_palette(name = 'brewer_red'))) + NoLegend() + NoAxes() +
  theme(plot.title = element_blank())
ggsave('./ALL/figures/UMAP/20221006_DNTT.png', width = 4, height = 4, dpi = 600, plot = p)
p = FeaturePlot(ALL, 'MS4A1', cols = c('grey', BuenColors::jdb_palette(name = 'brewer_red'))) + NoLegend() + NoAxes() +
  theme(plot.title = element_blank())
ggsave('./ALL/figures/UMAP/20221006_MS4A1.png', width = 4, height = 4, dpi = 600, plot = p)
