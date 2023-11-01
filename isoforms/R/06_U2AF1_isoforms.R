# compare U2AF1 isoform expression in AML3003

library(ComplexHeatmap)
library(ggplot2)
library(nanoranger.R)
library(Seurat)
library(SummarizedExperiment)

source("./R/celltypes.R")

AML3003.14 <- readRDS("./data/AML/objects/202208_AML3003_14.rds")

transcript.count.1 <- as.data.frame(data.table::fread("./data/isoforms/AML3003.1/transcript_count.csv.gz"))
rownames(transcript.count.1) <- transcript.count.1$transcript_id
colnames(transcript.count.1) <- paste0("AML3003.1_", colnames(transcript.count.1), "-1")
transcript.count.1 <- transcript.count.1[, intersect(colnames(AML3003.14), colnames(transcript.count.1))]

transcript.count.4 <- as.data.frame(data.table::fread("./data/isoforms/AML3003.4/transcript_count.csv.gz"))
rownames(transcript.count.4) <- transcript.count.4$transcript_id
colnames(transcript.count.4) <- paste0("AML3003.4_", colnames(transcript.count.4), "-1")
transcript.count.4 <- transcript.count.4[, intersect(colnames(AML3003.14), colnames(transcript.count.4))]

transcript.count <-
  cbind(
    transcript.count.1[intersect(rownames(transcript.count.1), rownames(transcript.count.4)), ],
    transcript.count.4[intersect(rownames(transcript.count.1), rownames(transcript.count.4)), ]
  )

AML3003.14 <- subset(AML3003.14, cells = colnames(transcript.count))
AML3003.14[["isoforms"]] <- CreateAssayObject(transcript.count)
DefaultAssay(AML3003.14) <- "isoforms"

AML3003.14 <- NormalizeData(AML3003.14)
AML3003.14 <- ScaleData(AML3003.14)

# U2AF1 isoforms
# gene ID: ENSG00000160201.12
# ENST00000291552.9 - variant a
# ENST00000380276.6 - variant b
# ENST00000464750.5 - variant c

boo <- as.data.frame(t(transcript.count))
# boo = boo[-which(rownames(boo) == 'gene_id'),]
# boo[boo == 0] = 0.1
boo[which(as.numeric(boo$ENST00000291552.9) == 0), "ENST00000291552.9"] <- 0.1
boo[which(as.numeric(boo$ENST00000380276.6) == 0), "ENST00000380276.6"] <- 0.1
p <- ggplot(
  boo, # [which(as.numeric(boo$ENST00000291552.9) != 0 | as.numeric(boo$ENST00000380276.6) != 0),],
  aes(x = as.numeric(`ENST00000291552.9`), y = as.numeric(`ENST00000380276.6`))
) +
  ggrastr::rasterize(geom_point(size = 0.5), dpi = 600) +
  scale_x_sqrt("U2AF1 variant a", limits = c(0, 200), breaks = c(0, 10, 50, 100, 200)) +
  scale_y_sqrt("U2AF1 variant b", limits = c(0, 200), breaks = c(0, 10, 50, 100, 200)) +
  geom_smooth(method = "loess", color = "firebrick") +
  theme_classic() +
  theme(
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./isoforms/figures/AML3003/plots/20230104_AML3003_14_U2AF1_variant_ab.svg", width = 2, height = 2, plot = p)

p <- FeaturePlot(AML3003.14, reduction = "ref.umap", features = "ENST00000291552.9", order = T, raster = T, raster.dpi = c(600, 600)) +
  scale_color_gradientn(colours = c("grey90", BuenColors::jdb_palette(name = "solar_rojos")))
ggsave("./isoforms/figures/AML3003/UMAP/20230104_AML3003_14_U2AF1_variant_a.svg", width = 4.5, height = 4, dpi = 600, plot = p)
ggsave("./isoforms/figures/AML3003/UMAP/20230104_AML3003_14_U2AF1_variant_a.png",
  width = 4, height = 4, dpi = 600,
  plot = p + NoAxes() + NoLegend() + theme(plot.title = element_blank())
)
q <- FeaturePlot(AML3003.14, reduction = "ref.umap", features = "ENST00000380276.6", order = T, raster = T, raster.dpi = c(600, 600)) +
  scale_color_gradientn(colours = c("grey90", BuenColors::jdb_palette(name = "solar_rojos")))
ggsave("./isoforms/figures/AML3003/UMAP/20230104_AML3003_14_U2AF1_variant_b.svg", width = 4.5, height = 4, dpi = 600, plot = q)
ggsave("./isoforms/figures/AML3003/UMAP/20230104_AML3003_14_U2AF1_variant_b.png",
  width = 4, height = 4, dpi = 600,
  plot = q + NoAxes() + NoLegend() + theme(plot.title = element_blank())
)
r <- FeaturePlot(AML3003.14, reduction = "ref.umap", features = "ENST00000464750.5", order = T, raster = T, raster.dpi = c(600, 600)) +
  scale_color_gradientn(colours = c("grey90", BuenColors::jdb_palette(name = "solar_rojos")))
ggsave("./isoforms/figures/AML3003/UMAP/20230104_AML3003_14_U2AF1_variant_c.svg", width = 4.5, height = 4, dpi = 600, plot = r)
ggsave("./isoforms/figures/AML3003/UMAP/20230104_AML3003_14_U2AF1_variant_c.png",
  width = 4, height = 4, dpi = 600,
  plot = r + NoAxes() + NoLegend() + theme(plot.title = element_blank())
)

DefaultAssay(AML3003.14) <- "isoforms"

boo <- as.data.frame(t(as.matrix(GetAssayData(AML3003.14, assay = "isoforms")[c("ENST00000291552.9", "ENST00000380276.6", "ENST00000464750.5"), ])))
boo$predicted.celltype <- AML3003.14$predicted.celltype[rownames(boo)]
# boo$sample = AML3003.14$orig.ident
boo <- boo %>%
  group_by(predicted.celltype) %>%
  summarize(
    U2AF1a = mean(ENST00000291552.9),
    U2AF1b = mean(ENST00000380276.6),
    U2AF1c = mean(ENST00000464750.5)
  )
boo <- as.data.frame(boo)
rownames(boo) <- boo$predicted.celltype
boo <- boo[, -1]
col_fun <- circlize::colorRamp2(breaks = seq(0, 4, 4 / 8), colors = BuenColors::jdb_palette(name = "solar_rojos"))
celltypes <- c(myeloid.clusters, TNK.clusters)
celltypes <- celltypes[which(celltypes != "Treg")]
svglite::svglite("./isoforms/figures/AML3003/heatmaps/20230104_AML3003_14_U2AF1_variants.svg", width = 2.3, height = 2.8)
Heatmap(boo[celltypes, ],
  cluster_rows = F, cluster_columns = F, col = col_fun, border = T, row_names_side = "left",
  column_labels = c("variant a", "variant b", "variant c"), row_split = c(rep("myeloid", 8), rep("TNK", 11)),
  row_title_gp = gpar(fontsize = 0), row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8)
)
dev.off()
