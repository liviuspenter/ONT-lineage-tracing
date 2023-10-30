# reanalysis data by van Galen et al., Cell 2019

library(dplyr)
library(Seurat)

source("./R/celltypes.R")

files <- list.files("./data/AML/van_Galen/", pattern = "*.dem.txt.gz")

so.merged <- NULL

for (f in files) {
  message(f)
  identifier <- stringr::str_split_fixed(f, pattern = "_", n = 2)[, 2]
  identifier <- gsub(identifier, pattern = ".dem.txt.gz", replacement = "")
  patient <- stringr::str_split_fixed(identifier, pattern = "-", n = 2)[, 1]
  timepoint <- stringr::str_split_fixed(identifier, pattern = "-", n = 2)[, 2]

  mat <- as.data.frame(data.table::fread(paste0("./data/AML/van_Galen/", f)))
  rownames(mat) <- mat$Gene
  mat <- mat[, -1]
  so <- CreateSeuratObject(mat, project = identifier, min.cells = 3, min.features = 200)
  so$patient <- patient
  so$timepoint <- timepoint
  if (is.null(so.merged)) {
    so.merged <- so
  } else {
    so.merged <- merge(so.merged, so)
  }
}

bm <- readRDS(file = "./data/reference/bm.reference.rds")
bm[["spca.annoy.neighbors"]] <- LoadAnnoyIndex(
  object = bm[["spca.annoy.neighbors"]],
  file = "./data/reference/reftmp.idx"
)

so <- so.merged
rm(so.merged)
so <- NormalizeData(so)
so <- FindVariableFeatures(so)
so <- ScaleData(so)
so <- RunPCA(so)
so <- RunUMAP(so, dims = 1:30)
so <- FindNeighbors(so)
so <- FindClusters(so, resolution = 0.3)

anchors <- FindTransferAnchors(
  reference = bm, query = so, k.filter = NA, reference.reduction = "spca",
  reference.neighbors = "spca.annoy.neighbors", dims = 1:50
)
AML <- MapQuery(
  anchorset = anchors, so, reference = bm,
  refdata = list(celltype = "celltype.l2", predicted_ADT = "ADT"),
  reference.reduction = "spca", reduction.model = "wnn.umap"
)

saveRDS("./data/AML/objects/20221003_van_Galen.rds", object = AML)
AML <- readRDS("./data/AML/objects/20221003_van_Galen.rds")

### read mutations
files <- list.files("./data/AML/van_Galen/", pattern = "*.anno.txt.gz")

mutation.df <- data.frame()
for (f in files) {
  message(f)
  identifier <- stringr::str_split_fixed(f, pattern = "_", n = 2)[, 2]
  identifier <- gsub(identifier, pattern = ".anno.txt.gz", replacement = "")
  patient <- stringr::str_split_fixed(identifier, pattern = "-", n = 2)[, 1]
  timepoint <- stringr::str_split_fixed(identifier, pattern = "-", n = 2)[, 2]

  mat <- as.data.frame(data.table::fread(paste0("./data/AML/van_Galen/", f)))
  mat <- mat[which(mat$MutTranscripts != "" | mat$WtTranscripts != ""), ]

  if (nrow(mat) > 0) {
    rownames(mat) <- mat$Cell
    mat$patient <- patient
    mat$timepoint <- timepoint
    mat$identifier <- identifier
    mat <- mat[, c("Cell", "MutTranscripts", "WtTranscripts", "patient", "timepoint", "identifier")]
    mat$MutTranscripts[is.na(mat$MutTranscripts)] <- ""
    mat$mutated <- ifelse(mat$MutTranscripts != "", "mutated", "wildtype")


    mutation.df <- rbind(mutation.df, mat)
  }
}
mutation.df$predicted.celltype <- AML$predicted.celltype[mutation.df$Cell]

### read nanopore data
files <- list.files("./data/AML/van_Galen/", pattern = "*.nanopore.txt.gz")

for (f in files) {
  message(f)
  identifier <- stringr::str_split_fixed(f, pattern = "_", n = 2)[, 2]
  identifier <- gsub(identifier, pattern = ".anno.txt.gz", replacement = "")
  patient <- stringr::str_split_fixed(identifier, pattern = "-", n = 2)[, 1]
  timepoint <- stringr::str_split_fixed(identifier, pattern = "-", n = 2)[, 2]

  mat <- as.data.frame(data.table::fread(paste0("./data/AML/van_Galen/", f)))
  mat <- mat[which(mat$NanoporeTranscripts != ""), ]
  mat$predicted.celltype <- AML$predicted.celltype[mat$Cell]
  mat$mutated <- "wildtype"
  mat$mutated[which(grepl("RUNX1-RUNX1T1", mat$NanoporeTranscripts))] <- "mutated"
  mat$mutated[which(grepl("FLT3.ITD/", mat$NanoporeTranscripts))] <- "mutated"
  mat$mutated[which(grepl("TP53.Q144P/", mat$NanoporeTranscripts))] <- "mutated"
  mat$mutated[which(grepl("TP53.P152R/", mat$NanoporeTranscripts))] <- "mutated"
  colnames(mat) <- c("Cell", "MutTranscripts", "predicted.celltype", "mutated")
  rownames(mat) <- mat$Cell
  mat$patient <- patient
  mat$timepoint <- timepoint
  mat$identifier <- identifier
  mat$WtTranscripts <- NA

  mutation.df <- rbind(mutation.df, mat)
}

### plot mutations

AML$genotyped <- "no"
AML$genotyped[mutation.df$Cell] <- "yes"
boo <- subset(AML, genotyped == "yes")
DimPlot(boo, reduction = "ref.umap", cells.highlight = list(
  "mutated" = mutation.df$Cell[which(mutation.df$mutated == "mutated")],
  "wildtype" = mutation.df$Cell[which(mutation.df$mutated == "wildtype")]
)) +
  scale_color_manual(values = c("wildtype" = "black", "mutated" = "red")) +
  NoLegend() +
  NoAxes()
ggsave("./AML/figures/van_Galen/UMAP/20221006_AML_mutations_van_Galen.png", width = 4, height = 4, dpi = 600)

### statistics
mutation.df.statistics <- mutation.df %>%
  filter(timepoint == "D0") %>%
  filter(predicted.celltype %in% c(myeloid.clusters, TNK.clusters, "Prog_Mk", "Prog_RBC")) %>%
  group_by(predicted.celltype, identifier) %>%
  summarize(
    mut.cells = length(which(mutated == "mutated")),
    wt.cells = length(which(mutated == "wildtype")),
    mut.freq = length(which(mutated == "mutated")) / length(predicted.celltype)
  ) %>%
  filter((mut.cells + wt.cells) > 2)
mutation.df.statistics$predicted.celltype <- factor(mutation.df.statistics$predicted.celltype, levels = names(AML.combined.colors))
ggplot(
  mutation.df.statistics,
  aes(x = predicted.celltype, y = 100 * mut.freq, color = predicted.celltype)
) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, color = "black", width = 0.5, size = 0.25) +
  geom_jitter(width = 0.2, size = 0.5) +
  scale_y_continuous("% mutated") +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text("Arial", size = 10, angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank()
  )
ggsave("./AML/figures/van_Galen/plots/20221006_AML_mutations_frequency_van_Galen.svg", width = 2.3, height = 2)

ggplot(
  mutation.df.statistics,
  aes(x = predicted.celltype, y = mut.cells + wt.cells, color = predicted.celltype)
) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, color = "black", width = 0.5, size = 0.25) +
  geom_jitter(width = 0.2, size = 0.5) +
  scale_y_log10("cells genotyped") +
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text("Arial", size = 10, angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank()
  )
ggsave("./AML/figures/van_Galen/plots/20221006_AML_cells_genotyped_van_Galen.svg", width = 2.3, height = 2)

View(mutation.df %>% filter(timepoint == "D0") %>% group_by(predicted.celltype, identifier) %>% summarize(
  mut.cells = length(which(mutated == "mutated")),
  wt.cells = length(which(mutated == "wildtype")),
  mut.freq = length(which(mutated == "mutated")) / length(predicted.celltype)
))

View(mutation.df %>% filter(predicted.celltype == "Prog_RBC") %>% group_by(identifier) %>% summarize(
  mut.cells = length(which(mutated == "mutated")),
  wt.cells = length(which(mutated == "wildtype")),
  mut.freq = length(which(mutated == "mutated")) / length(predicted.celltype)
))
