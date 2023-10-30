# perform analyses of mixing experiment K562 with Kasumi-1

# beside the K562 and Kasumi-1 cells, there appears to be a third population 'X' in the data
# as the origin of these cells is unclear, they were filtered out for the purpose of this analysis

# the mixing steps are as follows:

# Pool116_3: K562:Kasumi-1 1:100
# Pool116_4: K562:Kasumi-1 1:15
# Pool116_5: K562:Kasumi-1 1.5:1
# Pool116_6: K562:Kasumi-1 100:1

library(dplyr)
library(ggplot2)
library(Seurat)

cell.line.colors <- c("Kasumi1" = "firebrick", "K562" = "#ffcc00")

# create Seurat object
seurat.data <- Read10X(data.dir = paste0("./data/benchmarking/Pool116_3/filtered_feature_bc_matrix/"))
so <- CreateSeuratObject(counts = seurat.data, project = "K562.1", min.cells = 3, min.features = 200)
so <- RenameCells(so, add.cell.id = "K562.1")
seurat.data <- Read10X(data.dir = paste0("./data/benchmarking/Pool116_4/filtered_feature_bc_matrix/"))
so.2 <- CreateSeuratObject(counts = seurat.data, project = "K562.10", min.cells = 3, min.features = 200)
so.2 <- RenameCells(so.2, add.cell.id = "K562.10")
seurat.data <- Read10X(data.dir = paste0("./data/benchmarking/Pool116_5/filtered_feature_bc_matrix/"))
so.3 <- CreateSeuratObject(counts = seurat.data, project = "K562.50", min.cells = 3, min.features = 200)
so.3 <- RenameCells(so.3, add.cell.id = "K562.50")
seurat.data <- Read10X(data.dir = paste0("./data/benchmarking/Pool116_6/filtered_feature_bc_matrix/"))
so.4 <- CreateSeuratObject(counts = seurat.data, project = "K562.99", min.cells = 3, min.features = 200)
so.4 <- RenameCells(so.4, add.cell.id = "K562.99")

so <- merge(so, list(so.2, so.3, so.4))
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
so <- FindVariableFeatures(so)
so <- NormalizeData(so)
so <- ScaleData(so)
so <- RunPCA(so)
so <- RunUMAP(so, dims = 1:10)
so <- FindNeighbors(so)
so <- FindClusters(so, resolution = 0.1)
saveRDS(so, file = "./data/benchmarking/objects/20220815_mixing_experiment_1_10_50_99.rds")
mixing.so <- readRDS("./data/benchmarking/objects/20220815_mixing_experiment_1_10_50_99.rds")

# annotate cell lines by gene expression clusters
mixing.so$cell.line.expr <- "K562"
mixing.so$cell.line.expr[which(mixing.so$seurat_clusters %in% c(1, 4))] <- "Kasumi1"
mixing.so$cell.line.expr[which(mixing.so$seurat_clusters %in% c(2))] <- "X"

# read nanoranger output for TP53 mutations
mixing.data.1 <- as.data.frame(nanoranger.R::extract_mutation(
  BC.data.file = "./data/benchmarking/K562_Kasumi1_1/K562_Kasumi1_pileup_TP53.csv.gz",
  ALT = "T", REF = "C"
))
mixing.data.1$bc <- paste0("K562.1_", mixing.data.1$bc, "-1")
rownames(mixing.data.1) <- mixing.data.1$bc

mixing.data.10 <- as.data.frame(nanoranger.R::extract_mutation(
  BC.data.file = "./data/benchmarking/K562_Kasumi1_10/K562_Kasumi1_pileup_TP53.csv.gz",
  ALT = "T", REF = "C"
))
mixing.data.10$bc <- paste0("K562.10_", mixing.data.10$bc, "-1")
rownames(mixing.data.10) <- mixing.data.10$bc

mixing.data.50 <- as.data.frame(nanoranger.R::extract_mutation(
  BC.data.file = "./data/benchmarking/K562_Kasumi1_50/K562_Kasumi1_pileup_TP53.csv.gz",
  ALT = "T", REF = "C"
))
mixing.data.50$bc <- paste0("K562.50_", mixing.data.50$bc, "-1")
rownames(mixing.data.50) <- mixing.data.50$bc

mixing.data.99 <- as.data.frame(nanoranger.R::extract_mutation(
  BC.data.file = "./data/benchmarking/K562_Kasumi1_99/K562_Kasumi1_pileup_TP53.csv.gz",
  ALT = "T", REF = "C"
))
mixing.data.99$bc <- paste0("K562.99_", mixing.data.99$bc, "-1")
rownames(mixing.data.99) <- mixing.data.99$bc

mixing.data.1.raw <- data.table::fread("./data/benchmarking/K562_Kasumi1_1/K562_Kasumi1_pileup_TP53.csv.gz")
mixing.data.10.raw <- data.table::fread("./data/benchmarking/K562_Kasumi1_10/K562_Kasumi1_pileup_TP53.csv.gz")
mixing.data.50.raw <- data.table::fread("./data/benchmarking/K562_Kasumi1_50/K562_Kasumi1_pileup_TP53.csv.gz")
mixing.data.99.raw <- data.table::fread("./data/benchmarking/K562_Kasumi1_99/K562_Kasumi1_pileup_TP53.csv.gz")
mixing.data.1.count <- as.data.frame(mixing.data.1.raw %>% group_by(bc) %>% summarize(count = length(bc)))
rownames(mixing.data.1.count) <- paste0("K562.1_", mixing.data.1.count$bc, "-1")
mixing.data.1.count$mutated <- mixing.data.1[rownames(mixing.data.1.count), "mutated"]
mixing.data.10.count <- as.data.frame(mixing.data.10.raw %>% group_by(bc) %>% summarize(count = length(bc)))
rownames(mixing.data.10.count) <- paste0("K562.10_", mixing.data.10.count$bc, "-1")
mixing.data.10.count$mutated <- mixing.data.10[rownames(mixing.data.10.count), "mutated"]
mixing.data.50.count <- as.data.frame(mixing.data.50.raw %>% group_by(bc) %>% summarize(count = length(bc)))
rownames(mixing.data.50.count) <- paste0("K562.50_", mixing.data.50.count$bc, "-1")
mixing.data.50.count$mutated <- mixing.data.50[rownames(mixing.data.50.count), "mutated"]
mixing.data.99.count <- as.data.frame(mixing.data.99.raw %>% group_by(bc) %>% summarize(count = length(bc)))
rownames(mixing.data.99.count) <- paste0("K562.99_", mixing.data.99.count$bc, "-1")
mixing.data.99.count$mutated <- mixing.data.99[rownames(mixing.data.99.count), "mutated"]

mixing.data.1.reduced <- mixing.data.1[intersect(mixing.data.1$bc, colnames(mixing.so)), ]
mixing.data.1.reduced <- cbind(mixing.data.1.reduced, mixing.so@reductions$umap@cell.embeddings[mixing.data.1.reduced$bc, ])
mixing.data.10.reduced <- mixing.data.10[intersect(mixing.data.10$bc, colnames(mixing.so)), ]
mixing.data.10.reduced <- cbind(mixing.data.10.reduced, mixing.so@reductions$umap@cell.embeddings[mixing.data.10.reduced$bc, ])
mixing.data.50.reduced <- mixing.data.50[intersect(mixing.data.50$bc, colnames(mixing.so)), ]
mixing.data.50.reduced <- cbind(mixing.data.50.reduced, mixing.so@reductions$umap@cell.embeddings[mixing.data.50.reduced$bc, ])
mixing.data.99.reduced <- mixing.data.99[intersect(mixing.data.99$bc, colnames(mixing.so)), ]
mixing.data.99.reduced <- cbind(mixing.data.99.reduced, mixing.so@reductions$umap@cell.embeddings[mixing.data.99.reduced$bc, ])

mixing.data.1.reduced$count <- mixing.data.1.count[rownames(mixing.data.1.reduced), "count"]
mixing.data.10.reduced$count <- mixing.data.10.count[rownames(mixing.data.10.reduced), "count"]
mixing.data.50.reduced$count <- mixing.data.50.count[rownames(mixing.data.50.reduced), "count"]
mixing.data.99.reduced$count <- mixing.data.99.count[rownames(mixing.data.99.reduced), "count"]

# read souporcell output and annotate cell lines by genetic information
clusters.1 <- data.table::fread("./data/benchmarking/Pool116_3/souporcell/clusters.tsv")
clusters.1$barcode <- paste0("K562.1_", clusters.1$barcode)
clusters.1$cell.line <- "unassigned"
clusters.1$cell.line[which(clusters.1$assignment == "0")] <- "X"
clusters.1$cell.line[which(clusters.1$assignment == "1")] <- "K562"
clusters.1$cell.line[which(clusters.1$assignment == "2")] <- "Kasumi1"

clusters.10 <- data.table::fread("./data/benchmarking/Pool116_4/souporcell/clusters.tsv")
clusters.10$barcode <- paste0("K562.10_", clusters.10$barcode)
clusters.10$cell.line <- "unassigned"
clusters.10$cell.line[which(clusters.10$assignment == "0")] <- "K562"
clusters.10$cell.line[which(clusters.10$assignment == "1")] <- "X"
clusters.10$cell.line[which(clusters.10$assignment == "2")] <- "Kasumi1"

clusters.50 <- data.table::fread("./data/benchmarking/Pool116_5/souporcell/clusters.tsv")
clusters.50$barcode <- paste0("K562.50_", clusters.50$barcode)
clusters.50$cell.line <- "unassigned"
clusters.50$cell.line[which(clusters.50$assignment == "0")] <- "X"
clusters.50$cell.line[which(clusters.50$assignment == "1")] <- "Kasumi1"
clusters.50$cell.line[which(clusters.50$assignment == "2")] <- "K562"

# souporcell doesn't work with the condition of 1:99 - use gene expression annotation
clusters.99 <- data.table::fread("./data/benchmarking/Pool116_6/souporcell/clusters.tsv")
clusters.99$barcode <- paste0("K562.99_", clusters.99$barcode)
clusters.99$cell.line <- "unassigned"
clusters.99$cell.line[which(clusters.99$assignment == "0")] <- "K562"
clusters.99$cell.line[which(clusters.99$assignment == "1")] <- "Kasumi1"
clusters.99$cell.line <- mixing.so$cell.line.expr[clusters.99$barcode]

mixing.so$cell.line <- "unassigned"
mixing.so$cell.line[clusters.1$barcode] <- clusters.1$cell.line
mixing.so$cell.line[clusters.10$barcode] <- clusters.10$cell.line
mixing.so$cell.line[clusters.50$barcode] <- clusters.50$cell.line
mixing.so$cell.line[clusters.99$barcode] <- clusters.99$cell.line

clusters.1$mutated <- mixing.data.1[clusters.1$barcode, "mutated"]
clusters.10$mutated <- mixing.data.10[clusters.10$barcode, "mutated"]
clusters.50$mutated <- mixing.data.50[clusters.50$barcode, "mutated"]
clusters.99$mutated <- mixing.data.99[clusters.99$barcode, "mutated"]

mixing.data.1.count[clusters.1$barcode, "cell.line"] <- clusters.1$cell.line
mixing.data.10.count[clusters.10$barcode, "cell.line"] <- clusters.10$cell.line
mixing.data.50.count[clusters.50$barcode, "cell.line"] <- clusters.50$cell.line
# mixing.data.99.count[clusters.99$barcode, 'cell.line'] = clusters.99$cell.line

mixing.data.99.count[colnames(mixing.so)[which(mixing.so$orig.ident == "K562.99")], "cell.line"] <- mixing.so$cell.line.expr[which(mixing.so$orig.ident == "K562.99")]

### plot knee plots for each cell line
K562.count.1 <- mixing.data.1.count[which(mixing.data.1.count$cell.line == "K562"), ]
K562.count.1 <- K562.count.1[order(K562.count.1$count, decreasing = T), ]
K562.count.1$rank <- seq(1, nrow(K562.count.1))
Kasumi1.count.1 <- mixing.data.1.count[which(mixing.data.1.count$cell.line == "Kasumi1"), ]
Kasumi1.count.1 <- Kasumi1.count.1[order(Kasumi1.count.1$count, decreasing = T), ]
Kasumi1.count.1$rank <- seq(1, nrow(Kasumi1.count.1))

ggplot() +
  geom_line(data = K562.count.1[which(!is.na(K562.count.1$count)), ], aes(x = rank, y = count), color = cell.line.colors["K562"]) +
  geom_point(data = K562.count.1[which(!is.na(K562.count.1$count) & !is.na(K562.count.1$mutated)), ], aes(x = rank, y = count, color = mutated), size = 1) +
  geom_line(data = Kasumi1.count.1[which(!is.na(Kasumi1.count.1$count)), ], aes(x = rank, y = count), color = cell.line.colors["Kasumi1"]) +
  geom_point(data = Kasumi1.count.1[which(!is.na(Kasumi1.count.1$count) & !is.na(Kasumi1.count.1$mutated)), ], aes(x = rank, y = count, color = mutated), size = 1) +
  geom_hline(yintercept = 5) +
  scale_color_manual(values = c("mutated" = "magenta", "wildtype" = "cyan")) +
  scale_x_log10("cell rank") +
  scale_y_log10("read count") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/mixing/20220811_K562_Kasumi1_1_knee_plot.svg", width = 2, height = 2)

K562.count.10 <- mixing.data.10.count[which(mixing.data.10.count$cell.line == "K562"), ]
K562.count.10 <- K562.count.10[order(K562.count.10$count, decreasing = T), ]
K562.count.10$rank <- seq(1, nrow(K562.count.10))
Kasumi1.count.10 <- mixing.data.10.count[which(mixing.data.10.count$cell.line == "Kasumi1"), ]
Kasumi1.count.10 <- Kasumi1.count.10[order(Kasumi1.count.10$count, decreasing = T), ]
Kasumi1.count.10$rank <- seq(1, nrow(Kasumi1.count.10))

ggplot() +
  geom_line(data = K562.count.10[which(!is.na(K562.count.10$count)), ], aes(x = rank, y = count), color = cell.line.colors["K562"]) +
  geom_point(data = K562.count.10[which(!is.na(K562.count.10$count) & !is.na(K562.count.10$mutated)), ], aes(x = rank, y = count, color = mutated), size = 1) +
  geom_line(data = Kasumi1.count.10[which(!is.na(Kasumi1.count.10$count)), ], aes(x = rank, y = count), color = cell.line.colors["Kasumi1"]) +
  geom_point(data = Kasumi1.count.10[which(!is.na(Kasumi1.count.10$count) & !is.na(Kasumi1.count.10$mutated)), ], aes(x = rank, y = count, color = mutated), size = 1) +
  geom_hline(yintercept = 5) +
  scale_color_manual(values = c("mutated" = "magenta", "wildtype" = "cyan")) +
  scale_x_log10("cell rank") +
  scale_y_log10("read count") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/mixing/20220811_K562_Kasumi1_10_knee_plot.svg", width = 2, height = 2)

K562.count.50 <- mixing.data.50.count[which(mixing.data.50.count$cell.line == "K562"), ]
K562.count.50 <- K562.count.50[order(K562.count.50$count, decreasing = T), ]
K562.count.50$rank <- seq(1, nrow(K562.count.50))
Kasumi1.count.50 <- mixing.data.50.count[which(mixing.data.50.count$cell.line == "Kasumi1"), ]
Kasumi1.count.50 <- Kasumi1.count.50[order(Kasumi1.count.50$count, decreasing = T), ]
Kasumi1.count.50$rank <- seq(1, nrow(Kasumi1.count.50))

ggplot() +
  geom_line(data = K562.count.50[which(!is.na(K562.count.50$count)), ], aes(x = rank, y = count), color = cell.line.colors["K562"]) +
  geom_point(data = K562.count.50[which(!is.na(K562.count.50$count) & !is.na(K562.count.50$mutated)), ], aes(x = rank, y = count, color = mutated), size = 1) +
  geom_line(data = Kasumi1.count.50[which(!is.na(Kasumi1.count.50$count)), ], aes(x = rank, y = count), color = cell.line.colors["Kasumi1"]) +
  geom_point(data = Kasumi1.count.50[which(!is.na(Kasumi1.count.50$count) & !is.na(Kasumi1.count.50$mutated)), ], aes(x = rank, y = count, color = mutated), size = 1) +
  geom_hline(yintercept = 5) +
  scale_color_manual(values = c("mutated" = "magenta", "wildtype" = "cyan")) +
  scale_x_log10("cell rank") +
  scale_y_log10("read count") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/mixing//20220811_K562_Kasumi1_50_knee_plot.svg", width = 2, height = 2)


K562.count.99 <- mixing.data.99.count[which(mixing.data.99.count$cell.line == "K562"), ]
K562.count.99 <- K562.count.99[order(K562.count.99$count, decreasing = T), ]
K562.count.99$rank <- seq(1, nrow(K562.count.99))
Kasumi1.count.99 <- mixing.data.99.count[which(mixing.data.99.count$cell.line == "Kasumi1"), ]
Kasumi1.count.99 <- Kasumi1.count.99[order(Kasumi1.count.99$count, decreasing = T), ]
Kasumi1.count.99$rank <- seq(1, nrow(Kasumi1.count.99))

ggplot() +
  geom_line(data = K562.count.99[which(!is.na(K562.count.99$count)), ], aes(x = rank, y = count), color = cell.line.colors["K562"]) +
  geom_point(data = K562.count.99[which(!is.na(K562.count.99$count) & !is.na(K562.count.99$mutated)), ], aes(x = rank, y = count, color = mutated), size = 1) +
  geom_line(data = Kasumi1.count.99[which(!is.na(Kasumi1.count.99$count)), ], aes(x = rank, y = count), color = cell.line.colors["Kasumi1"]) +
  geom_point(data = Kasumi1.count.99[which(!is.na(Kasumi1.count.99$count) & !is.na(Kasumi1.count.99$mutated)), ], aes(x = rank, y = count, color = mutated), size = 1) +
  geom_hline(yintercept = 5) +
  scale_color_manual(values = c("mutated" = "magenta", "wildtype" = "cyan")) +
  scale_x_log10("cell rank") +
  scale_y_log10("read count") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/mixing//20220811_K562_Kasumi1_99_knee_plot.svg", width = 2, height = 2)

### plots for statistics
mutation.statistics <- bind_rows(
  clusters.1 %>% filter(barcode %in% colnames(mixing.so)) %>%
    group_by(cell.line) %>% summarize(
      mut = length(which(mutated == "mutated")),
      wt = length(which(mutated == "wildtype")),
      cells = length(cell.line),
      sample = "K562.1"
    ),
  clusters.10 %>% filter(barcode %in% colnames(mixing.so)) %>%
    group_by(cell.line) %>% summarize(
      mut = length(which(mutated == "mutated")),
      wt = length(which(mutated == "wildtype")),
      cells = length(cell.line),
      sample = "K562.10"
    ),
  clusters.50 %>% filter(barcode %in% colnames(mixing.so)) %>%
    group_by(cell.line) %>% summarize(
      mut = length(which(mutated == "mutated")),
      wt = length(which(mutated == "wildtype")),
      cells = length(cell.line),
      sample = "K562.50"
    ),
  clusters.99 %>% filter(barcode %in% colnames(mixing.so)) %>%
    group_by(cell.line) %>% summarize(
      mut = length(which(mutated == "mutated")),
      wt = length(which(mutated == "wildtype")),
      cells = length(cell.line),
      sample = "K562.99"
    )
)

mutation.statistics$genotyped <- mutation.statistics$mut + mutation.statistics$wt
mutation.statistics$genotyped.freq <- mutation.statistics$genotyped / mutation.statistics$cells
mutation.statistics$mut.freq <- mutation.statistics$mut / mutation.statistics$genotyped
mutation.statistics$wt.freq <- mutation.statistics$wt / mutation.statistics$genotyped

ggplot(mutation.statistics[which(mutation.statistics$cell.line %in% c("K562", "Kasumi1")), ], aes(x = genotyped, y = 100 * wt.freq)) +
  geom_point(aes(color = cell.line)) +
  scale_color_manual(values = cell.line.colors) +
  scale_x_continuous("cells genotyped") +
  scale_y_continuous("% TP53wt") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/mixing//20220811_K562_Kasumi1_wildtype_freq.svg", width = 1.5, height = 1.5)

ggplot(mutation.statistics[which(mutation.statistics$cell.line %in% c("K562", "Kasumi1")), ], aes(x = genotyped, y = 100 * mut.freq)) +
  geom_point(aes(color = cell.line)) +
  scale_color_manual(values = cell.line.colors) +
  scale_x_continuous("cells genotyped") +
  scale_y_continuous("% TP53R248G") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/mixing/20220811_K562_Kasumi1_mutated_freq.svg", width = 1.5, height = 1.5)

ggplot(mutation.statistics[which(mutation.statistics$cell.line %in% c("K562", "Kasumi1")), ], aes(x = sample, y = 100 * genotyped.freq)) +
  geom_col(aes(fill = cell.line), position = "dodge") +
  scale_fill_manual(values = cell.line.colors) +
  scale_x_discrete(labels = c("1% K562", "6% K562", "67% K562", "99% K562")) +
  scale_y_continuous("% genotyped") +
  theme_classic() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./benchmarking/figures/mixing/20220811_K562_Kasumi1_genotyping_freq.svg", width = 1.5, height = 2)

### insert Q136fs
TP53.insert.1 <- nanoranger.R::extract_indel("./data/benchmarking/K562_Kasumi1_1/K562_Kasumi1_pileup_TP53_2.csv.gz", REF = "T", CONSENSUS = 1)
TP53.insert.1$bc <- paste0("K562.1_", TP53.insert.1$bc, "-1")

TP53.insert.10 <- nanoranger.R::extract_indel("./data/benchmarking/K562_Kasumi1_10/K562_Kasumi1_pileup_TP53_2.csv.gz", REF = "T", CONSENSUS = 1)
TP53.insert.10$bc <- paste0("K562.10_", TP53.insert.10$bc, "-1")

TP53.insert.50 <- nanoranger.R::extract_indel("./data/benchmarking/K562_Kasumi1_50/K562_Kasumi1_pileup_TP53_2.csv.gz", REF = "T", CONSENSUS = 1)
TP53.insert.50$bc <- paste0("K562.50_", TP53.insert.50$bc, "-1")

TP53.insert.99 <- nanoranger.R::extract_indel("./data/benchmarking/K562_Kasumi1_99/K562_Kasumi1_pileup_TP53_2.csv.gz", REF = "T", CONSENSUS = 1)
TP53.insert.99$bc <- paste0("K562.99_", TP53.insert.99$bc, "-1")

TP53.insert.counts.1 <- as.data.frame(data.table::fread("./data/benchmarking/K562_Kasumi1_1/K562_Kasumi1_pileup_TP53_2.csv.gz") %>%
  filter(paste0("K562.1_", bc, "-1") %in% colnames(mixing.so)) %>%
  group_by(bc) %>%
  summarize(count = length(bc)))
TP53.insert.counts.1$bc <- paste0("K562.1_", TP53.insert.counts.1$bc, "-1")
rownames(TP53.insert.counts.1) <- TP53.insert.counts.1$bc
TP53.insert.counts.1[clusters.1$barcode, "cell.line"] <- clusters.1$cell.line
TP53.insert.counts.1[TP53.insert.1$bc, "mutated"] <- TP53.insert.1$mutated

TP53.insert.counts.1 <- TP53.insert.counts.1 %>%
  filter(!is.na(count)) %>%
  filter(cell.line %in% c("K562", "Kasumi1")) %>%
  group_by(cell.line) %>%
  arrange(-count, .by_group = T) %>%
  mutate(rank = row_number())

ggplot() +
  geom_line(data = TP53.insert.counts.1, aes(x = rank, y = count, color = cell.line)) +
  geom_point(data = TP53.insert.counts.1[which(TP53.insert.counts.1$count > 4), ], aes(x = rank, y = count, color = mutated), size = 1) +
  geom_hline(yintercept = 5) +
  scale_color_manual(values = c("mutated" = "magenta", "wildtype" = "cyan", cell.line.colors)) +
  scale_x_log10("cell rank") +
  scale_y_log10("read count") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/mixing/20220815_K562_Kasumi1_1_knee_plot_TP53_insert.svg", width = 2, height = 2)

###
TP53.insert.counts.10 <- as.data.frame(data.table::fread("./data/benchmarking/K562_Kasumi1_10/K562_Kasumi1_pileup_TP53_2.csv.gz") %>%
  filter(paste0("K562.10_", bc, "-1") %in% colnames(mixing.so)) %>%
  group_by(bc) %>%
  summarize(count = length(bc)))
TP53.insert.counts.10$bc <- paste0("K562.10_", TP53.insert.counts.10$bc, "-1")
rownames(TP53.insert.counts.10) <- TP53.insert.counts.10$bc
TP53.insert.counts.10[clusters.10$barcode, "cell.line"] <- clusters.10$cell.line
TP53.insert.counts.10[TP53.insert.10$bc, "mutated"] <- TP53.insert.10$mutated

TP53.insert.counts.10 <- TP53.insert.counts.10 %>%
  filter(!is.na(count)) %>%
  filter(cell.line %in% c("K562", "Kasumi1")) %>%
  group_by(cell.line) %>%
  arrange(-count, .by_group = T) %>%
  mutate(rank = row_number())

ggplot() +
  geom_line(data = TP53.insert.counts.10, aes(x = rank, y = count, color = cell.line)) +
  geom_point(data = TP53.insert.counts.10[which(TP53.insert.counts.10$count > 4), ], aes(x = rank, y = count, color = mutated), size = 1) +
  geom_hline(yintercept = 5) +
  scale_color_manual(values = c("mutated" = "magenta", "wildtype" = "cyan", cell.line.colors)) +
  scale_x_log10("cell rank") +
  scale_y_log10("read count") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/mixing/20220815_K562_Kasumi1_10_knee_plot_TP53_insert.svg", width = 2, height = 2)


###
TP53.insert.counts.50 <- as.data.frame(data.table::fread("./data/benchmarking/K562_Kasumi1_50/K562_Kasumi1_pileup_TP53_2.csv.gz") %>%
  filter(paste0("K562.50_", bc, "-1") %in% colnames(mixing.so)) %>%
  group_by(bc) %>%
  summarize(count = length(bc)))
TP53.insert.counts.50$bc <- paste0("K562.50_", TP53.insert.counts.50$bc, "-1")
rownames(TP53.insert.counts.50) <- TP53.insert.counts.50$bc
TP53.insert.counts.50[clusters.50$barcode, "cell.line"] <- clusters.50$cell.line
TP53.insert.counts.50[TP53.insert.50$bc, "mutated"] <- TP53.insert.50$mutated

TP53.insert.counts.50 <- TP53.insert.counts.50 %>%
  filter(!is.na(count)) %>%
  filter(cell.line %in% c("K562", "Kasumi1")) %>%
  group_by(cell.line) %>%
  arrange(-count, .by_group = T) %>%
  mutate(rank = row_number())

ggplot() +
  geom_line(data = TP53.insert.counts.50, aes(x = rank, y = count, color = cell.line)) +
  geom_point(data = TP53.insert.counts.50[which(TP53.insert.counts.50$count > 4), ], aes(x = rank, y = count, color = mutated), size = 1) +
  geom_hline(yintercept = 5) +
  scale_color_manual(values = c("mutated" = "magenta", "wildtype" = "cyan", cell.line.colors)) +
  scale_x_log10("cell rank") +
  scale_y_log10("read count") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/mixing//20220815_K562_Kasumi1_50_knee_plot_TP53_insert.svg", width = 2, height = 2)


###
TP53.insert.counts.99 <- as.data.frame(data.table::fread("./data/benchmarking/K562_Kasumi1_99/K562_Kasumi1_pileup_TP53_2.csv.gz") %>%
  filter(paste0("K562.99_", bc, "-1") %in% colnames(mixing.so)) %>%
  group_by(bc) %>%
  summarize(count = length(bc)))
TP53.insert.counts.99$bc <- paste0("K562.99_", TP53.insert.counts.99$bc, "-1")
rownames(TP53.insert.counts.99) <- TP53.insert.counts.99$bc
TP53.insert.counts.99[clusters.99$barcode, "cell.line"] <- clusters.99$cell.line
TP53.insert.counts.99[TP53.insert.99$bc, "mutated"] <- TP53.insert.99$mutated

TP53.insert.counts.99 <- TP53.insert.counts.99 %>%
  filter(!is.na(count)) %>%
  filter(cell.line %in% c("K562", "Kasumi1")) %>%
  group_by(cell.line) %>%
  arrange(-count, .by_group = T) %>%
  mutate(rank = row_number())

ggplot() +
  geom_line(data = TP53.insert.counts.99, aes(x = rank, y = count, color = cell.line)) +
  geom_point(data = TP53.insert.counts.99[which(TP53.insert.counts.99$count > 4), ], aes(x = rank, y = count, color = mutated), size = 1) +
  geom_hline(yintercept = 5) +
  scale_color_manual(values = c("mutated" = "magenta", "wildtype" = "cyan", cell.line.colors)) +
  scale_x_log10("cell rank") +
  scale_y_log10("read count") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/mixing/20220815_K562_Kasumi1_99_knee_plot_TP53_insert.svg", width = 2, height = 2)


mutation.statistics.insert <- bind_rows(
  TP53.insert.counts.1 %>% group_by(cell.line) %>% summarize(
    mut = length(which(mutated == "mutated")),
    wt = length(which(mutated == "wildtype")),
    cells = length(which(mixing.so$orig.ident == "K562.1" &
      mixing.so$cell.line == cell.line)),
    sample = "K562.1"
  ),
  TP53.insert.counts.10 %>% group_by(cell.line) %>% summarize(
    mut = length(which(mutated == "mutated")),
    wt = length(which(mutated == "wildtype")),
    cells = length(which(mixing.so$orig.ident == "K562.10" &
      mixing.so$cell.line == cell.line)),
    sample = "K562.10"
  ),
  TP53.insert.counts.50 %>% group_by(cell.line) %>% summarize(
    mut = length(which(mutated == "mutated")),
    wt = length(which(mutated == "wildtype")),
    cells = length(which(mixing.so$orig.ident == "K562.50" &
      mixing.so$cell.line == cell.line)),
    sample = "K562.50"
  ),
  TP53.insert.counts.99 %>% group_by(cell.line) %>% summarize(
    mut = length(which(mutated == "mutated")),
    wt = length(which(mutated == "wildtype")),
    cells = length(which(mixing.so$orig.ident == "K562.99" &
      mixing.so$cell.line.expr == cell.line)),
    sample = "K562.99"
  )
)

mutation.statistics.insert$genotyped <- mutation.statistics.insert$mut + mutation.statistics.insert$wt
mutation.statistics.insert$genotyped.freq <- mutation.statistics.insert$genotyped / mutation.statistics.insert$cells
mutation.statistics.insert$mut.freq <- mutation.statistics.insert$mut / mutation.statistics.insert$genotyped
mutation.statistics.insert$wt.freq <- mutation.statistics.insert$wt / mutation.statistics.insert$genotyped

ggplot(mutation.statistics.insert, aes(x = genotyped, y = 100 * wt.freq)) +
  geom_point(aes(color = cell.line)) +
  scale_color_manual(values = cell.line.colors) +
  scale_x_continuous("cells genotyped") +
  scale_y_continuous("% TP53wt") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/20220811_K562_Kasumi1_TP53_insert_wildtype_freq.svg", width = 1.5, height = 1.5)

ggplot(mutation.statistics.insert, aes(x = genotyped, y = 100 * mut.freq)) +
  geom_point(aes(color = cell.line)) +
  scale_color_manual(values = cell.line.colors) +
  scale_x_continuous("cells genotyped") +
  scale_y_continuous("% TP53Q136fs") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/mixing/20220811_K562_Kasumi1_TP53_insert_mutated_freq.svg", width = 1.5, height = 1.5)

ggplot(mutation.statistics.insert, aes(x = sample, y = 100 * genotyped.freq)) +
  geom_col(aes(fill = cell.line), position = "dodge") +
  scale_fill_manual(values = cell.line.colors) +
  scale_x_discrete(labels = c("1% K562", "6% K562", "67% K562", "99% K562")) +
  scale_y_continuous("% genotyped") +
  theme_classic() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./benchmarking/figures/mixing/20220811_K562_Kasumi1_TP53_insert_genotyping_freq.svg", width = 1.5, height = 2)

VlnPlot(subset(mixing.so, cell.line %in% c("K562", "Kasumi1")), group.by = "cell.line", features = "TP53", cols = cell.line.colors) +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 10, color = "black", "Arial", hjust = 1, vjust = 0.5),
    legend.position = "none"
  )
ggsave("./benchmarking/figures/mixing/20220817_K562_Kasumi1_TP53_expression.svg", width = 1.1, height = 2)

# non-zero values for TP53
TP53.data <- GetAssayData(subset(mixing.so, cell.line %in% c("K562", "Kasumi1")), slot = "counts", assay = "RNA")["TP53", ]
table(TP53.data[colnames(mixing.so)[which(mixing.so$cell.line == "K562")]])
table(TP53.data[colnames(mixing.so)[which(mixing.so$cell.line == "Kasumi1")]])

mixing.data.count <- bind_rows(mixing.data.1.count, mixing.data.10.count, mixing.data.50.count)
mixing.data.count

# downsampling experiment for TP53
downsample.df <- data.frame()
for (reads in seq(10000, 500000, 10000)) {
  mixing.data.downsampled <- as.data.frame(extract_mutation(
    BC.data.file = "./data/Mehdi/K562_Kasumi1_50/K562_Kasumi1_pileup_TP53.csv.gz",
    ALT = "T", REF = "C", downsample = reads
  ))
  mixing.data.downsampled$bc <- paste0("K562.50_", mixing.data.downsampled$bc, "-1")
  rownames(mixing.data.downsampled) <- mixing.data.downsampled$bc
  mixing.data.downsampled$cell.line <- mixing.so$cell.line[mixing.data.downsampled$bc]
  downsample.df <- rbind(downsample.df, data.frame(
    reads = reads,
    K562 = length(which(mixing.data.downsampled$cell.line == "K562")),
    Kasumi1 = length(which(mixing.data.downsampled$cell.line == "Kasumi1"))
  ))
}
write.csv2(downsample.df, "./data/benchmarking/K562_Kasumi1_50/20220817_downsample_experiment.csv", quote = F)
downsample.df <- read.csv2("./data/benchmarking/K562_Kasumi1_50/20220817_downsample_experiment.csv", row.names = 1) %>%
  as.data.frame()

ggplot(reshape2::melt(downsample.df, id.vars = "reads"), aes(x = reads, y = value, color = variable)) +
  geom_line(aes(group = variable)) +
  scale_x_continuous("1000x reads", breaks = c(100000, 200000, 300000, 400000, 500000), labels = c("100", "200", "300", "400", "500")) +
  scale_y_continuous("cells genotyped") +
  scale_color_manual(values = cell.line.colors) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    legend.position = "none"
  )
ggsave("./benchmarking/figures/mixing/20220817_K562_Kasumi1_TP53_downsampling.svg", width = 2, height = 1.5)

### BCRABL1 and RUNX1/RUNX1T1
BCRABL1.data <- nanoranger.R::extract_fusion_gene("./data/benchmarking/K562_Kasumi1_50/K562_Kasumi1_BCRABL1.csv.gz",
  WILDTYPE = "ABL1", FUSION = "BCR", downsample = NA
)
BCRABL1.data$bc <- paste0("K562.50_", BCRABL1.data$bc, "-1")
BCRABL1.data$cell.line <- mixing.so$cell.line[BCRABL1.data$bc]
BCRABL1.data <- BCRABL1.data[which(BCRABL1.data$cell.line %in% c("K562", "Kasumi1")), ]

BCRABL1.data <- merge(BCRABL1.data, data.table::fread("./data/benchmarking/K562_Kasumi1_50/fusion_20220810.csv.gz") %>%
  filter(gene %in% c("ABL1", "BCR")) %>%
  mutate(bc = paste0("K562.50_", bc, "-1")) %>%
  filter(bc %in% colnames(mixing.so)) %>%
  group_by(bc) %>%
  summarize(count = n()),
by = "bc"
)

BCRABL1.data <- BCRABL1.data %>%
  group_by(cell.line) %>%
  arrange(-count, .by_group = T) %>%
  mutate(rank = row_number())

ggplot(BCRABL1.data, aes(x = rank, y = count)) +
  geom_line(aes(color = cell.line, group = cell.line)) +
  geom_point(aes(color = mutated), size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("mutated" = "magenta", "wildtype" = "cyan", cell.line.colors)) +
  geom_hline(yintercept = 100) +
  theme_classic() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/mixing/20220816_K562_Kasumi1_50_BCRABL1_knee_plot.svg", width = 2, height = 2)

BCRABL1.statistics <- as.data.frame(BCRABL1.data %>%
  filter(count > 100) %>%
  group_by(cell.line) %>%
  summarize(
    genotyped = length(cell.line),
    mut = length(which(mutated == "mutated")),
    wt = length(which(mutated == "wildtype"))
  ))
BCRABL1.statistics$cells <- apply(BCRABL1.statistics, MARGIN = 1, FUN = function(x) {
  length(which(mixing.so$orig.ident == "K562.50" & mixing.so$cell.line == x["cell.line"]))
})
BCRABL1.statistics$genotyped.freq <- BCRABL1.statistics$genotyped / BCRABL1.statistics$cells
BCRABL1.statistics$mut.freq <- BCRABL1.statistics$mut / BCRABL1.statistics$genotyped
BCRABL1.statistics$target <- "BCRABL1"

AML1ETO.data <- nanoranger.R::extract_fusion_gene("./data/benchmarking/K562_Kasumi1_50/K562_Kasumi1_AML_ETO.csv", WILDTYPE = "RUNX1T1", FUSION = "RUNX1", downsample = NA)
AML1ETO.data$bc <- paste0("K562.50_", AML1ETO.data$bc, "-1")
AML1ETO.data$cell.line <- mixing.so$cell.line[AML1ETO.data$bc]
AML1ETO.data <- AML1ETO.data[which(AML1ETO.data$cell.line %in% c("K562", "Kasumi1")), ]

AML1ETO.data <- merge(AML1ETO.data, data.table::fread("./data/benchmarking/K562_Kasumi1_50/fusion_20220810.csv.gz") %>%
  filter(gene %in% c("RUNX1", "RUNX1T1")) %>%
  mutate(bc = paste0("K562.50_", bc, "-1")) %>%
  filter(bc %in% colnames(mixing.so)) %>%
  group_by(bc) %>%
  summarize(count = n()),
by = "bc"
)

AML1ETO.data <- AML1ETO.data %>%
  group_by(cell.line) %>%
  arrange(-count, .by_group = T) %>%
  mutate(rank = row_number())

ggplot(AML1ETO.data, aes(x = rank, y = count)) +
  geom_line(aes(color = cell.line, group = cell.line)) +
  geom_point(aes(color = mutated), size = 1) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("mutated" = "magenta", "wildtype" = "cyan", cell.line.colors)) +
  geom_hline(yintercept = 30) +
  theme_classic() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/mixing/20220816_K562_Kasumi1_50_AML1ETO_knee_plot.svg", width = 2, height = 2)

AML1ETO.statistics <- as.data.frame(AML1ETO.data %>%
  filter(count > 100) %>%
  group_by(cell.line) %>%
  summarize(
    genotyped = length(cell.line),
    mut = length(which(mutated == "mutated")),
    wt = length(which(mutated == "wildtype"))
  ))
AML1ETO.statistics$cells <- apply(AML1ETO.statistics, MARGIN = 1, FUN = function(x) {
  length(which(mixing.so$orig.ident == "K562.50" & mixing.so$cell.line == x["cell.line"]))
})
AML1ETO.statistics$genotyped.freq <- AML1ETO.statistics$genotyped / AML1ETO.statistics$cells
AML1ETO.statistics$mut.freq <- AML1ETO.statistics$mut / AML1ETO.statistics$genotyped
AML1ETO.statistics$target <- "AML1ETO"

ggplot(rbind(BCRABL1.statistics, AML1ETO.statistics), aes(x = target, y = 100 * genotyped.freq, fill = cell.line)) +
  geom_col(aes(fill = cell.line), position = "dodge") +
  scale_fill_manual(values = cell.line.colors) +
  scale_x_discrete(labels = c("BCR::ABL1", "RUNX1::RUNX1T1"), limits = c("BCRABL1", "AML1ETO")) +
  scale_y_continuous("% genotyped") +
  theme_classic() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
ggsave("./benchmarking/figures/mixing/20220817_K562_Kasumi1_genotyping_freq_fusion.svg", width = 1.1, height = 2)

ggplot(rbind(BCRABL1.statistics, AML1ETO.statistics), aes(x = genotyped, y = 100 * mut.freq, color = cell.line, label = target)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(label.size = 0, size = 2) +
  scale_color_manual(values = cell.line.colors) +
  scale_x_continuous(" cells genotyped") +
  scale_y_continuous("% fusion") +
  theme_classic() +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 8, color = "black"),
    axis.title = element_text("Arial", size = 8, color = "black")
  )
ggsave("./benchmarking/figures/mixing/20220817_K562_Kasumi1_mut_freq_fusion.svg", width = 1.5, height = 1.5)



### read length versus duplication
BC.data <- data.table::fread("./data/benchmarking/K562_Kasumi1_50/fusion_20220810.csv.gz")
BC.data$bc <- paste0(BC.data$bc, "-1")

BC.list <- BC.data %>%
  group_by(bc, umi) %>%
  summarize(n = n()) %>%
  filter(n > 4)

# run starcode and identify UMI clusters
UMIs.collapsed <- data.frame()
for (bc in unique(BC.list$bc)) {
  message(bc)
  filehandle.in <- tempfile()
  filehandle.out <- tempfile()
  write.table(x = BC.list[which(BC.list$bc == bc), c("umi", "n")], file = filehandle.in, sep = "\t", row.names = F, quote = F, col.names = F)

  system(paste0(STARCODE, " -d 3 -i ", filehandle.in, " -o ", filehandle.out, " --print-clusters"))
  starcode.output <- as.data.frame(read.csv2(filehandle.out, sep = "\t", header = F))
  colnames(starcode.output) <- c("umi", "n", "umi.non.collapsed")
  starcode.output$bc <- bc

  UMIs.collapsed <- rbind(UMIs.collapsed, starcode.output)

  unlink(filehandle.in)
  unlink(filehandle.out)
}

BC.data.condensed <- BC.data %>%
  group_by(gene, bc, umi) %>%
  summarize(
    n = length(gene),
    qlen.median = median(qlen)
  )
BC.data.condensed$bc.umi <- paste0(BC.data.condensed$bc, ".", BC.data.condensed$umi)
UMIs.collapsed$bc.umi <- paste0(UMIs.collapsed$bc, ".", UMIs.collapsed$umi)
rownames(UMIs.collapsed) <- UMIs.collapsed$bc.umi

BC.data.condensed <- BC.data.condensed[which(BC.data.condensed$bc.umi %in% UMIs.collapsed$bc.umi), ]
BC.data.condensed$count <- UMIs.collapsed[BC.data.condensed$bc.umi, "n"]

p <- ggplot() +
  ggrastr::rasterize(ggpointdensity::geom_pointdensity(
    data = BC.data.condensed[which(BC.data.condensed$gene == "BCR"), ],
    aes(x = -qlen.median, y = count), size = 0.5
  ), dpi = 600) +
  viridis::scale_color_viridis() +
  scale_x_continuous("distance from primer", limits = c(-1500, 0), breaks = c(-1500, -1000, -500, 0), labels = c("1500", "1000", "500", "0")) +
  scale_y_log10("reads per transcript") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/mixing/20220821_BCRABL1_transcript_length_duplication_rate.svg", width = 2.5, height = 1.5, plot = p)

p <- ggplot() +
  ggrastr::rasterize(ggpointdensity::geom_pointdensity(
    data = BC.data.condensed[which(BC.data.condensed$gene == "BCR"), ],
    aes(x = -qlen.median, y = count), size = 0.5
  ), dpi = 600) +
  viridis::scale_color_viridis() +
  scale_x_continuous("distance from primer", limits = c(-1500, 0), breaks = c(-1500, -1000, -500, 0), labels = c("1500", "1000", "500", "0")) +
  scale_y_log10("reads per transcript") +
  theme_classic() +
  theme( # legend.position = 'none',
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/mixing/20220821_BCRABL1_transcript_length_duplication_rate_with_legend.svg", width = 2.5, height = 1.5, plot = p)


q <- ggplot() +
  ggrastr::rasterize(ggpointdensity::geom_pointdensity(
    data = BC.data.condensed[which(BC.data.condensed$gene == "RUNX1T1"), ],
    aes(x = -qlen.median, y = count), size = 0.5
  ), dpi = 600) +
  viridis::scale_color_viridis() +
  scale_x_continuous("distance from primer", limits = c(-500, 0)) +
  scale_y_log10("reads per transcript") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/mixing/20220821_AML1ETO_transcript_length_duplication_rate.svg", width = 2.5, height = 2, plot = q)
