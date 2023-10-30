# compare mismatches and indels between Illumina and ONT sequencing data
# identify sequencing errors according to length of homopolymers

library(dplyr)
library(ggplot2)
library(ggrepel)

# sequencing error per read for Illumina and ONT data of same TCR library
Illumina.data <- data.table::fread("./data/TCR/benchmarking/ILL_CDR3_per_read_stats.csv.gz")
Illumina.data$mismatch.freq <- Illumina.data$mismatch / Illumina.data$reflen
Illumina.data$indel.freq <- Illumina.data$indel / Illumina.data$reflen
ONT.data <- data.table::fread("./data/TCR/benchmarking/ONT_CDR3_per_read_stats_new.csv.gz")
ONT.data$mismatch.freq <- ONT.data$mismatch / ONT.data$reflen
ONT.data$indel.freq <- ONT.data$indel / ONT.data$reflen

Illumina.mean.data <- data.table::fread("./data/TCR/benchmarking/ILL_CDR3_mean_stats.csv.gz")
ONT.mean.data <- data.table::fread("./data/TCR/benchmarking/ONT_CDR3_mean_stats_new.csv.gz")
CDR3.list <- data.table::fread("./data/TCR/benchmarking/good_filt_tcr3.csv.gz") %>% as.data.frame()
rownames(CDR3.list) <- CDR3.list$cdr3_reads

# statistics
Illumina.mean.data$mismatch.rank <- rank(Illumina.mean.data$mismatch, ties.method = "random")
Illumina.mean.data$indel.rank <- rank(Illumina.mean.data$indel, ties.method = "random")
Illumina.mean.data$CDR3.nt <- CDR3.list[Illumina.mean.data$clone, "cdr3_nt"]
Illumina.mean.data$reads <- CDR3.list[Illumina.mean.data$clone, "reads"]

ONT.mean.data$mismatch.rank <- rank(ONT.mean.data$mismatch, ties.method = "random")
ONT.mean.data$indel.rank <- rank(ONT.mean.data$indel, ties.method = "random")
ONT.mean.data$CDR3.nt <- CDR3.list[ONT.mean.data$clone, "cdr3_nt"]
ONT.mean.data$reads <- CDR3.list[ONT.mean.data$clone, "cloneCount"]

summary(ONT.mean.data$mismatch)
summary(ONT.mean.data$indel)

summary(Illumina.mean.data$mismatch)
summary(Illumina.mean.data$indel)

t.test(ONT.mean.data$mismatch, Illumina.mean.data$mismatch)
t.test(ONT.mean.data$indel, Illumina.mean.data$indel)

# plot from different angles
ggplot() +
  geom_histogram(
    data = Illumina.data, aes(x = mismatch, y = 100 * after_stat(count / sum(count))),
    binwidth = 1, fill = "grey90"
  ) +
  geom_histogram(
    data = ONT.data, aes(x = mismatch, y = 100 * after_stat(count / sum(count))),
    binwidth = 1, fill = "orange", alpha = 0.5
  ) +
  scale_x_continuous("mismatches per read") +
  scale_y_continuous("% reads") +
  theme_classic() +
  theme(
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./TCR/figures/plots/20230813_TCR3_mismatch_read.svg", width = 1.5, height = 2)

ggplot() +
  geom_histogram(
    data = Illumina.data, aes(x = 100 * mismatch.freq, y = 100 * after_stat(count / sum(count))),
    binwidth = 1, fill = "grey90"
  ) +
  geom_histogram(
    data = ONT.data, aes(x = 100 * mismatch.freq, y = 100 * after_stat(count / sum(count))),
    binwidth = 1, fill = "orange", alpha = 0.5
  ) +
  scale_x_continuous("% mismatches") +
  scale_y_continuous("% reads") +
  theme_classic() +
  theme(
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./TCR/figures/plots/20230813_TCR3_mismatch_read_relative.svg", width = 1.5, height = 2)


df <- merge(Illumina.mean.data, ONT.mean.data, by = "clone")

p <- ggplot(df, aes(x = mismatch.x, y = mismatch.y)) +
  ggrastr::rasterize(ggpointdensity::geom_pointdensity(size = 0.5), dpi = 600) +
  scale_x_continuous("% mismatch Illumina") +
  scale_y_continuous("% mismatch ONT") +
  scale_color_viridis_c() +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./TCR/figures/plots/20230813_TCR_Illumina_ONT_mismatch_rate.svg", width = 1.5, height = 2.5, plot = p)

p <- ggplot(df, aes(x = indel.x, y = indel.y)) +
  ggrastr::rasterize(ggpointdensity::geom_pointdensity(size = 0.5), dpi = 600) +
  scale_x_continuous("% indel Illumina") +
  scale_y_continuous("% indel ONT") +
  scale_color_viridis_c() +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./TCR/figures/plots/20230813_TCR_Illumina_ONT_indel_rate.svg", width = 1.5, height = 2.5, plot = p)



ggplot() +
  geom_histogram(
    data = Illumina.data, aes(x = indel, y = 100 * after_stat(count / sum(count))),
    binwidth = 1, fill = "grey90"
  ) +
  geom_histogram(
    data = ONT.data, aes(x = indel, y = 100 * after_stat(count / sum(count))),
    binwidth = 1, fill = "orange", alpha = 0.5
  ) +
  scale_x_continuous("indels per read") +
  scale_y_continuous("% reads") +
  theme_classic() +
  theme(
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./TCR/figures/plots/20230813_TCR3_indel_read.svg", width = 1.5, height = 2)

ggplot() +
  geom_histogram(
    data = Illumina.data, aes(x = 100 * indel.freq, y = 100 * after_stat(count / sum(count))),
    binwidth = 1, fill = "grey90"
  ) +
  geom_histogram(
    data = ONT.data, aes(x = 100 * indel.freq, y = 100 * after_stat(count / sum(count))),
    binwidth = 1, fill = "orange", alpha = 0.5
  ) +
  scale_x_continuous("% indels") +
  scale_y_continuous("% reads") +
  theme_classic() +
  theme(
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./TCR/figures/plots/20230813_TCR3_indel_read_relative.svg", width = 1.5, height = 2)

# correlate sequencing error with length of homopolymers

### package for finding homopolymers
# remotes::install_github("MarioniLab/sarlacc")
library(sarlacc)

# detect homopolymers
ONT.mean.data$homopolymer.G <- sapply(ONT.mean.data$CDR3.nt, FUN = function(x) {
  boo <- homopolymerFinder(DNAStringSet(x))[[1]]
  col <- which(mcols(boo)[, "base"] == "G")
  max(width(boo)[col])
})

ONT.mean.data$homopolymer.C <- sapply(ONT.mean.data$CDR3.nt, FUN = function(x) {
  boo <- homopolymerFinder(DNAStringSet(x))[[1]]
  col <- which(mcols(boo)[, "base"] == "C")
  max(width(boo)[col])
})

ONT.mean.data$homopolymer.A <- sapply(ONT.mean.data$CDR3.nt, FUN = function(x) {
  boo <- homopolymerFinder(DNAStringSet(x))[[1]]
  col <- which(mcols(boo)[, "base"] == "A")
  max(width(boo)[col])
})

ONT.mean.data$homopolymer.T <- sapply(ONT.mean.data$CDR3.nt, FUN = function(x) {
  boo <- homopolymerFinder(DNAStringSet(x))[[1]]
  col <- which(mcols(boo)[, "base"] == "T")
  max(width(boo)[col])
})

ONT.mean.data[ONT.mean.data == "-Inf"] <- 0

ONT.mean.data$homopolymer <- apply(ONT.mean.data, MARGIN = 1, FUN = function(x) {
  max(c(x["homopolymer.G"], x["homopolymer.C"], x["homopolymer.A"], x["homopolymer.T"]))
})



Illumina.mean.data$homopolymer.G <- sapply(Illumina.mean.data$CDR3.nt, FUN = function(x) {
  boo <- homopolymerFinder(DNAStringSet(x))[[1]]
  col <- which(mcols(boo)[, "base"] == "G")
  max(width(boo)[col])
})

Illumina.mean.data$homopolymer.C <- sapply(Illumina.mean.data$CDR3.nt, FUN = function(x) {
  boo <- homopolymerFinder(DNAStringSet(x))[[1]]
  col <- which(mcols(boo)[, "base"] == "C")
  max(width(boo)[col])
})

Illumina.mean.data$homopolymer.A <- sapply(Illumina.mean.data$CDR3.nt, FUN = function(x) {
  boo <- homopolymerFinder(DNAStringSet(x))[[1]]
  col <- which(mcols(boo)[, "base"] == "A")
  max(width(boo)[col])
})

Illumina.mean.data$homopolymer.T <- sapply(Illumina.mean.data$CDR3.nt, FUN = function(x) {
  boo <- homopolymerFinder(DNAStringSet(x))[[1]]
  col <- which(mcols(boo)[, "base"] == "T")
  max(width(boo)[col])
})

Illumina.mean.data[Illumina.mean.data == "-Inf"] <- 0

Illumina.mean.data$homopolymer <- apply(Illumina.mean.data, MARGIN = 1, FUN = function(x) {
  max(c(x["homopolymer.G"], x["homopolymer.C"], x["homopolymer.A"], x["homopolymer.T"]))
})

# plot statistics

ggplot() +
  geom_point(data = ONT.mean.data, aes(x = indel.rank, y = indel, color = as.numeric(homopolymer.G)), size = 0.5) +
  geom_label_repel(
    data = ONT.mean.data[which(ONT.mean.data$homopolymer.G > 9), ],
    aes(x = indel.rank, y = indel, label = stringr::str_split_fixed(clone, pattern = "_", n = 4)[, 1]),
    size = 2, label.size = 0, min.segment.length = 0.1
  ) +
  scale_x_continuous("ONT rank") +
  scale_y_continuous("% indel") +
  scale_color_gradientn("G homopolymer", colours = BuenColors::jdb_palette(name = "brewer_red")) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./TCR/figures/plots/20230813_TCR3_ONT_indel_rank_homopolymer_G.svg", width = 2.5, height = 2.5)

ggplot() +
  geom_point(data = Illumina.mean.data, aes(x = indel.rank, y = indel, color = as.numeric(homopolymer.G)), size = 0.5) +
  geom_label_repel(
    data = Illumina.mean.data[which(Illumina.mean.data$homopolymer.G > 9), ],
    aes(x = indel.rank, y = indel, label = stringr::str_split_fixed(clone, pattern = "_", n = 4)[, 1]),
    size = 2, label.size = 0, min.segment.length = 0.1
  ) +
  scale_x_continuous("Illumina rank") +
  scale_y_continuous("% indel") +
  scale_color_gradientn("G homopolymer", colours = BuenColors::jdb_palette(name = "brewer_red")) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./TCR/figures/plots/20230813_TCR3_Illumina_indel_rank_homopolymer_G.svg", width = 2.5, height = 2.5)



ggplot(ONT.mean.data, aes(x = mismatch.rank, y = mismatch, color = as.numeric(homopolymer.G))) +
  geom_point(size = 0.5) +
  scale_x_continuous("ONT rank") +
  scale_y_continuous("% indel") +
  scale_color_gradientn("G homopolymer", colours = BuenColors::jdb_palette(name = "brewer_red")) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./TCR/figures/plots/20230813_TCR3_ONT_indel_rank_homopolymer_G.svg", width = 2.5, height = 2.5)

ggplot(Illumina.mean.data, aes(x = mismatch.rank, y = mismatch, color = as.numeric(homopolymer.G))) +
  geom_point(size = 0.5) +
  scale_x_continuous("Illumina rank") +
  scale_y_continuous("% indel") +
  scale_color_gradientn("G homopolymer", colours = BuenColors::jdb_palette(name = "brewer_red")) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./TCR/figures/plots/20230813_TCR3_Illumina_indel_rank_homopolymer_G.svg", width = 2.5, height = 2.5)


ggplot() +
  ggrastr::rasterize(geom_point(
    data = Illumina.mean.data,
    aes(x = homopolymer.G, y = indel), size = 0.5, color = "grey"
  ), dpi = 600) +
  ggrastr::rasterize(geom_point(
    data = ONT.mean.data,
    aes(x = homopolymer.G, y = indel), size = 0.5, color = "orange", alpha = 0.5
  ), dpi = 600) +
  geom_smooth(data = Illumina.mean.data, aes(x = homopolymer.G, y = indel), color = "grey", fill = "lightblue") +
  geom_smooth(data = ONT.mean.data, aes(x = homopolymer.G, y = indel), color = "orange", fill = "lightblue") +
  scale_x_continuous("length guanine homopolymer") +
  scale_y_continuous("% indels") +
  theme_classic() +
  theme(
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text = element_text("Arial", size = 10, color = "black")
  )
ggsave("./TCR/figures/plots/20230813_TCR3_indel_rate_length_homopolymer_G.svg", width = 2, height = 2)
