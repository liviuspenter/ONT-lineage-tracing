# visualize nanoranger deconcatenation of transcripts

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)

libraries <- unique(stringr::str_split_fixed(list.files("./data/benchmarking/read_names/"), pattern = "\\.read_names", n = 2)[, 1])

stats.df <- data.frame(
  library = as.character(),
  reads = as.numeric(),
  transcripts = as.numeric()
)

reads.df <- data.frame()

for (l in libraries) {
  message(l)
  read.names <- data.table::fread(paste0("./data/benchmarking/read_names/", l, ".read_names.genome.gz"), header = F)
  # read.names = read.names[1:1000,]
  read.statistics <- as.data.frame(read.names %>% group_by(V1) %>% summarize(n = n()) %>% group_by(n) %>% summarize(transcripts = n()))
  read.statistics$transcripts.freq <- read.statistics$transcripts / sum(read.statistics$transcripts)
  read.statistics$library <- l
  reads.df <- rbind(reads.df, read.statistics)
  stats.df <- rbind(stats.df, data.frame(library = l, reads = length(unique(read.names$V1)), transcripts = nrow(read.names)))
}

write.csv2(reads.df, file = "./data/benchmarking/read_names/20230117_transcripts.csv", quote = F)
write.csv2(stats.df, file = "./data/benchmarking/read_names/20230117_reads.csv", quote = F)

reads.df <- data.table::fread("./data/benchmarking/read_names/20230117_transcripts.csv")
stats.df <- data.table::fread("./data/benchmarking/read_names/20230117_reads.csv")

ggplot(reads.df, aes(x = n, y = transcripts, color = library)) +
  geom_line(aes(group = library)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = BuenColors::jdb_palette(name = "corona")) +
  scale_x_log10("extracted transcripts per read", breaks = c(1, 2, 5, 10, 100, 500)) +
  scale_y_log10("sequenced reads", breaks = c(10, 1000, 100000, 5000000), labels = c("10", "1000", "100000", "5000000")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/deconcatenation/20230117_transcripts_per_read.svg", width = 2, height = 2)

boo <- as.data.frame(reads.df %>% group_by(library) %>% summarize(
  read1 = transcripts[which(n == 1)],
  read2 = sum(transcripts[which(n != 1)])
))

boo$read1.freq <- boo$read1 / (boo$read1 + boo$read2)
boo$read2.freq <- boo$read2 / (boo$read1 + boo$read2)

ggplot(boo %>% tidyr::pivot_longer(cols = c("read1.freq", "read2.freq"), names_to = "readtype"), aes(x = readtype, y = 100 * value)) +
  stat_summary(geom = "crossbar", fun = median, fun.min = median, fun.max = median, size = 0.5, width = 0.5) +
  geom_jitter(width = 0.1, size = 0.5, aes(color = library)) +
  scale_color_manual(values = BuenColors::jdb_palette(name = "corona")) +
  scale_x_discrete("extracted transcripts\nper read", labels = c("single", "multiple")) +
  scale_y_continuous("% sequenced reads") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/deconcatenation/20230117_transcripts_per_read_stats.svg", width = 1.5, height = 2)

ggplot(stats.df, aes(x = reads, y = transcripts, color = library)) +
  geom_abline(slope = 1) +
  geom_point() +
  scale_x_continuous("million obtained reads", breaks = c(1000000, 3000000, 5000000, 7000000), labels = c("1", "3", "5", "7"), limits = c(0, 8500000)) +
  scale_y_continuous("million extracted transcripts", breaks = c(1000000, 3000000, 5000000, 7000000), labels = c("1", "3", "5", "7"), limits = c(0, 8500000)) +
  scale_color_manual(values = BuenColors::jdb_palette(name = "corona")) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black")
  )
ggsave("./benchmarking/figures/deconcatenation/20230117_transcripts_per_read_scatter.svg", width = 2, height = 2)

stats.df$ratio <- stats.df$transcripts / stats.df$reads
