# visualize impact of nanoranger deconcatenation on read length

library(ComplexHeatmap)
library(dplyr)
library(ggplot2)

libraries = unique(stringr::str_split_fixed(list.files('./data/benchmarking/read_length/'), pattern = '\\.read_lengths', n=2)[,1])

stats.df = data.frame(library = as.character(),
                      deconcat = as.numeric(),
                      original = as.numeric())

reads.df = data.frame()

for (l in libraries) {
  message(l)
  deconcat = data.table::fread(paste0('./data/benchmarking/read_length/', l, '.read_lengths.deconcat.gz'))
  original = data.table::fread(paste0('./data/benchmarking/read_length/', l, '.read_lengths.original.gz'))
  genome = data.table::fread(paste0('./data/benchmarking/read_length/', l, '.read_lengths.genome.gz'))
  
  stats.df = rbind(stats.df, data.frame(library = l, 
                                        deconcat = nrow(deconcat),
                                        original = nrow(original),
                                        genome = nrow(genome)))
  
  read.stats = function(df, condition=NA) {
    df.reads = as.data.frame(df %>% mutate(bin = cut(V1, breaks = seq(0,3000,100), labels=F)))
    df.reads$bin = factor(df.reads$bin, levels = seq(1,30))
    df.reads = as.data.frame(df.reads %>% group_by(bin, .drop=F) %>% summarize(reads = length(bin)) %>% t())
    colnames(df.reads) = df.reads[1,]
    df.reads = df.reads[-1,]
    df.reads[1,] = sapply(as.numeric(df.reads[1,]), FUN = function(x) {x / sum(as.numeric(df.reads[1,]))}) 
    df.reads = df.reads[,seq(1,30)] # ignore longer reads
    df.reads$library = l
    df.reads$condition = condition
    
    df.reads
  }

  deconcat.reads = read.stats(deconcat, 'deconcat')
  original.reads = read.stats(original, 'original')
  genome.reads = read.stats(original, 'genome')

  reads.df = rbind(reads.df, deconcat.reads)
  reads.df = rbind(reads.df, original.reads)
  reads.df = rbind(reads.df, genome.reads)
}
reads.df.orig = reads.df

reads.df = reads.df[which(reads.df$condition %in% c('original', 'deconcat')),]
reads.df[,seq(1,30)] = apply(reads.df[,seq(1,30)], 2, function(x) as.numeric(x))
reads.df$library.name = paste0('', rep(seq(1,nrow(reads.df)/2), each=2))

col_fun = circlize::colorRamp2(breaks = seq(0,0.1,0.1/8), colors = BuenColors::jdb_palette(name = 'solar_rojos'))
svglite::svglite('./benchmarking/figures/coverage/20230116_coverage_original_deconcat.svg', width = 3, height = 2)
Heatmap(reads.df[,seq(1,30)], cluster_rows = F, cluster_columns = F, col = col_fun, border=T,
        row_split = factor(reads.df$condition, levels = c('original', 'deconcat')),
        row_labels = reads.df$library.name, row_names_side = 'left',
        column_labels = as.character(seq(100,3000,100)),
        row_names_gp = gpar(fontsize=8), column_names_gp = gpar(fontsize=8))
dev.off()

reads.df = reads.df[which(reads.df$condition %in% c('original', 'genome')),]
reads.df[,seq(1,30)] = apply(reads.df[,seq(1,30)], 2, function(x) as.numeric(x))
reads.df$library.name = paste0('', rep(seq(1,nrow(reads.df)/2), each=2))

col_fun = circlize::colorRamp2(breaks = seq(0,0.1,0.1/8), colors = BuenColors::jdb_palette(name = 'solar_rojos'))
svglite::svglite('./benchmarking/figures/coverage/20230116_coverage_original_genome.svg', width = 3, height = 2)
Heatmap(reads.df[,seq(1,30)], cluster_rows = F, cluster_columns = F, col = col_fun, border=T,
        row_split = factor(reads.df$condition, levels = c('original', 'genome')),
        row_labels = reads.df$library.name, row_names_side = 'left',
        column_labels = as.character(seq(100,3000,100)),
        row_names_gp = gpar(fontsize=8), column_names_gp = gpar(fontsize=8))
dev.off()