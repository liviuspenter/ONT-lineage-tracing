# demonstrate impact of clean-up PCR1 on pseudobulk coverage via depletion of TSO artifacts

library(ggplot2)
library(dplyr)
library(rtracklayer)

# DLI: no PCR1
# AML1019.1: PCR1 performed

coverage.1 = as.data.frame(data.table::fread('./data/benchmarking/bulk_coverage/20220125_coverage_AML_genotyping_no_cleanup'))
coverage.1$sample = 'DLI'
coverage.2 = as.data.frame(data.table::fread('./data/benchmarking/bulk_coverage/20220214_coverage_AML_genotyping'))
coverage.2$sample = 'AML1019.1'
coverage = rbind(coverage.1, coverage.2)

colnames(coverage) = c('chromosome', 'start', 'end', 'gene','position', 'reads', 'sample')
coverage$abs.position = coverage$start + coverage$position - 1
coverage$exon = NA

genes = import('./data/benchmarking/bulk_coverage/AML_genotyping_chromsomes.bed', format="bed")
primers = import('./data/benchmarking/bulk_coverage/AML_genotyping_nested_primers_no_chromosomes2.bed', format = 'bed')
mutations = import('./data/benchmarking/bulk_coverage/AML_mutations_hg38_no_chromosomes.bed', format = 'bed')
exomes = import('./data/benchmarking/bulk_coverage/Twist_Exome_Target_hg38.bed', format = 'bed')

for (gene in c('RUNX1', 'TP53')) {
  exones = ranges(intersect(genes[which(genes$name == gene)], exomes, ignore.strand = T))
  
  index = c()
  for (exon in seq_along(exones)) {
    index = c(index, which(coverage$abs.position %in% seq(start(exones)[exon], end(exones)[exon])))
    coverage$exon[which(coverage$abs.position %in% seq(start(exones)[exon], end(exones)[exon]))] = exon
  }
  
  coverage.sub = coverage[index,]
  coverage.sub = coverage.sub %>% group_by(sample) %>% mutate(x = row_number())
  #coverage.sub$x = seq(1, nrow(coverage.sub))
  exon.df = data.frame(start = as.numeric(), end = as.numeric(), label = as.character(), label.start = as.numeric())
  for (exon in seq_along(exones)) {
    exon.df = rbind(exon.df, data.frame(start = min(coverage.sub$x[which(coverage.sub$exon == exon)]),
                                        end = max(coverage.sub$x[which(coverage.sub$exon == exon)]),
                                        label = as.character(exon),
                                        label.start = min(coverage.sub$abs.position[which(coverage.sub$exon == exon)])))
  }
  exon.df$y = rep(c(1.1*max(coverage.sub$reads.rel), 1.3*max(coverage.sub$reads.rel), 1.5*max(coverage.sub$reads.rel)), length.out = nrow(exon.df))
  exon.df$mid = 0.5*(exon.df$end + exon.df$start)
  
  primer.df = data.frame(start = start(primers[which(grepl(gene, primers$name))]),
                         end = end(primers[which(grepl(gene, primers$name))]))
  
  mutations.df = data.frame(start = start(mutations[which(grepl(gene, mutations$name))]),
                            end = end(mutations[which(grepl(gene, mutations$name))]))
  
  coverage.sub = coverage.sub %>% group_by(sample) %>% mutate(reads.rel = reads / max(reads))
  
  p=ggplot() + geom_line(data=coverage.sub, aes(x=x, y=100*reads.rel, color=sample)) + 
    # primers
    geom_segment(data = primer.df, x = coverage.sub$x[which(coverage.sub$abs.position %in% primer.df$start & coverage.sub$sample == 'DLI')], 
                 y = 100*rep(1.15*max(coverage.sub$reads.rel), nrow(primer.df)), 
                 xend = coverage.sub$x[which(coverage.sub$abs.position %in% primer.df$end & coverage.sub$sample == 'DLI')], 
                 yend = 100*rep(1.15*max(coverage.sub$reads.rel), nrow(primer.df)), 
                 color='blue', size=5) + 
    # mutations
    geom_segment(data = mutations.df, x = coverage.sub$x[which(coverage.sub$abs.position %in% mutations.df$start & coverage.sub$sample == 'DLI')], 
                 y = 100*rep(1.15*max(coverage.sub$reads.rel), nrow(mutations.df)), 
                 xend = coverage.sub$x[which(coverage.sub$abs.position %in% mutations.df$end & coverage.sub$sample == 'DLI')]+1, 
                 yend = 100*rep(1.15*max(coverage.sub$reads.rel), nrow(mutations.df)), 
                 color='red', size=5) +
    scale_x_continuous('position on transcript') +
    scale_y_continuous('% coverage ') + 
    scale_color_manual(values = c('DLI' = 'grey', 'AML1019.1' = 'black')) + 
    ggtitle(gene) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.title = element_text(size=10, color='black'),
          axis.text = element_text(size=10, color='black'),
          plot.title = element_text('Arial', size=10, color='black', hjust=0.5, face = 'bold.italic'))
  ggsave(paste0('./benchmarking/figures/coverage/20220822_coverage_', gene, '.svg'), width = 2, height = 1.5, plot = p)
}
