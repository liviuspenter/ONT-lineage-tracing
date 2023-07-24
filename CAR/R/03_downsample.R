# downsample experiment of CAR reads and identified cells

library(dplyr)
library(ggplot2)
library(Seurat)

technology.colors = c('Both' = '#283350', 'ONT' = '#f93800', 'Illumina' = '#ffb500')

so = readRDS(file='./data/CAR/objects/20220819_CAR_so.rds')

### downsample
downsample.df = data.frame()
for (downsample in seq(10000, 2000000, 10000)) {
  message(downsample)
  CAR.data.ONT = as.data.frame(data.table::fread('./data/CAR/97_6_CART_CD28.csv.gz') %>% 
                                 filter(gene == 'CARTmod') %>%
                                 sample_n(downsample) %>% 
                                 group_by(bc) %>% 
                                 summarize(CARTmod = length(gene)))
  CAR.data.ONT$bc = paste0(CAR.data.ONT$bc, '-1')
  rownames(CAR.data.ONT) = CAR.data.ONT$bc
  CAR.data.ONT = CAR.data.ONT[intersect(rownames(CAR.data.ONT), colnames(so)),]
  downsample.df = rbind(downsample.df, data.frame(downsample = downsample, cells = length(which(CAR.data.ONT$CARTmod > 4))))
}

# plot downsample experiment
ggplot(downsample.df, aes(x=downsample, y=cells)) + 
  geom_line(color=technology.colors['ONT']) +
  scale_x_continuous('1000x reads ONT', breaks = c(500000, 1000000, 1500000, 2000000), labels = c('500', '1000', '1500', '2000')) + 
  scale_y_continuous('cells') + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./CAR/figures/plots/20220819_ONT_downsample.svg', width = 2.5, height = 2)
