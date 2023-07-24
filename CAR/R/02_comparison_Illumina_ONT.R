# comparison of Illumina and ONT CAR data

library(dplyr)
library(ggplot2)
library(Seurat)

technology.colors = c('Both' = '#283350', 'ONT' = '#f93800', 'Illumina' = '#ffb500')

so = readRDS(file='./data/CAR/objects/20220819_CAR_so.rds')

# read CAR ONT data
CAR.data.ONT = as.data.frame(data.table::fread('./data/CAR/97_6_CART_CD28.csv.gz') %>% 
                               filter(gene == 'CARTmod') %>% 
                               group_by(bc) %>% 
                               summarize(CARTmod = length(gene)))
CAR.data.ONT$bc = paste0(CAR.data.ONT$bc, '-1')
rownames(CAR.data.ONT) = CAR.data.ONT$bc
CAR.data.ONT = CAR.data.ONT[intersect(rownames(CAR.data.ONT), colnames(so)),]

# read Illumina data
CAR.data.Illumina = as.data.frame(data.table::fread('./data/CAR/97_6_CART_illumina.csv'))
colnames(CAR.data.Illumina)[1] = 'bc'
rownames(CAR.data.Illumina) = CAR.data.Illumina$bc
CAR.data.Illumina = CAR.data.Illumina[intersect(CAR.data.Illumina$bc, colnames(so)),]

# join data
CAR.data.combined = merge(CAR.data.ONT[,c('bc', 'CARTmod')], CAR.data.Illumina[,c('bc', 'CART_count')], by = 'bc', all=T)
CAR.data.combined[is.na(CAR.data.combined)] = 0
CAR.data.combined$CARTmod.plot = CAR.data.combined$CARTmod
CAR.data.combined$CART_count.plot = CAR.data.combined$CART_count
CAR.data.combined$CARTmod.plot[which(CAR.data.combined$CARTmod == 0)] = jitter(rep(0.1, length(which(CAR.data.combined$CARTmod == 0))), amount = 0.01)
CAR.data.combined$CART_count.plot[which(CAR.data.combined$CART_count == 0)] = jitter(rep(0.1, length(which(CAR.data.combined$CART_count == 0))), amount = 0.01)
CAR.data.combined = CAR.data.combined %>% arrange(-CARTmod) %>% mutate(ONT.rank = row_number())
CAR.data.combined = CAR.data.combined %>% arrange(-CART_count) %>% mutate(Illumina.rank = row_number())

# plot ONT versus Illumina knee plots
ggplot() + 
  geom_line(data=CAR.data.combined[which(CAR.data.combined$CARTmod != 0),], aes(x=ONT.rank, y=CARTmod), color=technology.colors['ONT']) + 
  geom_line(data=CAR.data.combined[which(CAR.data.combined$CART_count != 0),], aes(x=Illumina.rank, y=CART_count), color=technology.colors['Illumina']) + 
  #geom_line(data=CAR.data.combined[which(CAR.data.combined$CD28 != 0),], aes(x=CD28.rank, y=CD28), color='black') + 
  #geom_point(data=CAR.data.combined[which(CAR.data.combined$CD28.only == 'yes'),], aes(x=CD28.rank, y=CD28), color='green') + 
  scale_x_log10('cell rank') + 
  scale_y_log10('reads') + 
  geom_hline(yintercept = 5) +
  theme_classic() +
  theme(legend.position = 'none', 
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./CAR/figures/plots/20220819_knee_plot_ONT_Illumina.svg', width = 2, height = 2)

CAR.data.combined$technology.detected = 'Both'
CAR.data.combined$technology.detected[which(CAR.data.combined$CARTmod != 0 & CAR.data.combined$CART_count == 0)] = 'ONT'
CAR.data.combined$technology.detected[which(CAR.data.combined$CARTmod == 0 & CAR.data.combined$CART_count != 0)] = 'Illumina'

# plot ONT versus Illumina reads per cell barcode
ggplot() + 
  ggrastr::rasterize(geom_point(data=CAR.data.combined[which(CAR.data.combined$CARTmod != 0 | CAR.data.combined$CART_count != 0 ),], 
                                aes(x=CARTmod.plot, y=CART_count.plot, color=technology.detected), size=0.5), dpi=600) + 
  scale_x_log10('Reads ONT', breaks = c(0.1, 1, 10, 100), labels = c('ND', '1', '10', '100')) +
  scale_y_log10('Reads Illumina', breaks = c(0.1, 1, 10, 100), labels = c('ND', '1', '10', '100')) + 
  scale_color_manual(values = technology.colors) + 
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./CAR/figures/plots/20220819_reads_ONT_Illumina.svg', width = 2, height = 2)

length(which(CAR.data.combined$CARTmod == 0 & CAR.data.combined$CART_count != 0))
length(which(CAR.data.combined$CARTmod != 0 & CAR.data.combined$CART_count == 0))
length(which(CAR.data.combined$CARTmod == 0 & CAR.data.combined$CART_count == 0))
length(which(CAR.data.combined$CARTmod != 0 & CAR.data.combined$CART_count != 0))

# plot overlapping and exclusive cellsbetween ONT and Illumina
ggplot(data.frame(condition = c('Both', 'ONT', 'Illumina'), 
                  n = c(length(which(CAR.data.combined$CARTmod > 4 & CAR.data.combined$CART_count > 4)),
                        length(which(CAR.data.combined$CARTmod > 4 & CAR.data.combined$CART_count == 0)),
                        length(which(CAR.data.combined$CARTmod == 0 & CAR.data.combined$CART_count > 4)))), 
       aes(x=condition, y=n, fill=condition)) +
  geom_col() +
  scale_fill_manual(values = technology.colors) + 
  scale_x_discrete('technology') + 
  scale_y_continuous('cells with CAR detected') + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
ggsave('./CAR/figures/plots/20220819_cells_ONT_Illumina.svg', width = 1.3, height = 2)
