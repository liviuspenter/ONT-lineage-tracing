# analyse CD247 amplicon

library(ComplexHeatmap)
library(ggplot2)
library(Seurat)

technology.colors = c('Both' = '#283350', 'ONT' = '#f93800', 'Illumina' = '#ffb500', 'ONT2' = 'firebrick')

so = readRDS('./data/CAR/objects/20220819_CAR_so.rds')

CAR.data.ONT = as.data.frame(data.table::fread('./data/CAR/97_6_CART_CD28.csv.gz') %>% filter(gene == 'CARTmod') %>% group_by(bc) %>% summarize(CARTmod = length(gene)))
CAR.data.ONT$bc = paste0(CAR.data.ONT$bc, '-1')
rownames(CAR.data.ONT) = CAR.data.ONT$bc
CAR.data.ONT = CAR.data.ONT[intersect(rownames(CAR.data.ONT), colnames(so)),]

CAR.data.ONT.2 = as.data.frame(data.table::fread('./data/CAR/97_6_CART_CD247.csv.gz') %>% filter(gene == 'CARTmod') %>% group_by(bc) %>% summarize(CARTmod = length(gene)))
CAR.data.ONT.2$bc = paste0(CAR.data.ONT.2$bc, '-1')
rownames(CAR.data.ONT.2) = CAR.data.ONT.2$bc
CAR.data.ONT.2 = CAR.data.ONT.2[intersect(rownames(CAR.data.ONT.2), colnames(so)),]

CAR.data.Illumina = as.data.frame(data.table::fread('./data/CAR/97_6_CART_illumina.csv'))
colnames(CAR.data.Illumina)[1] = 'bc'
rownames(CAR.data.Illumina) = CAR.data.Illumina$bc
CAR.data.Illumina = CAR.data.Illumina[intersect(CAR.data.Illumina$bc, colnames(so)),]

CAR.data.combined = merge(CAR.data.ONT[,c('bc', 'CARTmod')], CAR.data.Illumina[,c('bc', 'CART_count')], by = 'bc', all=T)
CAR.data.combined = merge(CAR.data.combined, CAR.data.ONT.2[,c('bc', 'CARTmod')], by = 'bc', all=T)
CAR.data.combined[is.na(CAR.data.combined)] = 0
colnames(CAR.data.combined) = c('bc', 'CARTmod.CD28', 'CART_count', 'CARTmod.CD247')
CAR.data.combined$CARTmod.CD28.plot = CAR.data.combined$CARTmod.CD28
CAR.data.combined$CARTmod.CD247.plot = CAR.data.combined$CARTmod.CD247
CAR.data.combined$CART_count.plot = CAR.data.combined$CART_count
CAR.data.combined$CARTmod.CD28.plot[which(CAR.data.combined$CARTmod.CD28 == 0)] = jitter(rep(0.1, length(which(CAR.data.combined$CARTmod.CD28 == 0))), amount = 0.01)
CAR.data.combined$CARTmod.CD247.plot[which(CAR.data.combined$CARTmod.CD247 == 0)] = jitter(rep(0.1, length(which(CAR.data.combined$CARTmod.CD247 == 0))), amount = 0.01)
CAR.data.combined$CART_count.plot[which(CAR.data.combined$CART_count == 0)] = jitter(rep(0.1, length(which(CAR.data.combined$CART_count == 0))), amount = 0.01)

CAR.data.combined = CAR.data.combined %>% arrange(-CARTmod.CD28) %>% mutate(ONT.rank.CD28 = row_number())
CAR.data.combined = CAR.data.combined %>% arrange(-CARTmod.CD247) %>% mutate(ONT.rank.CD247 = row_number())
CAR.data.combined = CAR.data.combined %>% arrange(-CART_count) %>% mutate(Illumina.rank = row_number())

ggplot() + 
  geom_line(data=CAR.data.combined[which(CAR.data.combined$CARTmod.CD28 != 0),], aes(x=ONT.rank.CD28, y=CARTmod.CD28), color=technology.colors['ONT']) + 
  geom_line(data=CAR.data.combined[which(CAR.data.combined$CARTmod.CD247 != 0),], aes(x=ONT.rank.CD247, y=CARTmod.CD247), color=technology.colors['ONT2']) + 
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
ggsave('./CAR/figures/plots/20221007_knee_plot.svg', width = 2, height = 2)

CAR.data.combined$technology.detected = 'Both'
CAR.data.combined$technology.detected[which(CAR.data.combined$CARTmod.CD247 != 0 & CAR.data.combined$CART_count == 0)] = 'ONT'
CAR.data.combined$technology.detected[which(CAR.data.combined$CARTmod.CD247 == 0 & CAR.data.combined$CART_count != 0)] = 'Illumina'

ggplot() + 
  ggrastr::rasterize(geom_point(data=CAR.data.combined[which(CAR.data.combined$CARTmod.CD247 != 0 | CAR.data.combined$CART_count != 0 ),], 
                                aes(x=CARTmod.CD247.plot, y=CART_count.plot, color=technology.detected), size=0.5), dpi=600) + 
  scale_x_log10('Reads ONT', breaks = c(0.1, 1, 10, 100), labels = c('ND', '1', '10', '100')) +
  scale_y_log10('Reads Illumina', breaks = c(0.1, 1, 10, 100), labels = c('ND', '1', '10', '100')) + 
  scale_color_manual(values = technology.colors) + 
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'))
ggsave('./CAR/figures/plots/20221007_reads_ONT_Illumina_CD247.svg', width = 2, height = 2)

lt = list(CD28 = CAR.data.combined$bc[which(CAR.data.combined$CARTmod.CD28 != 0)],
          CD247 = CAR.data.combined$bc[which(CAR.data.combined$CARTmod.CD247 != 0)],
          Illumina = CAR.data.combined$bc[which(CAR.data.combined$CART_count != 0)])
m = make_comb_mat(lt)
svglite::svglite('./CAR/figures/plots/20221007_upset_plot.svg', width = 4, height = 3)
UpSet(m, comb_order = order(-comb_size(m)), 
      pt_size = unit(2, "mm"), lwd = 1, 
      comb_col = c(technology.colors['ONT'], technology.colors['ONT'], technology.colors['ONT'], 
                                                         technology.colors['ONT2'], technology.colors['ONT'], technology.colors['ONT2'],
                                                         technology.colors['Illumina']), 
      bg_col = c("grey90", "white"), bg_pt_col = "grey", top_annotation = upset_top_annotation(m, add_numbers=T, angle=90))
dev.off()