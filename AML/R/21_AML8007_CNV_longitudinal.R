# analysis longitudinal CNV changes in AML8007

library(ComplexHeatmap)
library(dplyr)
library(nanoranger.R)
library(ggplot2)
library(Seurat)

AML8007 = readRDS('./data/AML/objects/20220929_AML8007_all.rds')

cnv.data = read.csv2('./data/AML/numbat/AML8007.1/cnv_calls.csv', sep = '\t')
cnv.data$bc = paste0('AML8007.1_', cnv.data$cell)
cnv.data$p_cnv = as.numeric(cnv.data$p_cnv)
cnv.data = cnv.data %>% select(bc, seg, p_cnv) %>% filter(seg %in% c('1a', '3a', '5c')) %>% tidyr::pivot_wider(names_from = 'seg', values_from = 'p_cnv')
colnames(cnv.data) = c('bc', '1a', '3a', '5b')

cnv.data.2 = read.csv2('./data/AML/numbat/AML8007.2/cnv_calls.csv', sep = '\t')
cnv.data.2$bc = paste0('AML8007.2_', cnv.data.2$cell)
cnv.data.2$p_cnv = as.numeric(cnv.data.2$p_cnv)
cnv.data.2 = cnv.data.2 %>% select(bc, seg, p_cnv) %>% filter(seg %in% c('1a', '3a', '5c')) %>% tidyr::pivot_wider(names_from = 'seg', values_from = 'p_cnv')
colnames(cnv.data.2) = c('bc', '1a', '3a', '5b')

cnv.data.3 = read.csv2('./data/AML/numbat/AML8007.3/cnv_calls.csv', sep = '\t')
cnv.data.3$bc = paste0('AML8007.3_', cnv.data.3$cell)
cnv.data.3$p_cnv = as.numeric(cnv.data.3$p_cnv)
cnv.data.3 = cnv.data.3 %>% select(bc, seg, p_cnv) %>% filter(seg %in% c('1a', '3a', '5b')) %>% tidyr::pivot_wider(names_from = 'seg', values_from = 'p_cnv')

cnv.data.4 = read.csv2('./data/AML/numbat/AML8007.4/cnv_calls.csv', sep = '\t')
cnv.data.4$bc = paste0('AML8007.4_', cnv.data.4$cell)
cnv.data.4$p_cnv = as.numeric(cnv.data.4$p_cnv)
cnv.data.4 = cnv.data.4 %>% select(bc, seg, p_cnv) %>% filter(seg %in% c('1a', '3a', '5b')) %>% tidyr::pivot_wider(names_from = 'seg', values_from = 'p_cnv')

cnv.data = as.data.frame(bind_rows(cnv.data, cnv.data.2, cnv.data.3, cnv.data.4))
cnv.data$amp1 = ifelse(cnv.data$`1a` > 0.8, 'mutated', 'wildtype')
cnv.data$del3 = ifelse(cnv.data$`3a` > 0.8, 'mutated', 'wildtype') 
cnv.data$del5 = ifelse(cnv.data$`5b` > 0.8, 'mutated', 'wildtype')

cnv.data$predicted.celltype = AML8007$predicted.celltype[cnv.data$bc]
cnv.data$sample = stringr::str_split_fixed(cnv.data$bc, pattern = '_', n=2)[,1]
cnv.statistics = cnv.data %>% group_by(sample, predicted.celltype) %>% summarize(amp1.freq = length(which(amp1 == 'mutated')) / length(predicted.celltype),
                                                                                 del3.freq = length(which(del3 == 'mutated')) / length(predicted.celltype),
                                                                                 del5.freq = length(which(del5 == 'mutated')) / length(predicted.celltype),
                                                                                 predicted.celltype = predicted.celltype)
cnv.statistics$predicted.celltype = factor(cnv.statistics$predicted.celltype, levels = names(AML.combined.colors))

ggplot(data=cnv.statistics[which(cnv.statistics$predicted.celltype %in% c('HSC', 'LMPP', 'GMP', 'Prog_RBC')),], 
       aes(x=sample, y=100*amp1.freq, fill=predicted.celltype)) + geom_col(position = 'dodge') +
  scale_x_discrete(labels = c('Screening', 'Lead-in', 'C1', 'marrow CR')) +
  scale_y_continuous('% amp(1p)') + 
  scale_fill_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML8007/plots/20220930_8007_amp1.svg', width = 1.5, height = 2)


ggplot(data=cnv.statistics[which(cnv.statistics$predicted.celltype %in% c('HSC', 'LMPP', 'GMP', 'Prog_Mk','Prog_RBC')),], 
       aes(x=sample, y=100*amp1.freq, color=predicted.celltype)) + 
  geom_line(aes(group=predicted.celltype)) +
  scale_x_discrete(labels = c('Screening', 'Lead-in', 'C1', 'marrow CR')) +
  scale_y_continuous('% amp(1p)', limits = c(0,100)) + 
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML8007/plots/20220930_8007_amp1.svg', width = 1.2, height = 2)


ggplot(data=cnv.statistics[which(cnv.statistics$predicted.celltype %in% c('HSC', 'LMPP', 'GMP', 'Prog_RBC')),], 
       aes(x=sample, y=100*del3.freq, fill=predicted.celltype)) + geom_col(position = 'dodge') +
  scale_x_discrete(labels = c('Screening', 'Lead-in', 'C1', 'marrow CR')) +
  scale_y_continuous('% del(3p)') + 
  scale_fill_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML8007/plots/20220930_8007_del3.svg', width = 1.5, height = 2)


ggplot(data=cnv.statistics[which(cnv.statistics$predicted.celltype %in% c('HSC', 'LMPP', 'GMP', 'Prog_Mk','Prog_RBC')),], 
       aes(x=sample, y=100*del3.freq, color=predicted.celltype)) + 
  geom_line(aes(group=predicted.celltype)) +
  scale_x_discrete(labels = c('Screening', 'Lead-in', 'C1', 'marrow CR')) +
  scale_y_continuous('% del(3p)', limits = c(0,100)) + 
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML8007/plots/20220930_8007_del3.svg', width = 1.2, height = 2)


ggplot(data=cnv.statistics[which(cnv.statistics$predicted.celltype %in% c('HSC', 'LMPP', 'GMP', 'Prog_RBC')),], 
       aes(x=sample, y=100*del5.freq, fill=predicted.celltype)) + geom_col(position = 'dodge') +
  scale_x_discrete(labels = c('Screening', 'Lead-in', 'C1', 'marrow CR')) +
  scale_y_continuous('% del(5q)') + 
  scale_fill_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML8007/plots/20220930_8007_del5.svg', width = 1.5, height = 2)



ggplot(data=cnv.statistics[which(cnv.statistics$predicted.celltype %in% c('HSC', 'LMPP', 'GMP', 'Prog_Mk','Prog_RBC')),], 
       aes(x=sample, y=100*del5.freq, color=predicted.celltype)) + 
  geom_line(aes(group=predicted.celltype)) +
  scale_x_discrete(labels = c('Screening', 'Lead-in', 'C1', 'marrow CR')) +
  scale_y_continuous('% del(5q)', limits = c(0,100)) + 
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML8007/plots/20220930_8007_del5.svg', width = 1.2, height = 2)


Idents(AML8007) = 'predicted.celltype'
boo = as.data.frame(prop.table(table(Idents(AML8007), AML8007$orig.ident), margin = 2))

ggplot(boo[which(boo$Var1 %in% c('HSC', 'LMPP', 'GMP', 'Prog_Mk','Prog_RBC')),], 
       aes(x=Var2, y=100*Freq, color=Var1)) + 
  geom_point(size=0.5) + 
  geom_line(aes(group=Var1)) +
  scale_x_discrete(labels = c('Screening', 'Lead-in', 'C1', 'marrow CR')) +
  scale_y_continuous('% cells') + 
  scale_color_manual(values = AML.combined.colors) +
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/AML8007/plots/20220930_8007_celltype_kinetics.svg', width = 1.2, height = 2)

