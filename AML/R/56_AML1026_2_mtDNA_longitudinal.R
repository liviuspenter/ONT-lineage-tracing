# track mtDNA mutation longitudinally in AML1026

library(dplyr)
library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(nanoranger.R)
library(Seurat)
library(SummarizedExperiment)

# functionality from MAESTER (https://github.com/petervangalen/MAESTER-2021)
source('./R/210215_FunctionsGeneral.R')

AML.all = readRDS('./data/AML/objects/AML.all.rds')

maegtk.output.1 = readRDS('./data/AML/mtDNA/AML1026.1.maegatk/Pool96_15.rds')
af.dm.1 = data.matrix(computeAFMutMatrix(maegtk.output.1))*100
rownames(af.dm.1) = gsub(rownames(af.dm.1), pattern = '_', replacement = '')
colnames(af.dm.1) = paste0('AML1026.1_', colnames(af.dm.1))
maegtk.output.2 = readRDS('./data/AML/mtDNA/AML1026.2.maegatk/Pool96_16.rds')
af.dm.2 = data.matrix(computeAFMutMatrix(maegtk.output.2))*100
rownames(af.dm.2) = gsub(rownames(af.dm.2), pattern = '_', replacement = '')
colnames(af.dm.2) = paste0('AML1026.2_', colnames(af.dm.2))
maegtk.output.3 = readRDS('./data/AML/mtDNA/AML1026.3.maegatk/Pool96_17.rds')
af.dm.3 = data.matrix(computeAFMutMatrix(maegtk.output.3))*100
rownames(af.dm.3) = gsub(rownames(af.dm.3), pattern = '_', replacement = '')
colnames(af.dm.3) = paste0('AML1026.3_', colnames(af.dm.3))
maegtk.output.4 = readRDS('./data/AML/mtDNA/AML1026.4.maegatk/Pool96_18.rds')
af.dm.4 = data.matrix(computeAFMutMatrix(maegtk.output.4))*100
rownames(af.dm.4) = gsub(rownames(af.dm.4), pattern = '_', replacement = '')
colnames(af.dm.4) = paste0('AML1026.4_', colnames(af.dm.4))

df = as.data.frame(t(bind_cols(af.dm.1[c('10685G>A','15615G>A'),intersect(colnames(AML.all)[which(AML.all$orig.ident %in% c('AML1026.1'))], colnames(af.dm.1))],
           af.dm.2[c('10685G>A','15615G>A'),intersect(colnames(AML.all)[which(AML.all$orig.ident %in% c('AML1026.2'))], colnames(af.dm.2))],
           af.dm.3[c('10685G>A','15615G>A'),intersect(colnames(AML.all)[which(AML.all$orig.ident %in% c('AML1026.3'))], colnames(af.dm.3))],
           af.dm.4[c('10685G>A','15615G>A'),intersect(colnames(AML.all)[which(AML.all$orig.ident %in% c('AML1026.4'))], colnames(af.dm.4))])))
colnames(df) = c('10685G>A','15615G>A')
df$predicted.celltype = AML.all$predicted.celltype[rownames(df)]
df$sample = stringr::str_split_fixed(rownames(df), pattern = '_', n=2)[,1]

boo = df %>% group_by(sample, predicted.celltype) %>% summarize(mt.10685 = length(which(`10685G>A` != 0)),
                                                                mt.10685.freq = length(which(`10685G>A` != 0)) / length(predicted.celltype),
                                                                mt.15615 = length(which(`15615G>A` != 0)),
                                                                mt.15615.freq = length(which(`15615G>A` != 0)) / length(predicted.celltype))

Idents(AML.all) = 'predicted.celltype'
AML.all.df = as.data.frame(prop.table(table(Idents(AML.all), AML.all$orig.ident), margin = 2))

ggplot(data=AML.all.df[which(grepl('AML1026', AML.all.df$Var2) & AML.all.df$Var1 %in% c('HSC', 'LMPP', 'GMP')),], 
       aes(x=Var2, y=100*Freq)) + 
  geom_line(aes(group=Var1, color=Var1)) + 
  scale_x_discrete(labels = c('Screening', 'EOLN', 'C1', 'Relapse')) + 
  scale_y_continuous('% cells') + 
  scale_color_manual(values = AML.combined.colors) + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/mtDNA/plots/20220805_AML1026_cells_kinetics.svg', width = 1.2, height = 1.5)

boo = boo[which(boo$predicted.celltype %in% c('HSC', 'LMPP', 'GMP', 'Prog_Mk', 'Prog_RBC')),]
boo$predicted.celltype = factor(boo$predicted.celltype, levels = c('HSC', 'LMPP', 'GMP', 'Prog_Mk', 'Prog_RBC'))
ggplot() +
  geom_bar(data=boo[which(boo$predicted.celltype %in% c('HSC')),], 
           aes(x=sample, y=100*mt.10685.freq, fill=predicted.celltype), stat = 'identity', position = 'identity') + 
  geom_bar(data=boo[which(boo$predicted.celltype %in% c('LMPP')),], 
           aes(x=sample, y=100*mt.10685.freq, fill=predicted.celltype), stat = 'identity', position = 'identity') + 
  geom_bar(data=boo[which(boo$predicted.celltype %in% c('GMP')),], 
           aes(x=sample, y=100*mt.10685.freq, fill=predicted.celltype), stat = 'identity', position = 'identity') + 
  #geom_bar(data=boo[which(boo$predicted.celltype %in% c('Prog_Mk')),], 
  #         aes(x=sample, y=100*mt.10685.freq, fill=predicted.celltype), stat = 'identity', position = 'identity') + 
  #geom_bar(data=boo[which(boo$predicted.celltype %in% c('Prog_RBC')),], 
  #         aes(x=sample, y=100*mt.10685.freq, fill=predicted.celltype), stat = 'identity', position = 'identity') + 
  scale_fill_manual(values = AML.combined.colors, limits =  c('GMP', 'HSC','LMPP', 'Prog_Mk', 'Prog_RBC')) + 
  scale_x_discrete(labels = c('Screening', 'EOLN', 'C1', 'Relapse')) + 
  scale_y_continuous('% 10685G>A') + 
  theme_classic() + 
  theme(legend.position = 'none',
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
ggsave('./AML/figures/mtDNA/plots/20220805_AML1026_10685G>A_kinetics.svg', width = 1.2, height = 1.5)
