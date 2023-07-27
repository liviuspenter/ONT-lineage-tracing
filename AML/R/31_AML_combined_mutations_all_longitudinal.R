# somatic mutations before treatment with decitabine/ipilimumab and at response

library(ComplexHeatmap)
library(dplyr)
library(nanoranger.R)
library(ggplot2)
library(Seurat)

AML1002 = readRDS('./data/AML/objects/20220911_AML1002.rds')
AML1007 = readRDS('./data/AML/objects/20220825_AML1007.rds')
AML3003 = readRDS('./data/AML/objects/202208_AML3003_14.rds')
AML3005 = readRDS('./data/AML/objects/20220718_AML3005_13.rds')
AML8007 = readRDS('./data/AML/objects/20220911_AML8007.rds')

AML.longitudinal = data.frame(bc = c(colnames(AML1002), colnames(AML1007), colnames(AML3003), colnames(AML3005), colnames(AML8007)),
                              predicted.celltype = c(AML1002$predicted.celltype, 
                                                     AML1007$predicted.celltype, 
                                                     AML3003$predicted.celltype, 
                                                     AML3005$predicted.celltype,
                                                     AML8007$predicted.celltype),
                              UMAP = rbind(AML1002@reductions$ref.umap@cell.embeddings,
                                           AML1007@reductions$ref.umap@cell.embeddings,
                                           AML3003@reductions$ref.umap@cell.embeddings,
                                           AML3005@reductions$ref.umap@cell.embeddings,
                                           AML8007@reductions$ref.umap@cell.embeddings))

files = list.files('./data/AML/mutations/', pattern = 'AML*')

mutations.df = data.frame()
for (f in files) {
  message(f)
  boo = as.data.frame(read.csv2(file = paste0('./data/AML/mutations/',f), sep = '\t'))
  mutations.df = rbind(mutations.df, boo[,c('bc', 'alt', 'ref', 'mutated', 'vaf')])
}
mutations.df$sample = stringr::str_split_fixed(mutations.df$bc, pattern = '_', n=2)[,1]

mutations.df = mutations.df[which(mutations.df$bc %in% AML.longitudinal$bc),]
mutations.df$predicted.celltype = AML.longitudinal[mutations.df$bc, 'predicted.celltype']
mutations.df.condensed = mutations.df %>% group_by(bc) %>% summarize(mutation = ifelse('mutated' %in% mutated, 'mutated', 'wildtype'),
                                                                     predicted.celltype = predicted.celltype)
mutations.df.condensed$patient = stringr::str_split_fixed(mutations.df.condensed$bc, pattern = '\\.', n=2)[,1]
mutations.df.condensed$sample = stringr::str_split_fixed(mutations.df.condensed$bc, pattern = '_', n=2)[,1]
mutations.df.condensed$timepoint = 'Screening'
mutations.df.condensed$timepoint[which(mutations.df.condensed$sample %in% c('AML1002.3','AML1007.3', 'AML3003.4', 'AML3005.3', 'AML8007.4'))] = 'CR'

### 
boo = mutations.df.condensed %>% group_by(predicted.celltype, sample) %>%
  filter(sample %in% c('AML1002.1', 'AML1002.3','AML1007.1', 'AML1007.3', 'AML3003.1', 'AML3003.4', 'AML3005.1', 'AML3005.3', 'AML8007.1', 'AML8007.4')) %>%
  summarize(mut.count = length(which(mutation == 'mutated')),
            cells = length(predicted.celltype),
            mut.freq = length(which(mutation == 'mutated')) / length(predicted.celltype)) %>% ungroup() %>%
  tidyr::complete(predicted.celltype, sample)

boo$timepoint = 'Screening'
boo$timepoint[which(boo$sample %in% c('AML1002.3','AML1007.3', 'AML3003.4', 'AML3005.3', 'AML8007.4'))] = 'CR'
boo$timepoint = factor(boo$timepoint, levels = c('Screening', 'CR'))
boo$predicted.celltype.timepoint = paste0(boo$predicted.celltype, '.', boo$timepoint)
boo$patient = stringr::str_split_fixed(boo$sample, pattern = '\\.', n=2)[,1]
boo$predicted.celltype.patient = paste0(boo$predicted.celltype, '.', boo$patient)

ggplot(boo[which(boo$predicted.celltype %in% c('HSC', 'LMPP', 'GMP', 'Prog_Mk', 'Prog_RBC')),], 
       aes(x=predicted.celltype.timepoint, y=100*mut.freq, color=predicted.celltype)) + 
  scale_color_manual(values = AML.combined.colors) + 
  scale_x_discrete(limits = c('HSC.Screening', 'HSC.CR', 'LMPP.Screening', 'LMPP.CR', 'GMP.Screening', 'GMP.CR', 
                              'Prog_Mk.Screening', 'Prog_Mk.CR', 'Prog_RBC.Screening', 'Prog_RBC.CR')) + 
  scale_y_continuous('% mutated') +
  geom_line(aes(group=predicted.celltype.patient), color='grey') +
  geom_point() +
  stat_summary(geom = 'crossbar', fun = median, fun.min = 'median', fun.max = 'median', width=0.5, size=0.5, color='black') +
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_text('Arial', size=10, color='black'),
        axis.title = element_text('Arial', size=10, color='black'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
ggsave('./AML/figures/combined/plots/20221006_AML_mutations_longitudinal.svg', width = 3, height = 2)

mutations.df.condensed$UMAP1 = AML.longitudinal[mutations.df.condensed$bc, 'UMAP.refUMAP_1']
mutations.df.condensed$UMAP2 = AML.longitudinal[mutations.df.condensed$bc, 'UMAP.refUMAP_2']

ggplot(mutations.df.condensed[which(mutations.df.condensed$timepoint == 'Screening'),], aes(x=UMAP1, y=UMAP2, color=mutation)) + 
  geom_point() + 
  scale_color_manual(values = c('mutated' = 'red', 'wildtype' = 'black')) +
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
ggsave('./AML/figures/combined/UMAP/20221006_AML_mutations_screening.png', width = 4, height = 4, dpi = 600)

ggplot(mutations.df.condensed[which(mutations.df.condensed$timepoint == 'CR'),], aes(x=UMAP1, y=UMAP2, color=mutation)) + 
  geom_point() + 
  scale_color_manual(values = c('mutated' = 'red', 'wildtype' = 'black')) +
  theme_classic() + 
  theme(legend.position = 'none', 
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank())
ggsave('./AML/figures/combined/UMAP/20221006_AML_mutations_response.png', width = 4, height = 4, dpi = 600)
