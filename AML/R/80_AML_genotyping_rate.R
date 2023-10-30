# genotyping rate across all AML cases

files <- list.files("./data/AML/mutations/", pattern = "AML*|FLT3*")

AML.combined <- readRDS("./data/AML/objects/AML.all.rds")
FLT3.so <- readRDS("./data/AML/objects/20230523_FLT3-ITD.rds")

# gather mutations
mutations.df <- data.frame()
for (f in files) {
  boo <- as.data.frame(read.csv2(file = paste0("./data/AML/mutations/", f), sep = "\t"))
  boo$mutation <- stringr::str_split_fixed(f, pattern = "_", n = 2)[, 2]
  mutations.df <- rbind(mutations.df, boo[, c("bc", "alt", "ref", "mutated", "vaf", "mutation")])
}
mutations.df$sample <- stringr::str_split_fixed(mutations.df$bc, pattern = "_", n = 2)[, 1]
mutations.df$mutation <- stringr::str_split_fixed(mutations.df$mutation, pattern = "\\.", n = 2)[, 1]

mutations.df$bc <- gsub(mutations.df$bc, pattern = "AML1002.3", replacement = "AML1002\\.4")

stats <-
  mutations.df %>%
  filter(bc %in% c(colnames(AML.combined), colnames(FLT3.so))) %>%
  group_by(sample, mutation) %>%
  tally()

ggplot(stats, aes(x = reorder(mutation, n), y = n, group = sample)) +
  geom_jitter(width = 0.1, size = 0.5) +
  scale_y_sqrt("genotyped cells",
    breaks = c(10, 100, 500, 1000, 2500, 5000, 7500),
    labels = c("10", "100", "500", "1000", "2500", "5000", "7500"), limits = c(1, 6000)
  ) +
  theme_classic() +
  theme(
    axis.text = element_text("Arial", size = 10, color = "black"),
    axis.title = element_text("Arial", size = 10, color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )
ggsave("./AML/figures/combined/plots/20230807_genotyping_rate.svg", width = 2.5, height = 2.5)
