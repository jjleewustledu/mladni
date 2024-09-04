# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(anticlust))
sh(library(colormap))
sh(library(ggalluvial))
sh(library(ggplot2))
sh(library(ggseg))
sh(library(gtools))
sh(library(R.matlab))
sh(library(stringr))
sh(library(this.path))
sh(library(tidyverse))

sh(library(alluvial))
sh(library(ggfittext))
sh(library(ggrepel))

# === Set WD ===========

setwd("/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG/baseline_cn")

# ==== Read data =========

T <- read.csv('NMFHierarchies_build_table_for_ggalluvial.csv')
head(as.data.frame(T), n = 12)
is_alluvia_form(as.data.frame(T), axes = 1:6, silent = TRUE)

# ==== Lodes (long) format =========

T_lodes <- to_lodes_form(as.data.frame(T),
                         axes = 1:6,
                         id = "alluvium")
head(T_lodes, n = 50)
is_lodes_form(T_lodes, key = x, value = stratum, id = alluvium, silent = TRUE)

ggplot(data = T_lodes,
       aes(x = x,
           stratum = stratum,
           alluvium = alluvium,
           fill = Tag,
           label = after_stat(stratum))) +
  xlab("Pattern Spaces") +
  geom_flow(stat = "alluvium", lode.guidance = "frontback") +
  geom_stratum() +
  theme(legend.position = "bottom") +
  ggtitle("NMF Pattern Hierarchies in Voxels")

# ==== Alluvial (wide) format =========

T_wide <- data.frame(T)
head(T_wide)

ggplot(data = T_wide,
       aes(axis1 = Span.2, axis2 = Span.8, axis3 = Span.10, axis4 = Span.12, axis5 = Span.14, axis6 = Span.24,
           y = freq)) +
  scale_x_discrete(limits = c("2 Patterns", "8 Patterns", "10 Patterns", "12 Patterns", "14 Patterns", "24 Patterns"), expand = c(.05, .05)) +
  xlab("Pattern Spaces") +
  geom_alluvium(aes(fill = Tag)) +
  geom_stratum(alpha = .25, width = 1/8) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  ggtitle("NMF Pattern Hierarchies in Voxels",
          "stratified by spanning spaces and patterns of interest")

# ==== More suggestions from Tom Earnest ===========
#
#
# assignments <- as.data.frame(assignments.mat)
# rank.cols <- as.character(ranks)
# colnames(assignments) <- rank.cols
# assignments$Region <- regions$label
# 
# 
# ==== Create plot data ===========
# 
# plot.alluvial <- function(ranks) {
#   plot.data <- assignments %>%
#     pivot_longer(all_of(rank.cols), names_to = 'Rank', values_to = 'Component') %>%
#     mutate(Rank = as.numeric(Rank),
#            Component = factor(Component)) %>%
#     filter(Rank %in% ranks)
#   
#   
#   ggplot(plot.data, aes(x = Rank, stratum = Component, alluvium = Region,
#                         fill = Component, label = Component)) +
#     geom_flow(stat = "alluvium", lode.guidance = "frontback", aes.flow='forward') +
#     geom_stratum() +
#     theme_light()
# }
# 
# plot.alluvial(c(2, 6, 8))
# ggsave('alluvial_2_6_8.png', width=8, height=4)
# 
# plot.alluvial(c(2, 6, 8, 12, 16, 20))
# ggsave('alluvial_2_6_8_12_16_20.png', width=8, height=4)
