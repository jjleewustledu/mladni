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

# sh(library(alluvial))
# sh(library(ggfittext))
sh(library(ggrepel))

# === Functions ===========

hue_values <- function(values) {
  the_scaled_values <- scaled_values(values)
  N_colors <- length(unique(the_scaled_values))
  hue_color_palette <- hue_pal()
  colors <- hue_color_palette(N_colors)
  return(
    colors[ceiling(the_scaled_values * (N_colors - 1)) + 1]
  )
}

# === Set WD ===========

setwd("/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG/baseline_cn")

# ==== Read data =========

# df <- read.csv('NMFHierarchies_build_table_for_ggalluvial.csv')
df <- read.csv('NMFHierarchies_build_table_for_ggalluvial_tiny.csv')

# pivot the data to long format,
# there is one row for each voxel * each rank assignment
# 67019 voxels * 6 ranks shown in the plot = 4012114 rows
df.long <- df %>%
  mutate(Voxel=1:n()) %>%
  pivot_longer(starts_with("N_Patterns"), names_to = 'N_Patterns', values_to = 'Assignment') %>%
  mutate(Rank = as.numeric(str_extract(N_Patterns, '\\d+')),
         Assignment = factor(Assignment))

ggplot(df.long, aes(x = Rank, stratum = Assignment, alluvium=Voxel,
                    fill = Assignment, label = Assignment)) +
  geom_flow(stat = "flow", lode.guidance = "frontback", aes.flow='forward') +
  geom_stratum(aes()) +
  xlab("Number of Patterns in the Model Space") +
  ylab("Number of Voxels") +
  theme_light() +
  theme(legend.position = "bottom")

