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

sh(library(RColorBrewer))
sh(library(viridis))
sh(library(scales))

# === Functions ===========

scaled_values <- function(values, min_value = 0) {
  # returns values [0, 1] or [min_value, 1]
  
  min_ <- min(values)
  max_ <- max(values)
  scaled_values <- (values - min_) / (max_ - min_)
  scaled_values[scaled_values < min_value] <- min_value
  return(scaled_values)
}

viridis_values <- function(values) {
  the_scaled_values <- scaled_values(values)
  N_colors <- length(unique(the_scaled_values))
  colors <- viridis(N_colors)
  return(
    colors[ceiling(the_scaled_values * (N_colors - 1)) + 1]
  )
}

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

df <- read.csv('NMFHierarchies_build_table_for_ggalluvial2.csv')
# df <- read.csv('NMFHierarchies_build_table_for_ggalluvial_tiny.csv')

# pivot the data to long format,
# there is one row for each voxel * each rank assignment
# 67019 voxels * 7 ranks shown in the plot = 469133 rows
df.long <- df %>%
  mutate(Voxel=1:n()) %>%
  pivot_longer(starts_with("N_Patterns"), names_to = 'N_Patterns', values_to = 'Assignment') %>%
  mutate(Rank = as.numeric(str_extract(N_Patterns, '\\d+')),
         Assignment = factor(Assignment))

ggplot(df.long, aes(x = Rank, stratum = Assignment, alluvium=Voxel,
                    fill = Assignment, label = Assignment)) +
  geom_flow(alpha = 0.5, stat = "flow", aes.flow='forward') + # lode.guidance = "frontback", 
  geom_stratum() +
  xlab("Number of patterns in model space") +
  ylab("Number of voxels") +
  theme_classic(base_size = 24) +
  theme(legend.position = "right") + 
  scale_x_continuous(breaks=c(seq(2,24,by=2))) +
  scale_fill_manual(
    name="argmax clusters", 
    values=hue_values(1:24),
    labels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24")
  )

# See also
# https://github.com/PennLINC/multiscale/blob/b773efebace69f04ffc839ab59f66528b3406714/scripts/analyses/archive/BwRsqCentricOverview.Rmd
# lines 904 - 1023

