
# === Imports ==========

sh <- suppressPackageStartupMessages

sh(library(ADNIMERGE))
sh(library(anticlust))
sh(library(dplyr))
sh(library(lubridate))
sh(library(tidyr))
sh(library(this.path))

# === Set WD ==========

setwd(this.dir())

# === Load data for NMF ==========

PATH <- 'main_data.csv'
df <- read.csv(PATH) %>%
  mutate(DateTau = as_datetime(ymd(DateTau))) %>%
  filter(Group == 'TrainingBaseline') 

# === Load helper func ==========

source("../../../rscripts/create_nmf_input.R")

# === prep NMF ==========

features = colnames(df)[grepl('CTX_.*_SUVR', colnames(df), perl=T)]
outdata = 'adni_tau_nmf_matrix.csv'
outfeatures = 'adni_tau_nmf_regions.csv'
outids = 'adni_tau_nmf_ids.csv'

X <- create.nmf.input(df, features=features, id.cols = c("RID", "DateTau"),
                      outdata = outdata, outfeatures = outfeatures,
                      outids = outids)

# === Create reproducibility splits - single =============

# This is the original method

# continuous.vars <- df[, c("Age", "CorticalTauAverage")]
# categorical.vars <- df[, c("Gender", "CDRBinned")]
# v <- anticlustering(continuous.vars, 2, categories=categorical.vars)
# df$Split <- v
# 
# # verify that it worked
# g <- group_by(df, Split) %>%
#   summarise(Age=mean(Age),
#             Tau=mean(CorticalTauAverage),
#             Males=sum(Gender == 'Male'),
#             CDRlo=sum(CDRBinned=='0.0', na.rm=T),
#             CDRmid=sum(CDRBinned=='0.5', na.rm=T),
#             CDRhi=sum(CDRBinned=='1.0+', na.rm=T))
# 
# # save
# repA <- df[df$Split == 1, ]
# repB <- df[df$Split == 2, ]
# 
# outdata = 'adni_tau_nmf_matrix_repA.csv'
# outfeatures = NULL
# outids = 'adni_tau_nmf_ids_repA.csv'
# 
# XA <- create.nmf.input(repA, features=features, id.cols = c("RID", "DateTau"),
#                        outdata = outdata, outfeatures = outfeatures,
#                        outids = outids)
# 
# outdata = 'adni_tau_nmf_matrix_repB.csv'
# outfeatures = NULL
# outids = 'adni_tau_nmf_ids_repB.csv'
# XB <- create.nmf.input(repB, features=features, id.cols = c("RID", "DateTau"),
#                        outdata = outdata, outfeatures = outfeatures,
#                        outids = outids)

# === Create reproducibility splits - repeated =============

# number of repeats
N <- 50

splits <- matrix(data=NA, nrow=nrow(df), ncol=N)
continuous.vars <- df[, c("Age", "CorticalTauAverage")]
categorical.vars <- df[, c("Gender", "CDRBinned")]

for (i in 1:N) {
  v <- anticlustering(continuous.vars, 2, categories=categorical.vars)
  splits[, i] <- v
}

splits <- as.data.frame(splits)
write.table(splits, 'reproducibility_split_indices.csv', row.names = F, col.names = F, sep=',')