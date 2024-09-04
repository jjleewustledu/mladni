#install.packages("anticlust")
library(anticlust)
#install.packages("dplyr")
library(dplyr)
#install.packages("rapportools")
library(rapportools)
#install.packages("icesTAF")
library(icesTAF)
setwd("/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG")

subjects = read.csv("anticlust_ad_repeat.csv")
continuous.vars <- subjects[, c("Age")]
categorical.vars <- subjects[, c("Sex", "ApoE4")]

N <- 50
for (i in 1:N) {
  v <- anticlustering(continuous.vars, 2, categories=categorical.vars) # make 2 clusters
  subjects$Split <- v
  selectA = subjects$Split == 1 & !is.empty(c("Filelist"))
  repA = subjects[selectA, c("Filelist")]
  selectB = subjects$Split == 2 & !is.empty(c("Filelist"))
  repB = subjects[selectB, c("Filelist")]
  
  mkdir(paste0("baseline_ad_repA", i))
  mkdir(paste0("baseline_ad_repB", i))
  write.table(repA, paste0("baseline_ad_repA", i, "/nifti_files.csv"), row.names = F, col.names = F, sep = ',', quote = F)
  write.table(repB, paste0("baseline_ad_repB", i, "/nifti_files.csv"), row.names = F, col.names = F, sep = ',', quote = F)
}


