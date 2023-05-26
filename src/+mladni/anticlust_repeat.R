install.packages("anticlust")
library(anticlust)
install.packages("dplyr")
library(dplyr)
install.packages("rapportools")
library(rapportools)
setwd("/Volumes/PrecunealSSD/Singularity/ADNI/bids/derivatives")

clist <- c("cdr_0p5_aneg_emci", "cdr_0p5_aneg_lmci", 
           "cdr_0p5_apos_emci", "cdr_0p5_apos_lmci", 
           "cdr_0p5_aneg", "cdr_0p5_apos", "cdr_0p5_anan")

# "cn", "preclinical", 
# "cdr_0p5_aneg_emci", "cdr_0p5_aneg_lmci", 
# "cdr_0p5_apos_emci", "cdr_0p5_apos_lmci", 
# "cdr_0p5_aneg", "cdr_0p5_apos", "cdr_0p5_anan",
# "cdr_gt_0p5_aneg",  "cdr_gt_0p5_apos"

for (j in clist) {
  
  subjects = read.csv(paste0("table_",j,".csv"))
  continuous.vars <- subjects[, c("Age")]
  categorical.vars <- subjects[, c("Sex")]
  
  N <- 20
  for (i in 1:N) {
    v <- anticlustering(continuous.vars, 2, categories=categorical.vars)
    subjects$Split <- v
    selectA = subjects$Split == 1 & !is.empty(c("Filename"))
    repA = subjects[selectA, c("Filename")]
    selectB = subjects$Split == 2 & !is.empty(c("Filename"))
    repB = subjects[selectB, c("Filename")]
    write.table(repA, paste0("table_",j,"_repA", i, ".csv"), row.names = F, col.names = F, sep = ',', quote = F)
    write.table(repB, paste0("table_",j,"_repB", i, ".csv"), row.names = F, col.names = F, sep = ',', quote = F)
  }
}


