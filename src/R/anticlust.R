install.packages("anticlust")
library(anticlust)
install.packages("dplyr")
library(dplyr)
setwd("/Volumes/PrecunealSSD/Singularity/ADNI/bids/derivatives")
subjects = read.csv("table_cdr_0p5_anan.csv")
continuous.vars <- subjects[, c("Age")]
categorical.vars <- subjects[, c("Sex")]
v <- anticlustering(continuous.vars, 2, categories=categorical.vars)
subjects$Split <- v
repA = subjects[subjects$Split == 1, ]
repB = subjects[subjects$Split == 2, ]
write.table(repA, "table_cdr_0p5_anan_repA.csv", row.names = F, col.names = T, sep = ',', quote = F)
write.table(repB, "table_cdr_0p5_anan_repB.csv", row.names = F, col.names = T, sep = ',', quote = F)
