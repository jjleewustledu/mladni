#####################################################################################
#                Analysis code used in "Patterns of Neurodegeneration"              #
#                                    2023 Feb 28                                    #
#####################################################################################

require(mgcv)
library(visreg)
library(ggplot2)
library(reshape2)
library(corrplot)
library(voxel)
library(ppcor)
library(readr)
library(stats)
library(dplyr) # provides filter
library(pracma) # Matlab emulations

# add R functions to run models in individual networks and extract relevant statistics
source("/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG/PehlivanovaEtAllScripts-master/Rfunctions/regDataGam.R") 
source("/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG/PehlivanovaEtAllScripts-master/Rfunctions/regDataGamMultOut.R")

# load data
setwd("/Volumes/PrecunealSSD/Singularity/ADNI/NMF_FDG")
u <- read_csv("mladni_FDGDemographics_table_covariates_on_cn_assembled.csv")
u$MergeAgeSqr <- u$MergeAge^2 
u$AmyloidStatus <- as.numeric(factor(u$AmyloidStatus))
u$MergePtGender <- factor(u$MergePtGender)
u$MergePtGenderNum <- as.numeric(u$MergePtGender)
u$MergePtEthCat <- as.numeric(factor(u$MergePtEthCat))
u$MergePtRacCat <- as.numeric(factor(u$MergePtRacCat))
u$MergeApoE4 <- as.numeric(factor(u$MergeApoE4))
u$MergeDx <- factor(u$MergeDx)
u$MergeMmse <- as.numeric(factor(u$MergeMmse))
u$Phase <- as.numeric(factor(u$Phase))
u$CDGLOBAL <- as.numeric(factor(u$CDGLOBAL))
u$SITEID <- as.numeric(factor(u$SITEID))

# descriptive statistics for cdrsb
mean(u$MergeCdrsbBl)
sd(u$MergeCdrsbBl)

# associations with age and sex in GAM models
gmAgeSex01 <- gam(MergeCdrsbBl ~ s(MergeAge) + MergePtGender, data=u)
print(summary(gmAgeSex01))
gmAgeSex02 <- gam(MergeCdrsbBl ~ s(MergeAge) + MergePtGender + s(MergeAge, by=MergePtGender), data=u)
print(summary(gmAgeSex02))

# association with biomarkers, SES, cognition and other CD, ADNI confounders : pearson estimate       p.value statistic    n gp  Method
pcor.test(u$MergeCdrsbBl,u$AmyloidStatus,u[,c("MergeAge","MergePtGenderNum","MergeAgeSqr")],method="pearson") # 1 0.3074011 1.453391e-32  12.18604 1428  3 pearson
pcor.test(u$MergeCdrsbBl,u$MergeApoE4,u[,c("MergeAge","MergePtGenderNum","MergeAgeSqr")],method="pearson") # 1 0.2589835 2.848611e-23  10.11465 1428  3 pearson
pcor.test(u$MergeCdrsbBl,u$ICV,u[,c("MergeAge","MergePtGenderNum","MergeAgeSqr")],method="pearson") # not sig.

pcor.test(u$MergeCdrsbBl,u$MergePtEducat,u[,c("MergeAge","MergePtGenderNum","MergeAgeSqr")],method="pearson") # 1 -0.1489131 1.62389e-08  -5.68074 1428  3 pearson
pcor.test(u$MergeCdrsbBl,u$MergePtEthCat,u[,c("MergeAge","MergePtGenderNum","MergeAgeSqr")],method="pearson") # not sig.
pcor.test(u$MergeCdrsbBl,u$MergePtRacCat,u[,c("MergeAge","MergePtGenderNum","MergeAgeSqr")],method="pearson") # not sig.

pcor.test(u$MergeCdrsbBl,u$MergeMmse,u[,c("MergeAge","MergePtGenderNum","MergeAgeSqr")],method="pearson") # 1 -0.6597282 9.035285e-179 -33.11586 1428  3 pearson
pcor.test(u$MergeCdrsbBl,u$CDMEMORY,u[,c("MergeAge","MergePtGenderNum","MergeAgeSqr")],method="pearson") # 1 0.8198335       0  54.01014 1428  3 pearson
pcor.test(u$MergeCdrsbBl,u$CDORIENT,u[,c("MergeAge","MergePtGenderNum","MergeAgeSqr")],method="pearson") # 1 0.8180875       0  53.66159 1428  3 pearson
pcor.test(u$MergeCdrsbBl,u$CDJUDGE,u[,c("MergeAge","MergePtGenderNum","MergeAgeSqr")],method="pearson") # 1 0.8217448       0  54.39695 1428  3 pearson
pcor.test(u$MergeCdrsbBl,u$CDCOMMUN,u[,c("MergeAge","MergePtGenderNum","MergeAgeSqr")],method="pearson") # 1 0.8149814       0  53.05257 1428  3 pearson
pcor.test(u$MergeCdrsbBl,u$CDHOME,u[,c("MergeAge","MergePtGenderNum","MergeAgeSqr")],method="pearson") # 1 0.8068609       0   51.5233 1428  3 pearson
pcor.test(u$MergeCdrsbBl,u$CDCARE,u[,c("MergeAge","MergePtGenderNum","MergeAgeSqr")],method="pearson") # 1 0.4638428 6.287052e-77  19.75058 1428  3 pearson
pcor.test(u$MergeCdrsbBl,u$CDGLOBAL,u[,c("MergeAge","MergePtGenderNum","MergeAgeSqr")],method="pearson") # 1 0.7959614 1.647533e-312  49.60093 1428  3 pearson

pcor.test(u$MergeCdrsbBl,u$Phase,u[,c("MergeAge","MergePtGenderNum","MergeAgeSqr")],method="pearson") # 1 -0.1282658 1.18724e-06 -4.878827 1428  3 pearson
pcor.test(u$MergeCdrsbBl,u$SITEID,u[,c("MergeAge","MergePtGenderNum","MergeAgeSqr")],method="pearson") # not sig.

### main analyses of NMF networks
# create vector of model statements
Ncomp <- 28
nmf28Names <-names(u)[26:53] # NMF names
ms_<-paste0(nmf28Names,"~MergeCdrsbBl+s(MergeAge)+MergePtGender+AmyloidStatus+MergeApoE4+MergePtEducat+Phase")
#print(ms_)
ms2_<-paste0(nmf28Names,"~s(MergeCdrsbBl)+s(MergeAge)+MergePtGender+AmyloidStatus+MergeApoE4+MergePtEducat+Phase")
#print(ms2_)

# running basic models for association between NMF networks and cdrsb
modOut<-regDataGam(ms_,u,nmf28Names, 1)
allComp<-rownames(modOut)
sigComp<-rownames(modOut[which(modOut$tval < -11),]) # find top quartile
print(sigComp) # best components by tval:  5 8 9 11 15 19 21 24
modOut2<-regDataGam(ms2_,u,nmf28Names, 2)
allComp2<-rownames(modOut2)
sigComp2<-rownames(modOut[which(modOut2$F > 61),])
print(sigComp2) # best components by F:  2 9 11 13 18 19 20 21

# partial correlations between cdrsb and components 
for (idx in 1:32) {
  print(nmf28Names[idx])
  print(pcor.test(u$MergeCdrsbBl,u[,nmf28Names[idx]],u[,c("MergeAge","MergeAgeSqr","MergePtGender","AmyloidStatus","MergeApoE4","MergePtEducat","Phase")],method="pearson"))
} # all p < 2e-6, not discriminating

# correlations between age and components 
for (idx in 1:32) {
  print(nmf28Names[idx])
  print(pcor.test(u$MergeAge,u[,nmf28Names[idx]],u[,c("MergeCdrsbBl","MergePtGender","AmyloidStatus","MergeApoE4","MergePtEducat","Phase")],method="pearson"))
} # not sig.:  9 17 21 29

### age effects 
# interaction models in individual NMF networks
modInteraction<-paste0(nmf28Names,"~ MergeAge*MergeCdrsbBl + MergeCdrsbBl + MergeAge + MergeAgeSqr + MergePtGender")
pvals = NA
for (i in 1:length(modInteraction)){
  fprintf("======================= %s =======================\n", nmf28Names[i])
  foo<-lm(as.formula(modInteraction[i]),data=u)
  tab<-summary(foo)
  print(tab)
  pvals[i]<-tab$coefficients[5,4]
} # interaction of cdrsb*age significant at 1 3 4 5 6 7 8 9 10 11 13 14 15 16 18 19 20 21 23 24 25 26 27 30 32
  # not significant at 2 12 17 22 28 29 31
median(pvals)
min(pvals)
max(pvals)

### sex effects 
# interaction models in individual NMF networks
modInteraction<-paste0(nmf28Names,"~ MergePtGender*MergeCdrsbBl + MergeCdrsbBl + MergeAge + MergeAgeSqr + MergePtGender")
pvals = NA
for (i in 1:length(modInteraction)){
  fprintf("======================= %s =======================\n", nmf28Names[i])
  foo<-lm(as.formula(modInteraction[i]),data=u)
  tab<-summary(foo)
  print(tab)
  pvals[i]<-tab$coefficients[5,4]
} # interaction of cdrsb*age significant at 1 4 21 22 23 24 25 
# not significant otherwise
median(pvals)
min(pvals)
max(pvals)

### visualizations of component effects (scatterplots in Figure 4)
visreg_comp <- function(u, idx) {
  names(u)[25+idx] <- 'Component'
  gm.comp <- gam(Component ~ s(MergeCdrsbBl) + s(MergeAge) + MergePtGender + AmyloidStatus + MergeApoE4 + MergePtEducat + Phase, method="REML", data=u)
  #ylab_ <- sprintf("Component %i Weighted FDG", idx)
  #xlab_ <- "CDR sum of Boxes"
  visreg(gm.comp, 'MergeCdrsbBl', overlay=T,ylab="",xlab="",cex.axis=2.5,cex.lab=1,cex.main=1,points=list(cex=1))
} # 5 8 9 11 15 19 21 24
visreg_comp(u, 5)
visreg_comp(u, 8)
visreg_comp(u, 9)
visreg_comp(u, 11)
visreg_comp(u, 15)
visreg_comp(u, 19)
visreg_comp(u, 21)
visreg_comp(u, 24)

### interaction plots (Figure 5)
# dataset with cdrsb quartiles
plot_quarts_cdrsb <- function(u, idx) {
  MergeAge <- u$MergeAge
  AmyloidStatus <- u$AmyloidStatus
  MergeApoE4 <- u$MergeApoE4
  MergePtEducat <- u$MergePtEducat
  Phase <- u$Phase
  
  names(u)[25+idx] <- 'Component'
  v <- u %>% mutate(MergeCdrsbBlQ = ntile(u$MergeCdrsbBl, 4))
  v$MergeCdrsbBlQ <- sprintf('Q%i', v$MergeCdrsbBlQ)
  dataQuart<-v[v$MergeCdrsbBlQ %in% c('Q1','Q4'),]
  gm1 <- gam(Component ~ s(MergeAge) + MergeCdrsbBlQ + MergePtGender + AmyloidStatus + MergeApoE4 + MergePtEducat + Phase, data=dataQuart, method="REML")
  plot1 <- plotGAM(gm1, "MergeAge", "MergeCdrsbBlQ", orderedAsFactor = T, rawOrFitted = "raw", plotCI = T)
  
  print(plot1 + theme_bw() + ylab(sprintf("", idx)) + xlab("") +
    ggtitle("") + scale_size(range=c(8,20)) +
    theme(plot.title = element_text(hjust = 0.5,size=22,face="bold"), panel.background = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.title=element_text(size=24), legend.text = element_text(size = 24),
          legend.title = element_text(size=24), axis.text.x = element_text(size=24), axis.text.y = element_text(size=24)) +
    scale_colour_discrete(name="CDRSB"))
} # 5 8 9 11 15 19 21 24
plot_quarts_cdrsb(u, 5)
plot_quarts_cdrsb(u, 8)
plot_quarts_cdrsb(u, 9)
plot_quarts_cdrsb(u, 11)
plot_quarts_cdrsb(u, 15)
plot_quarts_cdrsb(u, 19)
plot_quarts_cdrsb(u, 21)
plot_quarts_cdrsb(u, 24)

# dataset with age quartiles
plot_quarts_age <- function(u, idx) {
  MergeCdrsbBl <- u$MergeCdrsbBl
  AmyloidStatus <- u$AmyloidStatus
  MergeApoE4 <- u$MergeApoE4
  MergePtEducat <- u$MergePtEducat
  Phase <- u$Phase
  
  names(u)[25+idx] <- 'Component'
  v <- u %>% mutate(MergeAgeQ = ntile(u$MergeAge, 4))
  v$MergeAgeQ <- sprintf('Q%i', v$MergeAgeQ)
  dataQuart<-v[v$MergeAgeQ %in% c('Q1','Q4'),]
  gm1 <- gam(Component ~ s(MergeCdrsbBl) + MergeAgeQ + MergePtGender + AmyloidStatus + MergeApoE4 + MergePtEducat + Phase, data=dataQuart, method="REML")
  plot1 <- plotGAM(gm1, "MergeCdrsbBl", "MergeAgeQ", orderedAsFactor = T, rawOrFitted = "raw", plotCI = T)

  ylab_ <- "" # sprintf("FDG Pattern %i", idx)
  xlab_ <- "" # "CDR sum of boxes"
  print(plot1 + theme_bw() + ylab(ylab_) + xlab(xlab_) +
          ggtitle("") + scale_size(range=c(8,20)) +
          theme(plot.title = element_text(hjust = 0.5,size=24,face="bold"), panel.background = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.title=element_text(size=24), legend.text = element_text(size = 24),
                legend.title = element_text(size=24), axis.text.x = element_text(size=24), axis.text.y = element_text(size=24)) +
          scale_colour_discrete(name="Age"))
} # 5 8 9 11 15 19 21 24
plot_quarts_age(u, 5)
plot_quarts_age(u, 8)
plot_quarts_age(u, 9)
plot_quarts_age(u, 11)
plot_quarts_age(u, 15)
plot_quarts_age(u, 19)
plot_quarts_age(u, 21)
plot_quarts_age(u, 24)

### sensitivity analyses
# adding Medu as a covariate 
#data_medu <- data_k[!is.na(data_k$meduCnbGo1),]
#pcor.test(data_medu$logk,data_medu$meduCnbGo1,data_medu[,c("ageAtScan","sex.o","ageSq")],method="pearson")

ms1 <-paste0(nmf28Names,"~s(MergeCdrsbBl)+s(MergeAge)+MergePtGender+AmyloidStatus+MergeApoE4+MergePtEducat+Phase")
ms1_model<-regDataGamMultOut(ms1,u,nmf28Names,2,1)
print(ms1_model)
ms2 <-paste0(nmf28Names,"~s(MergeAge)+s(MergeCdrsbBl)+MergePtGender+AmyloidStatus+MergeApoE4+MergePtEducat+Phase")
ms2_model<-regDataGamMultOut(ms2,u,nmf28Names,2,1)
print(ms2_model)
ms3 <-paste0(nmf28Names,"~MergeCdrsbBl+MergeAge+MergePtGender+AmyloidStatus+MergeApoE4+MergePtEducat+Phase")
ms3_model<-regDataGamMultOut(ms3,u,nmf28Names,1,1)
print(ms3_model)
ms4 <-paste0(nmf28Names,"~MergeAge+MergeCdrsbBl+MergePtGender+AmyloidStatus+MergeApoE4+MergePtEducat+Phase")
ms4_model<-regDataGamMultOut(ms3,u,nmf28Names,1,1)
print(ms4_model)
# components sorted by tval:  

#partial correlation for table 1, t-value best for 7
pcor.test(u$MergeCdrsbBl,u$Components_7,u[,c("MergeAge","MergePtGender","AmyloidStatus","MergeApoE4","MergePtEducat","Phase")],method="spearman")





### more ideas from Marieta #########################################################

# adding general cognitive factor as covariate 
#ms9 <-paste0(nmf28Names,"~logk+Overall_Accuracy+sex.o+s(ageAtScan)")
#ms9_model<-regDataGamMultOut(ms9,data_k,nmf28Names,1,2)

# sensitivity analysis on subjects not using medication
#data_k_noMed <-data_k[which(data_k$medsComb==0),] # subsetting data with subjects NOT using medications
#modOutNoMeds<-regDataGam(ms,data_k_noMed,nmf28Names,1) # w/o the noise component
#medSigComp<-rownames(modOutNoMeds[which(modOutNoMeds$fdr.p<.05),])

### multivariate prediction with NMF components 
# baseline covariates
covs1<-"sex.o+ageAtScan"
covs2<-"sex.o+ageAtScan+F1_Exec_Comp_Res_Accuracy+F3_Memory_Accuracy"
nmfcovs<-"Nmf20C1+Nmf20C2+Nmf20C3+Nmf20C4+Nmf20C5+Nmf20C6+Nmf20C7+Nmf20C8+Nmf20C9+Nmf20C10+Nmf20C11+Nmf20C12+Nmf20C13+Nmf20C14+Nmf20C15+Nmf20C16+Nmf20C18+Nmf20C19+Nmf20C20+"
covsMedu<-"sex.o+ageAtScan+meduCnbGo1"

# linear models
lmBase1<-lm(as.formula(paste0("logk~",covs1)),data=data_k)	          # model with demographic covariates
lmBase2<-lm(as.formula(paste0("logk~",covs2)),data=data_k)	          # model with demographics + cognitive variables
lmBase3<-lm(as.formula(paste0("logk~",covs1)),data=data_medu)             # model with demographic covariates in dataset with maternal education
lmBase4<-lm(as.formula(paste0("logk~",covsMedu)),data=data_medu)          # model with demographics + maternal in dataset with maternal education
lmCt1<-lm(as.formula(paste0("logk~",nmfcovs,covs1)),data=data_k)          # model with demographics + NMF networks
lmCt2<-lm(as.formula(paste0("logk~",nmfcovs,covs2)),data=data_k)          # model with demographics + cognitive variables + NMF networks
lmCtMedu<-lm(as.formula(paste0("logk~",nmfcovs,covsMedu)),data=data_medu) # model with demographics + maternal education + NMF networks

# model with just demographics
cor.test(predict(lmBase1), data_k$logk)
# model with just demographics + cognition
cor.test(predict(lmBase2),data_k$logk) 
# CT model against baseline of just age and sex
anova(lmCt1,lmBase1)
cor.test(predict(lmCt1),data_k$logk)
# CT model against baseline of just age and sex + cognition
anova(lmCt2,lmBase2)
cor.test(predict(lmCt2),data_k$logk)
# model with just demographics + medu
cor.test(predict(lmBase4),data_k$logk[-which(is.na(data_k$meduCnbGo1))])
# CT model against baseline of just age and sex + medu
anova(lmCtMedu,lmBase4)
cor.test(predict(lmCtMedu),data_medu$logk)

### plotting multivariate prediction (Figure 7)
# plotting values predicted from model with NMF networks + demographics against actual log k values
gamCovs1<-"sex.o+s(ageAtScan)"
nmfcovs<-"Nmf20C1+Nmf20C2+Nmf20C3+Nmf20C4+Nmf20C5+Nmf20C6+Nmf20C7+Nmf20C8+Nmf20C9+Nmf20C10+Nmf20C11+Nmf20C12+Nmf20C13+Nmf20C14+Nmf20C15+Nmf20C16+Nmf20C18+Nmf20C19+Nmf20C20+"
nmfModel<-gam(as.formula(paste0("logk~",nmfcovs,gamCovs1)),data=data_k,method="REML")
data_k$nmfPred<-predict(nmfModel)

plotNmfDemoModel <- lm(logk~nmfPred,data=data_k)
cor.test(data_k$logk,data_k$nmfPred,method="pearson")

# plotting
par("las"=1)
visreg(plotNmfDemoModel, 'nmfPred',xlab="Predicted Discount Rate (Log k)",ylab="Actual Discount Rate (Log k)",cex.axis=1.5,cex.lab=1.5,cex.main=1.8,points=list(cex=0.9),line.par=list(col="blue"))
legend(x='bottomright', legend='r = 0.33, p < 0.0001',bty = "n",cex=1.235)

