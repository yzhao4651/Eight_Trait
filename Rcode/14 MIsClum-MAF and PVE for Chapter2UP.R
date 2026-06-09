#################################################################################
########Get MAF and PVE for all of traits      ##################################
#################################################################################
###import all the dataset with MAF 0.05
###ECMML method with CV
###ECMML method with CV
getwd()
setwd("/Users/yonglizhao/Documents/Eight_Traits/")

EPS <- read.csv("GAPIT0.05Result/ECMML0.05.PS.csv",row.names=1)
str(EPS)
levels(as.factor(EPS$Trait.name))

#chang maf to MAF
colnames(EPS)[colnames(EPS)=="maf"] <- "MAF"
str(EPS)

###method with GLM, MLM, MLMM with CV
GPS <- read.csv("GAPIT0.05Result/GLMMLMMLMM0.05.PS.csv",row.names=1)
colnames(GPS)[colnames(GPS)=="maf"] <- "MAF"
str(GPS)
levels(as.factor(GPS$Trait.name))

#####method SUPER with PS, the results from SUPER method without PS will not be included because of two many false positive results without PS in the Super method
CMLM <- read.csv("GAPIT0.05Result/CMLMSUPER0.05.PS.PVE.2.csv",row.names=1)
colnames(CMLM)[colnames(CMLM)=="maf"] <- "MAF"
str(CMLM)
levels(as.factor(CMLM$Trait.name))
#####method SUPER
MLM <- read.csv("GAPIT0.05Result/MLMSUPER0.05.PS.PVE.a.2.csv",row.names=1)
colnames(MLM)[colnames(MLM)=="maf"] <- "MAF"
str(MLM)
levels(as.factor(MLM$Trait.name))

####rMLMM with flowering with CV
rMFPS <- read.csv("GAPIT0.05Result/rMLMMflo.PS.csv",row.names=1)
str(rMFPS)
levels(as.factor(rMFPS$Trait.name))
####rMLMM with flowering with NO CV
rMFPSNO <- read.csv("GAPIT0.05Result/rMLMMflo.PSNO.csv",row.names=1)
str(rMFPSNO)

###mrMLM with Culm trait with CV
rMCPS <- read.csv("GAPIT0.05Result/rMLMMCulm.PS.csv",row.names=1)
str(rMCPS)

levels(as.factor(rMCPS$Trait.name))

###mrMLM with Culm trait with NO CV
rMCPSNO <- read.csv("GAPIT0.05Result/rMLMMCulm.PSNO.csv",row.names=1)
str(rMCPSNO)

levels(as.factor(rMCPSNO$Trait.name))

###rrBLUP with Culm trait with NO CV
rrBLUP <- read.csv("GAPIT0.05Result/rrBLUPPCA.PVE.a.csv",row.names=1)
str(rrBLUP)

levels(as.factor(rrBLUP$Trait.name))

###FarmCPU without CV
FarmPSNO <- read.csv("GAPIT0.05Result/FarmCPU0.05.PVE.a.1.csv",row.names=1)
colnames(FarmPSNO)[colnames(FarmPSNO)=="maf"] <- "MAF"
str(FarmPSNO)

levels(as.factor(FarmPSNO$Trait.name))
####FarmCPU without CV
FarmPS <- read.csv("GAPIT0.05Result/FarmCPUCV0.05.PVE.a.1.csv",row.names=1)
colnames(FarmPS)[colnames(FarmPS)=="maf"] <- "MAF"
str(FarmPS)
levels(as.factor(FarmPS$Trait.name))


###combineing all of them together 
alltrait.method <- plyr::ldply(list(EPS,GPS,rMFPS,rMFPSNO,rMCPS,rMCPSNO,rrBLUP,FarmPS,FarmPSNO,CMLM,MLM))
str(alltrait.method)
###
###change the P.value to with zero decime 
alltrait.method$P.value.1 <- as.numeric(formatC(alltrait.method$P.value, format = "e", digits = 0))
alltrait.method$MAF.1<- as.numeric(formatC(alltrait.method$MAF, format = "e", digits = 3))
alltrait.method$PVE.1 <- as.numeric(formatC(alltrait.method$PVE, format = "e", digits = 3))
##write out the total traits 
write.csv(alltrait.method, file="data/alltraitsGWASresultupdated2.csv")

###select all the flowering trait

###import the all traits 
getwd()
setwd("../..")
all <- read.csv("data/alltraitsGWASresultupdated2.csv")
levels(all$Trait.name)
###import the all results and then selete the traits for Chapter 2

all$Trait.name <- as.factor(all$Trait.name)
levels(all$Trait.name)
levels(all$Trait.name)[levels(all$Trait.name)=="CmDW"] <- "CmDW_V"
levels(all$Trait.name)[levels(all$Trait.name)=="CmDW.V"] <- "CmDW_V_g"
levels(all$Trait.name)[levels(all$Trait.name)=="CmN."] <- "CmN"

levels(all$Trait.name)[levels(all$Trait.name)=="OWA.18"] <- "OWA18"
levels(all$Trait.name)[levels(all$Trait.name)=="OWA.19"] <- "OWA19"
#levels(all$Trait.name)[levels(all$Trait.name)=="CmN."] <- "CmN"
levels(all$Ind.SNP)[levels(all$Ind.SNP)=="+36088im"] <- "124+36088im"
#levels(all$Ind.SNP)[levels(all$Ind.SNP)=="124+36088im"] <- "124+36088im"

levels(as.factor(all$Trait.name))

###this one for eight traits
###this one for eight traits
###this one for eight traits
str(all)
levels(all$Trait.name)
all.chapter2 <- all[all$Trait.name=="Bcirc_cm"| all$Trait.name=="CCirc_cm"|
                      all$Trait.name=="CmD_LI_mm"| all$Trait.name=="CmDW_g"| 
                      all$Trait.name=="CmD_BI_mm"| all$Trait.name=="Cml_cm"|
                      all$Trait.name=="CmN"| all$Trait.name=="Yld",]
all.chapter2 <-  droplevels(all.chapter2)

###checking 
str(all.chapter2)
all.chapter2$P.value.1 <- as.character(formatC(all.chapter2$P.value.1, format = "e", digits = 0))
str(all.chapter2)
all.chapter2 <- droplevels(all.chapter2)
levels(as.factor(all.chapter2$Trait.name))

##write out the dataset chapter2
write.csv(all.chapter2, file="data/all.chapter2.5.csv",row.names = F)


### trying to get the Max of the MAF and PVE for each of trait

i <- "Bcirc_cm"
culm.MAF.PVE.Max <- list()
for (i in levels(all.chapter2$Trait.name)){
  culmall <- all.chapter2[all.chapter2$Trait.name==i,]
  culmall <- culmall[order(culmall$SNP, -(as.numeric(culmall$MAF.1))),] #sort by id and reverse of abs(value)
  culmall.MAF.2 <- culmall[ !duplicated(culmall$SNP), ] 
  culmall <- culmall[order(culmall$SNP, -(as.numeric(culmall$PVE.1))), ] #sort by id and reverse of abs(value)
  culmall.PVE.2 <- culmall[ !duplicated(culmall$SNP), ] 
  culm.MAF.PVE.Max [[i]]<- merge(x=culmall.MAF.2,y=culmall.PVE.2, by="SNP" )
}
allculm.MAF.PVE.max <- do.call(rbind,culm.MAF.PVE.Max)
str(allculm.MAF.PVE.max)

####write the data out set 

write.csv(allculm.MAF.PVE.max,file="data/allculm.MAF.PVE.Max.2.csv")
