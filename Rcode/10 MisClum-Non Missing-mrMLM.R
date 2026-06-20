############# this is prepare the dataset for mrMLMM software
############# this is prepare the dataset for mrMLMM software

####C.124.2560
####C.124.2560
myGMmrMLMM <- read.csv("data/myGMC.124.2560.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDC.125.2562.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMC.124.2560.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myY1<- read.csv("data/myYC.124.2560.csv")
Yldmakingup <- read.csv("data/alltraitsyldmakeup.csv",row.names=1)
Yldmakingup <- Yldmakingup[,c(1,4:6,3)]
str(Yldmakingup)
names(Yldmakingup)
Yld <- merge(myY1,Yldmakingup, by="Taxa", all.x =TRUE)
str(Yld)
myYmrMlMM <- Yld
names(myYmrMlMM)
##the trait names
## [1] "Taxa"      "CmDW_g"    "Cml_cm"    "CmD_BI_mm" "CmD_LI_mm" "CmN."      "Bcirc_cm"  "Yld_kg"    "SDW_kg"    "CCirc_cm"  "Lg" "GS"       
##[13] "FD"        "SRD"       "ADD"       "Yld"       "AIL"       "CmDW.V"    "Yld_1"  
## corresponding name for the rMLMM software
## [1] "<phenotype>"  "Trait1"  "Trait2"   "Trait3"    "Trait4"   "Trait5"     "Trait6"    "Trait7"    "Trait8"    "Trait9"    "Trait10"  "Trait11"    
##[13] "Trait12"  "Trait13"  "Trait14"  "Trait15"     "Trait16"     "Trait17"     "Trait18"  
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
names(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/DataC124/myYmrMlMMC.124.2560_1.csv", row.names = FALSE, na = "NA")

##this one for without population structure
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMC.124.2560.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMC.124.2560_1.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:18,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "~/Documents/R-corde for miscanthus project/Miscanthus/mrMLMM2/C124PSNO")

##this one for with population structure
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMC.124.2560.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMC.124.2560_1.csv",
      fileKin=NULL,filePS=filePS="~/Documents/R-corde for miscanthus project/Miscanthus/mrMLMM2/DataC125/mymrMLMMQC.125.124.2560.1.csv",
      Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:18,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "~/Documents/R-corde for miscanthus project/Miscanthus/mrMLMM2/C124PS")


####C.116.3293
####C.116.3293
myGMmrMLMM <- read.csv("data/myGMC.116.3293.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDC.116.3293.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMC.116.3293.csv", row.names = FALSE, na = "NA")

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)

myY1<- read.csv("data/myYC.116.3293.csv")
Yldmakingup <- read.csv("data/alltraitsyldmakeup.csv",row.names=1)
Yldmakingup <- Yldmakingup[,c(1,4:6,3)]
str(Yldmakingup)
names(Yldmakingup)
Yld <- merge(myY1,Yldmakingup, by="Taxa", all.x =TRUE)
str(Yld)
myYmrMlMM <- Yld
names(myYmrMlMM)
##the trait names
## [1] "Taxa"      "CmDW_g"    "Cml_cm"    "CmD_BI_mm" "CmD_LI_mm" "CmN."      "Bcirc_cm"  "Yld_kg"    "SDW_kg"    "CCirc_cm"  "Lg" "GS"       
##[13] "FD"        "SRD"       "ADD"       "Yld"       "AIL"       "CmDW.V"    "Yld_1"  
## corresponding name for the rMLMM software
## [1] "<phenotype>"  "Trait1"  "Trait2"  "Trait3"  "Trait4"   "Trait5"   "Trait6"  "Trait7"  "Trait8"  "Trait9"  "Trait10"  "Trait11"    
##[13] "Trait12"     "Trait13"     "Trait14"     "Trait15"     "Trait16"     "Trait17"     "Trait18"  
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/DataC116/myYmrMlMMC.116.3293_1.csv", row.names = FALSE, na = "NA")
str(myYmrMlMM)
library("mrMLM")
mrMLM(fileGen="~/Documents/R-corde for miscanthus project/Miscanthus/mrMLMM2/DataC116/subgenomrMLMMC.116.3293.csv",
      filePhe="~/Documents/R-corde for miscanthus project/Miscanthus/mrMLMM2/DataC116/myYmrMlMMC.116.3293_1.csv.csv",
      
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:18,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=TRUE,
      Plotformat ="jpeg",Resolution="Low", dir= "~/Documents/R-corde for miscanthus project/Miscanthus/mrMLMM2/C116PSNO")

### With PS
library("mrMLM")
mrMLM(fileGen="~/Documents/R-corde for miscanthus project/Miscanthus/mrMLMM2/DataC116/subgenomrMLMMculm.csv",
      filePhe="~/Documents/R-corde for miscanthus project/Miscanthus/mrMLMM2/DataC116/myYmrMlMMC.116.3293_1.csv.csv",
      fileKin=NULL,filePS="~/Documents/R-corde for miscanthus project/Miscanthus/mrMLMM2/DataC116/mymrMLMMQC.116.3293.1.csv",
      PopStrType= "Q",fileCov=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:18,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=TRUE,
      Plotformat ="jpeg",Resolution="Low", dir= "~/Documents/R-corde for miscanthus project/Miscanthus/mrMLMM2/C116PS")


####C.106.4202
####C.106.4202
myGMmrMLMM <- read.csv("data/myGMC.106.4202.csv")
colnames(myGMmrMLMM)[which(names(myGMmrMLMM) == "Name")] <- "rn"
###subset the matched Geotype (So this one will contain genotype for code 1)
subgeno <- plyr::join(data.frame(myGMmrMLMM),data.frame(allsnp2sub),by="rn")
###step5 import the phenotype myY for GAPIT to select the matched individual from genotype from step 4
myGDmrMLMM <- read.csv("data/myGDC.106.4202.csv",row.names=1)
subgenotran <- data.frame(t(myGDmrMLMM))-1
rn <- rownames(subgenotran)
rownames(subgenotran) <- NULL
subgenotran <- cbind(rn,subgenotran)
subgenomrMLMM <- plyr::join(subgeno,subgenotran,by="rn")
#str(subgenomrMLMM)
###change the name in order to fit the software requirment
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "rn")] <- "rs#"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Chromosome")] <- "chrom"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "Position")] <- "pos"
colnames(subgenomrMLMM)[which(names(subgenomrMLMM) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
#str(subgenomrMLMM)
write.csv(subgenomrMLMM, file = "mrMLMM2/subgenomrMLMMC.106.4202.csv", row.names = FALSE, na = "NA")
str(subgenomrMLMM)

###step6 chagne the name of the phenotype and also seperate them with less missing values
###change the name of pheotype in order to fit the software requirment
##change the name of individuals in the phenotype data to the same of individuals in the genotypes
#myYmrMlMM$Taxa <- gsub("-", "\\.", myYmrMlMM$Taxa)
myY1 <- read.csv("data/myYC.106.4202.csv")
Yldmakingup <- read.csv("data/alltraitsyldmakeup.csv",row.names=1)
Yldmakingup <- Yldmakingup[,c(1,4:6,3)]
str(Yldmakingup)
names(Yldmakingup)
Yld <- merge(myY1,Yldmakingup, by="Taxa", all.x =TRUE)
str(Yld)
myYmrMlMM <- Yld
names(myYmrMlMM)
##the trait names
## [1] "Taxa"      "CmDW_g"    "Cml_cm"    "CmD_BI_mm" "CmD_LI_mm" "CmN."      "Bcirc_cm"  "Yld_kg"    "SDW_kg"    "CCirc_cm"  "Lg" "GS"       
##[13] "FD"        "SRD"       "ADD"       "Yld"       "AIL"       "CmDW.V"    "Yld_1"  
## corresponding name for the rMLMM software
## [1] "<phenotype>"  "Trait1"  "Trait2"  "Trait3"  "Trait4"   "Trait5"   "Trait6"  "Trait7"  "Trait8"  "Trait9"  "Trait10"  "Trait11"    
##[13] "Trait12"     "Trait13"     "Trait14"     "Trait15"     "Trait16"     "Trait17"     "Trait18"  
myYmrMlMM$Taxa <- make.names(myYmrMlMM$Taxa)
###check the name of Taxa
#chaning the name of the Taxa to phenotype in order to fit the requirment 
colnames(myYmrMlMM)[which(names(myYmrMlMM) == "Taxa")] <- "<phenotype>"
names(myYmrMlMM)
##depending on the missing value and then separate the data to three data set: 
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Trait",i-1,sep="")
  }  
  return(x)
}
myYmrMlMM <- colRename(myYmrMlMM)
write.csv(myYmrMlMM, file = "mrMLMM2/DataC106/myYmrMlMMC.106.4202_1.csv", row.names = FALSE, na = "NA")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMC.106.4202.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMC.106.4202_1.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:18,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "~/Documents/R-corde for miscanthus project/Miscanthus/mrMLMM2/C106PSNO")

library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\subgenomrMLMMC.106.4202.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMC.106.4202_1.csv",
      fileKin=NULL,filePS="~/Documents/R-corde for miscanthus project/Miscanthus/mrMLMM2/DataC116/mymrMLMMQC.106.4202.1.csv",Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:18,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=FALSE,
      Plotformat ="jpeg",Resolution="Low", dir= "~/Documents/R-corde for miscanthus project/Miscanthus/mrMLMM2/C106PS")
