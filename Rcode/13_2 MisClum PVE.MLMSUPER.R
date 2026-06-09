
############################################################################
###Get the PVE for the data from the method of SUPER + MLM with PS##########
############################################################################

### this one for MLMSUPER
myYSRD <- read.csv("data/myYimputedSNP19SRD.csv")
myYOWA <- read.csv("data/myYimputedSNP19OWA2.csv")
myY1 <- read.csv("data/myYimputedSNP19.csv")
myY1 <- myY1[,-c(2,18)]
myY <- Reduce(function(x, y) merge(x, y, all.x=TRUE), list(myY1,myYOWA,myYSRD))
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
flo.1 <- read.csv("GAPIT0.05Result/SUPER0.05.PS.csv")
str(flo.1)
levels(flo.1$Ind.SNP)
flo <- flo.1[flo.1$Ind.SNP=="124+36088im" & flo.1$Method=="MLM+SUPER",]
str(flo)
names(flo)
levels(flo$Trait.name)
flo <- droplevels(flo)
levels(flo$Ind.SNP)
levels(as.factor(flo$Method))
levels(as.factor(flo$Trait.name))
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="fprin"),]
flo <- flo[!(flo$Trait.name=="OWA"),]
levels()
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FLOSUPPER.7.15.116im")
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FLOSUPPER.8.3.116im")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/MLMSUPER0.05.36088imUP1")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}


PVE.fun <- function(sample,myY,myGD){
  if (nrow(sample)==1){
    ftdGDSNPs <- data.frame(myGD$Taxa,myGD[names(myGD) %in% sample$SNP])
    colnames(ftdGDSNPs)[colnames(ftdGDSNPs)=="myGD.Taxa"] <- "Taxa"
    flowerSNPs <- merge(myY, ftdGDSNPs,by="Taxa")
    a <- ncol(myY)+1
    lm1 <-lm(flowerSNPs[[levels(sample$Trait.name)]] ~ flowerSNPs[[a]])
    af <- anova(lm1)
    afss <- af$"Sum Sq"
    PctExp <- afss/sum(afss)*100
    PVE3 <- data.frame(cbind(colnames(ftdGDSNPs[2]),PctExp[[1]]))
    colnames(PVE3) <- c("SNP","PVE")
  } else {
    ftdGDSNPs <- data.frame(myGD$Taxa,myGD[names(myGD) %in% sample$SNP])
    colnames(ftdGDSNPs)[colnames(ftdGDSNPs)=="myGD.Taxa"] <- "Taxa"
    flowerSNPs <- merge(myY, ftdGDSNPs,by="Taxa")
    number=1
    a <- ncol(myY)+1
    b <- a+nrow(sample)-1
    out_variable = colnames(flowerSNPs[a:b])
    outcome <- matrix(NA, nrow=ncol(flowerSNPs[a:b]),
                      ncol = 1, dimnames = list(out_variable, "PctExp"))
    for(j in a:b){
      lm1 <-lm(flowerSNPs[[levels(sample$Trait.name)]] ~ flowerSNPs[[j]])
      af <- anova(lm1)
      afss <- af$"Sum Sq"
      PctExp <- afss/sum(afss)*100
      outcome[number] <-PctExp[[1]]
      number=number+1
      PVE3 <- data.frame(outcome)
      colnames(PVE3) <-  c("PVE")
      SNP <- rownames(PVE3)
      rownames(PVE3) <- NULL
      PVE3 <- cbind(SNP,PVE3)
    }
  }
  return(PVE3)
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("MLMSUPER0.05.36088imUP1/*.", extension, sep = ""))

#myY <- read.csv("data/myYimputedSNP19.csv")
myGD <- read.csv("data/myGDimputedSNP19.csv")

mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
MLMCUPERim <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(MLMCUPERim$filename)

MLMCUPERim$filename <- gsub('MLMSUPER0.05.36088imUP1/Trait.', '', MLMCUPERim$filename)
MLMCUPERim$filename <- gsub('.csv', '', MLMCUPERim$filename)

levels(as.factor(MLMCUPERim$filename))

MLMCUPERim.1 <- unite(MLMCUPERim,"SNP.T",SNP,filename,sep="/")
flo.2<- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

MLMCUPERim.flo <- merge(flo.2,MLMCUPERim.1,by="SNP.T")
MLMCUPERim.flo.S <- separate(MLMCUPERim.flo, SNP.T, c("SNP","Trait.name"),sep="/")
write.csv(MLMCUPERim.flo.S,file="GAPIT0.05Result/MLMCUPERim.flo.S1.csv")
levels(as.factor(MLMCUPERim.flo.S$Trait.name))


###for Culm 106
###for Culm 106
###for Culm 106
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYC.106.4202.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDC.106.4202.csv")
flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="C106+4202" & flo.1$Method=="MLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/MLMSUPER0.05.C106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("MLMSUPER0.05.C106/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
MLMSUPER0.05.C106 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(MLMSUPER0.05.C106$filename)

MLMSUPER0.05.C106$filename <- gsub('MLMSUPER0.05.C106/Trait.', '', MLMSUPER0.05.C106$filename)
MLMSUPER0.05.C106$filename <- gsub('.csv', '', MLMSUPER0.05.C106$filename)

levels(as.factor(MLMSUPER0.05.C106$filename))

MLMSUPER0.05.C106.1 <- unite(MLMSUPER0.05.C106,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

MLMSUPER0.05.C106.flo <- merge(flo.2,MLMSUPER0.05.C106.1,by="SNP.T")
MLMSUPER0.05.C106.flo.S <- separate(MLMSUPER0.05.C106.flo, SNP.T, c("SNP","Trait.name"),sep="/")
MLMSUPER0.05.C106.flo.S <- MLMSUPER0.05.C106.flo.S[-18]
write.csv(MLMSUPER0.05.C106.flo.S,file="GAPIT0.05Result/MLMSUPER0.05.C106.flo.S.csv")

###for Culm 116
###for Culm 116
###for Culm 116
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYC.116.3293.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDC.116.3293.csv")
flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="C116+3293" & flo.1$Method=="MLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- droplevels(flo)
levels(as.factor(flo$Trait.name))
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/MLMSUPER0.05.C116")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("MLMSUPER0.05.C116/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
MLMSUPER0.05.C116 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(MLMSUPER0.05.C116$filename)

MLMSUPER0.05.C116$filename <- gsub('MLMSUPER0.05.C116/Trait.', '', MLMSUPER0.05.C116$filename)
MLMSUPER0.05.C116$filename <- gsub('.csv', '', MLMSUPER0.05.C116$filename)

levels(as.factor(MLMSUPER0.05.C116$filename))

MLMSUPER0.05.C116.1 <- unite(MLMSUPER0.05.C116,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

MLMSUPER0.05.C116.flo <- merge(flo.2,MLMSUPER0.05.C116.1,by="SNP.T")
MLMSUPER0.05.C116.flo.S <- separate(MLMSUPER0.05.C116.flo, SNP.T, c("SNP","Trait.name"),sep="/")

MLMSUPER0.05.C116.flo.S <- MLMSUPER0.05.C116.flo.S[-18]
write.csv(MLMSUPER0.05.C116.flo.S,file="GAPIT0.05Result/MLMSUPER0.05.C116.flo.S.csv")

###for Culm 124
###for Culm 124
###for Culm 124
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYC.124.2560.csv")
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDC.124.2560.csv")

#flo.1 <- read.csv("Result.8.3/MLMCMLMSUPER0.05.Adj.P.PS.csv")
flo.1 <- read.csv("GAPIT0.05Result/SUPER0.05.PS.csv")
str(flo.1)
flo <- flo.1[flo.1$Ind.SNP=="C124+2560" & flo.1$Method=="MLM+SUPER",]
#flo <- flo.1[flo.1$Ind.SNP=="124.2560" & flo.1$Method=="MLM+SUPER",]
flo <- flo[!(flo$Trait.name=="Surv"),]


setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/MLMSUPER0.05.C124")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("MLMSUPER0.05.C124/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
MLMSUPER0.05.C124 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(MLMSUPER0.05.C124$filename)

MLMSUPER0.05.C124$filename <- gsub('MLMSUPER0.05.C124/Trait.', '', MLMSUPER0.05.C124$filename)
MLMSUPER0.05.C124$filename <- gsub('.csv', '', MLMSUPER0.05.C124$filename)

levels(as.factor(MLMSUPER0.05.C124$filename))

MLMSUPER0.05.C124.1 <- unite(MLMSUPER0.05.C124,"SNP.T",SNP,filename, sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

MLMSUPER0.05.C124.flo <- merge(flo.2,MLMSUPER0.05.C124.1,by="SNP.T")
MLMSUPER0.05.C124.flo.S <- separate(MLMSUPER0.05.C124.flo, SNP.T, c("SNP","Trait.name"),sep="/")
write.csv(MLMSUPER0.05.C124.flo.S,file="GAPIT0.05Result/MLMSUPER0.05.C124.flo.S.csv")


###combine all of the result
MLMSUPER0.05.PVE.a.2 <- do.call("rbind",list(MLMSUPER0.05.C106.flo.S,MLMSUPER0.05.C116.flo.S,MLMSUPER0.05.C124.flo.S, MLMCUPERim.flo.S))

write.csv(MLMSUPER0.05.PVE.a.2,file="GAPIT0.05Result/MLMSUPER0.05.PS.PVE.a.2.csv")
###combine 
#MLMCMLMSUPER0.05.PVE.a <- do.call("rbind", list(MLMSUPER0.05.PVE.a, CMLMSUPER0.05.PVE.a))

#write.csv(MLMCMLMSUPER0.05.PVE.a,file="GAPIT0.05Result/AllMLMCMLMSUPER0.05.PS.PVE.a.csv")

