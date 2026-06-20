##############################################################################
##########Calculating PVE for the result from FarmCPU#########################
##############################################################################

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myYSRD <- read.csv("data/myYimputedSNP19SRD.csv")
myYOWA <- read.csv("data/myYimputedSNP19OWA2.csv")
myY1 <- read.csv("data/myYimputedSNP19.csv")
myY1 <- myY1[,-c(2,18)]
myY2 <- Reduce(function(x, y) merge(x, y, all.x=TRUE), list(myY1,myYOWA,myYSRD))

Yldmakingup <- read.csv("data/alltraitsyldmakeup.csv",row.names=1)
Yldmakingup <- Yldmakingup[,c(1,4:6,3)]
str(Yldmakingup)
myY <- merge(myY2,Yldmakingup, by="Taxa", all.x=TRUE)
str(myY)
colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
#flo.1 <- read.csv("GAPIT005all/FarmCPU0.05.Adj.P.PS.csv")
flo.1 <- read.csv("GAPIT005Result/FarmCPU0.05.Adj.P.PS.csv")
str(flo.1)
levels(flo.1$Ind.SNP)
levels(flo.1$Trait.name)
flo <- flo.1[flo.1$Ind.SNP=="124+36088im" & flo.1$Method=="FarmCPU",]
str(flo)
names(flo)
flo <- droplevels(flo)
levels(flo$Ind.SNP)
levels(as.factor(flo$Method))
levels(as.factor(flo$Trait.name))
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- flo[!(flo$Trait.name=="Surv"),]
flo <- flo[!(flo$Trait.name=="OWA"),]
flo <- flo[!(flo$Trait.name=="mc"),]
names(flo)
names(myY)
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- droplevels(flo)
levels(flo$Trait.name)

#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FLOSUPPER.7.15.116im")
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FLOSUPPER.7.25.116im")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.36088imup")
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
fileNames <- Sys.glob(paste("FarmCPUCV0.05.36088imup/*.", extension, sep = ""))

#myY <- read.csv("data/myYimputedSNP19.csv")
myGD <- read.csv("data/myGDimputedSNP19.csv")

mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCVim <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(as.factor(FarmCPUCVim$filename))

FarmCPUCVim$filename <- gsub('FarmCPUCV0.05.36088imup/Trait.', '', FarmCPUCVim$filename)
FarmCPUCVim$filename <- gsub('.csv', '', FarmCPUCVim$filename)

levels(as.factor(FarmCPUCVim$filename))

FarmCPUCVim.1 <- unite(FarmCPUCVim,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")
flo.2$SNP.T %in% FarmCPUCVim.1$SNP.T
FarmCPUCVim.flo <- merge(flo.2,FarmCPUCVim.1,by="SNP.T")
FarmCPUCVim.flo.S <- separate(FarmCPUCVim.flo, SNP.T, c("SNP","Trait.name"),sep="/")
levels(as.factor(FarmCPUCVim.flo.S$Trait.name))
str(FarmCPUCVim.flo.S)

###for Culm 106
###for Culm 106
###for Culm 106
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYC.106.4202_up.csv")
names(myY)
#Yldmakingup <- read.csv("data/alltraitsyldmakeup.csv",row.names=1)
#str(Yldmakingup)
#myY <- merge(myY1,Yldmakingup, by="Taxa", all.x=TRUE)
str(myY)

colnames(myY)[colnames(myY)=="CmDW.V"] <- "CmDW"
#colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
#names(myY)
#colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDC.106.4202_up.csv")

#flo.1 <- read.csv("GAPIT005all/FarmCPU0.05.Adj.P.PS.csv")
flo.1 <- read.csv("GAPIT005Result/FarmCPU0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="C106+4202" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.C106")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.C106/*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV0.05.C106 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV0.05.C106$filename)

FarmCPUCV0.05.C106$filename <- gsub('FarmCPUCV0.05.C106/Trait.', '', FarmCPUCV0.05.C106$filename)
FarmCPUCV0.05.C106$filename <- gsub('.csv', '', FarmCPUCV0.05.C106$filename)

levels(as.factor(FarmCPUCV0.05.C106$filename))

FarmCPUCV0.05.C106.1 <- unite(FarmCPUCV0.05.C106,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV0.05.C106.flo <- merge(flo.2,FarmCPUCV0.05.C106.1,by="SNP.T")
FarmCPUCV0.05.C106.flo.S <- separate(FarmCPUCV0.05.C106.flo, SNP.T, c("SNP","Trait.name"),sep="/")


###for Culm 116
###for Culm 116
###for Culm 116
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYC.116.3293_up.csv")
#Yldmakingup <- read.csv("data/alltraitsyldmakeup.csv",row.names=1)
#str(Yldmakingup)
#myY <- merge(myY1,Yldmakingup, by="Taxa", all.x=TRUE)
str(myY)
names(myY)
colnames(myY)[colnames(myY)=="CmDW.V"] <- "CmDW"
#colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
#names(myY)
#colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
names(myY)
myGD <- read.csv("data/myGDC.116.3293_up.csv")

#flo.1 <- read.csv("GAPIT005all/FarmCPU0.05.Adj.P.PS.csv")
flo.1 <- read.csv("GAPIT005Result/FarmCPU0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="C116+3293" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.C116")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.C116/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV0.05.C116 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV0.05.C116$filename)

FarmCPUCV0.05.C116$filename <- gsub('FarmCPUCV0.05.C116/Trait.', '', FarmCPUCV0.05.C116$filename)
FarmCPUCV0.05.C116$filename <- gsub('.csv', '', FarmCPUCV0.05.C116$filename)

levels(as.factor(FarmCPUCV0.05.C116$filename))

FarmCPUCV0.05.C116.1 <- unite(FarmCPUCV0.05.C116,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV0.05.C116.flo <- merge(flo.2,FarmCPUCV0.05.C116.1,by="SNP.T")
FarmCPUCV0.05.C116.flo.S <- separate(FarmCPUCV0.05.C116.flo, SNP.T, c("SNP","Trait.name"),sep="/")

###for Culm 124
###for Culm 124
###for Culm 124
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
myY <- read.csv("data/myYC.124.2560_up.csv")

#Yldmakingup <- read.csv("data/alltraitsyldmakeup.csv",row.names=1)
#str(Yldmakingup)
#myY <- merge(myY1,Yldmakingup, by="Taxa", all.x=TRUE)
str(myY)
colnames(myY)[colnames(myY)=="CmDW.V"] <- "CmDW"
#colnames(myY)[colnames(myY)=="HD_50."] <- "HD_50"
names(myY)
#colnames(myY)[colnames(myY)=="FD_50."] <- "FD_50"
colnames(myY)[colnames(myY)=="CmN."] <- "CmN"
colnames(myY)[colnames(myY)=="TFN."] <- "TFN"
#names(myY)
myGD <- read.csv("data/myGDC.124.2560_up.csv")
#flo.1 <- read.csv("GAPIT005all/FarmCPU0.05.Adj.P.PS.csv")
flo.1 <- read.csv("GAPIT005Result/FarmCPU0.05.Adj.P.PS.csv")
flo <- flo.1[flo.1$Ind.SNP=="C124+2560" & flo.1$Method=="FarmCPU",]
flo <- flo[!(flo$Trait.name=="Surv"),]
colnames(flo)[colnames(flo)=="Name"] <- "SNP"
flo <- droplevels(flo)
levels(flo$Trait.name)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/FarmCPUCV0.05.C124")
for (i in levels(flo$Trait.name)){
  Trait <- flo[which(flo$Trait.name== i),]
  if(nrow(Trait) >0) {
    write.csv(Trait,file= paste('Trait', i, 'csv', sep = '.'))
  }
}

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
extension <- "csv"
fileNames <- Sys.glob(paste("FarmCPUCV0.05.C124/*.", extension, sep = ""))
source("Function/PVE.fun.R")
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample.PVE <- PVE.fun(sample,myY,myGD)
  mzList[[i]] = data.frame(sample.PVE, filename = rep(fileNames[i], length(nrow(sample))))
}

#resultPVE1 <- plyr::ldply(mzList, data.frame)
FarmCPUCV0.05.C124 <- plyr::rbind.fill(mzList)
#resultPVE = do.call("rbind", mzList)
##combine all of the result with flo

levels(FarmCPUCV0.05.C124$filename)

FarmCPUCV0.05.C124$filename <- gsub('FarmCPUCV0.05.C124/Trait.', '', FarmCPUCV0.05.C124$filename)
FarmCPUCV0.05.C124$filename <- gsub('.csv', '', FarmCPUCV0.05.C124$filename)

levels(as.factor(FarmCPUCV0.05.C124$filename))

FarmCPUCV0.05.C124.1 <- unite(FarmCPUCV0.05.C124,"SNP.T",SNP,filename,sep="/")
flo.2 <- unite(flo,"SNP.T",SNP,Trait.name,sep="/")

FarmCPUCV0.05.C124.flo <- merge(flo.2,FarmCPUCV0.05.C124.1,by="SNP.T")
FarmCPUCV0.05.C124.flo.S <- separate(FarmCPUCV0.05.C124.flo, SNP.T, c("SNP","Trait.name"),sep="/")

###combine all of the result
FarmCPUCV0.05.PVE.a <- do.call("rbind",list(FarmCPUCVim.flo.S,FarmCPUCV0.05.C106.flo.S,FarmCPUCV0.05.C116.flo.S,FarmCPUCV0.05.C124.flo.S))

#write.csv(FarmCPUCV0.05.PVE.a,file="GAPIT005all/FarmCPUCV0.05.PVE.a.1.csv")
write.csv(FarmCPUCV0.05.PVE.a,file="GAPIT005Result/FarmCPUCV0.05.PVE.a.1.csv")








