#######################################################################################
######Get all of the mean and variance of cross-validation for all of the traits#######
#######################################################################################
###Mean → estimate overall prediction performance######################################
###Variance → measure how stable/reliable that performance is#########################

###import all the .csv file from different file 
###emerging all the data together 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPNEW/")#set up the folder for the results to put in. 
#setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GPNEW/")
extension <- "csv"
fileNames <- Sys.glob(paste("GP*/Rsquare*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i], row.names=1)
  require(pastecs)
  options(scipen=100)
  options(digits=4)
  stat.desc(sample)
  mz1=stat.desc(sample)
  names <- rownames(mz1)
  rownames(mz1) <- NULL
  data <- cbind(names,mz1)
  mzList[[i]] = data.frame(data, filename = rep(fileNames[i], length(mz1)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
GP.Summ<- plyr::ldply(mzList, data.frame)
str(GP.Summ)
levels(GP.Summ$filename)
GP.Summ$filename <- gsub('.csv', '', GP.Summ$filename)
GP.Summ$filename <- gsub('Rsquare ', '', GP.Summ$filename)
GP.Summ$filename <- gsub('GP/', 'GPIM/', GP.Summ$filename)
GP.Summ$filename <- gsub('GP', '', GP.Summ$filename)
str(GP.Summ)
levels(as.factor(GP.Summ$filename))
##separated by the"/"
library(tidyverse)
GP.Summ.all <- separate(GP.Summ, filename,c("Ind.SNP","Trait"),"/")
GP.Summ.all$Trait <- gsub('GR', '', GP.Summ.all$Trait)
str(GP.Summ.all)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
write.csv(GP.Summ.all, file= "data/GPtotal3.Summ.all.csv")
levels(as.factor(GP.Summ.all$Trait))
levels(as.factor(GP.Summ.all$Ind.SNP))


GP.all.1 <- GP.Summ.all
levels(as.factor(GP.all.1$Trait))
levels(as.factor(GP.all.1$Ind.SNP))
levels(as.factor(GP.all.1$names))
GP.all.1 <- droplevels(GP.all.1)
GP.all.1$Ind.SNP <- as.factor(GP.all.1$Ind.SNP)

levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RR0"] <- "F116im/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo106R0"] <- "F106/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo116R0"] <- "F116/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo96R0"] <- "F96/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA102R0"] <- "OWA102/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA112R0"] <- "OWA112/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA122R0"] <- "OWA122/NOPC"
#levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA135R"] <- "OWA135/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC106R0"] <- "C106/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC116R0"] <- "C116/NOPC"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC124R0"] <- "C124/NOPC"

levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RR1"] <- "F116im/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo106R1"] <- "F106/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo116R1"] <- "F116/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo96R1"] <- "F96/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA102R1"] <- "OWA102/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA112R1"] <- "OWA112/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA122R1"] <- "OWA122/PC1"
#levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA135R"] <- "OWA135/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC106R1"] <- "C106/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC116R1"] <- "C116/PC1"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC124R1"] <- "C124/PC1"

levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RR2"] <- "F116im/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo106R2"] <- "F106/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo116R2"] <- "F116/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo96R2"] <- "F96/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA102R2"] <- "OWA102/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA112R2"] <- "OWA112/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA122R2"] <- "OWA122/PC1~2"
#levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA135R"] <- "OWA135/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC106R2"] <- "C106/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC116R2"] <- "C116/PC1~2"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC124R2"] <- "C124/PC1~2"

levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RR3"] <- "F116im/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo106R3"] <- "F106/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo116R3"] <- "F116/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="flo96R3"] <- "F96/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA102R3"] <- "OWA102/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA112R3"] <- "OWA112/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA122R3"] <- "OWA122/PC1~3"
#levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="OWA135R"] <- "OWA135/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC106R3"] <- "C106/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC116R3"] <- "C116/PC1~3"
levels(GP.all.1$Ind.SNP)[levels(GP.all.1$Ind.SNP)=="RC124R3"] <- "C124/PC1~3"


GP.all.2 <- separate(GP.all.1, Ind.SNP,c("Ind.SNP","PC"),"/")

levels(as.factor(GP.all.2$Trait))
levels(as.factor(GP.all.2$Ind.SNP))
levels(as.factor(GP.all.2$names))
levels(as.factor(GP.all.2$PC))
GP.all.2 <- droplevels(GP.all.2)

GP.all.2$Trait <- as.factor(GP.all.2$Trait)
write.csv(GP.all.2,file="data/GP.all.2.csv")


