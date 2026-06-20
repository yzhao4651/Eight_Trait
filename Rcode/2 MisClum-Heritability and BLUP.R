################# this one is to get the heritability and 
###import the data
###import the data
###check the data format
#str(normadata)
#this one include 
normadata <- read.csv("data/traits1718normalited_Yld.csv",na.strings = c("",".","NA"),row.names=1)
###check the data format
str(normadata)
###change the format of the several variables 
normadata$GS <- as.numeric(as.character(normadata$GS))
normadata$SRD <- as.numeric(as.character(normadata$SRD))
normadata$Entry = as.factor(normadata$Entry)
normadata$Rep=as.factor(normadata$Rep)
normadata$Year=as.factor(normadata$Year)
###check the data format
str(normadata)
names(normadata)
###get the heritablity of all the traits using the function in the Function

source("Function/Heritability.OWA.R")

Herit <- Heritability(2,17,normadata)

write.csv(Herit, file = "data/HeritablityChapter2.csv", na = ".")

###get the heritablity of all the traits using the function in the Function
getwd()
source("/Users/yonglizhao/Documents/Eight_Traits/Rcodes updates/Function/BLUP.R")
#install.packages("lme4")
library(lme4)

library(Matrix)
ranefvalueall <- ranefvalue(2,17,normadata)

colnames(ranefvalueall)[which(names(ranefvalueall)=="CmDW")] <- "CmDW_g"
colnames(ranefvalueall)[which(names(ranefvalueall)=="Cml")] <- "Cml_cm"
colnames(ranefvalueall)[which(names(ranefvalueall)=="CmD_BI")] <- "CmD_BI_mm"
colnames(ranefvalueall)[which(names(ranefvalueall)=="CmD_LI")] <- "CmD_LI_mm"
colnames(ranefvalueall)[which(names(ranefvalueall)=="CmN")] <- "CmN."
colnames(ranefvalueall)[colnames(ranefvalueall)=="Bcirc"] <- "Bcirc_cm"
colnames(ranefvalueall)[colnames(ranefvalueall)=="CCirc"] <- "CCirc_cm"

colnames(ranefvalueall)[colnames(ranefvalueall)=="SDW_kg"] <- "SDW_kg"

colnames(ranefvalueall)[colnames(ranefvalueall)=="yldcheck"] <- "Yld_1"


write.csv(ranefvalueall, file = "data/ranefvalueChapter2.csv", row.names = T, na = ".")




