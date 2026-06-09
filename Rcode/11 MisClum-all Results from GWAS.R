#################################################################
###Organizing results from all software           ###############
#################################################################

#################################################################
###These results do not separate from other traits###############
#################################################################
###check ECMMLCV and ECMMLNOCV
##Result from ECMML method 
##Result from ECMML method 
##Result from ECMML method 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GAPIT005/ECMMLCV0.05/")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/ECMMLCV0.05/")
extension <- "csv"
fileNames <- Sys.glob(paste("*ECMMLCV0.05/*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  sample <- adj_P_function(sample, 4, 4)
  mz.idx = which(sample$Adj.P.P.value <  0.05)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
ECMMLCVa <- plyr::ldply(mzList, data.frame)
str(ECMMLCVa)
levels(ECMMLCVa$filename)
ECMMLCVa$filename <- gsub('.GWAS.Results.csv', '', ECMMLCVa$filename)
ECMMLCVa$filename <- gsub('ECMMLCV0.05/GAPIT.', '.PS_Y.GAPIT.', ECMMLCVa$filename)
#ECMMLCVa$filename <- gsub('ECMMLCV0.05/', 'PS_Y.', ECMMLCVa$filename)
#ECMMLCVa$filename <- gsub('/GAPIT.', '.GAPIT.', ECMMLCVa$filename)
levels(as.factor(ECMMLCVa$filename))
##install packages to get function unit
#install.packages(c("tidyr", "devtools"))
library(tidyr)
ECMMLCVa.1 <- separate(ECMMLCVa, filename, c("Ind.SNP","PS","Software","Method","Trait.name"), "\\.")
str(ECMMLCVa.1)
levels(as.factor(ECMMLCVa$nobs))
ECMMLCVa.1$Ind.SNP <- gsub("All", '+36088im', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("C124", 'C124+2560', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("C106", 'C106+4202', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("C116", 'C116+3293', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("F106", 'F106+4185', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("F116", 'F116+3098', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("F96", 'F96+4814', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("O102", 'O102+4322', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("O112", 'O112+3450', ECMMLCVa.1$Ind.SNP)
ECMMLCVa.1$Ind.SNP <- gsub("O122", 'O122+2646', ECMMLCVa.1$Ind.SNP)

str(ECMMLCVa.1)
ECMMLCVa.1$PVE <- (ECMMLCVa.1$Rsquare.of.Model.with.SNP-ECMMLCVa.1$Rsquare.of.Model.without.SNP)*100
#ECMMLCVa.2 <- ECMMLCVa.1[,c(1:15)]
#str(ECMMLCVa.2)
str(ECMMLCVa.1)
levels(as.factor(ECMMLCVa.1$Ind.SNP))
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
write.csv(ECMMLCVa.1, file="GAPIT0.05Result/ECMML0.05.PS.csv")


##ECMML without CV 
##ECMML without CV
##ECMML without CV
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GAPIT005/ECMML0.05/")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/ECMML0.05/")
extension <- "csv"
fileNames <- Sys.glob(paste("*ECMML0.05/*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  sample <- adj_P_function(sample, 4, 4)
  mz.idx = which(sample$Adj.P.P.value <  0.05)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
ECMMLNOCVa <-  plyr::ldply(mzList, data.frame)
str(ECMMLNOCVa)
levels(ECMMLNOCVa$filename)
ECMMLNOCVa$filename <- gsub('.GWAS.Results.csv', '', ECMMLNOCVa$filename)
ECMMLNOCVa$filename <- gsub('ECMML0.05/GAPIT.', '.PS_N.GAPIT.', ECMMLNOCVa$filename)
levels(as.factor(ECMMLNOCVa$filename))
ECMMLNOCVa.1 <- separate(ECMMLNOCVa, filename, c("Ind.SNP","PS","Software","Method","Trait.name"), "\\.")
str(ECMMLNOCVa.1)
levels(as.factor(ECMMLNOCVa$nobs))
ECMMLNOCVa.1$Ind.SNP <- gsub("All", '124+36088im', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("C124", 'C124+2560', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("C106", 'C106+4202', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("C116", 'C116+3293', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("F106", 'F106+4185', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("F116", 'F116+3098', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("F96", 'F96+4814', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("O102", 'O102+4322', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("O112", 'O112+3450', ECMMLNOCVa.1$Ind.SNP)
ECMMLNOCVa.1$Ind.SNP <- gsub("O122", 'O122+2646', ECMMLNOCVa.1$Ind.SNP)

str(ECMMLNOCVa.1)
ECMMLNOCVa.1$PVE <- (ECMMLNOCVa.1$Rsquare.of.Model.with.SNP-ECMMLNOCVa.1$Rsquare.of.Model.without.SNP)*100
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
write.csv(ECMMLNOCVa.1, file="GAPIT0.05Result/ECMML0.05.PSNO.csv")

##MLMSUPPER and CMLMMSUPPER with CV 
##MLMSUPPER and CMLMMSUPPER with CV 
##MLMSUPPER and CMLMMSUPPER with CV 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GAPIT005/SUPERCV0.05/")
#setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/SUPERCV0.05/")
#setwd("GAPIT005/MLMSUPPER/")
getwd()
setwd("..")
getwd()
setwd("/Users/yonglizhao/Documents/Eight_Traits/GAPIT005/SUPERCV0.05/")
extension <- "csv"
fileNames <- Sys.glob(paste("*MLMSUPPERCV0.05/*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  sample <- adj_P_function(sample, 4, 4)
  mz.idx = which(sample$Adj.P.P.value <  0.05)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
SUPERCVa <-  plyr::ldply(mzList, data.frame)
str(SUPERCVa)
SUPERCVa$filename <- as.factor(SUPERCVa$filename)
levels(as.factor(SUPERCVa$filename))

#SUPERCVa$filename <- gsub('GAPIT005/', '', SUPERCVa$filename)
SUPERCVa$filename <- gsub('.GWAS.Results.csv', '', SUPERCVa$filename)
SUPERCVa$filename <- gsub('CMLMSUPPERCV0.05/GAPIT.SUPER', '.CMLM+SUPER.PS_Y.GAPIT', SUPERCVa$filename)
SUPERCVa$filename <- gsub('MLMSUPPERCV0.05/GAPIT.SUPER', '.MLM+SUPER.PS_Y.GAPIT', SUPERCVa$filename)
#SUPERCVa$filename <- gsub('/GAPIT.', '.GAPIT.', SUPERCVa$filename)
levels(as.factor(SUPERCVa$filename))
str(SUPERCVa.1)
library(tidyr)

SUPERCVa.1 <- separate(SUPERCVa, filename, c("Ind.SNP","Method","PS","Software","Trait.name"), "\\.")
levels(as.factor(SUPERCVa.1$Ind.SNP))
levels(as.factor(SUPERCVa.1$Method))
levels(as.factor(SUPERCVa.1$PS))
levels(as.factor(SUPERCVa.1$Software))
levels(as.factor(SUPERCVa.1$Trait.name))

SUPERCVa.1$Trait.name <- as.factor(SUPERCVa.1$Trait.name)
levels(SUPERCVa.1$Trait.name)[levels(SUPERCVa.1$Trait.name)=="CmDW"] <- "CmDW.V"

levels(as.factor(SUPERCVa.1$Trait.name))

SUPERCVa.1$Ind.SNP <- gsub("C124", 'C124+2560', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("C106", 'C106+4202', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("C116", 'C116+3293', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("F106", 'F106+4185', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("F116", 'F116+3098', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("F96", 'F96+4814', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("O102", 'O102+4322', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("O112", 'O112+3450', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("O122", 'O122+2646', SUPERCVa.1$Ind.SNP)
SUPERCVa.1$Ind.SNP <- gsub("All", "124+36088im", SUPERCVa.1$Ind.SNP)

levels(as.factor(SUPERCVa.1$Ind.SNP))
str(SUPERCVa.1)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")

write.csv(SUPERCVa.1, file="GAPIT0.05Result/SUPER0.05.PS.csv")


##MLMSUPPER and CMLMSUPPER NO CV 
##MLMSUPPER and CMLMSUPPER NO CV 
##MLMSUPPER and CMLMSUPPER NO CV 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GAPIT005/SUPER0.05/")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/SUPER0.05/")
#setwd("GAPIT005/MLMSUPPER/")
extension <- "csv"
fileNames <- Sys.glob(paste("*MLMSUPPER0.05/*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  sample <- adj_P_function(sample, 4, 4)
  mz.idx = which(sample$Adj.P.P.value <  0.05)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
SUPERNOCVa <- plyr::ldply(mzList, data.frame)
#MLMCMLMSUPPERa <- do.call("rbind", mzList)
str(SUPERNOCVa)
levels(SUPERNOCVa$filename)
str(SUPERNOCVa)
levels(SUPERNOCVa$filename)

SUPERNOCVa$filename <- gsub('.GWAS.Results.csv', '', SUPERNOCVa$filename)
SUPERNOCVa$filename <- gsub('CMLMSUPPER0.05/GAPIT.SUPER', '.CMLM+SUPER.PS_N.GAPIT', SUPERNOCVa$filename)
SUPERNOCVa$filename <- gsub('MLMSUPPER0.05/GAPIT.SUPER', '.MLM+SUPER.PS_N.GAPIT', SUPERNOCVa$filename)

levels(as.factor(SUPERNOCVa$filename))
str(SUPERNOCVa.1)
SUPERNOCVa.1 <- separate(SUPERNOCVa, filename, c("Ind.SNP","Method","PS","Software","Trait.name"), "\\.")
levels(as.factor(SUPERNOCVa.1$nobs))
levels(as.factor(SUPERNOCVa.1$Ind.SNP))
str(SUPERNOCVa.1)
SUPERNOCVa.1$Trait.name <- as.factor(SUPERNOCVa.1$Trait.name)
levels(SUPERNOCVa.1$Trait.name)[levels(SUPERNOCVa.1$Trait.name)=="CmDW"] <- "CmDW.V"

levels(as.factor(SUPERCVa.1$Trait.name))

SUPERNOCVa.1$Ind.SNP <- gsub("C124", 'C124+2560', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("C106", 'C106+4202', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("C116", 'C116+3293', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("F106", 'F106+4185', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("F116", 'F116+3098', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("F96", 'F96+4814', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("O102", 'O102+4322', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("O112", 'O112+3450', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("O122", 'O122+2646', SUPERNOCVa.1$Ind.SNP)
SUPERNOCVa.1$Ind.SNP <- gsub("All", "124+36088im", SUPERNOCVa.1$Ind.SNP)

levels(as.factor(SUPERNOCVa.1$Ind.SNP))
str(SUPERNOCVa.1)
#SUPERNOCVa.1$PVE <- (SUPERNOCVa.1$Rsquare.of.Model.with.SNP-SUPERNOCVa.1$Rsquare.of.Model.without.SNP)*100
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
write.csv(SUPERNOCVa.1, file="GAPIT0.05Result/SUPER0.05.PSNO.csv")


###GLM and MLM and MLMM with CV
###GLM and MLM and MLMM with CV
###GLM and MLM and MLMM with CV

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GAPIT005/GLMMLMMLMMCV0.05/")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/GLMMLMMLMMCV0.05/")
extension <- "csv"
fileNames <- Sys.glob(paste("*GLMMLMMLMMCV0.05/GAPIT.*LM.*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  sample <- adj_P_function(sample, 4, 4)
  mz.idx = which(sample$Adj.P.P.value <  0.05)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
GLMMLMMLMMCVa <- plyr::ldply(mzList, data.frame)
str(GLMMLMMLMMCVa)
levels(GLMMLMMLMMCVa$filename)
levels(as.factor(GLMMLMMLMMCVa$nobs))

#GLMMLMMLMMCVa$filename <- gsub('GAPIT005/', '', GLMMLMMLMMCVa$filename)
GLMMLMMLMMCVa$filename <- gsub('.GWAS.Results.csv', '', GLMMLMMLMMCVa$filename)
GLMMLMMLMMCVa$filename <- gsub('GLMMLMMLMMCV0.05/', '.PS_Y.', GLMMLMMLMMCVa$filename)
#GLMMLMMLMMCVa$filename <- gsub('/GAPIT.', '.GAPIT.', GLMMLMMLMMCVa$filename)
levels(as.factor(GLMMLMMLMMCVa$filename))
str(GLMMLMMLMMCVa)
GLMMLMMLMMCVa.1 <- separate(GLMMLMMLMMCVa, filename, c("Ind.SNP","PS","Software","Method","Trait.name"), "\\.")
levels(as.factor(GLMMLMMLMMCVa.1$nobs))
levels(as.factor(GLMMLMMLMMCVa.1$Ind.SNP))


GLMMLMMLMMCVa.1$Ind.SNP <- gsub("C124", 'C124+2560', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("C106", 'C106+4202', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("C116", 'C116+3293', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("F106", 'F106+4185', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("F116", 'F116+3098', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("F96", 'F96+4814', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("O102", 'O102+4322', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("O112", 'O112+3450', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("O122", 'O122+2646', GLMMLMMLMMCVa.1$Ind.SNP)
GLMMLMMLMMCVa.1$Ind.SNP <- gsub("All", "124+36088im", GLMMLMMLMMCVa.1$Ind.SNP)

levels(as.factor(GLMMLMMLMMCVa.1$Ind.SNP))
str(GLMMLMMLMMCVa.1)
GLMMLMMLMMCVa.1$PVE <- (GLMMLMMLMMCVa.1$Rsquare.of.Model.with.SNP-GLMMLMMLMMCVa.1$Rsquare.of.Model.without.SNP)*100
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
write.csv(GLMMLMMLMMCVa.1, file="GAPIT0.05Result/GLMMLMMLMM0.05.PS.csv")

###GLM and MLM, MLMMwithout CV
###GLM and MLM, MLMMwithout CV
###GLM and MLM, MLMMwithout CV

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GAPIT005/GLMMLMMLMM0.05/")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/GLMMLMMLMM0.05/")
extension <- "csv"
fileNames <- Sys.glob(paste("*GLMMLMMLMM0.05/GAPIT.*LM.*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  sample <- adj_P_function(sample, 4, 4)
  mz.idx = which(sample$Adj.P.P.value <  0.05)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
GLMMLMMLMMNOCVa <- plyr::ldply(mzList, data.frame)
str(GLMMLMMLMMNOCVa)
levels(GLMMLMMLMMNOCVa$filename)

#GLMMLMMLMMNOCVa$filename <- gsub('GAPIT005/', '', GLMMLMMLMMNOCVa$filename)
GLMMLMMLMMNOCVa$filename <- gsub('.GWAS.Results.csv', '', GLMMLMMLMMNOCVa$filename)
GLMMLMMLMMNOCVa$filename <- gsub('GLMMLMMLMM0.05/', '.PS_N.', GLMMLMMLMMNOCVa$filename)
#GLMMLMMLMMNOCVa$filename <- gsub('/GAPIT.', '.GAPIT.', GLMMLMMLMMNOCVa$filename)
levels(as.factor(GLMMLMMLMMNOCVa$filename))
str(GLMMLMMLMMNOCVa)
GLMMLMMLMMNOCVa.1 <- separate(GLMMLMMLMMNOCVa, filename, c("Ind.SNP","PS","Software","Method","Trait.name"), "\\.")
levels(as.factor(GLMMLMMLMMNOCVa.1$nobs))
levels(as.factor(GLMMLMMLMMNOCVa.1$Ind.SNP))


GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("C124", 'C124+2560', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("C106", 'C106+4202', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("C116", 'C116+3293', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("F106", 'F106+4185', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("F116", 'F116+3098', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("F96", 'F96+4814', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("O102", 'O102+4322', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("O112", 'O112+3450', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("O122", 'O122+2646', GLMMLMMLMMNOCVa.1$Ind.SNP)
GLMMLMMLMMNOCVa.1$Ind.SNP <- gsub("All", "124+36088im", GLMMLMMLMMNOCVa.1$Ind.SNP)

levels(as.factor(GLMMLMMLMMNOCVa.1$Ind.SNP))
str(GLMMLMMLMMNOCVa.1)
GLMMLMMLMMNOCVa.1$PVE <- (GLMMLMMLMMNOCVa.1$Rsquare.of.Model.with.SNP-GLMMLMMLMMNOCVa.1$Rsquare.of.Model.without.SNP)*100
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
write.csv(GLMMLMMLMMNOCVa.1, file="GAPIT0.05Result/GLMMLMMLMM0.05.PSNO.csv")


##FarmCPU with CV
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GAPIT005/FarmCPUCV0.05")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/FarmCPUCV0.05")
extension <- "csv"
fileNames <- Sys.glob(paste("*FarmCPUCV0.05/FarmCPU.*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  sample <- adj_P_function(sample, 4, 4)
  mz.idx = which(sample$Adj.P.P.value < 0.05 )
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}

#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
FarmCPUCVa <- plyr::ldply(mzList, data.frame)
str(FarmCPUCVa)
levels(FarmCPUCVa$filename)
levels(FarmCPUCVa$nobs)

#FarmCPUCVa$filename <- gsub('GAPIT005/', '', FarmCPUCVa$filename)
FarmCPUCVa$filename <- gsub('.GWAS.Results.csv', '', FarmCPUCVa$filename)
FarmCPUCVa$filename <- gsub('FarmCPUCV0.05', '.PS_Y', FarmCPUCVa$filename)
FarmCPUCVa$filename <- gsub('/FarmCPU.', '.FarmCPU.FarmCPU.', FarmCPUCVa$filename)
levels(as.factor(FarmCPUCVa$filename))
str(FarmCPUCVa)
FarmCPUCVa.1 <- separate(FarmCPUCVa, filename, c("Ind.SNP","PS","Software","Method","Trait.name"), "\\.")
levels(as.factor(FarmCPUCVa.1$bs))
levels(as.factor(FarmCPUCVa.1$Ind.SNP))


FarmCPUCVa.1$Ind.SNP <- gsub("C124", 'C124+2560', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("C106", 'C106+4202', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("C116", 'C116+3293', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("F106", 'F106+4185', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("F116", 'F116+3098', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("F96", 'F96+4814', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("O102", 'O102+4322', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("O112", 'O112+3450', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("O122", 'O122+2646', FarmCPUCVa.1$Ind.SNP)
FarmCPUCVa.1$Ind.SNP <- gsub("All", "124+36088im", FarmCPUCVa.1$Ind.SNP)

levels(as.factor(FarmCPUCVa.1$Ind.SNP))
str(FarmCPUCVa.1)
#FarmCPUCVa.1$PVE <- (FarmCPUCVa.1$Rsquare.of.Model.with.SNP-FarmCPUCVa.1$Rsquare.of.Model.without.SNP)*100
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
write.csv(FarmCPUCVa.1,file="GAPIT0.05Result/FarmCPU0.05.Adj.P.PS.csv")

#write.csv(FarmCPUNOCVa.1,file="GAPIT0.05Result/FarmCPU0.05.PSNO.csv")

##FarmCPU NO CV 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GAPIT005/FarmCPU0.05")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/FarmCPU0.05")
extension <- "csv"
fileNames <- Sys.glob(paste("*FarmCPU0.05/FarmCPU.*Results.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  sample <- adj_P_function(sample, 4, 4)
  mz.idx = which(sample$Adj.P.P.value < 0.05 )
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
FarmCPUNOCVa <- plyr::ldply(mzList, data.frame)
str(FarmCPUNOCVa)
levels(FarmCPUNOCVa$filename)
levels(as.factor(FarmCPUNOCVa$nobs))

FarmCPUNOCVa$filename <- gsub('.GWAS.Results.csv', '', FarmCPUNOCVa$filename)
levels(as.factor(FarmCPUNOCVa$filename))
FarmCPUNOCVa$filename <- gsub('FarmCPU0.05', '.PS_N', FarmCPUNOCVa$filename)
levels(as.factor(FarmCPUNOCVa$filename))
FarmCPUNOCVa$filename <- gsub('/FarmCPU.', '.FarmCPU.FarmCPU.', FarmCPUNOCVa$filename)
levels(as.factor(FarmCPUNOCVa$filename))
str(FarmCPUNOCVa)
FarmCPUNOCVa.1 <- separate(FarmCPUNOCVa, filename, c("Ind.SNP","PS","Software","Method","Trait.name"), "\\.")
levels(as.factor(FarmCPUNOCVa.1$nobs))
levels(as.factor(FarmCPUNOCVa.1$Ind.SNP))


FarmCPUNOCVa.1$Ind.SNP <- gsub("C124", 'C124+2560', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("C106", 'C106+4202', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("C116", 'C116+3293', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("F106", 'F106+4185', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("F116", 'F116+3098', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("F96", 'F96+4814', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("O102", 'O102+4322', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("O112", 'O112+3450', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("O122", 'O122+2646', FarmCPUNOCVa.1$Ind.SNP)
FarmCPUNOCVa.1$Ind.SNP <- gsub("All", "124+36088im", FarmCPUNOCVa.1$Ind.SNP)

levels(as.factor(FarmCPUNOCVa.1$Ind.SNP))
str(FarmCPUNOCVa.1)
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
write.csv(FarmCPUNOCVa.1,file="GAPIT0.05Result/FarmCPU0.05.PSNO.csv")

#################################################################
##########These MLMM results only for culm traits###############
#################################################################

###rmlmm Clum CV with imputed SNPs
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/IMCulmimPS/*Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  names(sample)[names(sample)=="Var_phen.total."] <- "Var_phen..total."
  mz.idx = which(sample$X..log10.P..!= Inf)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMCulmCVim <- do.call("rbind", mzList)
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait2"] <- "CmDW_g"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait3"] <- "Cml_cm"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait4"] <- "CmD_BI_mm"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait5"] <- "CmD_LI_mm"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait6"] <- "CmN." 
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait7"] <- "Bcirc_cm"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait8"] <- "Yld_kg"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait9"] <- "SDW_kg"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait10"] <- "CCirc_cm"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait39"] <- "Yld"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait40"] <- "AIL"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait41"] <- "CmDW.V"
#levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait17"] <- "SRD"
levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait42"] <- "Yld_1"
rMLMMCulmCVim <- droplevels(rMLMMCulmCVim)
str(rMLMMCulmCVim)

###Clum CV with non missing
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/C*PS/*Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  names(sample)[names(sample)=="Var_phen.total."] <- "Var_phen..total."
  colnames(sample)[colnames(sample)=="Marker.position..bp."] <- "Marker.position..bp."
  if (nrow(sample)>0){
    mz.idx = which(sample$X..log10.P..!= Inf)
    mz1 = sample[mz.idx, ]
    mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
  }
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMCVCulmnm <- plyr::ldply(mzList, data.frame)
str(rMLMMCVCulmnm)

levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait1"] <- "CmDW_g"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait2"] <- "Cml_cm"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait3"] <- "CmD_BI_mm"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait4"] <- "CmD_LI_mm"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait5"] <- "CmN." 
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait6"] <- "Bcirc_cm"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait7"] <- "Yld_kg"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait8"] <- "SDW_kg"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait9"] <- "CCirc_cm"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait10"] <- "Lg"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait11"] <- "GS"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait12"] <- "FD"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait13"] <- "SRD"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait14"] <- "ADD"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait15"] <- "Yld"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait16"] <- "AIL"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait17"] <- "CmDW.V"
levels(rMLMMCVCulmnm$Trait.name)[levels(rMLMMCVCulmnm$Trait.name)=="Trait18"] <- "Yld_1" 
str(rMLMMCVCulmnm)

###combina all of the Results together 

rMLMMClumCV <- plyr::ldply(list(rMLMMCulmCVim,rMLMMCVCulmnm,rMLMMCVOWAim,rMLMMCVOWAnm,rMLMMCVSRDim), data.frame)

str(rMLMMClumCV)
levels(rMLMMClumCV$filename)
levels(as.factor(rMLMMClumCV$nobs))

rMLMMClumCV$filename <- gsub('mrMLMM2', 'mrMLMM/2', rMLMMClumCV$filename)
rMLMMClumCV$filename <- gsub('PS', '/PS_Y', rMLMMClumCV$filename)
#rMLMMClumCV$filename <- gsub('Resultfloall1', 'F116+36088im', rMLMMClumCV$filename)
#rMLMMClumCV$filename <- gsub('/GAPIT.', '.GAPIT.', rMLMMClumCV$filename)
levels(as.factor(rMLMMClumCV$filename))
str(rMLMMClumCV)
rMLMMClumCV.1 <- separate(rMLMMClumCV, filename, c("Software","2","Ind.SNP","PS","discard"), "/")
levels(as.factor(rMLMMClumCV.1$nobs))
levels(as.factor(rMLMMClumCV.1$Ind.SNP))

rMLMMClumCV.1$Ind.SNP <- gsub("C124", 'C124+2560', rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("C106", 'C106+4202', rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("C116", 'C116+3293', rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("F106", 'F106+4185', rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("F116", 'F116+3098', rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("F96", 'F96+4814', rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("ResultOWA2all", "124+36088im", rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("SRDim", "124+36088im", rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("IMCulmim", "124+36088im", rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("O102", "O102+4322", rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("O112", "O112+3450", rMLMMClumCV.1$Ind.SNP)
rMLMMClumCV.1$Ind.SNP <- gsub("O122", "O122+2646", rMLMMClumCV.1$Ind.SNP)



levels(as.factor(rMLMMClumCV.1$Ind.SNP))
str(rMLMMClumCV.1)
colnames(rMLMMClumCV.1)[colnames(rMLMMClumCV.1)=="RS."] <- "SNP"
colnames(rMLMMClumCV.1)[colnames(rMLMMClumCV.1)=="Marker.position..bp."] <- "Position"
rMLMMClumCV.1$P.value <- 10^-(rMLMMClumCV.1$X..log10.P..)
colnames(rMLMMClumCV.1)[colnames(rMLMMClumCV.1)=="r2...."] <- "PVE"
str(rMLMMClumCV.1)
levels(as.factor(rMLMMClumCV.1$Ind.SNP))
#rMLMMClumCV.1$PVE <- (rMLMMClumCV.1$Rsquare.of.Model.with.SNP-rMLMMClumCV.1$Rsquare.of.Model.without.SNP)*100
rMLMMClumCV.2 <- rMLMMClumCV.1[,c(4:6,2:3,7:8,10:15,17:18,20)]
#rMLMMClumCV.2 <- rMLMMClumCV.1[,c(4:6,2:3,7:8,10:15,17:18,21)]
str(rMLMMClumCV.2)
levels(as.factor(rMLMMClumCV.2$PS))
write.csv(rMLMMClumCV.2,file="GAPIT0.05Result/rMLMMCulm.PS.csv")


###rmlmm Clum NO CV with imputed SNPs
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/IMCulmimPSNO/*Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  names(sample)[names(sample)=="Var_phen.total."] <- "Var_phen..total."
  colnames(sample)[colnames(sample)=="Marker.Position..bp."] <- "Marker.position..bp."
  mz.idx = which(sample$X..log10.P..!= Inf)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMCulmNOCVim <- do.call("rbind", mzList)
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait2"] <- "CmDW_g"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait3"] <- "Cml_cm"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait4"] <- "CmD_BI_mm"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait5"] <- "CmD_LI_mm"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait6"] <- "CmN." 
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait7"] <- "Bcirc_cm"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait8"] <- "Yld_kg"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait9"] <- "SDW_kg"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait10"] <- "CCirc_cm"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait11"] <- "Lg"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait12"] <- "GS"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait16"] <- "FD"
#levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait17"] <- "SRD"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait18"] <- "ADD"

levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait39"] <- "Yld"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait40"] <- "AIL"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait41"] <- "CmDW.V"
#levels(rMLMMCulmCVim$Trait.name)[levels(rMLMMCulmCVim$Trait.name)=="Trait17"] <- "SRD"
levels(rMLMMCulmNOCVim$Trait.name)[levels(rMLMMCulmNOCVim$Trait.name)=="Trait42"] <- "Yld_1"
rMLMMCulmNOCVim <- droplevels(rMLMMCulmNOCVim)
str(rMLMMCulmNOCVim)


###Clum NO CV non missing
extension <- "csv"
fileNames <- Sys.glob(paste("mrMLMM2/C*PSNO/*Final result.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i])
  names(sample)[names(sample)=="Var_phen.total."] <- "Var_phen..total."
  colnames(sample)[colnames(sample)=="Marker.Position..bp."] <- "Marker.position..bp."
  mz.idx = which(sample$X..log10.P..!= Inf)
  mz1 = sample[mz.idx, ]
  mzList[[i]] = data.frame(mz1, filename = rep(fileNames[i], length(mz.idx)))
}
#resultImputedSNPCV <- plyr::ldply(mzList, data.frame)
#resultImputedSNPCV <- plyr::rbind.fill(mzList)
rMLMMNOCVCulmnm <- plyr::ldply(mzList, data.frame)
str(rMLMMNOCVCulmnm)
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait1"] <- "CmDW_g"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait2"] <- "Cml_cm"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait3"] <- "CmD_BI_mm"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait4"] <- "CmD_LI_mm"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait5"] <- "CmN." 
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait6"] <- "Bcirc_cm"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait7"] <- "Yld_kg"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait8"] <- "SDW_kg"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait9"] <- "CCirc_cm"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait10"] <- "Lg"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait11"] <- "GS"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait12"] <- "FD"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait13"] <- "SRD"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait14"] <- "ADD"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait15"] <- "Yld"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait16"] <- "AIL"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait17"] <- "CmDW.V"
levels(rMLMMNOCVCulmnm$Trait.name)[levels(rMLMMNOCVCulmnm$Trait.name)=="Trait18"] <- "Yld_1" 
str(rMLMMNOCVCulmnm)

###Combine all of the results NO CV
rMLMMClumNOCV <- plyr::ldply(list(rMLMMCulmNOCVim,rMLMMNOCVCulmnm,rMLMMOWANOCVim,rMLMMOWANOCVnm,rMLMMSRDNOCVim), data.frame)
str(rMLMMClumNOCV)
levels(rMLMMClumNOCV$filename)
levels(as.factor(rMLMMClumNOCV$nobs))

rMLMMClumNOCV$filename <- gsub('mrMLMM2', 'mrMLMM/2', rMLMMClumNOCV$filename)
rMLMMClumNOCV$filename <- gsub('PSNO', '/PS_N', rMLMMClumNOCV$filename)
#rMLMMClumNOCV$filename <- gsub('Resultfloall1', 'F116+36088im', rMLMMClumNOCV$filename)
#rMLMMClumNOCV$filename <- gsub('/GAPIT.', '.GAPIT.', rMLMMClumNOCV$filename)
levels(as.factor(rMLMMClumNOCV$filename))
str(rMLMMClumNOCV)
rMLMMClumNOCV.1 <- separate(rMLMMClumNOCV, filename, c("Software","2","Ind.SNP","PS","discard"), "/")
levels(as.factor(rMLMMClumNOCV.1$nobs))
levels(as.factor(rMLMMClumNOCV.1$Ind.SNP))


rMLMMClumNOCV.1$Ind.SNP <- gsub("C124", 'C124+2560', rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("C106", 'C106+4202', rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("C116", 'C116+3293', rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("F106", 'F106+4185', rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("F116", 'F116+3098', rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("F96", 'F96+4814', rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("ResultOWA2all", "124+36088im", rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("SRDim", "124+36088im", rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("IMCulmim", "124+36088im", rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("O102", "O102+4322", rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("O112", "O112+3450", rMLMMClumNOCV.1$Ind.SNP)
rMLMMClumNOCV.1$Ind.SNP <- gsub("O122", "O122+2646", rMLMMClumNOCV.1$Ind.SNP)


levels(as.factor(rMLMMClumNOCV.1$Ind.SNP))
str(rMLMMClumNOCV.1)
colnames(rMLMMClumNOCV.1)[colnames(rMLMMClumNOCV.1)=="RS."] <- "SNP"
colnames(rMLMMClumNOCV.1)[colnames(rMLMMClumNOCV.1)=="Marker.position..bp."] <- "Position"
rMLMMClumNOCV.1$P.value <- 10^-(rMLMMClumNOCV.1$X..log10.P..)
colnames(rMLMMClumNOCV.1)[colnames(rMLMMClumNOCV.1)=="r2...."] <- "PVE"
str(rMLMMClumNOCV.1)
levels(as.factor(rMLMMClumNOCV.1$Ind.SNP))
levels(as.factor(rMLMMClumNOCV.1$PS))
#rMLMMClumNOCV.1$PVE <- (rMLMMClumNOCV.1$Rsquare.of.Model.with.SNP-rMLMMClumNOCV.1$Rsquare.of.Model.without.SNP)*100
rMLMMClumNOCV.2 <- rMLMMClumNOCV.1[,c(4:6,2:3,7:8,10:15,17:18,20)]
str(rMLMMClumNOCV.2)
write.csv(rMLMMClumNOCV.2,file="GAPIT0.05Result/rMLMMCulm.PSNO.csv")



###result from rrBLUP
###result from rrBLUP
####import result from RRBLUP all of imputed SNPs
rrBLUPimall <- read.csv("rrBLUPim/rrBLUPgwasResultsqq.csv",row.names = 1)
str(rrBLUPimall)
names(rrBLUPimall)
rrBLUPimall.1 <- rrBLUPimall[, -c(4,20,42,58)]

str(rrBLUPimall.1)
rrBLUPimall.2 <- rrBLUPimall.1[rowSums(rrBLUPimall.1[4:39])!=0,]
names(rrBLUPimall.2 )

rrBLUPimall.3 <- rrBLUPimall.2[, c(1:3,40:75)]
str(rrBLUPimall.3)
str(rrBLUPimall.3)
names(rrBLUPimall.3)

adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
  }
  return(data)
}
plot.p.adjusted <- adj_P_function(rrBLUPimall.3, 4,39)
str(plot.p.adjusted)
data <- plot.p.adjusted
rrblupList = list()
for(i in 4:39){
  rrBLUP <- data[,c(1:3,i)][data[[i+36]] < 0.05, ]
  if(nrow(rrBLUP)>=1){
    rrblupList[[i]] <- data.frame(rrBLUP, Trait.name = rep(colnames(data[i]), nrow(rrBLUP)))
    colnames(rrblupList[[i]])[colnames(rrblupList[[i]])== colnames(data[i]) ] <- "P.value"
  }
}
rrblupSNPim <- plyr::ldply(rrblupList, data.frame)
rrblupSNPim$filename  <- paste("rrBLUP124im+36088/")
str(rrblupSNPim)
####import result from RRBLUP for flower 106. 116.122

extension <- "csv"
fileNames <- Sys.glob(paste("rrBLUPF*/rrBLUP_GWAS_results.flo*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample <- sample[rowSums(sample[4:26])!=0,]
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  plot.p.adjusted <- adj_P_function(sample, 27, 49)
  data <- plot.p.adjusted[,c(1:3,27:72)]
  rrblupList = list()
  for(j in 4:26){
    rrBLUP <- data[,c(1:3,j)][data[[j+23]] < 0.05, ]
    if(nrow(rrBLUP)>=1){
      rrblupList[[j]] <- data.frame(rrBLUP, Trait.name = rep(colnames(data[j]), nrow(rrBLUP)))
      colnames(rrblupList[[j]])[colnames(rrblupList[[j]])== colnames(data[j]) ] <- "P.value"
    }
  }
  rrblupSNPfnon <- plyr::ldply(rrblupList, data.frame)
  #rrblupSNPfnon <- plyr::rbind.fill(rrblupList)
  mzList[[i]] = data.frame(rrblupSNPfnon, filename = rep(fileNames[i], nrow(rrblupSNPfnon)))  
}
rrblupSNPfnm <- plyr::ldply(mzList, data.frame)



####import result from OWA
####import result from OWA
extension <- "csv"
#rrBLUPimall <- read.csv("rrBLUPOWA/rrBLUP_GWAS_results.*.csv",row.names = 1)
fileNames <- Sys.glob(paste("rrBLUPimOWA/rrBLUP_GWAS_results.*.", extension, sep = ""))
i <- 1
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample <- sample[rowSums(sample[4:5])!=0,]
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  plot.p.adjusted <- adj_P_function(sample, 6, 7)
  data <- plot.p.adjusted[,c(1:3,6:9)]
  rrblupList = list()
  for(j in 4:5){
    rrBLUP <- data[,c(1:3,j)][data[[j+2]] < 0.05, ]
    if(nrow(rrBLUP)>=1){
      rrblupList[[j]] <- data.frame(rrBLUP, Trait.name = rep(colnames(data[j]), nrow(rrBLUP)))
      colnames(rrblupList[[j]])[colnames(rrblupList[[j]])== colnames(data[j]) ] <- "P.value"
    }
  }
  rrblupSNPfnon <- plyr::ldply(rrblupList, data.frame)
  #rrblupSNPfnon <- plyr::rbind.fill(rrblupList)
  mzList[[i]] = data.frame(rrblupSNPfnon, filename = rep(fileNames[i], nrow(rrblupSNPfnon)))  
}
rrblupSNPOWAall <- plyr::ldply(mzList, data.frame)

####import result from OWA
####import result from OWA
extension <- "csv"
#rrBLUPimall <- read.csv("rrBLUPOWA/rrBLUP_GWAS_results.*.csv",row.names = 1)
fileNames <- Sys.glob(paste("rrBLUPO*/rrBLUP_GWAS_results.*.", extension, sep = ""))
i <- 1
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample <- sample[rowSums(sample[4:5])!=0,]
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  plot.p.adjusted <- adj_P_function(sample, 6, 7)
  data <- plot.p.adjusted[,c(1:3,6:9)]
  rrblupList = list()
  for(j in 4:5){
    rrBLUP <- data[,c(1:3,j)][data[[j+2]] < 0.05, ]
    if(nrow(rrBLUP)>=1){
      rrblupList[[j]] <- data.frame(rrBLUP, Trait.name = rep(colnames(data[j]), nrow(rrBLUP)))
      colnames(rrblupList[[j]])[colnames(rrblupList[[j]])== colnames(data[j]) ] <- "P.value"
    }
  }
  rrblupSNPfnon <- plyr::ldply(rrblupList, data.frame)
  #rrblupSNPfnon <- plyr::rbind.fill(rrblupList)
  mzList[[i]] = data.frame(rrblupSNPfnon, filename = rep(fileNames[i], nrow(rrblupSNPfnon)))  
}
rrblupSNPOWAnm <- plyr::ldply(mzList, data.frame)




####import result from rrblup
####import result from rrblup
extension <- "csv"
#rrBLUPimall <- read.csv("rrBLUPOWA/rrBLUP_GWAS_results.*.csv",row.names = 1)
fileNames <- Sys.glob(paste("rrBLUPSRD/rrBLUP_GWAS_results.*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample <- sample[rowSums(sample[4])!=0,]
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  plot.p.adjusted <- adj_P_function(sample, 5, 5)
  data <- plot.p.adjusted[,c(1:3,5:6)]
  rrblupList = list()
  for(j in 4:4){
    rrBLUP <- data[,c(1:3,j)][data[[j+1]] < 0.05, ]
    if(nrow(rrBLUP)>=1){
      rrblupList[[j]] <- data.frame(rrBLUP, Trait.name = rep(colnames(data[j]), nrow(rrBLUP)))
      colnames(rrblupList[[j]])[colnames(rrblupList[[j]])== colnames(data[j]) ] <- "P.value"
    }
  }
  rrblupSNPfnon <- plyr::ldply(rrblupList, data.frame)
  #rrblupSNPfnon <- plyr::rbind.fill(rrblupList)
  mzList[[i]] = data.frame(rrblupSNPfnon, filename = rep(fileNames[i], nrow(rrblupSNPfnon)))  
}
rrblupSNPSRD <- plyr::ldply(mzList, data.frame)

####import result from Culm
####import result from Culm
extension <- "csv"
#rrBLUPimall <- read.csv("rrBLUPOWA/rrBLUP_GWAS_results.*.csv",row.names = 1)
fileNames <- Sys.glob(paste("rrBLUPC*/rrBLUP_GWAS_results.*.", extension, sep = ""))
mzList = list()
for(i in 1:length(fileNames)){
  sample = read.csv(fileNames[i],row.names = 1)
  sample <- sample[rowSums(sample[,c(4:17)])!=0,]
  names(sample)
  adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
    for(i in start_var:end_var){
      data[ ,paste0("Adj.P.",colnames(data[i]))] <- p.adjust(data[[i]], method ="fdr", n = length(data[[i]]))
    }
    return(data)
  }
  plot.p.adjusted <- adj_P_function(sample, 18, 31)
  data <- plot.p.adjusted[,c(1:3,18:45)]
  names(data)
  rrblupList = list()
  for(j in 4:17){
    rrBLUP <- data[,c(1:3,j)][data[[j+14]] < 0.05, ]
    if(nrow(rrBLUP)>=1){
      rrblupList[[j]] <- data.frame(rrBLUP, Trait.name = rep(colnames(data[j]), nrow(rrBLUP)))
      colnames(rrblupList[[j]])[colnames(rrblupList[[j]])== colnames(data[j]) ] <- "P.value"
    }
  }
  rrblupSNPfnon <- plyr::ldply(rrblupList, data.frame)
  #rrblupSNPfnon <- plyr::rbind.fill(rrblupList)
  mzList[[i]] = data.frame(rrblupSNPfnon, filename = rep(fileNames[i], nrow(rrblupSNPfnon)))  
}
rrblupSNPClumnm <- plyr::ldply(mzList, data.frame)

####combine all of results from rrBLUP together 
rrBLUPalltrait <- plyr::ldply(list(rrblupSNPClumnm,rrblupSNPfnm,rrblupSNPim,rrblupSNPOWAall,rrblupSNPOWAnm,rrblupSNPSRD),data.frame)
str(rrBLUPalltrait)
levels(rrBLUPalltrait$Trait.name)
levels(as.factor(rrBLUPalltrait$filename))

rrBLUPalltrait.1 <- separate(rrBLUPalltrait, Trait.name, c("Method","Trait.name"), "\\.")
str(rrBLUPalltrait.1)

rrBLUPalltrait.1$filename <- gsub('rrBLUP', 'rrBLUP/', rrBLUPalltrait.1$filename)
levels(as.factor(rrBLUPalltrait.1$filename))

rrBLUPalltrait.2 <- separate(rrBLUPalltrait.1, filename, c("Software","Ind.SNP","discard","discard2"), "/")
str(rrBLUPalltrait.2)

levels(as.factor(rrBLUPalltrait.2$Ind.SNP))
rrBLUPalltrait.2$Ind.SNP <- gsub('C124', 'C124+2560', rrBLUPalltrait.2$Ind.SNP)
rrBLUPalltrait.2$Ind.SNP <- gsub('C116', 'C116+3293', rrBLUPalltrait.2$Ind.SNP)
rrBLUPalltrait.2$Ind.SNP <- gsub('C106', 'C106+4202', rrBLUPalltrait.2$Ind.SNP)
rrBLUPalltrait.2$Ind.SNP <- gsub('F106', 'F106+4185', rrBLUPalltrait.2$Ind.SNP)

rrBLUPalltrait.2$Ind.SNP <- gsub("O102", "O102+4322", rrBLUPalltrait.2$Ind.SNP)
rrBLUPalltrait.2$Ind.SNP <- gsub("O112", "O112+3450", rrBLUPalltrait.2$Ind.SNP)
rrBLUPalltrait.2$Ind.SNP <- gsub("O122", "O122+2646", rrBLUPalltrait.2$Ind.SNP)
#rrBLUPalltrait.2$Ind.SNP <- gsub("124+36088im", "C124im+36088", rrBLUPalltrait.2$Ind.SNP)
rrBLUPalltrait.2$Ind.SNP <- gsub("imOWA", "124im+36088", rrBLUPalltrait.2$Ind.SNP)

str(rrBLUPalltrait.2)
levels(as.factor(rrBLUPalltrait.2$Ind.SNP))
colnames(rrBLUPalltrait.2)[colnames(rrBLUPalltrait.2)=="Name"] <- "SNP"
rrBLUPalltrait.3 <- rrBLUPalltrait.2[1:8]
rrBLUPalltrait.3$PS <- paste("PCA=3")
str(rrBLUPalltrait.3)
write.csv(rrBLUPalltrait.3,file="GAPIT0.05Result/rrBLUPalltrait.PCA.csv")





