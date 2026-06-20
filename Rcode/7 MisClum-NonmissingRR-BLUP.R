####rrBLUP method
####import the genotype data with NON missing values SNPs 
####import the genotype data with NON missing values SNPs 
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")

###Culm116
###Culm116
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYC.116.3293_up.csv")
str(myY)
myGD <- read.csv("data/myGDC.116.3293_up.csv", row.names=1)
myGM <- read.csv("data/myGMC.116.3293_up.csv")
myQ<- read.csv("data/myQC.116.3293_up.csv", row.names=1)
str(myQ)
###tansform the genotype data set into rrblup data set
myGDt <- t(myGD)-1
#myGDtdataframe <- data.frame(myGDt)
#head(myGDtdataframe)
### change the row name into the first column 
myGDtdataframe <- data.frame(myGDt)
Name <- rownames(myGDtdataframe)
rownames(myGDtdataframe) <- NULL
myGDtdataframe <- cbind(Name,myGDtdataframe)
###combine the genotype with the chromosome position (dataset: myGM)
genorrblup <- merge(myGM,myGDtdataframe,by="Name")
###change the ID of the name of phenotype to same to the ID of the genotype ( like change "_" to "." )
head(myY)
as.character(myY$Taxa)
myY$Taxa <- gsub("-",".", myY$Taxa)
## myY$Taxa <- make.names(mY$Taxa) # make.names(mY$Taxa) would also work
#head(myY)
#str(myY)
###question: how can incorprated the population structure? 
###question: how can incorprated the population structure? 
###run rrblup
###run rrblup
#install.packages("rrBLUP")
#install.packages("qqman")
### this one is without Kin matrix, 
library(rrBLUP)
library(qqman)
myGDk <- as.matrix(myGD - 1)
dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
#dimnames(myGDk)[[1]] <- make.names(dimnames(myGDk)[[1]])
###get kinship matrix
kmatrix <- A.mat(myGDk,min.MAF=0.05, max.missing=NULL,impute.method="mean",tol=0.02,
                 n.core=1,shrink=FALSE,return.imputed=FALSE)
####do GWAS analysis again with Kmatrix
gwasResults.C116 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                         min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.C116, file = "rrBLUPC116/rrBLUP_GWAS_results.C116.csv")
str(gwasResults.C116)

###this one for qq plots
out_start=4
out_end=17
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPC116",
                      paste("C116 QQ plot of VAR_", names(gwasResults.C116[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.C116[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=17
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPC116",
                      paste("C116 Manhattan plots of _", names(gwasResults.C116[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.C116$Name,
                       CHR = gwasResults.C116$Chromosome,
                       BP = gwasResults.C116$Position,
                       P = 10 ^ -gwasResults.C116[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}

gwasResults.C116 <- read.csv("rrBLUPC116/rrBLUP_GWAS_results.C116.csv",row.names=1)
str(gwasResults.C116)
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("rrBLUP.",colnames(data[i]))] <- 10^-data[i]
  }
  return(data)
}
gwasResults.C116 <- adj_P_function(gwasResults.C116, 4, 17)

write.csv(gwasResults.C116, file = "rrBLUPC116/rrBLUP_GWAS_results.C116.csv")


###C.124.2560
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
myY <- read.csv("data/myYC.124.2560_up.csv")
str(myY)
myGD <- read.csv("data/myGDC.124.2560_up.csv", row.names=1)
myGM <- read.csv("data/myGMC.124.2560_up.csv")
myQ<- read.csv("data/myQC.124.2560_up.csv", row.names=1)
str(myQ)
###tansform the genotype data set into rrblup data set
myGDt <- t(myGD)-1
#myGDtdataframe <- data.frame(myGDt)
#head(myGDtdataframe)
### change the row name into the first column 
myGDtdataframe <- data.frame(myGDt)
Name <- rownames(myGDtdataframe)
rownames(myGDtdataframe) <- NULL
myGDtdataframe <- cbind(Name,myGDtdataframe)
###combine the genotype with the chromosome position (dataset: myGM)
genorrblup <- merge(myGM,myGDtdataframe,by="Name")
###change the ID of the name of phenotype to same to the ID of the genotype ( like change "_" to "." )
head(myY)
as.character(myY$Taxa)
myY$Taxa <- gsub("-",".", myY$Taxa)
## myY$Taxa <- make.names(mY$Taxa) # make.names(mY$Taxa) would also work
#head(myY)
#str(myY)
###question: how can incorprated the population structure? 
###question: how can incorprated the population structure? 
###run rrblup
###run rrblup
#install.packages("rrBLUP")
#install.packages("qqman")
### this one is without Kin matrix, 
library(rrBLUP)
library(qqman)
myGDk <- as.matrix(myGD - 1)
dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
#dimnames(myGDk)[[1]] <- make.names(dimnames(myGDk)[[1]])
###get kinship matrix
kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                 n.core=1,shrink=FALSE,return.imputed=FALSE)
####do GWAS analysis again with Kmatrix
gwasResults.C124 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                         min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.C124, file = "rrBLUPC124/rrBLUP_GWAS_results.C124.csv")
str(gwasResults.C124)

###this one for qq plots
out_start=4
out_end=17
for (i in out_start:out_end){
  mypath <- file.path("rrBLUPC124",
                      paste("C124 QQ plot of VAR_", names(gwasResults.C124[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.C124[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=17
for (i in out_start:out_end){
  mypath <- file.path("rrBLUPC124",
                      paste("C124 Manhattan plots of _", names(gwasResults.C124[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.C124$Name,
                       CHR = gwasResults.C124$Chromosome,
                       BP = gwasResults.C124$Position,
                       P = 10 ^ -gwasResults.C124[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}

gwasResults.C124 <- read.csv("rrBLUPC124/rrBLUP_GWAS_results.C124.csv",row.names=1)
str(gwasResults.C124)
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("rrBLUP.",colnames(data[i]))] <- 10^-data[i]
  }
  return(data)
}
gwasResults.C124 <- adj_P_function(gwasResults.C124, 4, 17)

write.csv(gwasResults.C124, file = "rrBLUPC124/rrBLUP_GWAS_results.C124.csv")


###Culm106.4202
myY <- read.csv("data/myYC.106.4202_up.csv")
str(myY)
myGD <- read.csv("data/myGDC.106.4202_up.csv", row.names=1)
myGM <- read.csv("data/myGMC.106.4202_up.csv")
myQ<- read.csv("data/myQC.106.4202_up.csv", row.names=1)
str(myQ)
###tansform the genotype data set into rrblup data set
myGDt <- t(myGD)-1
#myGDtdataframe <- data.frame(myGDt)
#head(myGDtdataframe)
### change the row name into the first column 
myGDtdataframe <- data.frame(myGDt)
Name <- rownames(myGDtdataframe)
rownames(myGDtdataframe) <- NULL
myGDtdataframe <- cbind(Name,myGDtdataframe)
###combine the genotype with the chromosome position (dataset: myGM)
genorrblup <- merge(myGM,myGDtdataframe,by="Name")
###change the ID of the name of phenotype to same to the ID of the genotype ( like change "_" to "." )
head(myY)
as.character(myY$Taxa)
myY$Taxa <- gsub("-",".", myY$Taxa)
## myY$Taxa <- make.names(mY$Taxa) # make.names(mY$Taxa) would also work
#head(myY)
#str(myY)
###question: how can incorprated the population structure? 
###question: how can incorprated the population structure? 
###run rrblup
###run rrblup
#install.packages("rrBLUP")
#install.packages("qqman")
### this one is without Kin matrix, 
library(rrBLUP)
library(qqman)
myGDk <- as.matrix(myGD - 1)
dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
#dimnames(myGDk)[[1]] <- make.names(dimnames(myGDk)[[1]])
###get kinship matrix
kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                 n.core=1,shrink=FALSE,return.imputed=FALSE)
####do GWAS analysis again with Kmatrix
gwasResults.C106 <- GWAS(myY, genorrblup, fixed=NULL, K=kmatrix, n.PC=3,
                         min.MAF=0.05, n.core=1, P3D=TRUE, plot=FALSE)
write.csv(gwasResults.C106, file = "rrBLUPC106/rrBLUP_GWAS_results.C106.csv")
str(gwasResults.C106)

###this one for qq plots
out_start=4
out_end=17
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPC106",
                      paste("C106 QQ plot of VAR_", names(gwasResults.C106[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  par(mar = c(5,4,1,1))
  qq(10 ^ -gwasResults.C106[[i]])
  dev.off()
}
###this one for Manhattan plots
out_start=4
out_end=17
for (i in out_start:out_end){
  mypath <- file.path("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/rrBLUPC106",
                      paste("C106 Manhattan plots of _", names(gwasResults.C106[i]), ".jpeg", sep = ""))
  jpeg(file=mypath) 
  require(qqman)
  manhattan(data.frame(SNP = gwasResults.C106$Name,
                       CHR = gwasResults.C106$Chromosome,
                       BP = gwasResults.C106$Position,
                       P = 10 ^ -gwasResults.C106[[i]]), genomewideline = FALSE, suggestiveline = FALSE)
  dev.off()
}

gwasResults.C106 <- read.csv("rrBLUPC106/rrBLUP_GWAS_results.C106.csv",row.names=1)
str(gwasResults.C106)
adj_P_function <- function(data,start_var, end_var,na.rm=TRUE){
  for(i in start_var:end_var){
    data[ ,paste0("rrBLUP.",colnames(data[i]))] <- 10^-data[i]
  }
  return(data)
}
gwasResults.C106 <- adj_P_function(gwasResults.C106, 4, 17)
write.csv(gwasResults.C106, file = "rrBLUPC106/rrBLUP_GWAS_results.C106.csv")
