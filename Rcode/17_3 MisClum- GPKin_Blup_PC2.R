#####################################################################################
####This script performs genomic prediction using SNP data and phenotypes (BLUPs) ###
####          It uses 10-fold cross-validation repeated 100 times                 ###
###              only using PC1,PC2                                              #### 
#####################################################################################

###import the dataset with phenotype 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
###import the genotype
myGD <- read.csv("data/myGDimputedSNP19.csv",row.names=1)
###change to matrix 
myGD.matrix <- as.matrix(myGD)
GD_comp <- prcomp(myGD.matrix)
###PC screen plot
pdf(paste("PC screen plot of fprinW", 2 ,".pdf",sep="")) 
par(mfrow=c(1,2))
plot(GD_comp)
plot(GD_comp,type="line", main=paste0("PC screen plot of fprinW"), cex.main=0.9)
abline(h=1,lty=3, col="red")
dev.off()
###plot the PC1 and PC1
#plot(GD_comp$x[, 1], GD_comp$x[, 2], main = "PCA", xlab = "PC1", ylab = "PC2")
###get pc1,2 and 3 for each individual
str(myGD)
GD<-cbind(data.frame(dimnames(myGD)[[1]]), data.frame(GD_comp$x[,1:3]))
str(GD)
names(GD)
colnames(GD)[colnames(GD)=="dimnames.myGD...1.."] <- "Taxa"
singlesiteBLUPs <- read.csv("data/myYimputedSNP19updated.csv", header = TRUE)
str(singlesiteBLUPs)
## combine the data set with PCA
Phen.PC3 <- plyr::join(singlesiteBLUPs,GD,by="Taxa")
str(Phen.PC3)
#install.packages("lme4","Matrix")
ranefvalue <- function(out_start,out_end,y){
  require(lme4)
  require(Matrix)
  out_nvar=out_end-out_start+1
  out_variable = colnames(y[out_start:out_end])
  for (i in out_start:out_end){
    outcome = colnames(y)[i]
    model <- lm(get(outcome)~1+PC1+PC2,
                na.action = na.exclude,
                data=y)
    beta <- residuals(model)
    y[ ,paste0("GR",colnames(y[i]))] <-  beta
  }
  return(y)
}

residual <- ranefvalue(2,39,Phen.PC3)
str(residual)
singlesiteBLUPs <- data.frame(residual[,-c(2:42)],row.names = 1)
str(singlesiteBLUPs)
###calculte the K matrix 
# trait to try out
thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPR2/")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPRR2/")
#looping for trait#
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:38){
  thisTrait <- traitlist[k,]
  
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))
  ###trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]
  myGDk <- as.matrix(myGD1 - 1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
  ###get kinship matrix
  library(rrBLUP)
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    
    scrambleInd <- sample(theseInd) # random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
    nPerRep <- floor(length(theseInd)/thisKfold)
    predValues <- numeric(length(theseInd)) # to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]
      pred <- scrambleInd[firstind:lastind]
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1])
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }
    
    cor(singlesiteBLUPs[theseInd,thisTrait], predValues, use = "complete") # R-squared value (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    Rsquare=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Rsquare)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Rsquare",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
}

###import the dataset with phenotype 
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
###import the genotype
myGD <- read.csv("data/myGDimputedSNP19OWA2.csv",row.names=1)
###change to matrix 
myGD.matrix <- as.matrix(myGD)
GD_comp <- prcomp(myGD.matrix)
###PC screen plot
pdf(paste("PC screen plot of fprinW", 2 ,".pdf",sep="")) 
par(mfrow=c(1,2))
plot(GD_comp)
plot(GD_comp,type="line", main=paste0("PC screen plot of fprinW"), cex.main=0.9)
abline(h=1,lty=3, col="red")
dev.off()
###plot the PC1 and PC1
#plot(GD_comp$x[, 1], GD_comp$x[, 2], main = "PCA", xlab = "PC1", ylab = "PC2")
###get pc1,2 and 3 for each individual
str(myGD)
GD<-cbind(data.frame(dimnames(myGD)[[1]]), data.frame(GD_comp$x[,1:3]))
str(GD)
names(GD)
colnames(GD)[colnames(GD)=="dimnames.myGD...1.."] <- "Taxa"
singlesiteBLUPs <- read.csv("data/myYimputedSNP19OWA2.csv", header = TRUE)
str(singlesiteBLUPs)
## combine the data set with PCA
Phen.PC3 <- plyr::join(singlesiteBLUPs,GD,by="Taxa")
str(Phen.PC3)
#install.packages("lme4","Matrix")
ranefvalue <- function(out_start,out_end,y){
  require(lme4)
  require(Matrix)
  out_nvar=out_end-out_start+1
  out_variable = colnames(y[out_start:out_end])
  for (i in out_start:out_end){
    outcome = colnames(y)[i]
    model <- lm(get(outcome)~1+PC1+PC2,
                na.action = na.exclude,
                data=y)
    beta <- residuals(model)
    y[ ,paste0("GR",colnames(y[i]))] <-  beta
  }
  return(y)
}

residual <- ranefvalue(2,3,Phen.PC3)
str(residual)
singlesiteBLUPs <- data.frame(residual[,-c(2:6)],row.names = 1)
str(singlesiteBLUPs)
###calculte the K matrix 
# trait to try out
thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPR2/")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPRR2/")
#looping for trait#
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:2){
  thisTrait <- traitlist[k,]
  
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))
  ###trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]
  myGDk <- as.matrix(myGD1 - 1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
  ###get kinship matrix
  library(rrBLUP)
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    
    scrambleInd <- sample(theseInd) # random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
    nPerRep <- floor(length(theseInd)/thisKfold)
    predValues <- numeric(length(theseInd)) # to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]
      pred <- scrambleInd[firstind:lastind]
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1])
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }
    
    cor(singlesiteBLUPs[theseInd,thisTrait], predValues, use = "complete") # R-squared value (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    Rsquare=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Rsquare)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Rsquare",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
}

###C.106.4202
###this is for culm combination nonmissing
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYC.106.4202.csv", header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDC.106.4202.csv",row.names=1)
myGD.matrix <- as.matrix(myGD)
GD_comp <- prcomp(myGD.matrix)
###PC screen plot
pdf(paste("PC screen plot of fprinW", 2 ,".pdf",sep="")) 
par(mfrow=c(1,2))
plot(GD_comp)
plot(GD_comp,type="line", main=paste0("PC screen plot of fprinW"), cex.main=0.9)
abline(h=1,lty=3, col="red")
dev.off()
###plot the PC1 and PC1
#plot(GD_comp$x[, 1], GD_comp$x[, 2], main = "PCA", xlab = "PC1", ylab = "PC2")
###get pc1,2 and 3 for each individual
str(myGD)
GD<-cbind(data.frame(dimnames(myGD)[[1]]), data.frame(GD_comp$x[,1:3]))
str(GD)
names(GD)
colnames(GD)[colnames(GD)=="dimnames.myGD...1.."] <- "Taxa"
## combine the data set with PCA
Phen.PC3 <- plyr::join(singlesiteBLUPs,GD,by="Taxa")
str(Phen.PC3)
#install.packages("lme4","Matrix")
ranefvalue <- function(out_start,out_end,y){
  require(lme4)
  require(Matrix)
  out_nvar=out_end-out_start+1
  out_variable = colnames(y[out_start:out_end])
  for (i in out_start:out_end){
    outcome = colnames(y)[i]
    model <- lm(get(outcome)~1+PC1+PC2,
                na.action = na.exclude,
                data=y)
    beta <- residuals(model)
    y[ ,paste0("GR",colnames(y[i]))] <-  beta
  }
  return(y)
}

residual <- ranefvalue(2,15,Phen.PC3)
str(residual)
singlesiteBLUPs <- data.frame(residual[,-c(2:18)],row.names = 1)
str(singlesiteBLUPs)
###calculte the K matrix 
# trait to try out
#thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPC106R2/")
#looping for trait#
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPRC106R2/")
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:14){
  thisTrait <- traitlist[k,]
  
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))
  ###trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]
  myGDk <- as.matrix(myGD1 - 1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
  ###get kinship matrix
  library(rrBLUP)
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    
    scrambleInd <- sample(theseInd) # random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
    nPerRep <- floor(length(theseInd)/thisKfold)
    predValues <- numeric(length(theseInd)) # to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]
      pred <- scrambleInd[firstind:lastind]
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1])
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }
    
    cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    Rsquare=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Rsquare)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Rsquare",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
} 

###C.116.3293
###C.116.3293
###this is for culm combination nonmissing
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYC.116.3293.csv", header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDC.116.3293.csv",row.names=1)
myGD.matrix <- as.matrix(myGD)
GD_comp <- prcomp(myGD.matrix)
###PC screen plot
pdf(paste("PC screen plot of fprinW", 2 ,".pdf",sep="")) 
par(mfrow=c(1,2))
plot(GD_comp)
plot(GD_comp,type="line", main=paste0("PC screen plot of fprinW"), cex.main=0.9)
abline(h=1,lty=3, col="red")
dev.off()
###plot the PC1 and PC1
#plot(GD_comp$x[, 1], GD_comp$x[, 2], main = "PCA", xlab = "PC1", ylab = "PC2")
###get pc1,2 and 3 for each individual
str(myGD)
GD<-cbind(data.frame(dimnames(myGD)[[1]]), data.frame(GD_comp$x[,1:3]))
str(GD)
names(GD)
colnames(GD)[colnames(GD)=="dimnames.myGD...1.."] <- "Taxa"
## combine the data set with PCA
Phen.PC3 <- plyr::join(singlesiteBLUPs,GD,by="Taxa")
str(Phen.PC3)
#install.packages("lme4","Matrix")
ranefvalue <- function(out_start,out_end,y){
  require(lme4)
  require(Matrix)
  out_nvar=out_end-out_start+1
  out_variable = colnames(y[out_start:out_end])
  for (i in out_start:out_end){
    outcome = colnames(y)[i]
    model <- lm(get(outcome)~1+PC1+PC2,
                na.action = na.exclude,
                data=y)
    beta <- residuals(model)
    y[ ,paste0("GR",colnames(y[i]))] <-  beta
  }
  return(y)
}

residual <- ranefvalue(2,15,Phen.PC3)
str(residual)
singlesiteBLUPs <- data.frame(residual[,-c(2:18)],row.names = 1)
str(singlesiteBLUPs)


###calculte the K matrix 
# trait to try out
#thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPC116R2/")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPRC116R2/")
#looping for trait#
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:14){
  thisTrait <- traitlist[k,]
  
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))
  ###trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]
  myGDk <- as.matrix(myGD1 - 1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
  ###get kinship matrix
  library(rrBLUP)
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    
    scrambleInd <- sample(theseInd) # random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
    nPerRep <- floor(length(theseInd)/thisKfold)
    predValues <- numeric(length(theseInd)) # to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]
      pred <- scrambleInd[firstind:lastind]
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1])
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }
    
    cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    Rsquare=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Rsquare)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  
  write.csv(output.cor,file=paste("Rsquare",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
} 

###C.124.2560
###C.124.2560
###this is for flo combination nonmissing
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYC.124.2560.csv", header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDC.124.2560.csv",row.names=1)
myGD.matrix <- as.matrix(myGD)
GD_comp <- prcomp(myGD.matrix)
###PC screen plot
pdf(paste("PC screen plot of fprinW", 2 ,".pdf",sep="")) 
par(mfrow=c(1,2))
plot(GD_comp)
plot(GD_comp,type="line", main=paste0("PC screen plot of fprinW"), cex.main=0.9)
abline(h=1,lty=3, col="red")
dev.off()
###plot the PC1 and PC1
#plot(GD_comp$x[, 1], GD_comp$x[, 2], main = "PCA", xlab = "PC1", ylab = "PC2")
###get pc1,2 and 3 for each individual
str(myGD)
GD<-cbind(data.frame(dimnames(myGD)[[1]]), data.frame(GD_comp$x[,1:3]))
str(GD)
names(GD)
colnames(GD)[colnames(GD)=="dimnames.myGD...1.."] <- "Taxa"
## combine the data set with PCA
Phen.PC3 <- plyr::join(singlesiteBLUPs,GD,by="Taxa")
str(Phen.PC3)
#install.packages("lme4","Matrix")
ranefvalue <- function(out_start,out_end,y){
  require(lme4)
  require(Matrix)
  out_nvar=out_end-out_start+1
  out_variable = colnames(y[out_start:out_end])
  for (i in out_start:out_end){
    outcome = colnames(y)[i]
    model <- lm(get(outcome)~1+PC1+PC2,
                na.action = na.exclude,
                data=y)
    beta <- residuals(model)
    y[ ,paste0("GR",colnames(y[i]))] <-  beta
  }
  return(y)
}
residual <- ranefvalue(2,15,Phen.PC3)
str(residual)
singlesiteBLUPs <- data.frame(residual[,-c(2:18)],row.names = 1)
str(singlesiteBLUPs)
###calculte the K matrix 
# trait to try out
#thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPC125R2/")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPRC124R2/")
#looping for trait#
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:14){
  thisTrait <- traitlist[k,]
  
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))
  ###trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]
  myGDk <- as.matrix(myGD1 - 1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])
  ###get kinship matrix
  library(rrBLUP)
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    
    scrambleInd <- sample(theseInd) # random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold # you can save some processing time by pre-specifying the length of a vector or list
    nPerRep <- floor(length(theseInd)/thisKfold)
    predValues <- numeric(length(theseInd)) # to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]
      pred <- scrambleInd[firstind:lastind]
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1])
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }
    
    cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    Rsquare=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Rsquare)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  
  write.csv(output.cor,file=paste("Rsquare",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
} 

