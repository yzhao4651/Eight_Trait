#####################################################################################
####This script performs genomic prediction using SNP data and phenotypes (BLUPs) ###
####          It uses 10-fold cross-validation repeated 100 times                 ###
###              without using any PC                                            #### 
#####################################################################################
#load("data/160324EMimputedSNP_Msi.RData")
#names(myA.EM.Msi) # $A is the kinship matrix
#singlesiteBLUPs <- read.csv("data/170106BLUPs_to_redo_for_GS_and_GWAS.csv", row.names = 1, header = TRUE)

###import the dataset with phenotype this one include all of the traits
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYimputedSNP19updated.csv",row.names = 1, header = TRUE)#phenotype data (each column = a trait)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDimputedSNP19.csv",row.names=1)#genotype matrix (SNP data)
###calculte the K matrix 
# trait to try out
thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GP/")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPPC0/")
#looping for trait#
library(rrBLUP)
traitlist<-as.matrix(colnames(singlesiteBLUPs))
for (k in 1:38){
  thisTrait <- traitlist[k,]
  # what individuals should be used for this trait
  theseInd <- which(!is.na(singlesiteBLUPs[[thisTrait]]))#Only individuals with non-missing phenotype values for the current trait are used.
  # trying to get the matrix for each trait
  myGD1 <- myGD[dimnames(myGD)[[1]]%in%dimnames(singlesiteBLUPs)[[1]][theseInd],]#Match genotype data to selected individuals
  myGDk <- as.matrix(myGD1 - 1) #Recode SNP values (e.g., 0/1/2 → -1/0/1)
  dimnames(myGDk)[[1]]<- gsub("-",".",dimnames(myGDk)[[1]])#change the name with "_" to "."
  ###get kinship matrix
  library(rrBLUP)
  #A kinship matrix (K matrix) is computed using SNP data, representing genetic relatedness among individuals.
  kmatrix <- A.mat(myGDk,min.MAF=0.05,max.missing=NULL,impute.method="mean",tol=0.02,
                   n.core=1,shrink=FALSE,return.imputed=FALSE)
  
  #kmatrix.df <- as.data.frame(kmatrix)
  dimnames(kmatrix)[[1]]<- gsub("\\.","-",dimnames(kmatrix)[[1]])
  
  output.cor<-NULL
  output.PredValues<-NULL#
  output.seednumber<-NULL#
  ##Each trait is evaluated 100 times with different random splits for stability
  for (j in 1:100){
    ## run genomic prediction (can loop to do this 100 or more times)
    ## 10-Fold Cross-Validation
    thisKfold <- 10
    seednumber<- sample(-1000000:1000000,1)#
    set.seed(seednumber)
    scrambleInd <- sample(theseInd) ##random individual order for this iteration
    genPredOutput <- list()
    length(genPredOutput) <- thisKfold ##you can save some processing time by pre-specifying the length of a vector or list
    ##Individuals are randomly shuffled and divided into 10 folds.
    nPerRep <- floor(length(theseInd)/thisKfold) ##floor function ensures each fold has an integer size
    predValues <- numeric(length(theseInd)) ##to store breeding values for each individual
    names(predValues) <- dimnames(singlesiteBLUPs)[[1]][theseInd]
    for(i in 1:thisKfold){
      # identify individuals for training and prediction sets
      firstind <- (i-1) * nPerRep + 1
      if(i == thisKfold){
        lastind <- length(theseInd)
      } else {
        lastind <- i * nPerRep
      }
      train <- scrambleInd[-(firstind:lastind)]#train: training set
      pred <- scrambleInd[firstind:lastind]#  pred: validation (test) set
      
      # phenotypes for training set only
      thisphen <- rep(NA, dim(singlesiteBLUPs)[1]) #Training individuals → phenotype available,# Testing individuals → set to NA
      thisphen[train] <- singlesiteBLUPs[train,thisTrait]
      thisphen <- thisphen[theseInd]
      genPredOutput[[i]] <- kin.blup(data = data.frame(pheno = thisphen, geno = row.names(singlesiteBLUPs)[theseInd]),
                                     geno = "geno", pheno = "pheno", GAUSS = FALSE,
                                     K = as.matrix(kmatrix))
      predValues[match(pred, theseInd)] <- genPredOutput[[i]]$g[match(pred, theseInd)]
    }

    cor(singlesiteBLUPs[theseInd,thisTrait], predValues, use = "complete") # Correlation between observed and predicted values (prediction accuracy)
    
    pdf(paste(thisTrait,j,"Obs.vs.Pred.pdf"))
    #
    plot(singlesiteBLUPs[theseInd,thisTrait], predValues)
    #
    dev.off()#
    # R-squared value (prediction accuracy) #Correlation between observed and predicted values:
    #Used as prediction accuracy
    Correlation=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) 
  
    output.cor<-c(output.cor,Correlation)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Correlation",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
}


####C.106.4202
####C.106.4202
###this is for flo combination nonmissing
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYC.106.4202_up.csv",row.names = 1, header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDC.106.4202_up.csv",row.names=1)
###calculte the K matrix 
# trait to try out
#thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPC106/")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPRC106R0/")
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
    Correlation=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Correlation)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Correlation",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
} 


####C.116.3293
####C.116.3293
###this is for flo combination nonmissing
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYC.116.3293_up.csv",row.names = 1, header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDC.116.3293_up.csv",row.names=1)

###calculte the K matrix 
# trait to try out
#thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPC116/")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPRC116R0/")
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
    Correlation=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Correlation)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Correlation",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
} 

#####C.124.2560
#####C.124.2560
###this is for flo combination nonmissing
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
singlesiteBLUPs <- read.csv("data/myYC.124.2560_up.csv",row.names = 1, header = TRUE)
str(singlesiteBLUPs)## checking the datasets
###import the genotype
myGD <- read.csv("data/myGDC.124.2560_up.csv",row.names=1)
###calculte the K matrix 
# trait to try out
#thisTrait <- "OWA" # could replace this with a for loop to go through all traits in file
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPC125/")
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/GPRC124R0/")
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
    Correlation=cor(singlesiteBLUPs[theseInd,thisTrait], predValues) # R-squared value (prediction accuracy)
    #
    output.cor<-c(output.cor,Correlation)
    output.PredValues=cbind(output.PredValues,predValues)#
    output.seednumber<-c(output.seednumber,seednumber)
    #
  }#
  write.csv(output.cor,file=paste("Correlation",thisTrait,".csv"))
  write.csv(output.PredValues,file=paste("predValues",thisTrait,".csv"))#
  write.csv(output.seednumber,file=paste("seednumber",thisTrait,".csv"))
} 


