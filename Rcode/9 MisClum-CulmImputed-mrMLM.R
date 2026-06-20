########################################################################
#######preparing imputed datasets and run mrMLM packages#############
########################################################################

###this one has the SNPS with imputed number not int
###this one has the SNPS with imputed number not int

####import all of the genotype data
setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
subgenomrMLMM <- read.csv("data/subgenomrMLMM.csv")
subgenomrMLMM1 <- data.frame(subgenomrMLMM[,c(1,5:157)],row.names = 1)
str(subgenomrMLMM1)
###seperate two parts
subgenomrMLMMtran <- data.frame(t(subgenomrMLMM1))

####import all of the phenotyp data
###Clum trait
mymrMLMMculm <- read.csv("mrMLMM2/myYmrMLMMculm4.csv")
str(mymrMLMMculm)
###get the same
genomymrMLMMculm <-subgenomrMLMMtran[match(mymrMLMMculm$X.phenotype.,rownames(subgenomrMLMMtran), nomatch=0),]
str(genomymrMLMMculm)
###select the MAF > 0.01
source("Function/SelectMAF-mrMLMM.R")
genomymrMLMMculm <- Select.MAF(genomymrMLMMculm)
str(genomymrMLMMculm)
##get the mateched SNP 
genomymrMLMMculm1 <-subgenomrMLMM[1:4][match(genomymrMLMMculm$rn,subgenomrMLMM$rs., nomatch=0),]
colnames(genomymrMLMMculm1)[which(names(genomymrMLMMculm1) == "rs.")] <- "rn"
genomymrMLMMculm <- plyr::join_all(list(genomymrMLMMculm1,genomymrMLMMculm),by="rn")
str(genomymrMLMMculm)
###change the name in order to fit the software requirment
colnames(genomymrMLMMculm)[which(names(genomymrMLMMculm) == "rn")] <- "rs#"
colnames(genomymrMLMMculm)[which(names(genomymrMLMMculm) == "genotype.for.code.1")] <- "genotype for code 1"
###write out the dataset
write.csv(genomymrMLMMculm, file = "data/subgenomrMLMMculm_up.csv", row.names = FALSE, na = "NA")

###import the myQ
myQ<- read.csv("/Users/yonglizhao/DOcuments/Eight_Traits/data/myQimputedSNP19_chapter2.csv")
myQ$Taxa <- make.names(myQ$Taxa)
str(myQ)
###import the phenotype dataset and then get the popluation sturcture matched with the phenotype dataset 
mymrMLMMClum <- read.csv("/Users/yonglizhao/DOcuments/Eight_Traits/mrMLMM2/myYmrMLMMculm4.csv")
str(mymrMLMMClum)
mymrMLMMQClum <- myQ[myQ$Taxa %in% mymrMLMMClum$X.phenotype., ]
str(mymrMLMMQClum)
colnames(mymrMLMMQClum)[which(names(mymrMLMMQClum) == "Taxa")] <- "<Trait>"
str(mymrMLMMQClum)
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste("Q",i-1,sep="")
  }  
  return(x)
}
mymrMLMMQClum <- colRename(mymrMLMMQClum)
write.csv(mymrMLMMQClum,file ="/Users/yonglizhao/DOcuments/Eight_Traits/data/mymrMLMMQculm.csv",row.names = FALSE)

###import 
mymrMLMMQClum <- read.csv("data/mymrMLMMQculm.csv",header=F)
str(mymrMLMMQClum)
##
colnames(mymrMLMMQClum)[which(names(mymrMLMMQClum) == "V1")] <- "<Covariate>"
colRename<-function(x){  
  for(i in 2:ncol(x)){
    colnames(x)[i] <- paste(" ")
  }  
  return(x)
}
mymrMLMMQClum <- colRename(mymrMLMMQClum)
write.csv(mymrMLMMQClum,file="mrMLMM2/mymrMLMMQculm.1.csv",row.names = FALSE)



###this is the function used for run the data 
###without PS
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\data\\subgenomrMLMMculm_up.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm4.csv",
      fileKin=NULL,filePS=NULL,Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=TRUE,
      Plotformat ="jpeg",Resolution="Low", dir= "C:/Users/Admin/Desktop/Miscanthus/Miscanthus/mrMLMM2/IMCulmimPSNO")

### With PS
library("mrMLM")
mrMLM(fileGen="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\data\\subgenomrMLMMculm_up.csv",
      filePhe="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\myYmrMlMMculm4.csv",
      fileKin=NULL,filePS="C:\\Users\\Admin\\Desktop\\Miscanthus\\Miscanthus\\mrMLMM2\\mymrMLMMQculm.1.csv",Genformat="Num",
      method=c("mrMLM","FASTmrMLM","FASTmrEMMA","pLARmEB","pKWmEB","ISIS EM-BLASSO"),
      Likelihood="REML",
      trait=1:14,
      SearchRadius=20,CriLOD=3,SelectVariable=50,
      Bootstrap=FALSE,DrawPlot=TRUE,
      Plotformat ="jpeg",Resolution="Low", dir= "C:/Users/Admin/Desktop/Miscanthus/Miscanthus/mrMLMM2/IMCulmimPS")