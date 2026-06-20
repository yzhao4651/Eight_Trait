####GAPIT 
###imported all imputed datasets for GAPIT
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYimputedSNP19_chapter2.csv")
myGD <- read.csv("data/myGDimputedSNP19_chapter2.csv")
myGM <- read.csv("data/myGMimputedSNP19_chapter2.csv")
myQ3<- read.csv("data/myQimputedSNP19_chapter2.csv", row.names = 1)

###imputed SNP for FarmCPU
###imputed SNP for FarmCPU
###imputed SNP for FarmCPU

  install.packages("bigmemory")
  install.packages("biganalytics")
  library("bigmemory")
  library("biganalytics")
  library("bigmemory")
  library("biganalytics")
  library("compiler") #this library is already installed in R
  source("http://zzlab.net/GAPIT/gapit_functions.txt")
  source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")

  str(myY)
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthusupdated result/FarmCPUCV0.05")
out_start=2
out_end=38
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    CV=myQ3,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}  

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthusupdated result/FarmCPU0.05")
out_start=2
out_end=38
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
} 


###Non-missing SNPs for FarmCPU
###Non-missing SNPs for FarmCPU
###Non-missing SNPs for FarmCPU

####ALL of R codes below using Bottom as Supper and changed different top. but without incorporating CV, 
####becasue i saw the example does not use the PCA, but CV in the manual 
##########Using MLM SUPPER with CV
install.packages("bigmemory")
install.packages("biganalytics")
library("bigmemory")
library("biganalytics")
library("bigmemory")
library("biganalytics")
library("compiler") #this library is already installed in R
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYC.106.4202_up.csv")
str(myY)
myGD <- read.csv("data/myGDC.106.4202_up.csv")
myGM <- read.csv("data/myGMC.106.4202_up.csv")
myQ3<- read.csv("data/myQC.106.4202_up.csv",row.names = 1)

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthusupdated result/C106FarmCPUCV0.05")
out_start=2
out_end=14
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    CV=myQ3,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}  

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthusupdated result/C106FarmCPU0.05")
out_start=2
out_end=14
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    method.bin="optimum",
    threshold.output=1,
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
} 


setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYC.116.3293_up.csv")
myGD <- read.csv("data/myGDC.116.3293_up.csv")
myGM <- read.csv("data/myGMC.116.3293_up.csv")
myQ3<- read.csv("data/myQC.116.3293_up.csv",row.names = 1)
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthusupdated result/C116FarmCPUCV0.05")
out_start=2
out_end=14
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    CV=myQ3,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}  

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthusupdated result/C116FarmCPU0.05")
out_start=2
out_end=14
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
} 


setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYC.124.2560_up.csv")
str(myY)
myGD <- read.csv("data/myGDC.124.2560_up.csv")
myGM <- read.csv("data/myGMC.124.2560_up.csv")
myQ3<- read.csv("data/myQC.124.2560_up.csv",row.names = 1)

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/C124FarmCPUCV0.05")
out_start=2
out_end=14
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    CV=myQ3,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
}  

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/C124FarmCPU0.05")
out_start=2
out_end=14
for (i in out_start:out_end){ 
  myFarmCPU <- FarmCPU(
    Y=myY[,c(1,i)],
    GD=myGD,
    GM=myGM,
    threshold.output=1,
    method.bin="optimum",
    MAF.calculate=TRUE,
    maf.threshold=0.05
  )  
} 


######################################################
####GAPIT ############################################
####GAPIT ############################################
######################################################

###imported all the data need for GAPIT

setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/")
myY <- read.csv("data/myYimputedSNP19_chapter2.csv")
#myKI <- read.csv("data/kmatrix.df.csv", head = FALSE)
myGD <- read.csv("data/myGDimputedSNP19_chapter2.csv")
myGM <- read.csv("data/myGMimputedSNP19_chapter2.csv")
myQ<- read.csv("data/myQimputedSNP19_chapter2.csv")

###install the packaged need for GAPIT
source("http://www.bioconductor.org/BiocManager.R")
install.packages("BiocManager")
BiocManager::install("multtest")
install.packages("gplots")
install.packages("LDheatmap")
install.packages("genetics")
install.packages("ape")
install.packages("EMMREML")
install.packages("scatterplot3d") #The downloaded link at: http://cran.r-project.org/package=scatterplot3d
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

###########Using CMLM SUPPER with CV
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT1/CMLMSUPPER")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/CMLMSUPPERCV0.05")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ,
  SNP.MAF =0.05,
  sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1
)
###########Using MLM SUPPER with CV
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT1/MLMSUPPER")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/MLMSUPPERCV0.05")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ,
  SNP.MAF =0.05,
  kinship.cluster=c("average"), kinship.group=c("Mean"), 
  sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1
)

####Enhanced Compression with CV
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT1/ECMMLCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/ECMMLCV0.05")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ,
  SNP.MAF =0.05,
  kinship.cluster=c("average"), kinship.group=c("Mean"), 
  group.from=100,
  group.to=36088, 
  group.by=10
)

###########"GLM","MLM","MLMM" with CV
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT1/GLMMLMMLMMFarmCPUCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/GLMMLMMLMMCV0.05")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ, 
  SNP.MAF =0.05,
  model=c("GLM","MLM","MLMM")
)




###this without CV
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/CMLMSUPPER0.05")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  SNP.MAF =0.05,
  sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1
)

setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT1/MLMSUPPER")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/MLMSUPPER0.05")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  SNP.MAF =0.05,
  sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1
)

####Enhanced Compression without CV
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT1/ECMMLNOCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/ECMML0.05")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  SNP.MAF =0.05,
  kinship.cluster=c("average"), kinship.group=c("Mean"), 
  group.from=100,
  group.to=36088, 
  group.by=10
)



###########"GLM","MLM","MLMM" without CV
setwd("~/Documents/R-corde for miscanthus project/Miscanthus/Result GAPIT1/GLMMLMMLMM0.05")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/GAPIT005/GLMMLMMLMM0.05")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  SNP.MAF =0.05,
  model=c("GLM","MLM","MLMM")
)

####GAPIT 
###imported all the data need for GAPIT
###install the packaged need for GAPIT
source("http://www.bioconductor.org/biocLite.R")
install.packages("BiocManager")
biocLite("multtest")
BiocManager::install("multtest")
install.packages("gplots")
install.packages("LDheatmap")
install.packages("genetics")
install.packages("ape")
install.packages("EMMREML")
install.packages("scatterplot3d") #The downloaded link at: http://cran.r-project.org/package=scatterplot3d
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(ape)
library(EMMREML)
library(compiler) #this library is already installed in R
library("scatterplot3d")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

####Basic Scenario with Model selection in order to find the optimal number of PCs for different traits. 
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYC.106.4202_up.csv")
myGD <- read.csv("data/myGDC.106.4202_up.csv")
myGM <- read.csv("data/myGMC.106.4202_up.csv")
myQ<- read.csv("data/myQC.106.4202_up.csv")
#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/ECMMLNOCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C106ECMMLNOCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
  group.from=100,
  group.to=5000, 
  group.by=10
)
####Enhanced Compression with CV
#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/ECMMLCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C106ECMMLCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ,
  kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
  group.from=100,
  group.to=5000, 
  group.by=10
)

###########"GLM","MLM","MLMM","FarmCPU" without CV
#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/GLMMLMMLMMFarmCPU")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C106GLMMLMMLMMFarmCPU")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  model=c("GLM","MLM","MLMM")
)
###########"GLM","MLM","MLMM","FarmCPU" with CV
#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/GLMMLMMLMMFarmCPUCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C106GLMMLMMLMMFarmCPUCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ, 
  model=c("GLM","MLM","MLMM","FarmCPU")
)


####ALL of R codes below using Bottom as Supper and changed different top. but without incorporating CV, 
####becasue i saw the example does not use the PCA, but CV in the manual 
##########Using MLM SUPPER with CV
#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/MLMSUPPER")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C106MLMSUPPER")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ,
  sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1
)

#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/MLMSUPPER")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C106CMLMSUPPER")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ,
  sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1
)

#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITflo106/MLMSUPPER")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C106CMLMSUPPERNOCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1
)


###########Using CMLM SUPPER with CV
###
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYC.116.3293_up.csv")
myGD <- read.csv("data/myGDC.116.3293_up.csv")
myGM <- read.csv("data/myGMC.116.3293_up.csv")
myQ<- read.csv("data/myQC.116.3293_up.csv")

####Basic Scenario with Model selection in order to find the optimal number of PCs for different traits. 
#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/ECMMLNOCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C116ECMMLNOCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
  group.from=100,
  group.to=5000, 
  group.by=10
)

####Enhanced Compression with CV
#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/ECMMLCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C116ECMMLCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ,
  kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
  group.from=100,
  group.to=5000, 
  group.by=10
)

###########"GLM","MLM","MLMM","FarmCPU" without CV
#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/GLMMLMMLMMFarmCPU")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C116GLMMLMMLMMFarmCPU")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  model=c("GLM","MLM","MLMM","FarmCPU")
)

###########"GLM","MLM","MLMM","FarmCPU" with CV
#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/GLMMLMMLMMFarmCPUCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C116GLMMLMMLMMFarmCPUCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ, 
  model=c("GLM","MLM","MLMM")
)

####ALL of R codes below using Bottom as Supper and changed different top. but without incorporating CV, 
####becasue i saw the example does not use the PCA, but CV in the manual 
##########Using MLM SUPPER with CV
#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/MLMSUPPER")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C116CMLMSUPPER")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ,
  sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1
)

#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITflo106/MLMSUPPER")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C116CMLMSUPPERNOCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1
)



###########Using CMLM SUPPER with CV
###
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus")
myY <- read.csv("data/myYC.124.2560_up.csv")
myGD <- read.csv("data/myGDC.124.2560_up.csv")
myGM <- read.csv("data/myGMC.124.2560_up.csv")
myQ<- read.csv("data/myQC.124.2560_up.csv")

####Basic Scenario with Model selection in order to find the optimal number of PCs for different traits. 
#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/ECMMLNOCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C124ECMMLNOCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
  group.from=100,
  group.to=5000, 
  group.by=10
)

####Enhanced Compression with CV
#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/ECMMLCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C124ECMMLCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ,
  kinship.cluster=c("average", "complete", "ward"), kinship.group=c("Mean", "Max"), 
  group.from=100,
  group.to=5000, 
  group.by=10
)

###########"GLM","MLM","MLMM","FarmCPU" without CV
#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/GLMMLMMLMMFarmCPU")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C124GLMMLMMLMMFarmCPU")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  model=c("GLM","MLM","MLMM","FarmCPU")
)

###########"GLM","MLM","MLMM","FarmCPU" with CV
#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/GLMMLMMLMMFarmCPUCV")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C124GLMMLMMLMMFarmCPUCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ, 
  model=c("GLM","MLM","MLMM")
)

####ALL of R codes below using Bottom as Supper and changed different top. but without incorporating CV, 
####becasue i saw the example does not use the PCA, but CV in the manual 
##########Using MLM SUPPER with CV
#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITCulm106/MLMSUPPER")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C124CMLMSUPPER")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  CV=myQ,
  sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1
)
#setwd("~/Documents/R-corde for miscanthus project/Miscanthus/GAPITflo106/MLMSUPPER")
setwd("C:/Users/Admin/Desktop/Miscanthus/Miscanthus/Result GAPIT1/C124CMLMSUPPERNOCV")
myGAPIT <- GAPIT(
  Y=myY,
  GD=myGD,
  GM=myGM,
  sangwich.top="CMLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  LD=0.1
)
