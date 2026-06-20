####this one is to obtain the imputed SNP, polulation instrcute for GAPIT and rrBLUP software
####this one is to obtain the imputed SNP, polulation instrcute for GAPIT and rrBLUP software

############step 1 emerge all the phenotype data with the index data
############step 1 emerge all the phenotype data with the index data
### data 1 import the index data with number matched Entry of phenotype and name of SNP 
Altable <- read.csv("data/Altable19.csv", header = TRUE, na="")
###check the data format
str(Altable)

###data2 import the  all traits 
alltraitsblup <- read.csv("data/ranefvalueChapter2.csv", stringsAsFactors = FALSE, header = TRUE,row.names = 1)

###check the data format 
str(alltraitsblup)
####combine the traits data with indexdata  
allblup <- plyr::join_all(list(Altable[1:2],alltraitsblup), by='Entry')
###check the data format 
str(allblup)
###rename of the column name of Acession to Taxa 
colnames(allblup)[colnames(allblup)=="Accession"] <- "Taxa"
###remove the colomn of Entry
allblup <- allblup[,-2]
allblup[allblup == "."] <- NA
allblup <- allblup[!(rowSums(is.na(allblup)) == 13), ]

###check the data format 
str(allblup)
###write out the data set 
write.csv(allblup, file = "data/alltraits_chapter2.csv",row.names = T)



############step 1 import phenotypic data:myY 
############step 1 import phenotypic data:myY
allblupClu <- read.csv("data/alltraits_chapter2.csv", na.strings = c("",".","NA"),row.names = 1)
str(allblupClu)

############# step2 
############# step2
############loading the genotype SNP data with missing value called datacomb3####
load("data/160322filteredSNPs.RData")

########## inspect the genotype data
#str(datacomb3) # 594 genotypes, 46,177 markers
#dimnames(datacomb3)[[1]] # this is how the genotypes are named in the SNP data

#############changing the rowname to the first columne###########################
source("Function/subsetNONmissingSNP.GAPIT.R")
SNP <- subsetNONmisingSNP.GAPIT(datacomb3,allblup)

##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 30, 
### Remove the row with 30% missing value (more than 30% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.3),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.106.4202 <- Select.MAF(SNPnmcol)

#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)

### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.106.4202))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.106.4202 <- getmyGD(myGM,myGD.106.4202)

myGM.106.4202 <- getmyGM(myGM,myGD.106.4202)

###check if they are shared the same name in the same order
identical(as.character(myGM.106.4202$Name), colnames(myGD.106.4202))
write.csv(myGM.106.4202, file = "data/myGMC.106.4202_up.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.106.4202 <- FirstColumn(myGD.106.4202)

write.csv(myGD.106.4202, file = "data/myGDC.106.4202_up.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
myY.106.4202 <- allblupClu[match(myGD.106.4202$Taxa,allblupClu$Taxa, nomatch=0),]
myY.106.4202 <- droplevels(myY.106.4202)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])

#### write out the phenotype with correct TAXA
write.csv(myY.106.4202, file = "data/myYC.106.4202_up.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myY.106.4202$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.106.4202 <- myQ[match(myY.106.4202$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.106.4202, file = "data/myQC.106.4202_up.csv",row.names = FALSE)

###0.351
###0.351
###0.351
##############step3 (myGD) select suitable SNP and individual for the next analysis 
##############step3 (myGD) select suitable SNP and individual for the next analysis 
########1)remove row with missing value more than 351, 
### Remove the row with 351% missing value (more than 351% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.336),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.116.3293 <- Select.MAF(SNPnmcol)

#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)#############
myGM <- read.csv("data/myGM.csv",row.names = 1)

### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.116.3293))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.116.3293 <- getmyGD(myGM,myGD.116.3293)

myGM.116.3293 <- getmyGM(myGM,myGD.116.3293)

###check if they are shared the same name in the same order
identical(as.character(myGM.116.3293$Name), colnames(myGD.116.3293))
write.csv(myGM.116.3293, file = "data/myGMC.116.3293_up.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.116.3293 <- FirstColumn(myGD.116.3293)

write.csv(myGD.116.3293, file = "data/myGDC.116.3293_up.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
myY.116.3293 <- allblupClu[match(myGD.116.3293$Taxa,allblupClu$Taxa, nomatch=0),]
myY.116.3293 <- droplevels(myY.116.3293)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])

#### write out the phenotype with correct TAXA
write.csv(myY.116.3293, file = "data/myYC.116.3293_up.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myY.116.3293$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.116.3293 <- myQ[match(myY.116.3293$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.116.3293, file = "data/myQC.116.3293_up.csv",row.names = FALSE)

###0.42
###0.42
###0.42
### Remove the row with 42% missing value (more than 42% NA)
SNPnmrow <- SNP[-which(rowMeans(is.na(SNP)) > 0.42),]
########2) remove column with missng more than 0
### Remove the column with 0 missing value (more than 0 NA)
SNPnmcol <- SNPnmrow[, -which(colMeans(is.na(SNPnmrow)) > 0)]
rownames(SNPnmcol)
SNPnmcol <- SNPnmcol[-c(121),]
rownames(SNPnmcol)
###remove one individauls 
rownames(SNP)

############select the MAF > 0.05
source("Function/SelectMAF-GAPIT.R")
myGD.124.2560 <- Select.MAF(SNPnmcol)
write.csv(myGD.124.2560, file = "data/myGDC.124.2560_up.csv",row.names = FALSE)

#############step4 myGM:Genetic mapping information
#############step4 myGM:Genetic mapping information
########## GM (Genetic mapping dataset)############# 
myGM <- read.csv("data/myGM.csv",row.names = 1)

### Function to get myGM 
# put the SNPs back in order by chromosome and position
identical(as.character(myGM$Name), colnames(myGD.124.2560))
source("Function/getmyGM.R")
source("Function/getmyGD.R")
myGD.124.2560 <- getmyGD(myGM,myGD.124.2560)
myGM.124.2560 <- getmyGM(myGM,myGD.124.2560)

###check if they are shared the same name in the same order
identical(as.character(myGM.124.2560$Name), colnames(myGD.124.2560))
write.csv(myGM.124.2560, file = "data/myGMC.124.2560_up.csv",row.names = FALSE)
source("Function/FirstColumn.R")
myGD.124.2560 <- FirstColumn(myGD.124.2560)
write.csv(myGD.124.2560, file = "data/myGDC.124.2560_up.csv",row.names = FALSE)


##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##############step5 myY: Phenotype data for GAPIT and FarmCPU 
##import the dataset with myY.124.2560 

myY.124.2560 <- allblupClu[match(myGD.124.2560$Taxa,allblupClu$Taxa, nomatch=0),]

myY.124.2560 <- droplevels(myY.124.2560)
######### order the Taxa
#myYorderg <- na.omit(myYorderg[order(myYorderg$Taxa),])

#### write out the phenotype with correct TAXA
write.csv(myY.124.2560, file = "data/myYC.124.2560_up.csv", row.names = FALSE, na = "")

##############step6 myQ: population structure
##############step6 myQ: population structure
########## Q value (population structure)###############################
#####import Q data that Lindsay sent to me 
myQ <- read.csv("data/myQ.csv", stringsAsFactors = FALSE,header = TRUE)
####change the Sample_name to Taxa
colnames(myQ)[colnames(myQ)=="Sample_name"] <- "Taxa"
####check if the name match from both data
myY.124.2560$Taxa %in% myQ$Taxa
###elect myQ dataset using the phenotype data
myQ.124.2560 <- myQ[match(myY.124.2560$Taxa, myQ$Taxa, nomatch=0),]
##check the format
#str(myQ.clum.107.3707)
####write out the dataset
write.csv(myQ.124.2560, file = "data/myQC.124.2560_up.csv",row.names = FALSE)