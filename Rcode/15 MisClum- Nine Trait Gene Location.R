###########################################################################
###This one used for searching the genes around assoicated SNPs############
###########################################################################

getwd()
setwd("/Users/yonglizhao/Documents/Eight_Traits")
Charpter2 <- read.csv("data/Chapter2.unique.updated.0227.22.csv")

###install gfortran 8.2 since Matrix could not not compatible for this version macOS 10.13

install.packages("Matrix")### choose not to update R 

###intall the packages if using the higher version of macOS 10.13
#install.packages("BiocManager")
#BiocManager::install("rtracklayer")
#BiocManager::install("GenomicFeatures")
#BiocManager::install("IRanges")
#BiocManager::install("AnnotationDbi")
#BiocManager::install("GenomicFeatures")
#BiocManager::install("GenomicRanges")
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("biocLite")
#BiocManager::install("Repitools")

# Then install this packages if the computer version is macOS 10.13
remotes::install_version("survival", version = "3.4-0")
remotes::install_version("gsmoothr", version = "0.1.7")
install.packages("Rsolnp", type = "binary")
BiocManager::install(c(
  "AnnotationDbi",
  "rtracklayer",
  "GenomicFeatures",
  "IRanges",
  "GenomicRanges",
  "Repitools"
), ask = FALSE, update = FALSE)


library(Matrix)
library(AnnotationDbi)
library(rtracklayer)
library(GenomicFeatures)
library(IRanges)
library(GenomicRanges)
library(Repitools)
library(survival)
####"Sbicolor_454_v3.1.1.gene.gff3" downloaded and then saved in the folder
txdb <- makeTxDbFromGFF("Sbicolor_454_v3.1.1.gene.gff3",
                        format = "gff3",
                        dataSource = "Phytosome 12",
                        organism = "Sorghum bicolor")

###this one is to get the gene that the SNP hit 
mySNPs <- Charpter2
#source("Rcodes/Function/Findneargenes.R")
Findneargenes <- function(mySNPs, search_radius){
  saved_genes <- list()
  length(saved_genes) <- nrow(mySNPs)
  names(saved_genes) <- as.character(mySNPs$SNP)
  chromnames <- paste("Chr", formatC(mySNPs$Chromosome, width = 2, flag = '0'), sep = "")
  for(i in 1:nrow(mySNPs)){
    gr <- GRanges(chromnames[i],
                  IRanges(mySNPs$Position[i] - search_radius,
                          mySNPs$Position[i] + search_radius))
    mygenes <- transcriptsByOverlaps(txdb, gr)
    saved_genes[[i]] <- mygenes$tx_name
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        Gene = unlist(saved_genes))
  return(gene_df)
}

gene_df.genes <- Findneargenes(Charpter2, 1e3)
str(gene_df.genes)
nlevels(as.factor(gene_df.genes$Gene))
#seqinfo(txdb)
#To save the txdb database for later and avoid having to recreate it every time we use it, we can use saveDb() and, later, loadDb()
#BiocManager::install("AnnotationDbi")
#library(AnnotationDbi)
#saveDb(txdb, 'txdb.Sbicolor_454_v3.1.1')
#Once the database is stored to disk, it can be reloaded later, almost instantly.
#library(AnnotationDbi)
#txdb = loadDb(file = 'txdb.Sbicolor_454_v3.1.1')
#genes(txdb)
#seqinfo(txdb)
###install.packags library(IRanges)
#BiocManager::install("IRanges")
### fucntion to get the neargenes 
#setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/")
### this one for flowering traits
#full_data.flo <- read.csv("data/flo-106.116.122.106im.csv")
#str(full_data.flo)
#levels(full_data.flo$SNP)
#nlevels(full_data.flo$SNP)
##this function need input dataset and search_radius 
#source("https://bioconductor.org/biocLite.R")
#biocLite("Repitools")
library(Repitools)
columns(txdb)

###this one is to get which chromosome that gene located that the SNP hit 
###this one is to get which chromosome that gene located that the SNP hit 
###this one is to get which chromosome that gene located that the SNP hit 
mySNPs <- Charpter2
#source("Function/Findneargenes.R")
Findneargenes <- function(mySNPs, search_radius){
  saved_genes <- list()
  length(saved_genes) <- nrow(mySNPs)
  names(saved_genes) <- as.character(mySNPs$SNP)
  chromnames <- paste("Chr", formatC(mySNPs$Chromosome, width = 2, flag = '0'), sep = "")
  for(i in 1:nrow(mySNPs)){
    gr <- GRanges(chromnames[i],
                  IRanges(mySNPs$Position[i] - search_radius,
                          mySNPs$Position[i] + search_radius))
    mygenes <- transcriptsByOverlaps(txdb, gr)
    t <- annoGR2DF(mygenes) 
    saved_genes[[i]] <- t$chr
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        chr = unlist(saved_genes))
  return(gene_df)
}


gene_df.chr <- Findneargenes(Charpter2, 1e3)
str(gene_df.chr)


###this one is to get that start point that gene located that the SNP hit 
###this one is to get that start point that gene located that the SNP hit 
###this one is to get that start point that gene located that the SNP hit 
Findneargenes <- function(mySNPs, search_radius){
  saved_genes <- list()
  length(saved_genes) <- nrow(mySNPs)
  names(saved_genes) <- as.character(mySNPs$SNP)
  chromnames <- paste("Chr", formatC(mySNPs$Chromosome, width = 2, flag = '0'), sep = "")
  for(i in 1:nrow(mySNPs)){
    gr <- GRanges(chromnames[i],
                  IRanges(mySNPs$Position[i] - search_radius,
                          mySNPs$Position[i] + search_radius))
    mygenes <- transcriptsByOverlaps(txdb, gr)
    
    t <- annoGR2DF(mygenes)
    saved_genes[[i]] <- t$start
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        start = unlist(saved_genes))
  return(gene_df)
}

gene_df.start <- Findneargenes(Charpter2, 1e3)
str(gene_df.start)

###this one is to get that end point that gene located that the SNP hit 
###this one is to get that end point that gene located that the SNP hit
###this one is to get that end point that gene located that the SNP hit
Findneargenes <- function(mySNPs, search_radius){
  saved_genes <- list()
  length(saved_genes) <- nrow(mySNPs)
  names(saved_genes) <- as.character(mySNPs$SNP)
  chromnames <- paste("Chr", formatC(mySNPs$Chromosome, width = 2, flag = '0'), sep = "")
  for(i in 1:nrow(mySNPs)){
    gr <- GRanges(chromnames[i],
                  IRanges(mySNPs$Position[i] - search_radius,
                          mySNPs$Position[i] + search_radius))
    mygenes <- transcriptsByOverlaps(txdb, gr)
    
    
    t <- annoGR2DF(mygenes)
    saved_genes[[i]] <- t$end
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        end = unlist(saved_genes))
  return(gene_df)
}

gene_df.end <- Findneargenes(Charpter2, 1e3)
str(gene_df.end)

###this one is to get that width of the gene located that the SNP hit
###this one is to get that width of the gene located that the SNP hit
###this one is to get that width of the gene located that the SNP hit
Findneargenes <- function(mySNPs, search_radius){
  saved_genes <- list()
  length(saved_genes) <- nrow(mySNPs)
  names(saved_genes) <- as.character(mySNPs$SNP)
  chromnames <- paste("Chr", formatC(mySNPs$Chromosome, width = 2, flag = '0'), sep = "")
  for(i in 1:nrow(mySNPs)){
    gr <- GRanges(chromnames[i],
                  IRanges(mySNPs$Position[i] - search_radius,
                          mySNPs$Position[i] + search_radius))
    mygenes <- transcriptsByOverlaps(txdb, gr)
    
    
    t <- annoGR2DF(mygenes)
    saved_genes[[i]] <- t$width
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        width = unlist(saved_genes))
  return(gene_df)
}

gene_df.width <- Findneargenes(Charpter2, 1e3)
str(gene_df.width)

###this one is to get that strand of the gene located that the SNP hit
Findneargenes <- function(mySNPs, search_radius){
  saved_genes <- list()
  length(saved_genes) <- nrow(mySNPs)
  names(saved_genes) <- as.character(mySNPs$SNP)
  chromnames <- paste("Chr", formatC(mySNPs$Chromosome, width = 2, flag = '0'), sep = "")
  for(i in 1:nrow(mySNPs)){
    gr <- GRanges(chromnames[i],
                  IRanges(mySNPs$Position[i] - search_radius,
                          mySNPs$Position[i] + search_radius))
    mygenes <- transcriptsByOverlaps(txdb, gr)
    
    
    t <- annoGR2DF(mygenes)
    saved_genes[[i]] <- t$strand
  }
  gene_df <- data.frame(SNP = rep(mySNPs$SNP, 
                                  times = sapply(saved_genes, length)),
                        strand = unlist(saved_genes))
  return(gene_df)
}

gene_df.strand <- Findneargenes(Charpter2, 1e3)
str(gene_df.strand)

###combine all of them together 
gene_df12.27 <- do.call("cbind",list(gene_df.genes,gene_df.chr, gene_df.start, gene_df.end, gene_df.width, gene_df.strand))
str(gene_df12.27)

gene_df12.27.a <- gene_df12.27[,c(1,2,4,6,8,10,12)]
gene_df12.27.a <- droplevels(gene_df12.27.a)

str(gene_df12.27.a)

nlevels(as.factor(gene_df12.27.a$SNP)) ##343
nlevels(as.factor(gene_df12.27.a$Gene))###576

write.csv(gene_df12.27.a,file="data/AHG.Chapter2UP.T.hitandnohit.csv")

##merge all of data with unique

Chapter2UP.T <- merge(Charpter2, gene_df12.27.a,by="SNP",all.x=TRUE)
Chapter2UP.T <- droplevels(Chapter2UP.T)
#nlevels(char2.maf.t.12.27$SNP)##1857
nlevels(as.factor(Chapter2UP.T$Gene))##576
nlevels(as.factor(Chapter2UP.T$SNP))###546
###get all of the no any hit SNP
Chapter2UP.T.nohhit <- Chapter2UP.T[is.na(Chapter2UP.T$Gene),]
Chapter2UP.T.nohhit <- droplevels(Chapter2UP.T.nohhit)
nlevels(as.factor(Chapter2UP.T.nohhit$SNP)) ###203

write.csv(Chapter2UP.T.nohhit,file="data/Chapter2UP.T.nohit.csv")

###all of them with hit gene
Chapter2UP.T.hit <- Chapter2UP.T[!(is.na(Chapter2UP.T$Gene)),]
str(Chapter2UP.T.hit)
Chapter2UP.T.hit$dis_start.positon <- Chapter2UP.T.hit$start - Chapter2UP.T.hit$Position
Chapter2UP.T.hit$dis_positon.end <- Chapter2UP.T.hit$Position- Chapter2UP.T.hit$end
Chapter2UP.T.hit <- droplevels(Chapter2UP.T.hit)
nlevels(as.factor(Chapter2UP.T.hit$SNP)) ###343
nlevels(as.factor(Chapter2UP.T.hit$Gene))###576

write.csv(Chapter2UP.T.hit,file="data/Chapter2UP.T.hit.csv")

###and then get gene anotation 

annotations <- read.delim("Sbicolor_454_v3_1_1_annotation_info.txt",
                          stringsAsFactors = FALSE)
head(annotations)
str(annotations)
###checking if the gene matched each other 
#Chapter2UP.T.hit$Gene %in% annotations$transcriptName
myrows <- match(Chapter2UP.T.hit$Gene, annotations$transcriptName)
allhitsgenes.char2 <- annotations[myrows,]
head(allhitsgenes.char2)
###combine two dataset in order to get the SNP with gene anotation in one file. 
AHG.char2.1 <- cbind(Chapter2UP.T.hit,allhitsgenes.char2)
colnames(AHG.char2.1)[colnames(AHG.char2.1)=="Gene"] <- "Gene"


install.packages("dplyr")
remotes::install_version("ragg", version = "1.2.5")
install.packages("tidyverse")
library(tidyverse)
library(dplyr)
##seperate the Gene with number and then delete the duplicated genes
AHG.char2.2 <- separate(AHG.char2.1, Gene, c("Gene1","Gene2","Gene.num"), "\\.") 
AHG.char2.3 <- unite(AHG.char2.2,"SNP.Gene",SNP,Gene1,Gene2,sep="/")
AHG.char2.4<- AHG.char2.3[!(duplicated(AHG.char2.3$SNP.Gene)),]
AHG.char2.5 <- separate(AHG.char2.4, SNP.Gene, c("SNP","Gene1","Gene2"), "/")
AHG.char2.5 <- unite(AHG.char2.5, Gene, c("Gene1","Gene2"), sep=".")
AHG.char2.5 <- droplevels(AHG.char2.5)
nlevels(as.factor(AHG.char2.5$SNP)) ###343
nlevels(as.factor(AHG.char2.5$Gene))###395

write.csv(AHG.char2.5,file="data/AHG.char2.12.27.anota.csv")
###select the Genes without studying: which means the X.pac is NA 
char2.hit.gene <- read.csv("data/AHG.char2.12.27.anota.csv")
str(char2.hit.gene)
nlevels(as.factor(char2.hit.gene$Gene))##395

###this one works even this is hit with sobic genes, but many of them has no any study  
char2.NOX.pacId <- char2.hit.gene[is.na(char2.hit.gene$X.pacId),]
char2.NOX.pacId <- droplevels(char2.NOX.pacId)
nlevels(as.factor(char2.NOX.pacId$SNP))###49
nlevels(as.factor(char2.NOX.pacId$Gene))###56
write.csv(char2.NOX.pacId,file="data/AHG.char2.NOX.pacId.csv")

## with some hit 
char2.X.pacId <- char2.hit.gene[!is.na(char2.hit.gene$X.pacId),]
char2.X.pacId <- droplevels(char2.X.pacId)
nlevels(as.factor(char2.X.pacId$SNP))###294
nlevels(as.factor(char2.X.pacId$Gene))###339
write.csv(char2.X.pacId,file="data/AHG.char2.X.pacId.csv")

###check how many hit with ara 
char2.NO.ara <- char2.X.pacId[char2.X.pacId$Best.hit.arabi.name=="",]
char2.NO.ara <- droplevels(char2.NO.ara)
names(char2.NO.ara)
nlevels(as.factor(char2.NO.ara$SNP))###41
nlevels(as.factor(char2.NO.ara$Gene))###45
write.csv(char2.NO.ara,file="data/AHG.char2.X.pacId.NO.ara.csv")


###check how many hit with ara 
char2.X.pacId.ara <- char2.X.pacId[!(char2.X.pacId$Best.hit.arabi.name==""),]
char2.X.pacId.ara <- droplevels(char2.X.pacId.ara)
names(char2.X.pacId.ara)
nlevels(as.factor(char2.X.pacId.ara$SNP))###269
nlevels(as.factor(char2.X.pacId.ara$Gene))###294
write.csv(char2.X.pacId.ara,file="data/AHG.char2.X.pacId.ara.csv")

###check the GO information
### and then using the data file to find the SNPs hit GO term 
char2.X.pacId.GO <- char2.X.pacId[!(char2.X.pacId$GO==""),]
char2.X.pacId.GO <- droplevels(char2.X.pacId.GO)
names(char2.X.pacId.GO)
nlevels(as.factor(char2.X.pacId.GO$SNP))###188
nlevels(as.factor(char2.X.pacId.GO$Gene))###199
write.csv(char2.X.pacId.GO,file="data/AHG.char2.X.pacI.GO.csv")

#### rice 
char2.X.pacId.rice <- char2.X.pacId[!(char2.X.pacId$Best.hit.rice.name==""),]
char2.X.pacId.rice <- droplevels(char2.X.pacId.rice)
names(char2.X.pacId.rice)
nlevels(as.factor(char2.X.pacId.rice$SNP))###283
nlevels(as.factor(char2.X.pacId.rice$Gene))###312
write.csv(char2.X.pacId.rice,file="data/AHG.char2.X.pacI.rice.csv")



