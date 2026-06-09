########################################################################################
###############Preparing Figures########################################################
########################################################################################

###############Eight trait article  for Figure 3,4,5 ###################################
###############Eight trait article  for Figure 3,4,5 ###################################
###############Eight trait article  for Figure 3,4,5 ###################################
###making the image with different methods and traits and ind.SNP
###three ways with method, ind.SNP, and traits
#install.packages("directlabels", "ggplot2")
library(directlabels)
library(ggplot2)
getwd()
#setwd("/User/yonglizhao/Documents/R-corde for miscanthus project/")
##Remove the SNPs from GLM and rrBLUP for CmN trait because they did not include the K matrix 
## draw the line for the trend of the number of the SNP-traits associations instead of using the bar
#install.packages("directlabels")
library(directlabels)
library(ggplot2)
library(directlabels)
library(gridExtra)
library(ggpubr)
###############Eight trait article  for Figure 3 #####################################
###############Eight trait article  for Figure 3 #####################################
###############Eight trait article  for Figure 3 #####################################
###the number of all Ind. SNP
###the number of all Ind. SNP
###the number of all Ind. SNP
filename3 <- read.csv("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/data/all.chapter2.5.csv",row.names = 1)
filename <- filename3
filename$Trait.name <- as.factor(filename$Trait.name)
levels(filename$Trait.name)

levels(filename$Ind.SNP)[levels(filename$Ind.SNP)=="143+36088im"] <- "124+36088im"
levels(filename$Trait.name)[levels(filename$Trait.name)=="Bcirc_cm"] <- "Bcirc"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CCirc_cm"] <- "CCirc"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmD_BI_mm"] <- "CmD_BI"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmD_LI_mm"] <- "CmD_LI"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmDW_g"] <- "CmDW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="Cml_cm"] <- "Cml"
levels(filename$Trait.name)[levels(filename$Trait.name)=="SDW_kg"] <- "SDW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="Yld_kg"] <- "Yld"

levels(as.factor(filename$Ind.SNP))

filename$Ind.SNP <- as.factor(filename$Ind.SNP)
levels(filename$Ind.SNP)[levels(filename$Ind.SNP)=="124+36088im"] <- "C124im"
levels(filename$Ind.SNP)[levels(filename$Ind.SNP)=="C106+4202"] <- "C106"
levels(filename$Ind.SNP)[levels(filename$Ind.SNP)=="C116+3293"] <- "C116"
levels(filename$Ind.SNP)[levels(filename$Ind.SNP)=="C124+2560"] <- "C124"

str(filename)
t.H.1.M <- as.data.frame(table(filename[,c(12,16)]))
t.H.1.MCM <- t.H.1.M
t.H.1.M $Ind.SNP_1 <- paste("Ind.SNP_1")
t.H.1.MCM <- t.H.1.M [t.H.1.M $Trait.name!="CmN",]
#t.H.1.MCM <- filename
levels(t.H.1.MCM$Trait.name)
t.H.1.MCM$Ind.SNP<- factor(t.H.1.MCM$Ind.SNP,levels = c("C124im", "C106", "C116", "C124"))
t.H.1.MCM$Trait.name <- factor(t.H.1.MCM$Trait.name,levels = c("CmDW","Cml","CmD_BI","CmD_LI","CmN","Bcirc","CCirc","Yld" ))
levels(t.H.1.MCM$Trait.name)


##this one for eight traits without CmN
##this one for eight traits without CmN
##this one for eight traits without CmN
Trait.x <- ggplot(data=t.H.1.MCM, aes(x=Ind.SNP, y=Freq, fill=Ind.SNP,group=Trait.name, colour=Trait.name)) +
  geom_bar(stat="identity",width=0.4, position=position_dodge(width=-0.8),
           linetype = 0)+ facet_grid(Trait.name~Ind.SNP_1,switch="y",space="free")+
  geom_text(aes(x=Ind.SNP, y=Freq,label=Freq),size=8,hjust = 0.5,vjust=0.5,color="black")+
  theme(plot.margin=margin(0.8,0.9,0.5,0.25,'cm'),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5,vjust=0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 15),
        axis.title=element_text(size = 15),
        axis.title.y=element_text(size=20,face = "bold"),
        plot.caption=element_text(size=15),
        strip.text.y = element_text(face="bold", size=15),
        strip.text.x=element_blank(),
        legend.position="None",
        axis.text.x=element_text(face="bold",size=17,angle=0), axis.text.y=element_text(face="bold",size=16))+ 
  labs(y="The number of QTNs")

Trait.x

##this one for CmN trait
##this one for CmN trait
##this one for CmN trait
t.H.1.MCM <- t.H.1.M [t.H.1.M $Trait.name=="CmN",]
#t.H.1.MCM <- filename
levels(t.H.1.MCM$Trait.name)
t.H.1.MCM$Ind.SNP<- factor(t.H.1.MCM$Ind.SNP,levels = c("C124im", "C106", "C116", "C124"))
t.H.1.MCM$Trait.name <- factor(t.H.1.MCM$Trait.name,levels = c("CmDW","Cml","CmD_BI","CmD_LI","CmN","Bcirc","CCirc","Yld" ))
levels(t.H.1.MCM$Trait.name)
Trait.x.2 <- ggplot(data=t.H.1.MCM, aes(x=Ind.SNP, y=Freq, fill=Ind.SNP,group=Trait.name, colour=Trait.name)) +
  geom_bar(stat="identity",width=0.4, position=position_dodge(width=-0.8),
           linetype = 0)+ facet_grid(Trait.name~Ind.SNP_1,switch="y")+
  geom_text(aes(x=Ind.SNP, y=Freq,label=Freq),size=8,hjust = 0.5,vjust=0.5,color="black")+
  theme(plot.margin=margin(0.8,0.9,0.5,0.25,'cm'),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5,vjust=0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 15),
        axis.title=element_text(size = 15),
        axis.title.y=element_blank(),
        plot.caption=element_text(size=15),
        strip.text.y = element_text(face="bold", size=15),
        strip.text.x=element_blank(),
        legend.position="None",
        axis.text.x=element_text(face="bold",size=16,angle=0), axis.text.y=element_text(face="bold",size=16))+ 
  labs(y="The number of QTNs")
Trait.x.2

library(gridExtra)
library(ggpubr)
#figall.5 <- ggarrange(Trait.x,Trait.x.2, labels = c("A", "B"), ncol = 1, nrow = 2)
library(gridExtra)
figall.6 <- grid.arrange(Trait.x,Trait.x.2, ncol = 2, widths = c(6, 4))

###############Eight trait article  for Figure 4 #####################################
###############Eight trait article  for Figure 4 #####################################
###############Eight trait article  for Figure 4 #####################################
###the SNPs-traits associations from Method
###the SNPs-traits associations from Method
library(ggplot2)
library(directlabels)

#filename <- read.csv("Chapter2UP/all.chapter2up.csv",row.names = 1)
#levels(filename$Ind.SNP)[levels(filename$Ind.SNP)=="143+36088im"] <- "124+36088im"

filename3 <- read.csv("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/data/all.chapter2.5.csv",row.names = 1)
filename <- filename3
filename$Trait.name <- as.factor(filename$Trait.name)
str(filename)
levels(filename$Trait.name)
levels(filename$Trait.name)[levels(filename$Trait.name)=="Bcirc_cm"] <- "Bcirc"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CCirc_cm"] <- "CCirc"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmD_BI_mm"] <- "CmD_BI"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmD_LI_mm"] <- "CmD_LI"
levels(filename$Trait.name)[levels(filename$Trait.name)=="CmDW_g"] <- "CmDW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="Cml_cm"] <- "Cml"
levels(filename$Trait.name)[levels(filename$Trait.name)=="SDW_kg"] <- "SDW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="Yld_kg"] <- "Yld"
#filename <- filename[!(filename$Method=="GLM"|filename=="rrBLUP"),]
#filename <- droplevels(filename)
t.H.1.M <- as.data.frame(table(filename[,c(15:16)]))
t.H.1.MCM <- t.H.1.M
write.csv(t.H.1.M,file="Chapter2UP/QTNacrossmethod.csv")
#t.H.1.M$Freq[t.H.1.M$Freq==0] <- NA

t.H.1.MCM$Method_1 <- paste("Method_total")
str(t.H.1.MCM)
t.H.1.MCM_1 <- t.H.1.MCM[t.H.1.MCM$Trait.name!="CmN",]
t.H.1.MCM_1$Method <- as.factor(t.H.1.MCM_1$Method)
t.H.1.MCM_1$Method <- factor(t.H.1.MCM_1$Method,levels = c("ISIS EM-BLASSO","pLARmEB","FASTmrMLM","FarmCPU","mrMLM","FASTmrEMMA","pKWmEB","CMLM+SUPER","rrBLUP","MLM+SUPER","GLM"))

t.H.1.MCM_1$Trait.name <- factor(t.H.1.MCM_1$Trait.name,levels = c("CmDW","Cml","CmD_BI","CmD_LI","CmN","Bcirc","CCirc","Yld" ))


##this one for eight traits without CmN
##this one for eight traits without CmN
##this one for eight traits without CmN
Method.x <- ggplot(data=t.H.1.MCM_1, aes(x=Method, y=Freq, fill=Method,group=Trait.name, colour=Trait.name)) +
  geom_bar(stat="identity",width=0.5, position=position_dodge(width=0.4),
           linetype = 0)+ facet_grid(Trait.name~Method_1,switch="y",space="free")+
  geom_text(aes(x=Method, y=Freq,label=Freq),size=6,hjust = 0.5,vjust=0.5,color="black")+
  scale_y_continuous(breaks=seq(0,50,by=5))+
  theme(plot.margin=margin(0.25,0.25,0.25,0.25,'cm'),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5,vjust=0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=20,face = "bold"),
        plot.caption=element_text(size=12),
        strip.text.y = element_text(face="bold", size=12),
        strip.text.x=element_blank(),
        legend.position="None",
        axis.text.x=element_text(face="bold",size=13,angle=90), axis.text.y=element_text(size=12))+ 
  labs(y="The number of QTNs")

Method.x

##this one for CmN trait
##this one for CmN trait
##this one for CmN trait

t.H.1.MCM_2 <- t.H.1.MCM[t.H.1.MCM$Trait.name=="CmN",]
levels(as.factor(t.H.1.MCM$Method))
t.H.1.MCM_2$Method <- as.factor(t.H.1.MCM_2$Method)
t.H.1.MCM_2$Method <- factor(t.H.1.MCM_2$Method,levels = c("ISIS EM-BLASSO","pLARmEB","FASTmrMLM","FarmCPU","mrMLM","FASTmrEMMA","pKWmEB","CMLM+SUPER","rrBLUP","MLM+SUPER","GLM"))
Method.x.2 <- ggplot(data=t.H.1.MCM_2, aes(x=Method, y=Freq, fill=Method,group=Trait.name, colour=Trait.name)) +
  geom_bar(stat="identity",width=0.5, position=position_dodge(width=0.4),
           linetype = 0)+ facet_grid(Trait.name~Method_1,switch="y",space="free")+
  geom_text(aes(x=Method, y=Freq,label=Freq),size=6,hjust = 0.5,vjust=0.5,color="black")+
  theme(plot.margin=margin(0.25,0.25,0.25,0.25,'cm'),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5,vjust=0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_blank(),
        plot.caption=element_text(size=12),
        strip.text.y = element_text(face="bold", size=12),
        strip.text.x=element_blank(),
        legend.position="None",
        axis.text.x=element_text(face="bold",size=13,angle=90), axis.text.y=element_text(size=12))+ 
  labs(y="The number of QTNs")
library(gridExtra)
library(ggpubr)
#figall.5 <- ggarrange(Trait.x,Trait.x.2, labels = c("A", "B"), ncol = 1, nrow = 2)
library(gridExtra)
figall.7 <- grid.arrange(Method.x,Method.x.2, ncol = 2, widths = c(5, 5))




###################Figure 5S #########################################
###################Figure 5S #########################################
###################Figure 5S #########################################

####the prediction accuracies.
####the prediction accuracies.
####the prediction accuracies.

###import all of the data 
GP.all.2 <- read.csv("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus/data/GP.all.2.csv",row.names = 1)

#GP.all.2 <- read.csv(file="data/GP.all.2.csv",row.names = 1)
str(GP.all.2)
levels(as.factor(GP.all.2$Trait))

###get all of the traits for the charpter2 

GP.all.culm <- GP.all.2[GP.all.2$Trait=="Bcirc_cm "| GP.all.2$Trait=="CCirc_cm "| GP.all.2$Trait=="CmD_BI_mm "| GP.all.2$Trait=="CmDW_g "|
                          GP.all.2$Trait=="Cml_cm " | GP.all.2$Trait=="CmN. " | GP.all.2$Trait=="Yld "| GP.all.2$Trait== "CmD_LI_mm ",]

###chasing the PC1 to 3 
GP.all.culm$Ind.SNP <- as.factor(GP.all.culm$Ind.SNP)
levels(GP.all.culm$Ind.SNP)

levels(GP.all.culm$Ind.SNP)[levels(GP.all.culm$Ind.SNP)=="F116im"] <- "C124im"

levels(GP.all.culm$Ind.SNP)

GP.all.c <- GP.all.culm[GP.all.culm$Ind.SNP=="C106" | GP.all.culm$Ind.SNP=="C116" | GP.all.culm$Ind.SNP=="C124" | GP.all.culm$Ind.SNP=="C124im",]

GP.all.c <- GP.all.culm[GP.all.culm$Ind.SNP=="C106" | GP.all.culm$Ind.SNP=="C116" | GP.all.culm$Ind.SNP=="C124" | GP.all.culm$Ind.SNP=="C124im",]

levels(as.factor(GP.all.c$Ind.SNP))


GP.all.c <- droplevels(GP.all.c)
str(GP.all.c)
levels(as.factor(GP.all.c$Trait))
levels(as.factor((GP.all.c$Ind.SNP)))

GP.all.c$Trait <- as.factor(GP.all.c$Trait)
levels(GP.all.c$Trait)[levels(GP.all.c$Trait)=="Bcirc_cm "] <- "Bcirc"
levels(GP.all.c$Trait)[levels(GP.all.c$Trait)=="CCirc_cm "] <- "CCirc"
levels(GP.all.c$Trait)[levels(GP.all.c$Trait)=="CmD_BI_mm "] <- "CmD_BI"
levels(GP.all.c$Trait)[levels(GP.all.c$Trait)=="CmD_LI_mm "] <- "CmD_LI"
levels(GP.all.c$Trait)[levels(GP.all.c$Trait)=="CmDW_g "] <- "CmDW"
levels(GP.all.c$Trait)[levels(GP.all.c$Trait)=="Cml_cm "] <- "Cml"
#levels(GP.all.c$Trait)[levels(GP.all.c$Trait)=="SDW_kg "] <- "SDW"
levels(GP.all.c$Trait)[levels(GP.all.c$Trait)=="Yld "] <- "Yld"
levels(GP.all.c$Trait)[levels(GP.all.c$Trait)=="CmN. "] <- "CmN"

levels(as.factor((GP.all.c$Trait)))
### plot the image with FT_1
GP.all.c.m <- GP.all.c[GP.all.c$names=="mean",]
str(GP.all.c.m)
write.csv(GP.all.c.m,file="Chapter2UP/GP.chapter2.gp.mean.csv")
GP.all.c.m <- droplevels(GP.all.c.m)
GP.all.c.m$Ind.SNP<- factor(GP.all.c.m$Ind.SNP,levels = c("C124im", "C124", "C116", "C106"))
str(GP.all.c.m)
library(ggplot2)
#install.packages("ggrepel")
library(ggrepel)

library(directlabels) ### this the package for the geom_dl function 
###this one is for PC were separated and Trait were group together 

levels(as.factor(GP.all.c.m$x))
GP.all.c.m$y <- ifelse(GP.all.c.m$x<0, 0.00, GP.all.c.m$x)
levels(as.factor(GP.all.c.m$y))

GP.all.c.m$Trait <- factor(GP.all.c.m$Trait,levels = c("CmDW","Cml","CmD_BI","CmD_LI","CmN","Bcirc","CCirc","Yld" ))
Culmall.Trait.group.1 <- ggplot(data = GP.all.c.m, aes(x=PC, y=y, fill=PC, group=PC, colour=PC)) +
  geom_bar(stat="identity",width=0.5, position=position_dodge(width=0.5),
           linetype = 0)+ 
  facet_grid(Trait~Ind.SNP,switch="y",space="free")+
  geom_text(aes(x=PC, y=y,label=round(y,2)),size=6, hjust = 0.5,vjust=0.5,colour="black")+
  theme(plot.margin=margin(0.25,0.25,0.25,0.25,'cm'),
        plot.title = element_text(size = 20, face = "bold", hjust = 0,vjust=0.5), plot.subtitle=element_text(size=20),axis.title.x=element_blank(),
        axis.text.y.left =element_text(size = 16),
        axis.text=element_text(size = 20),
        axis.title=element_text(size = 20),
        axis.title.y=element_text(size = 20,face = "bold"),
        plot.caption=element_text(size=20),
        strip.text.y = element_text(face="bold", size=15),
        strip.text.x=element_text(face="bold", size=17),
        legend.position="None",
        axis.text.x=element_text(face="bold",size=17,angle=20), axis.text.y=element_text(size=20))+
  labs(y="Mean of prediction accuracies",title="A")

###Figure 5 SA
###Figure 5 SA
###Figure 5 SA
Culmall.Trait.group.1 

####plot variance 

GP.all.c.s <- GP.all.c[GP.all.c$names=="std.dev",]
GP.all.c.s <- droplevels(GP.all.c.s)
GP.all.c.s$Ind.SNP<- factor(GP.all.c.s$Ind.SNP,levels = c("C124im", "C124", "C116", "C106"))
str(GP.all.c.s)
write.csv(GP.all.c.s,file="Chapter2UP/GP.chapter2.gp.std.csv")

###switch the x-axis to PC
GP.all.c.s$Trait <- factor(GP.all.c.s$Trait,levels = c("CmDW","Cml","CmD_BI","CmD_LI","CmN","Bcirc","CCirc","Yld" ))
Culmall.PCx.s.2 <- ggplot(data = GP.all.c.s, aes(x=PC, y=x, fill=PC, group=PC, colour=PC)) +
  geom_bar(stat="identity",width=0.5, position=position_dodge(width=0.4),
           linetype = 0)+ 
  facet_grid(Trait~Ind.SNP,switch="y",space="free")+
  geom_text(aes(x=PC, y=x,label=round(x,2)),size=6, hjust = 0.5,vjust=0.5,colour="black")+
  theme(plot.margin=margin(0.25,0.25,0.25,0.25,'cm'),
        plot.title = element_text(size = 20, face = "bold", hjust = 0,vjust=0.5), plot.subtitle=element_text(size=20),axis.title.x=element_blank(),
        axis.text.y.left =element_text(size = 16),
        axis.text=element_text(size = 20),
        axis.title=element_text(size = 20),
        axis.title.y=element_text(size = 20,face = "bold"),
        plot.caption=element_text(size=20),
        strip.text.y = element_text(face="bold", size=15),
        strip.text.x=element_text(face="bold", size=20),
        legend.position="None",
        axis.text.x=element_text(face="bold",size=17,angle=20), axis.text.y=element_text(size=20))+
  labs(y="Standard deviation of prediction accuracies",title="B")
####Figure 5S A
Culmall.PCx.s.2





