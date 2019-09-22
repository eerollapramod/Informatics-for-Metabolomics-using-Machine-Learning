rm(list=ls())
graphics.off() 

install.packages("matlab")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("mixOmics"))

install.packages("caret")
install.packages("tidyr")

require(matlab)
require(mixOmics)
require(caret)
library(tidyr)

source("common.r")

##########################################################################
# ############# Data Preperation and Exploratory Analysis ################
##########################################################################

# read data
enose <- read.table("Enose_data.csv", sep=",", header=TRUE)
colnames(enose)[1]<-"Samples"
bacterial_counts <- read.table("Bacterial_Counts.csv", sep=",", header=TRUE)
hplc <- read.table("HPLC_data.csv",sep = ",",header = T)

############# Merge data ###############
all.data <- common.match_csv_rows("Enose_data.csv","Bacterial_Counts.csv")
sensory <- read.table("Sensory_score.csv", sep=",", header=TRUE,row.names = 1)
colnames(sensory)[1]<-"sensory"
all.data.sensory <- merge(all.data, sensory, by="row.names")
rownames(all.data.sensory) <- all.data.sensory$Row.names
all.data.sensory<- all.data.sensory[,-1]

###################### PCA ########################
## PCA for Enose data for further analysis of count distribution according to Sensory scores
enoseForPCA <- read.table("Enose_data.csv", sep=",", header=TRUE,row.names = 1)
pcaForEnose<-pca(enoseForPCA,ncomp = 4,scale = T)
plotIndiv(pcaForEnose,ind.names=T, style="lattice",title = "PCA of counts by Enose data")
## PCA for HPLC data for further analysis of count distribution according to Sensory scores
hplcForPCA <- read.table("HPLC_data.csv",sep = ",",header = T,row.names = 1)
pcaForHPLC <- pca(hplcForPCA, ncomp = 4, scale = TRUE)
plotIndiv(pcaForHPLC,ind.names=T, style="lattice",title = "PCA of counts by HPLC data")

# Generating histograms for TVC counts during fresh, semi fresh and spoiled stages
par(mfrow=c(2,2))
hist(all.data.sensory[(all.data.sensory$sensory==1),]$TVC,breaks = 10,plot=T,xlab = "Fresh Counts",col = "green",main="TVC (Fresh)")
hist(all.data.sensory[(all.data.sensory$sensory==2),]$TVC,breaks = 10,plot=T,xlab = "Semi-Fresh Counts",col = "yellow",main="TVC (Semi-Fresh)")
hist(all.data.sensory[(all.data.sensory$sensory==3),]$TVC,breaks = 15,plot=T,xlab = "Spoiled Counts",col = "red",main="TVC (Spoiled)")
# Generating histograms for Pseudomonads counts during fresh, semi fresh and spoiled stages
par(mfrow=c(2,2))
hist(all.data.sensory[(all.data.sensory$sensory==1),]$Pseudomonads,plot=T,breaks = 10,xlab = "Fresh Counts",col = "green",main="Pseudomonads Counts (Fresh)")
hist(all.data.sensory[(all.data.sensory$sensory==2),]$Pseudomonads,plot=T,breaks = 10,xlab = "Semi-Fresh Counts",col = "yellow",main="Pseudomonads Counts (Semi-Fresh)")
hist(all.data.sensory[(all.data.sensory$sensory==3),]$Pseudomonads,plot=T,breaks = 15,xlab = "Spoiled Counts",col = "red",main="Pseudomonads Counts (Spoiled)")


### plotting bacterial counts under different temperatures during different time intervals ###
# seperating temperatre and time
seperated.counts <- separate(bacterial_counts[-c(1:2),],"Samples",into=c("temp","time"),sep="F")
# Total Viable Counts plot
par(mfrow=c(1,1))
plot(seperated.counts[1:11,2],seperated.counts[1:11,3],type="l",xlab="time (hours)",ylab="count",main="TVC Counts Under different TEMPs over Time",xlim = c(2,20), ylim = c(2,10),xaxt="n")
lines(seperated.counts[12:20,2],seperated.counts[12:20,3],col="yellow")
lines(seperated.counts[21:30,2],seperated.counts[21:30,3],col="green")
lines(seperated.counts[31:41,2],seperated.counts[31:41,3],col="blue")
lines(seperated.counts[42:50,2],seperated.counts[42:50,3],col="red")
legend(x="bottomright",y=0.92, c("0˚C","5˚C","10˚C","15˚C","20˚C"),c("black","yellow","green","blue","red"))
axis(1,seq(from=2,to=20,by=2))

# Pseudomonas plot
plot(seperated.counts[1:11,2],seperated.counts[1:11,4],type="l",xlab="time (hours)",ylab="count",main="Pseudomonas Counts Under different TEMPs over Time",xlim = c(2,20), ylim = c(2,10),xaxt="n")
lines(seperated.counts[12:20,2],seperated.counts[12:20,4],col="yellow")
lines(seperated.counts[21:30,2],seperated.counts[21:30,4],col="green")
lines(seperated.counts[31:41,2],seperated.counts[31:41,4],col="blue")
lines(seperated.counts[42:50,2],seperated.counts[42:50,4],col="red")
legend(x="bottomright",y=0.92, c("0˚C","5˚C","10˚C","15˚C","20˚C"),c("black","yellow","green","blue","red"))
axis(1,seq(from=2,to=20,by=2))

#########################################
############# scatter plots  ############
#########################################
merged.enose.counts <- merge(enose, bacterial_counts,by.x = "Samples")
merged.enose.counts <- data.frame(merged.enose.counts[,-1],row.names = merged.enose.counts[,1])

merged.hplc.counts <- merge(hplc,bacterial_counts,by.x = "Samples")
merged.hplc.counts <- data.frame(merged.hplc.counts[,-1],row.names = merged.hplc.counts[,1])

#∞∞∞∞∞∞ enose - TVC ∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞
par(mfrow=c(2,3))
for (i in 1:8){
  plot(merged.enose.counts[,i],merged.enose.counts[,9],col ="red",ylab = "Bacterial Count", main = "TVC by Enose data")
}
#∞∞∞∞∞∞∞∞ enose - pseudomonads ∞∞∞∞∞∞∞∞∞
par(mfrow=c(2,3))
for (i in 1:8){
  plot(merged.enose.counts[,i],merged.enose.counts[,10],col ="blue",ylab = "Bacterial Count", main = "Pseudomonad by Enose data")
}
#∞∞∞∞∞∞∞∞∞∞∞∞ HPLC -TVC ∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞
par(mfrow=c(2,3))
for (i in 1:19){
  plot(merged.hplc.counts[,i],merged.hplc.counts[,9],col ="green",ylab = "Bacterial Count", main = "TVC by HPLC data")
}
#∞∞∞∞∞∞∞∞∞∞∞ HPLC - pseudomonads ∞∞∞∞∞∞∞
par(mfrow=c(2,3))
for (i in 1:19){
  plot(merged.hplc.counts[,i],merged.hplc.counts[,9],col ="orange",ylab = "Bacterial Count", main = "Pseudomonads by HPLC data")
}
############## Boxplots ####################

# boxplots of TVC counts against Sensory scores
par(mfrow=c(1,1))
boxplot(all.data.sensory$TVC ~ all.data.sensory$sensory, main = "TVC Count by Sensory Score",ylab = "TVC Count",
        names=c("1(Fresh)","2(Semi-Fresh)","3(Spoiled)"),col = c("green","yellow","red"))
# boxplots of Paeudomonads against Sensory scores
boxplot(all.data.sensory$Pseudomonads ~ all.data.sensory$sensory, main = "Pseudomonads Count by Sensory Score",
        ylab = "Pseudomonads Count",names=c("1(Fresh)","2(Semi-Fresh)","3(Spoiled)"),col = c("green","yellow","red"))

################## END CODE ####################

