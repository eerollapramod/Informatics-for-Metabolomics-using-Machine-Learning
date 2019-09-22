rm(list=ls())
graphics.off()
# load required packages
install.packages("randomForest")
install.packages("caret")
install.packages("MASS")
install.packages('HH')

require(randomForest)
require(caret)
library(MASS)
require(HH)

source('common.R')

#######################################################################################
########## Regression models for TVC & Pseudomonas Counts by Enose data ##############
#######################################################################################
# Read and merge depedent and independent variable datasets
all.data <- common.match_csv_rows("Enose_data.csv","Bacterial_Counts.csv")

# Split training and test set for TVC count
set.seed(90)
train.index <- createDataPartition(all.data$Pseudomonads, p = .7,list = FALSE, times = 1)
train <- all.data[train.index,]
test <- all.data[-train.index,]

# Loading sensory (classification) data to determine classification 
sensory <- read.table('Sensory_score.csv', sep=',', header=TRUE,row.names=1)
all.data.sensory <- merge(all.data, sensory, by='row.names')
colnames(all.data.sensory)[12] <- "sensory"

#######################################################################
########### Regression models for TVC by Enose data ###################
#######################################################################

# Plotting TVC count distribution
hist(all.data$TVC, breaks=30, freq=F,col = "green", main="Bacterial count distribution (TVC)",
     xlab="Log10 of TVC/g", xlim=c(3, 10))

################### Linear Model for TVC count by Enose data #####################
lm.model.fit.TVC <- train(TVC ~ ., method='lmStepAIC', data=train)
predicted <- predict(lm.model.fit.TVC,test)

# calculate RMSE for lm and plot predicted vs actual+RMSE
RMSE.lm <- RMSE(test$TVC,predicted)
plot(predicted,test$TVC,xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue",
     main=paste("TVC linear Model RMSE:",round(RMSE.lm,digits = 2)))

# ploting confidence intervals
psd <- test[,9]
LM <- lm(predicted ~ psd) 
ci.plot(LM,main = "95% Confidence & Prediction Intervals for LM (TVC)")

################ k-nearest neighbours for TVC count by Enose data ################
knn.model.fit.TVC <- train(TVC ~ ., method='knn', data=train,tuneGrid=expand.grid(k=1:20))
predicted.knn <- predict(knn.model.fit.TVC,test)

# calculate RMSE for kNN and plot predicted vs actual+RMSE
RMSE.knn <- RMSE(test$TVC,predicted.knn)
plot(predicted.knn,test$TVC,xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("k-nearest neighbours RMSE:",round(RMSE.knn,digits = 2)))

# Plot visualizing the minimisation of the k parameter
plot(knn.model.fit.TVC$results$k, knn.model.fit.TVC$results$RMSE, xaxt="n",col = "blue",
     ylab="RMSE", xlab="k value", main="RMSE for kNN regression")
axis(1, at=1:20)

# ploting confidence intervals for k-nearest neighbours
KNN <- lm(predicted.knn ~ psd) 
ci.plot(KNN,main = "95% Confidence & Prediction Intervals for KNN (TVC)")
summary(knn.model.fit.TVC)
################### Random Forests for TVC count by Enose data #####################
RF.model.fit.TVC <- train(TVC ~ ., method='rf', train)
predicted.rf <- predict(RF.model.fit.TVC,test)

# calculate RMSE for kNN and plot predicted vs actual+RMSE
RMSE.rf <- RMSE(test$TVC,predicted.rf)
plot(predicted.rf,test$TVC,xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g", col="blue",
     main=paste("Random Forests RMSE:",round(RMSE.rf,digits = 2)))

# ploting confidence intervals for Random Forests
RF <- lm(predicted.rf ~ psd) 
ci.plot(RF,main = "95% Confidence & Prediction Intervals for RF (TVC)")

###########################################################################
####### Regression models for Pseudomonas count by enose data ############
###########################################################################
# Plotting Pseudomonas count distribution
hist(all.data$Pseudomonads, breaks=30, freq=F,col = "green", main="Bacterial count (Pseudomonas) distribution",
     xlab="Log10 of Pse/g", xlim=c(3, 10))

################### Linear Model for Pseudomonas count by Enose data #####################
lm.model.fit.Pse <- train(Pseudomonads ~ ., method='lmStepAIC', data=train)
predicted.lm.Pse <- predict(lm.model.fit.Pse,test)

# calculate RMSE for lm and plot predicted vs actual+RMSE
RMSE.lm.Pse <- RMSE(test$Pseudomonads,predicted.lm.Pse)
plot(predicted.lm.Pse,test$Pseudomonads,xlab="Predcted log10 Pse/g",col="blue",
     ylab="Actual log10 Pse/g", 
     main=paste("linear Model RMSE:",round(RMSE.lm.Pse,digits = 2)))

# ploting confidence intervals
psd.Pse <- test[,10]
LM.Pse <- lm(predicted.lm.Pse ~ psd.Pse) 
ci.plot(LM.Pse,main = "95% Confidence & Prediction Intervals for LM (Pse)")

################ k-nearest neighbours for TVC count by Enose data ################
knn.model.fit.Pse <- train(Pseudomonads ~ ., method='knn', data=train,tuneGrid=expand.grid(k=1:20))
predicted.knn.Pse <- predict(knn.model.fit.Pse,test)

# calculate RMSE for kNN and plot predicted vs actual+RMSE
RMSE.knn.Pse <- RMSE(test$Pseudomonads,predicted.knn.Pse)
plot(predicted.knn.Pse,test$Pseudomonads,xlab="Predcted log10 Pse/g",
     ylab="Actual log10 Pse/g", col = "blue",
     main=paste("k-nearest neighbours RMSE:",round(RMSE.knn.Pse,digits = 2)))

# Plot visualizing the minimisation of the k parameter
plot(knn.model.fit.Pse$results$k, knn.model.fit.Pse$results$RMSE, xaxt="n",col = "blue",
     ylab="RMSE", xlab="k value", main="RMSE for kNN regression")
axis(1, at=1:20)

# ploting confidence intervals for k-nearest neighbours
KNN.Pse <- lm(predicted.knn.Pse ~ psd.Pse) 
ci.plot(KNN.Pse,main = "95% Confidence & Prediction Intervals for KNN (Pse)")

################### Random Forests for TVC count by Enose data #####################
RF.model.fit.Pse <- train(Pseudomonads ~ ., method='rf', train)
predicted.rf.Pse <- predict(RF.model.fit.Pse,test)

# calculate RMSE for kNN and plot predicted vs actual+RMSE
RMSE.rf.Pse <- RMSE(test$Pseudomonads,predicted.rf.Pse)
plot(predicted.rf.Pse,test$Pseudomonads,xlab="Predcted log10 Pse/g",
     ylab="Actual log10 Pse/g", col = "blue",
     main=paste("Random Forests RMSE:",round(RMSE.rf.Pse,digits = 2)))

# ploting confidence intervals for Random Forests
RF.Pse <- lm(predicted.rf.Pse ~ psd.Pse) 
ci.plot(RF.Pse,main = "95% Confidence & Prediction Intervals for RF (Pse)")

#######################################################################################
########## Regression models for TVC & Pseudomonas Counts by HPLC data ###############
#######################################################################################
# clear workspace
rm(list=ls())
graphics.off()

# load required packages
require(randomForest)
library(caret)
library(MASS)
require(HH)
source('common.R')

# Read and merge depedent and independent variable datasets
all.data <- common.match_csv_rows("HPLC_data.csv","Bacterial_Counts.csv")

# Split training and test set for TVC count
set.seed(200)
train.index <- createDataPartition(all.data$TVC, p = .7,times = 1,list = F)
train <- all.data[train.index,]
test <- all.data[-train.index,]

######################################################################
########### Regression models for TVC by HPLC data ###################
######################################################################

# Plotting TVC count distribution
hist(all.data$TVC, breaks=30, freq=F, main="Bacterial count (TVC) distribution",col = "green",
     xlab="log10 of TVC/g", xlim=c(3, 10))

################### Linear Model for TVC count by HPLC data #####################
lm.model.fit.TVC <- lm(TVC ~ ., data=train)
predicted.lm.TVC <- predict(lm.model.fit.TVC,test)

# calculate RMSE for lm and plot predicted vs actual+RMSE
RMSE.lm.TVC <- RMSE(test$TVC,predicted.lm.TVC)
plot(predicted.lm.TVC,test$TVC,xlab="Predcted log10 TVC/g",col = "blue",
     ylab="Actual log10 TVC/g", 
     main=paste("linear Model RMSE:",round(RMSE.lm.TVC,digits = 2)))

# ploting confidence intervals
psd.TVC <- test[,20]
LM.TVC <- lm(predicted.lm.TVC ~ psd.TVC) 
ci.plot(LM.TVC,main = "95% Confidence & Prediction Intervals for LM (TVC)")

################ k-nearest neighbours for TVC count by HPLC data ################
knn.model.fit.TVC <- train(TVC ~ ., method='knn', data=train,tuneGrid=expand.grid(k=1:20))
predicted.knn.TVC <- predict(knn.model.fit.TVC,test)

# calculate RMSE for kNN and plot predicted vs actual+RMSE
RMSE.knn.TVC <- RMSE(test$TVC,predicted.knn.TVC)
plot(predicted.knn.TVC,test$TVC,xlab="Predcted log10 TVC/g",
     ylab="Actual log10 TVC/g", col = "blue",
     main=paste("k-nearest neighbours RMSE:",round(RMSE.knn.TVC,digits = 2)))

# Plot visualizing the minimisation of the k parameter
plot(knn.model.fit.TVC$results$k, knn.model.fit.TVC$results$RMSE, xaxt="n",col = "blue",
     ylab="RMSE", xlab="k value", main="RMSE for kNN regression")
axis(1, at=1:20)

# ploting confidence intervals for k-nearest neighbours
KNN.TVC <- lm(predicted.knn.TVC ~ psd.TVC) 
ci.plot(KNN.TVC,main = "95% Confidence & Prediction Intervals for KNN (TVC)")

################### Random Forests for TVC count by HPLC data #####################
RF.model.fit.TVC <- train(TVC ~ ., method='rf', train)
predicted.rf.TVC <- predict(RF.model.fit.TVC,test)

# calculate RMSE for kNN and plot predicted vs actual+RMSE
RMSE.rf.TVC <- RMSE(test$TVC,predicted.rf.TVC)
plot(predicted.rf.TVC,test$TVC,xlab="Predcted log10 TVC/g",ylab="Actual log10 TVC/g",col = "blue", 
     main=paste("Random Forests RMSE:",round(RMSE.rf.TVC,digits = 2)))

# ploting confidence intervals for Random Forests
RF.TVC <- lm(predicted.rf.TVC ~ psd.TVC) 
ci.plot(RF.TVC,main = "95% Confidence & Prediction Intervals for RF (TVC)")

###########################################################################
####### Regression models for Pseudomonas count by HPLC data ############
###########################################################################

# Plotting Pseudomonas count distribution
hist(all.data$Pseudomonads, breaks=30, freq=F,col = "green", main="Bacterial count (Pseudomonas) distribution",
     xlab="Log10 of Pse/g", xlim=c(3, 10))

################### Linear Model for Pseudomonas count by HPLC data #####################
lm.model.fit.Pse <- lm(Pseudomonads ~ ., data=train)
predicted.lm.Pse <- predict(lm.model.fit.Pse,test)

# calculate RMSE for lm and plot predicted vs actual+RMSE
RMSE.lm.Pse <- RMSE(test$Pseudomonads,predicted.lm.Pse)
plot(predicted.lm.Pse,test$Pseudomonads,xlab="Predcted log10 Pse/g",col = "blue",
     ylab="Actual log10 Pse/g", 
     main=paste("linear Model RMSE:",round(RMSE.lm.Pse,digits = 2)))

# ploting confidence intervals
psd.Pse <- test[,10]
LM.Pse <- lm(predicted.lm.Pse ~ psd.Pse) 
ci.plot(LM.Pse,main = "95% Confidence & Prediction Intervals for LM (Pse)")

################ k-nearest neighbours for TVC count by HPLC data ################
knn.model.fit.Pse <- train(Pseudomonads ~ ., method='knn', data=train,tuneGrid=expand.grid(k=1:20))
predicted.knn.Pse <- predict(knn.model.fit.Pse,test)

# calculate RMSE for kNN and plot predicted vs actual+RMSE
RMSE.knn.Pse <- RMSE(test$Pseudomonads,predicted.knn.Pse)
plot(predicted.knn.Pse,test$Pseudomonads,xlab="Predcted log10 Pse/g",
     ylab="Actual log10 Pse/g", col = "blue",
     main=paste("k-nearest neighbours RMSE:",round(RMSE.knn.Pse,digits = 2)))

# Plot visualizing the minimisation of the k parameter
plot(knn.model.fit.Pse$results$k, knn.model.fit.Pse$results$RMSE, xaxt="n",col = "blue",
     ylab="RMSE", xlab="k value", main="RMSE for kNN regression")
axis(1, at=1:20)

# ploting confidence intervals for k-nearest neighbours
KNN.Pse <- lm(predicted.knn.Pse ~ psd.Pse) 
ci.plot(KNN.Pse,main = "95% Confidence & Prediction Intervals for KNN (Pse)")

################### Random Forests for TVC count by HPLC data #####################
RF.model.fit.Pse <- train(Pseudomonads ~ ., method='rf', train)
predicted.rf.Pse <- predict(RF.model.fit.Pse,test)

# calculate RMSE for kNN and plot predicted vs actual+RMSE
RMSE.rf.Pse <- RMSE(test$Pseudomonads,predicted.rf.Pse)
plot(predicted.rf.Pse,test$Pseudomonads,xlab="Predcted log10 Pse/g",
     ylab="Actual log10 Pse/g", col = "blue",
     main=paste("Random Forests RMSE:",round(RMSE.rf.Pse,digits = 2)))

# ploting confidence intervals for Random Forests
RF.Pse <- lm(predicted.rf.Pse ~ psd.Pse) 
ci.plot(RF.Pse,main = "95% Confidence & Prediction Intervals for RF (Pse)")

##################### END CODE ######################

