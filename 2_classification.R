# clear workspace
rm(list=ls())
graphics.off() 

# instal/load packages/libraries 
install.packages("class")
install.packages("caret", dependencies = TRUE)
install.packages('e1071', dependencies=TRUE)
install.packages("kernlab")
install.packages("gmodels")
install.packages("tree")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("mixOmics"))

install.packages("pROC")

library(class)
library(caret)
library(kernlab)
library(gmodels)
library(tree)
library(mixOmics)
require(pROC)

source('common.R')

# Reading and merging data
merged <- common.match_csv_rows("Enose_data.csv", "Sensory_score.csv")
merged2 <- common.match_csv_rows("HPLC_data.csv","Sensory_score.csv")

########## creating Data Partitions ################
class.merged <- as.factor(merged[,9])
class.merged2 <- as.factor(merged2[,20])

set.seed(1600)
DP1 <- createDataPartition(class.merged, times = 1, p = 0.7)
DP2 <- createDataPartition(class.merged2, times = 1, p = 0.7)
# Training and testing sets for sensory scores by enose data 
Train <- merged[DP1$Resample1,-1]
trainclass <- class.merged[DP1$Resample1]
Test <- merged[-DP1$Resample1,-1]
testclass <- class.merged[-DP1$Resample1]
# Training and testing sets for sensory scores by HPLC data 
Train.hplc <- merged2[DP2$Resample1,-1]
trainclass.hplc <- class.merged2[DP2$Resample1]
Test.hplc <- merged2[-DP2$Resample1,-1]
testclass.hplc <- class.merged2[-DP2$Resample1]

########################################################################################
##### Objective 1 :  Classification models for Sensory scores using Enose data #########
########################################################################################

###### PCA for merged enose data & sensory scores #######
pca.enose <- pca(merged[,-ncol(merged)], ncomp=4, scale=TRUE)
plotIndiv(pca.enose, ind.names=as.character(merged[,ncol(merged)]),title = "PCA plot of Enose data",
          group=as.factor(merged[,ncol(merged)]), style="lattice")

############ Decision Trees model for sensory score by enose data ##############
dtrees.model = tree(trainclass ~ ., data=Train)
plot(dtrees.model)
text(dtrees.model)
summary(dtrees.model)
# predicting model using test set
predicted.dtrees <- predict(dtrees.model, Test, type="class")
# cross table
crosstable.dtrees <- CrossTable(testclass,predicted.dtrees,prop.r = F,prop.c = F,prop.t = F,prop.chisq = F)
# confusion matrix
CM.dtrees <- confusionMatrix(predicted.dtrees, testclass, positive="3")
print(CM.dtrees)
cat('Tree accuracy: ', CM.dtrees$overall[1])

####### ROC curve for Decision Trees model #######
m.roc.dtrees <- multiclass.roc(as.numeric(predicted.dtrees),as.numeric(testclass))
rocs.dtrees <- m.roc.dtrees$rocs
plot.roc(main="DT ROC for Sensory Scores by Enose data",rocs.dtrees[[1]])
lines(rocs.dtrees[[2]],col=2)
lines(rocs.dtrees[[3]],col=3)
legend("bottomright",y=NULL, legend = c("Class 1","Class 2","Class 1"), c("black","red","green"))
# pruning is not need due to tree is being too small to overfit

############ SVM model for sensory score by enose data #########
svm.model<- ksvm(trainclass ~ ., data=Train,C=1)
predicted.svm <- predict(svm.model,Test,type="response")
# Cross Table
crosstable.svm <- CrossTable(testclass,predicted.svm,prop.c = F,prop.r = F,prop.t = F,prop.chisq = F)
print(crosstable.svm)
# Confusion Matrix
CM.svm <- confusionMatrix(predicted.svm,testclass,positive = "3")
# print confusion matrix
CM.svm
# print SVM accuracy
cat('SVM accuracy: ', CM.svm$overall[1])

### Measuring success rate of each svm model over 100 iterations
iteration <- as.array(c(1:100))
accuracy.vs.sampling.hplc <- apply(iteration, 1, function(i) {
  results.hplc <- common.run(function(Train.hplc, trainclass.hplc, Test.hplc, testclass.hplc) {
    m2 = ksvm(trainclass.hplc ~ ., data=Train.hplc, C=1)
    return(m2)
  }, merged2)
  kernel.predicted.hplc <- predict(results.hplc$model, results.hplc$testSet, type="response")
  kernel.confusion.matrix.hplc <- confusionMatrix(kernel.predicted.hplc, results.hplc$testCl, positive="3")
  return(kernel.confusion.matrix.hplc$overall[1])
})

evolution.of.accuracy.mean.hplc = c()
for(i in 1:length(accuracy.vs.sampling.hplc)) {
  evolution.of.accuracy.mean.hplc <- c(evolution.of.accuracy.mean.hplc, mean(accuracy.vs.sampling.hplc[1:i]))
}
plot(iteration, evolution.of.accuracy.mean.hplc, type='l',main="SVM success rate (100 iterations)")

############ kNN model for sensory score by enose data #########
kNN.model <- knn(train = Train,test = Test, k = 3, cl = trainclass)
crosstable.knn <- CrossTable(testclass,kNN.model,prop.r = F,prop.c = F,prop.t = F,prop.chisq = F)
crosstable.knn
CM.knn <- confusionMatrix(kNN.model, testclass, positive="3")
CM.knn
# print SVM accuracy
cat('k-NN accuracy: ', CM.knn$overall[1])

# measuring success rate for knn model over 20 iterations
normalise<-function(x){
  return ((x-min(x))/(max(x)-min(x)))
}

norm.data <- apply(merged[,1:ncol(merged)-1], 2, normalise)
norm.data <- as.data.frame(norm.data)
norm.data$sensory <- merged[,ncol(merged)]

array.of.ks <- as.array(c(1:20))
accuracy.vs.k <- apply(array.of.ks, 1, function(k) {
  set.seed(1600)
  results <- common.run(function(trnS, trnC, tstS, tstC) { return(knn(trnS, tstS, trnC, k)) }, norm.data)
  cross.table2 <- CrossTable(results$testCl, results$model, prop.chisq=FALSE)
  CM.knn <- confusionMatrix(results$model, results$testCl, positive="3")
  return(CM.knn$overall[1])
})
plot(accuracy.vs.k, type='l',ylab="Accuracy of K",main="Effect of k over Accuracy")
##### Plotting Variable importance is not available svm classification model 

########################################################################################
##### Objective 1 :  Classification models for Sensory scores using HPLC data #########
########################################################################################

###### PCA #######
pca.hplc <- pca(merged2[,-ncol(merged2)], ncomp=4, scale=TRUE)
plotIndiv(pca.hplc, ind.names=as.character(merged2[,ncol(merged2)]),
          group=as.factor(merged2[,ncol(merged2)]), style="lattice",title = "PCA plot of HPLC data")

############ Decision Trees model for sensory score by HPLC data #########
dtrees.model.hplc = tree(trainclass.hplc ~ ., data=Train.hplc)
plot(dtrees.model.hplc)
text(dtrees.model.hplc)
summary(dtrees.model.hplc)
# predicting model using test set
predicted.dtrees.hplc <- predict(dtrees.model.hplc, Test.hplc, type="class")
# cross table
crosstable.dtrees.hplc <- CrossTable(testclass.hplc,predicted.dtrees.hplc,prop.r = F,prop.c = F,prop.t = F,prop.chisq = F)
# confusion matrix
CM.dtrees.hplc <- confusionMatrix(predicted.dtrees.hplc, testclass.hplc, positive="3")
cat('Tree accuracy: ', CM.dtrees.hplc$overall[1])

####### ROC curve for Decision Trees model  #######
m.roc.dtrees.hplc <- multiclass.roc(as.numeric(predicted.dtrees.hplc),as.numeric(testclass.hplc))
rocs.dtrees.hplc <- m.roc.dtrees.hplc$rocs
plot.roc(main="DT ROC for Sensory Scores by HPLC data",(rocs.dtrees.hplc[[1]]))
lines(rocs.dtrees.hplc[[2]],col=2)
lines(rocs.dtrees.hplc[[3]],col=3)
legend("bottomright",y=NULL, legend = c("Class 1","Class 2","Class 1"), c("black","red","green"))
# pruning is not need due to tree is being too small to overfit

############ SVM model for sensory score by HPLC data #########
svm.model.hplc<- ksvm(trainclass.hplc ~ ., data=Train.hplc,C=1)
predicted.svm.hplc <- predict(svm.model.hplc,Test.hplc,type="response")
# Cross Table
crosstable.svm.hplc <- CrossTable(testclass.hplc,predicted.svm.hplc,prop.c = F,prop.r = F,prop.t = F,prop.chisq = F)
# Confusion Matrix
CM.svm.hplc <- confusionMatrix(predicted.svm.hplc,testclass.hplc,positive = "3")
# print confusion matrix
CM.svm.hplc
# print SVM accuracy
cat('SVM accuracy: ', CM.svm.hplc$overall[1])

### Measuring success rate of each svm model over 100 iterations for sensory score by HPLC data
iteration <- as.array(c(1:100))
accuracy.vs.sampling <- apply(iteration, 1, function(i) {
  results <- common.run(function(Train, trainclass, Test, testclass) {
    m = ksvm(trainclass ~ ., data=Train, C=1)
    return(m)
  }, merged)
  kernel.predicted <- predict(results$model, results$testSet, type="response")
  kernel.confusion.matrix <- confusionMatrix(kernel.predicted, results$testCl, positive="3")
  return(kernel.confusion.matrix$overall[1])
})

evolution.of.accuracy.mean = c()
for(i in 1:length(accuracy.vs.sampling)) {
  evolution.of.accuracy.mean <- c(evolution.of.accuracy.mean, mean(accuracy.vs.sampling[1:i]))
}
plot(iteration, evolution.of.accuracy.mean, type='l',ylim=c(0.0,0.6),main="SVM success rate on HPLC data (100 iterations)")

############ kNN model for sensory score by HPLC data #########
kNN.model.hplc <- knn(train = Train.hplc,test = Test.hplc, k = 3, cl = trainclass.hplc)
crosstable.knn.hplc <- CrossTable(testclass.hplc,kNN.model.hplc,prop.r = F,prop.c = F,prop.t = F,prop.chisq = F)
crosstable.knn.hplc
CM.knn.hplc <- confusionMatrix(kNN.model.hplc, testclass.hplc, positive="3")
CM.knn.hplc

# measuring success rate for knn model over 100 iterations
normalise<-function(x){
  return ((x-min(x))/(max(x)-min(x)))
}

norm.data.hplc <- apply(merged2[,1:ncol(merged2)-1], 2, normalise)
norm.data.hplc <- as.data.frame(norm.data.hplc)
norm.data.hplc$sensory <- merged2[,ncol(merged2)]

array.of.ks <- as.array(c(1:20))
accuracy.vs.k <- apply(array.of.ks, 1, function(k) {
  set.seed(1600) # need to select same subsets when partitioning
  results <- common.run(function(trnS, trnC, tstS, tstC) { return(knn(trnS, tstS, trnC, k)) }, norm.data.hplc)
  cross.table2 <- CrossTable(results$testCl, results$model, prop.chisq=FALSE)
  CM.knn <- confusionMatrix(results$model, results$testCl, positive="3")
  return(CM.knn$overall[1])
})
plot(accuracy.vs.k, type='l',ylab="Accuracy of K",main="Effect of k on Accuracy (Sensory by HPLC) ")
##### Plotting Variable importance is not available svm classification model 

####################################
######### End of CODE ##############
####################################

