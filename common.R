library(caret)

# assumes lines can be matched based on row.names, which then become the data.frame rownames
common.match_csv_rows <- function (enose_file, sensory_file) {
  enose   <- read.table(enose_file, sep=",", header=TRUE, row.names=1)
  sensory <- read.table(sensory_file, header = TRUE, sep =",",row.names=1)
  # Matching rows from enose to rows from sensory
  merged <- merge(enose, sensory, by="row.names")
  rownames(merged) = merged[,1]
  return(as.data.frame(merged[,-1]))
}

# uses caret's createDataPartition to select specific rows from a vector of classes (eg sensory)
# which can then be used to split data into training and test sets
common.partition_data <- function (classes) {
  trainIndex <- createDataPartition(as.factor(classes), p = .7, 
                                    list = FALSE, 
                                    times = 1)
  return(trainIndex);
}

# applies the supplied function to the provided data.frame after partitioning it
# into a training set, a test set and relevant class vectors
# which are then returned along with the outcome of the function
# as a named list
# 
# Usage:
#
#   results.k3 <- common.run(function(trnS, trnC, tstS, tstC) { return(knn(trnS, tstS, trnC, 3)) }, some_data)
#
#   results.k3$model     # provides access to the knn() results
#   results.k3$trainSet  # provides access to the training set
#   results.k3$testSet   # provides access to the test set
#   results.k3$trainCl   # provides access to the training classes
#   results.k3$testCl    # provides access to the test classes
# 
common.run <- function (operation, data.frame) {
  trainIndex <- common.partition_data(data.frame[,ncol(data.frame)])

  trainSet <- data.frame[trainIndex,]
  testSet <- data.frame[-trainIndex,]

  trainCl <- trainSet[,ncol(trainSet)]
  trainCl <- as.factor(trainCl)
  testCl <- testSet[,ncol(trainSet)]
  testCl <- as.factor(testCl)

  trainSet <- trainSet[,1:(ncol(trainSet)-1)]
  testSet <- testSet[,1:(ncol(testSet)-1)]
  model <- operation(trainSet, trainCl, testSet, testCl)
  
  return(
    list(
      trainSet = trainSet,
      testSet = testSet,
      trainCl = trainCl,
      testCl = testCl,
      model = model 
    )
  )
}

