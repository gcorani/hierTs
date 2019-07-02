batchSynth <- function () {
  reps <- 250
  correl <- c(-0.9, -0.5, 0.01, 0.1, 0.3, 0.5, 0.9)
  n <- c(5, 10, 20, 100)
  
  
  for (currentCorr in correl) {  
    for (currentN in n) {
      start_time <- Sys.time()
      for (i in 1:reps){
        hierRec(dset="synthetic", seed=i, synth_n = currentN, synthCorrel = currentCorr)
      }
      end_time <- Sys.time()
      print (end_time - start_time)
    }
  }
}

setwd(paste(getwd(),"results", sep="/"))
fileList <- list.files()
dataFrame <- data.frame()
colnames(dataFrame) <- c("correl","sampleSize","mseMint/mseBayesDiag","mseMint/mseBayesCorr","mseBase/mseBayesCorr")
counter <- 1
for (currentFile in fileList){
  if (grepl("mseHierReconcsynthetic", currentFile)){
    currentData <- read_csv(currentFile)
    dataFrame[counter,"correl"] <- currentData$correlB1_U[1]
    dataFrame[counter,"sampleSize"] <- currentData$sampleSize[1]
    dataFrame[counter,"mseMint/mseBayesDiag"] <- median(currentData$mseMint / currentData$mseBayesDiag)
    dataFrame[counter,"mseMint/mseBayesCorr"] <- median(currentData$mseMint / currentData$mseBayesCorr)
    dataFrame[counter,"mseBase/mseBayesCorr"] <- median(currentData$mseBase / currentData$mseBayesCorr)
  }
  counter <- counter + 1
}
dataFrame <- na.omit(dataFrame)
