batchSynth <- function (correl, n) {
  source("hierRec.R")
  library(readr)
  library(ggplot2)
  reps <- 1000
  for (i in 1:reps){
    hierRec(dset="synthetic", seed=i, synth_n = n, synthCorrel = correl)
  }
  #we need to load the same file in which hierRec has written the results
  dset <- paste0("synthetic","_correl",correl,"_n",n)
  filename <- paste0("results/mseHierReconc",dset,".csv",sep="")
  currentData <- read_csv(filename)
  #plot the graph
  title <- paste0( "plot_n", as.character(currentData$sampleSize[1]), "_correlB1_U_", 
                   currentData$correlB1_U[1])
  ggplot(currentData, mapping = aes(correlB2_U, mseMint/mseBayesCorr)) +
    # geom_point(aes(color = class)) +
    geom_smooth() +
    ggtitle(title) + 
    geom_hline(yintercept=1)
  
  ggsave(paste0("results/",title,".pdf"))
}

setwd(paste(getwd(),"results", sep="/"))
fileList <- list.files()
dataFrame <- data.frame()
# colnames(dataFrame) <- c("correl","sampleSize","mseMint/mseBayesDiag","mseMint/mseBayesCorr","mseBase/mseBayesCorr")
counter <- 1
library(ggplot2)
library(readr)
for (currentFile in fileList){
  if (grepl("mseHierReconcsynthetic", currentFile)){
    currentData <- read_csv(currentFile)
    dataFrame[counter,"correl"] <- currentData$correlB1_U[1]
    dataFrame[counter,"sampleSize"] <- currentData$sampleSize[1]
    if (currentData$sampleSize[1]==100 & currentData$correlB1_U[1]==-0.01){
      browser()
    }
    dataFrame[counter,"mseMint/mseBayesDiag"] <- median(currentData$mseMint / currentData$mseBayesDiag)
    dataFrame[counter,"mseMint/mseBayesCorr"] <- median(currentData$mseMint / currentData$mseBayesCorr)
    dataFrame[counter,"mseBase/mseBayesCorr"] <- median(currentData$mseBase / currentData$mseBayesCorr)
    
    #plot the graph
    title <- paste ( "n=", as.character(currentData$sampleSize[1]), "  correlB1_U:",
                     currentData$correlB1_U[1])
    ggplot(currentData, mapping = aes(correlB2_U, log(mseMint/mseBayesCorr) )) +
      ggtitle(title) +  
      geom_hline(yintercept=0)  + geom_point() + geom_smooth(method = lm) + ylim(c(-0.1,0.1))
    
    # 
    # ggsave(paste(title,".pdf"))
    
  }
  counter <- counter + 1
}
dataFrame <- na.omit(dataFrame)
UU