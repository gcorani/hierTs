#!/opt/R/bin/Rscript
args <- commandArgs (TRUE)
n <- as.integer(args[1])

setwd("/homeb/corani/hierTsCode")
source("/homeb/corani/hierTsCode/hierRec.R")

library(readr)

# reps <- 2500
reps <- 20
for (i in 1:reps){
  hierRec(dset="syntheticLarge", seed=i, synth_n = n)
}
#we need to load the same file in which hierRec has written the results
dset <- paste0("largeSynthetic_n",n)
filename <- paste0("results/",dset,".csv",sep="")
currentData <- read_csv(filename)
currentData <- na.omit(currentData)


#summarize and save the results
dataFrame <- data.frame(
  currentData$fmethod[1],
  currentData$sampleSize[1],
  median(currentData$mseMintSample/currentData$mseBayesSample),
  median(currentData$mseCombMintShr/currentData$mseBayesGlasso),
  median(currentData$mseBayesDiag/currentData$mseBayesGlasso),
  median(currentData$mseCombMintShr/currentData$mseBase)
)

colnames(dataFrame) <- c("fmethod", "sampleSize",  "mintSample/BayesSample",
                         "MintShr/BayesGlasso", "BayesDiag/BayesGlass", "Mint/Base")

filename <- "results/summaryLargeSynthetic.csv"
writeNames <- TRUE
if(file.exists(filename)){
  writeNames <- FALSE
}
write.table(dataFrame, file=filename, append = TRUE, sep=",", row.names = FALSE, col.names = writeNames)

