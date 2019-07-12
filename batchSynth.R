#!/opt/R/bin/Rscript
args <- commandArgs (TRUE)
correl <- as.numeric(args[1])
n <- as.integer(args[2])

setwd("/homeb/corani/hierTsCode")
source("/homeb/corani/hierTsCode/hierRec.R")

library(readr)
library(ggplot2)
 reps <- 5000
 reps <- 20
for (i in 1:reps){
  hierRec(dset="synthetic", seed=i, synth_n = n, synthCorrel = correl)
}
#we need to load the same file in which hierRec has written the results
dset <- paste0("synthetic","_correl",correl,"_n",n)
filename <- paste0("results/",dset,".csv",sep="")
currentData <- read_csv(filename)
currentData <- na.omit(currentData)

#plot the graph
title <- paste0( "plot_n", as.character(currentData$sampleSize[1]), "_correlB1_U_", 
                 currentData$correlB1_U[1])

ggplot(currentData, mapping = aes(correlB2_U, log(mseMintSample/mseBayesSample) )) +
  ggtitle(title) +  
  geom_hline(yintercept=0)  + geom_point() + geom_smooth(method = lm) + ylim(c(-0.1,0.1)) 
ggsave(paste0("results/sampleCovar_",title,".pdf"))

ggplot(currentData, mapping = aes(correlB2_U, log(mseCombMintShr/mseBayesGlasso) )) +
  ggtitle(title) +  
  geom_hline(yintercept=0)  + geom_point() + geom_smooth(method = lm) + ylim(c(-0.1,0.1)) 
ggsave(paste0("results/shrCovar_",title,".pdf"))


#summarize and save the results
dataFrame <- data.frame(
  currentData$fmethod[1],
  currentData$sampleSize[1],
  currentData$correlB1_U[1],
  median(currentData$mseMintSample/currentData$mseBayesSample,  na.rm = TRUE),
  median(currentData$mseCombMintShr/currentData$mseBayesGlasso,  na.rm = TRUE),
  median(currentData$mseBayesDiag/currentData$mseBayesGlasso, na.rm = TRUE),
  median(currentData$mseCombMintShr/currentData$mseBase, na.rm = TRUE)
)

colnames(dataFrame) <- c("fmethod", "sampleSize", "correlB1_U",  "mintSample/BayesSample",
                         "MintShr/BayesGlasso", "BayesDiag/BayesGlass", "Mint/Base")

filename <- "results/summarySyntheticx.csv"
writeNames <- TRUE
if(file.exists(filename)){
  writeNames <- FALSE
}
write.table(dataFrame, file=filename, append = TRUE, sep=",", row.names = FALSE, col.names = writeNames)

