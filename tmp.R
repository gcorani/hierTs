for (i in 1:200) { 
  print(i)
  hierRec("synthetic", seed = i, synth_n = 5, synthCorrel = 0.01)
}


currentData <- read_csv(filename)
#plot the graph
title <- paste0( "plot_n", as.character(currentData$sampleSize[1]), "_correlB1_U_", 
                 currentData$correlB1_U[1])

library(ggplot2)
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
    median(currentData$mseMintSample/currentData$mseBayesSample),
    median(currentData$mseCombMintShr/currentData$mseBayesGlasso)
)

colnames(dataFrame) <- c("fmethod", "sampleSize", "correlB1_U",  "mintSample/BayesSample",
                         "MintShr/BayesGlasso")

filename <- "results/summary.csv"
writeNames <- TRUE
if(file.exists(filename)){
  writeNames <- FALSE
}
write.table(dataFrame, file=filename, append = TRUE, sep=",", row.names = FALSE, col.names = writeNames)
  
