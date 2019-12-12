parseHierResults_synthetic <- function (){
  #parse the results of hierarchical non-temporal reconciliation on synthetic dataa set
  #readt the mse, extract the proportion of favorable signs and the produces the boxplot
  library(readr)
  results <- read_csv(paste("results/mseHierReconcsynthetic.csv"))
  keys <- unique(cbind(results$fmethod, results$synth_n, results$synthCorrel, results$howManyBottom)) #because some experiements on the cluste?r are duplicated
  configs = dim(keys)[1]
  
  #init the lists of results
  summaryEachConfig <- data.frame(
    howManyBottom=rep(-1,configs),
    numExperiments=rep(-1,configs),
    sampleSize=rep(-1,configs),
    fmethod=rep("ets",configs),
    bottomCorrelation=rep(-1,configs),
    propDiagBeatBase=rep(-1,configs),
    propDiagBeatMint=rep(-1,configs),
    propCorrBeatBase=rep(-1,configs),
    propCorrBeatMint=rep(-1,configs),
    medianBaseBayes=rep(-1,configs),
    medianBaseBayesCorr=rep(-1,configs),
    medianMintBayes=rep(-1,configs),
    medianMintBayesCorr=rep(-1,configs),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:dim(keys)[1]){
    key=keys[i,]
    idx = ( results$fmethod == key[1] &
              results$synth_n == key[2] &
              results$synthCorrel == key[3] & 
              results$howManyBottom == key[4])
    subResults <- results[idx,]
    summaryEachConfig$fmethod[i] <- key[1]
    summaryEachConfig$sampleSize[i] <- key[2]
    summaryEachConfig$numExperiments[i] <- sum(idx)
    summaryEachConfig$bottomCorrelation[i] <- key[3]
    summaryEachConfig$howManyBottom[i] <- key[4]
    summaryEachConfig$propDiagBeatBase[i] <- mean(subResults$mseBayes<subResults$mseBase)
    summaryEachConfig$propCorrBeatBase[i] <- mean(subResults$mseBayesCorr<subResults$mseBase)
    summaryEachConfig$propDiagBeatMint[i] <- mean(subResults$mseBayes<subResults$mseCombMint)
    summaryEachConfig$propCorrBeatMint[i] <- mean(subResults$mseBayesCorr<subResults$mseCombMint)
    summaryEachConfig$medianBaseBayes[i] <- median(subResults$mseBase/subResults$mseBayes)
    summaryEachConfig$medianBaseBayesCorr[i] <- median(subResults$mseBase/subResults$mseBayesCorr)
    summaryEachConfig$medianMintBayes[i] <- median(subResults$mseCombMint/subResults$mseBayes)
    summaryEachConfig$medianMintBayesCorr[i] <- median(subResults$mseCombMint/subResults$mseBayesCorr)
  }
  
options(digits=2) 
filename = paste("results/summarySynthetic.csv")
write.table(summaryEachConfig,file=filename,sep=",",row.names = FALSE)
return (summaryEachConfig)
}
