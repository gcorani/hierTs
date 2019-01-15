parseHierResults_eachH <- function (dset){
  #parse the results of hierarchical non-temporal reconciliation
  #readt the mse, extract the proportion of favorable signs and the produces the boxplot
  library(readr)
  results <- read_csv(paste("results/mseHierReconc",dset,".csv",sep=""))
  results <- unique(results) #because some experiements on the cluster are duplicated
  fmethods <- unique(results$fmethod)
  dsets <- unique(results$dset)
  horizons <- unique(results$h)
  configs <- length(fmethods) * length(dsets) * length(horizons)
  
  #we need first to instantiate the data frame with placeholder values, and then we fill the correct values
  
  favorableProps <- data.frame(dset=rep(dsets[1],configs),
                               n=rep(horizons[1],configs),
                               h=rep(horizons[1],configs),
                               fmethod=rep(fmethods[1],configs),
                               propBeatBase=rep(-1,configs),
                               propBeatMint=rep(-1,configs),
                               medianBaseBayes=rep(-1,configs),
                               medianMintBayes=rep(-1,configs),
                               propCorrBeatBase=rep(-1,configs),
                               propCorrBeatMint=rep(-1,configs),
                               medianBaseBayesCorr=rep(-1,configs),
                               medianMintBayesCorr=rep(-1,configs),
                               stringsAsFactors = FALSE
  )
  
  counter <- 1
  for (dset in dsets){
    for (fmethod in fmethods){
      for (h in horizons){
        print(paste(fmethod,dset,h))
        favorableProps$dset[counter] <- dset
        favorableProps$fmethod[counter] <- fmethod
        favorableProps$h[counter] <- h
        idx = results$fmethod==fmethod & results$dset==dset & results$h==h
        favorableProps$n[counter] <- sum(idx)
        if (sum(idx)>0){
          subresults <- results[idx,]
          favorableProps$propBeatBase[counter] <- mean (subresults$mseBase>subresults$mseBayes)
          favorableProps$propBeatMint[counter] <- mean (subresults$mseCombMint>subresults$mseBayes)
          favorableProps$medianBaseBayes[counter] <- median(subresults$mseBase / subresults$mseBayes)
          favorableProps$medianMintBayes[counter] <- median(subresults$mseCombMint / subresults$mseBayes)
          
          favorableProps$propCorrBeatBase[counter] <- mean (subresults$mseBase>subresults$mseBayesCorr)
          favorableProps$propCorrBeatMint[counter] <- mean (subresults$mseCombMint>subresults$mseBayesCorr)
          favorableProps$medianBaseBayesCorr[counter]  <- median(subresults$mseBase / subresults$mseBayesCorr)
          favorableProps$medianMintBayesCorr[counter]  <- median(subresults$mseCombMint / subresults$mseBayesCorr)
          
          counter <- counter + 1
        }
      }
    }
  }
  options(digits=2) 
  filename = paste("results/summary",dset,"EachH.csv",sep="")
  write.table(favorableProps,file=filename,sep=",",row.names = FALSE)
  return (favorableProps)
}
