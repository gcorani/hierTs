parseHierResults_aggregatedH <- function (){
  #parse the results of hierarchical non-temporal reconciliation
  #aggregating over different h, instead of treating them separately
  #readt the mse, extract the proportion of favorable signs and the produces the boxplot
  library(readr)
  results <- read_csv("results/mseHierReconc.csv")
  fmethods <- unique(results$fmethod)
  dsets <- unique(results$dset)
  # horizons <- unique(results$h)
  configs <- length(fmethods) * length(dsets) 
  #the 7 fields are fmethod, dset, prop against each method
  # header <- c("dset","fmethod","h","propBeatBase","propBeatBu","propBeatComb",
  # "propBeatCombWls","propBeatMint")
  # favorableProps <- matrix(nrow = configs, ncol = length(header))
  
  #we need first to instantiate the data frame with placeholder values, and then we fill the correct values
  favorableProps <- data.frame(dset=rep(dsets[1],configs),
                               # h=rep(horizons[1],configs),
                               fmethod=rep(fmethods[1],configs),
                               propBeatBase=rep(-1,configs),
                               # propBeatBu=rep(-1,configs),
                               # propBeatComb=rep(-1,configs),
                               # propBeatCombWls=rep(-1,configs),
                               propBeatMint=rep(-1,configs),
                               stringsAsFactors = FALSE
  )
  
  counter <- 1
  for (dset in dsets){
    for (fmethod in fmethods){
      # for (h in horizons){
      print(paste(fmethod,dset))
      favorableProps$dset[counter] <- dset
      favorableProps$fmethod[counter] <- fmethod
      # favorableProps$h[counter] <- h
      idx = results$fmethod==fmethod & results$dset==dset #& results$h==h
      if (sum(idx)>0){
        subresults <- results[idx,]
        favorableProps$propBeatBase[counter] <- mean (subresults$mseBase>subresults$mseBayes)
        # favorableProps$propBeatBu[counter] <- mean (subresults$mseBu>subresults$mseBayes)
        # favorableProps$propBeatComb[counter] <- mean (subresults$mseComb>subresults$mseBayes)
        # favorableProps$propBeatCombWls[counter] <- mean (subresults$mseCombWls>subresults$mseBayes)
        favorableProps$propBeatMint[counter] <- mean (subresults$mseCombMint>subresults$mseBayes)
        
        #generate the bplot with ggplot2
        library(ggplot2)
        pdfname <- paste("results/GGPLOThier","_",dset,"_",fmethod,".pdf",sep = "")
        denom <- subresults$mseBase 
        resLenght <- length(subresults$mseBase)
        relMse <- rbind(
          # matrix(subresults$mseBu/denom), matrix(subresults$mseCombWls/denom),
                        matrix(subresults$mseCombMint/denom), matrix(subresults$mseBayes/denom))
        # label <-  factor(rbind(matrix(rep("Bu",resLenght)),matrix(rep("CombWls",resLenght)),
        #                        matrix(rep("Mint",resLenght)),matrix(rep("Bayes",resLenght))),
        #                  levels = c("Bu","CombWls","Mint","Bayes"))
        label <-  factor(rbind(matrix(rep("Mint",resLenght)),matrix(rep("Bayes",resLenght))),
                         levels = c("Mint","Bayes"))
        dataPlot <- as.data.frame(relMse)
        dataPlot$label <- label
        currentPlot <- ggplot(dataPlot, aes(x = label, y = log10(relMse))) + geom_boxplot()  +
          stat_boxplot(geom = "errorbar", width = 0.5) +  #draw the whiskers
          scale_x_discrete(name = "") +
          scale_y_continuous(name = "Log10 (MSE / MSE base) ")
        print(currentPlot)
        ggsave(pdfname, width = 4, height = 3)
        
        
        
        
        
        
        counter <- counter + 1
      }
    }
  }
  write.table(favorableProps,file="results/hierFavorableProps.csv",sep=",",row.names = FALSE)
  return (favorableProps)
}