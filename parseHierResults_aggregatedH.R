parseHierResults_aggregatedH <- function (dset){
  #parse the results of hierarchical non-temporal reconciliation
  #readt the mse, extract the proportion of favorable signs and the produces the boxplot
  library(readr)
  results <- read_csv(paste("results/mseHierReconc",dset,".csv",sep=""))
  fmethods <- unique(results$fmethod)
  dsets <- unique(results$dset)
  configs <- length(fmethods) * length(dsets) 
  #the 7 fields are fmethod, dset, prop against each method
  # header <- c("dset","fmethod","h","propBeatBase","propBeatBu","propBeatComb",
  # "propBeatCombWls","propBeatMint")
  # favorableProps <- matrix(nrow = configs, ncol = length(header))
  
  #we need first to instantiate the data frame with placeholder values, and then we fill the correct values
  favorableProps <- data.frame(dset=rep(dsets[1],configs),
                               fmethod=rep(fmethods[1],configs),
                               propBeatBase=rep(-1,configs),
                               propBeatMint=rep(-1,configs),
                               medianBaseBayes=rep(-1,configs),
                               medianBaseMint=rep(-1,configs),
                               stringsAsFactors = FALSE
  )
  
  counter <- 1
  for (dset in dsets){
    for (fmethod in fmethods){
      print(paste(fmethod,dset))
      favorableProps$dset[counter] <- dset
      favorableProps$fmethod[counter] <- fmethod
      idx = results$fmethod==fmethod & results$dset==dset 
      if (sum(idx)>0){
        subresults <- results[idx,]
        favorableProps$propBeatBase[counter] <- mean (subresults$mseBase>subresults$mseBayes)
        favorableProps$propBeatMint[counter] <- mean (subresults$mseCombMint>subresults$mseBayes)
        favorableProps$medianBaseBayes[counter] <- median(subresults$mseBase / subresults$mseBayes)
        favorableProps$medianBaseMint[counter] <- median(subresults$mseCombMint / subresults$mseBayes)
        
        #generate the bplot with ggplot2
        library(ggplot2)
        pdfname <- paste("results/GGPLOThier","_",dset,"_",fmethod,".pdf",sep = "")
        denom <- subresults$mseBase 
        resLenght <- length(subresults$mseBase)
        relMse <- rbind(matrix(subresults$mseCombMint/denom), matrix(subresults$mseBayes/denom))
        
        # relMse <- rbind(matrix(subresults$mseBu/denom), matrix(subresults$mseCombWls/denom),
        #                 matrix(subresults$mseCombMint/denom), matrix(subresults$mseBayes/denom))
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
  filename=paste("results/summary",dset,".csv",sep="")
  write.table(favorableProps,file=filename,sep=",",row.names = FALSE)
  return (favorableProps)
}