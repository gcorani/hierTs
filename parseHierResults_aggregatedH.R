parseHierResults_aggregatedH <- function (dset){
  #parse the results of hierarchical non-temporal reconciliation
  #readt the mse, extract the proportion of favorable signs and the produces the boxplot
  library(readr)
  results <- read_csv(paste("results/mseHierReconc",dset,".csv",sep=""))
  results <- unique(results) #because some experiements on the cluster are duplicated
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
      print(paste(fmethod,dset))
      favorableProps$dset[counter] <- dset
      favorableProps$fmethod[counter] <- fmethod
      idx = results$fmethod==fmethod & results$dset==dset 
      if (sum(idx)>0){
        subresults <- results[idx,]
        favorableProps$propBeatBase[counter] <- mean (subresults$mseBase>subresults$mseBayes)
        favorableProps$propBeatMint[counter] <- mean (subresults$mseCombMint>subresults$mseBayes)
        favorableProps$medianBaseBayes[counter] <- median(subresults$mseBase / subresults$mseBayes)
        favorableProps$medianMintBayes[counter] <- median(subresults$mseCombMint / subresults$mseBayes)
        
        favorableProps$propCorrBeatBase[counter] <- mean (subresults$mseBase>subresults$mseBayesCorr)
        favorableProps$propCorrBeatMint[counter] <- mean (subresults$mseCombMint>subresults$mseBayesCorr)
        favorableProps$medianBaseBayesCorr[counter]  <- median(subresults$mseBase / subresults$mseBayesCorr)
        favorableProps$medianMintBayesCorr[counter] <- median(subresults$mseCombMint / subresults$mseBayesCorr)
        
        #generate the bplot with ggplot2
        library(ggplot2)
        pdfname <- paste("results/plot","_",dset,"_",fmethod,".pdf",sep = "")
        denom <- subresults$mseBase 
        resLenght <- length(subresults$mseBase)
        #old code, 3 models
        # relMse <- rbind(matrix(subresults$mseCombMint/denom), matrix(subresults$mseBayes/denom), matrix(subresults$mseBayesCorr/denom))
        # label <-  factor(rbind(matrix(rep("Mint",resLenght)),matrix(rep("Bayes",resLenght)),matrix(rep("Bayes (corr)",resLenght))),
                         # levels = c("Mint","Bayes","Bayes (corr)"))
        #new code, 2 models (minT and Bayes corr)
        relMse <- rbind(matrix(subresults$mseCombMint/denom), matrix(subresults$mseBayesCorr/denom))
        label <-  factor(rbind(matrix(rep("Mint",resLenght)),matrix(rep("Bayes (corr)",resLenght))),
                         levels = c("Mint","Bayes (corr)"))
        
        dataPlot <- as.data.frame(relMse)
        dataPlot$label <- label
        currentPlot <- ggplot(dataPlot, aes(x = label, y = log10(relMse))) + geom_boxplot()  +
          stat_boxplot(geom = "errorbar", width = 0.5) +  #draw the whiskers
          scale_x_discrete(name = "") +
          scale_y_continuous(name = "Log10 ( mse / mse(base) ) ")
        
        scaling <- 1.8 #to avoid large outliers that make the boxplot unreadable
        if (dset=="tourism"){
          scaling<- 1.1  
        }
        else if (fmethod=="ets"){
          scaling<- 3 
        }
        
        
        ylim1 = boxplot.stats(log(dataPlot$V1))$stats[c(1, 5)]
        currentPlot = currentPlot + coord_cartesian(ylim = ylim1*scaling)  + geom_hline(yintercept = 0, color='darkblue')
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
