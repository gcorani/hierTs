parseHierResults <- function (){
  #parse the results of hierarchical non-temporal reconciliation
  #readt the mse, extract the proportion of favorable signs and the produces the boxplot
  library(readr)
  results <- read_csv("results/mseHierReconc.csv")
  fmethods <- unique(results$fmethod)
  dsets <- unique(results$dset)
  horizons <- unique(results$h)
  configs <- length(fmethods) * length(dsets) * length(horizons)
  #the 7 fields are fmethod, dset, prop against each method
  header <- c("dset","fmethod","h","propBeatBase","propBeatBu","propBeatComb",
              "propBeatCombWls","propBeatMint")
  favorableProps <- matrix(nrow = configs, ncol = length(header))
  colnames(favorableProps) <- header
  
  counter <- 1
  for (dset in dsets){
    for (fmethod in fmethods){
      for (h in horizons){
        print(paste(fmethod,dset,h))
        favorableProps[counter,"dset"] <- dset
        favorableProps[counter,"fmethod"] <- fmethod
        favorableProps[counter,"h"] <- h
        idx = results$fmethod==fmethod & results$dset==dset & results$h==h
        if (sum(idx)>0){
          subresults <- results[idx,]
          favorableProps[counter,"propBeatBase"] <- mean (subresults$mseBase>subresults$mseBayes)
          favorableProps[counter,"propBeatBu"] <- mean (subresults$mseBu>subresults$mseBayes)
          favorableProps[counter,"propBeatComb"] <- mean (subresults$mseComb>subresults$mseBayes)
          favorableProps[counter,"propBeatCombWls"] <- mean (subresults$mseCombWls>subresults$mseBayes)
          favorableProps[counter,"propBeatMint"] <- mean (subresults$mseCombMint>subresults$mseBayes)
          #generate the bplot
          # pdfname <- paste("results/hier","_",dset,"_",fmethod,"_h",h,".pdf",sep = "")
          # pdf(pdfname) 
          # denom <- subresults$mseBase 
          # a <-  cbind(subresults$mseBu/denom, subresults$mseComb/denom, subresults$mseCombWls/denom, subresults$mseCombMint/denom, subresults$mseBayes/denom)
          # boxplot(log10(a),names=c("bu","comb","combWLS","mint","Bayes"), outline=TRUE, ylab="Relative MSE (log10)")
          # dev.off()
          
          #generate the bplot with ggplot2
          library(ggplot2)
          pdfname <- paste("results/GGPLOThier","_",dset,"_",fmethod,"_h",h,".pdf",sep = "")
          denom <- subresults$mseBase 
          resLenght <- length(subresults$mseBase)
          # relMse <- rbind(matrix(subresults$mseBu/denom), matrix(subresults$mseComb/denom), matrix(subresults$mseCombWls/denom),
          # matrix(subresults$mseCombMint/denom), matrix(subresults$mseBayes/denom))
          relMse <- rbind(matrix(subresults$mseBu/denom), matrix(subresults$mseCombWls/denom),
                          matrix(subresults$mseCombMint/denom), matrix(subresults$mseBayes/denom))
          # label <-  factor(rbind(matrix(rep("Bu",resLenght)),matrix(rep("Comb",resLenght)),matrix(rep("CombWls",resLenght)),
          # matrix(rep("Mint",resLenght)),matrix(rep("Bayes",resLenght))),
          # levels = c("Bu","Comb","CombWls","Mint","Bayes"))
          label <-  factor(rbind(matrix(rep("Bu",resLenght)),matrix(rep("CombWls",resLenght)),
                                 matrix(rep("Mint",resLenght)),matrix(rep("Bayes",resLenght))),
                           levels = c("Bu","CombWls","Mint","Bayes"))
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
  }
  
  return (favorableProps)
}