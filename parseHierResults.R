#this helping function generates boxplot of the mae
  #generate the bplot with ggplot2
ggBplot <- function (subresults){
  library(ggplot2)
  pdfname <- paste("results/plot","_",dset,"_",fmethod,".pdf",sep = "")
  resLenght <- length(subresults$mseBase)
  mse <- rbind(matrix(subresults$mseCombMintShr), matrix(subresults$mseBayesShr), matrix(subresults$mseBase))
  label <-  factor(rbind(matrix(rep("MinT",resLenght)),matrix(rep("Bayes",resLenght)),
                         matrix(rep("Base",resLenght))),
                   levels = c("MinT","Bayes","Base"))
  
  dataPlot <- as.data.frame(mse)
  dataPlot$label <- label
  currentPlot <- ggplot(dataPlot, aes(x = label, y = mse)) + geom_boxplot()  +
    stat_boxplot(geom = "errorbar", width = 0.5) +  #draw the whiskers
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Mse") + 
    ggtitle(paste (dset, "(",fmethod,")"))
  
  scaling <- 1.8 #to avoid large outliers that make the boxplot unreadable
  if (dset=="tourism"){
    scaling<- 1.1  
  }
  else if (fmethod=="ets"){
    scaling<- 3 
  }
  
  
  ylim1 = boxplot.stats((dataPlot$V1))$stats[c(1, 5)]
  currentPlot = currentPlot + coord_cartesian(ylim = ylim1*scaling)  #+ geom_hline(yintercept = 0, color='darkblue', linetype="dashed")
  print(currentPlot)
  ggsave(pdfname, width = 4, height = 3)
}

parseHierResults <- function (dset){
  #parse the results of hierarchical non-temporal reconciliation
  #readt the mse, extract the proportion of favorable signs and the produces the boxplot
  library(readr)
  source('bayesianSignedRank.R')
  results <- read_csv(paste("results/mse_",dset,".csv",sep=""))
  results <- unique(results) #because some experiements on the cluster are duplicated
  fmethods <- unique(results$fmethod)
  horizons <- unique(results$h)
  #we realize that we are only interested in h up to 4
  horizons <- horizons[horizons<5]
  
  configs <- length(fmethods) * length(horizons) 
  
  #we need first to instantiate the data frame with placeholder values, and then we fill the correct values
  comparison <- data.frame(
    cases = rep(fmethods[1],configs),
    h = rep(1,configs),
    fmethod=rep(fmethods[1],configs),
    # medianBaseMint=rep(-1,configs),
    medianBaseBayesShr=rep(-1,configs),
    # medianBaseBayesGlasso=rep(-1,configs),
    medianMintBayesShr =rep(-1,configs),
    medianBaseMint =rep(-1,configs),
    # medianMintBayesGlasso =rep(-1,configs),
    # pValMedianMintBayesShr=rep(-1,configs),
    # pValMedianMintBayesGlasso=rep(-1,configs),
    stringsAsFactors = FALSE
  )
  
  aggrComparison <- comparison[1:length(fmethods),]
  
  
  #analysis for each h
  counter <- 1
  for (fmethod in fmethods){
    for (h in horizons){
      print(paste(fmethod,dset))
      comparison$fmethod[counter] <- fmethod
      idx = results$fmethod==fmethod & results$h==h 
      if (sum(idx)>0){
        subresults <- results[idx,]
        comparison$cases[counter] <- sum(idx)
        comparison$fmethod[counter] <- fmethod
        comparison$h[counter] <- h
        # comparison$medianBaseMint[counter] <- median(subresults$mseBase / subresults$mseCombMintShr)
        comparison$medianBaseBayesShr[counter] <- round ( median(subresults$mseBase / subresults$mseBayesShr), digits = 2)
        # comparison$medianBaseBayesGlasso[counter] <- median(subresults$mseBase / subresults$mseBayesGlasso)
        comparison$medianMintBayesShr[counter] <- round ( median(subresults$mseCombMintShr / subresults$mseBayesShr), digits = 2)
        comparison$medianBaseMint[counter] <- round ( median(subresults$mseBase / subresults$mseCombMintShr), digits = 2)
        # comparison$medianMintBayesGlasso[counter] <- median(subresults$mseCombMintShr / subresults$mseBayesGlasso)
        # comparison$pValMedianMintBayesShr[counter] <- wilcox.test(log(subresults$mseCombMintShr/ subresults$mseBayesShr),
                                                                  # alternative="less")$p.value
        # comparison$pValMedianMintBayesGlasso[counter] <- wilcox.test(log(subresults$mseCombMintShr / subresults$mseBayesGlasso), 
                                                                     # alternative="less")$p.value
      }
      counter <- counter + 1
    }
}
filename=paste("results/summaryEachH_",dset,".csv",sep="")
write.table(comparison,file=filename,sep=",",row.names = FALSE)


#creation of the aggregated results and related graphs
# analysis aggregated over h
 counter <- 1
for (fmethod in fmethods){
    aggrComparison$fmethod[counter] <- fmethod
    idx = results$fmethod==fmethod
    # if (sum(idx)>0){
    #   subresults <- results[idx,]
    #   aggrComparison$cases[counter] <- sum(idx)
    #   aggrComparison$fmethod[counter] <- fmethod
    #   aggrComparison$medianBaseMint[counter] <- median(subresults$mseBase / subresults$mseCombMintShr)
    #   aggrComparison$medianBaseBayesShr[counter] <- median(subresults$mseBase / subresults$mseBayesShr)
    #   aggrComparison$medianMintBayesShr[counter] <- median(subresults$mseCombMintShr / subresults$mseBayesShr)
    # }
    # counter <- counter + 1
    ggBplot(results[idx,])
  }
 
}
