#this helping function generates boxplot of the mae
  #generate the bplot with ggplot2
ggBplot <- function (subresults, dset){
  library(ggplot2)
  fmethod <- unique(subresults$fmethod)
  pdfname <- paste("results/plot","_",dset,"_",fmethod,".pdf",sep = "")
  resLenght <- length(subresults$mseBase)
  mse <- rbind(matrix(subresults$mseCombMintShr), matrix(subresults$mseBayesShr), matrix(subresults$mseBase))
  label <-  factor(rbind(matrix(rep("MinT",resLenght)),matrix(rep("Bayes",resLenght)),
                         matrix(rep("Base",resLenght))),
                   levels = c("MinT","Bayes","Base"))
  dataPlot <- as.data.frame(mse)
  dataPlot$label <- label
  currentPlot <- ggplot(dataPlot, aes(x = label, y = sqrt(mse))) + geom_boxplot()  +
    stat_boxplot(geom = "errorbar", width = 0.5) +  #draw the whiskers
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "rmse") + 
    ggtitle(paste (dset, paste0("(",fmethod,")")))
  
  scaling <- 1.2 #to avoid large outliers that make the boxplot unreadable
  if (dset=="tourism"){
    scaling <- 1.1
  }
  
  ylim1 = boxplot.stats(sqrt(dataPlot$V1))$stats[c(1, 5)]
  currentPlot = currentPlot + coord_cartesian(ylim = ylim1*scaling)  #+ geom_hline(yintercept = 0, color='darkblue', linetype="dashed")
  print(currentPlot)
  ggsave(pdfname, width = 4, height = 3)
}

parseHierResults <- function (dset){
  #parse the results of hierarchical non-temporal reconciliation
  #readt the mse, extract the proportion of favorable signs and the produces the boxplot
  library(readr)
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
    medianBaseBayesShr=rep(-1,configs),
    medianMintBayesShr =rep(-1,configs),
    medianBaseMint =rep(-1,configs),
    rmseMint=rep(-1,configs),
    rmseBayes=rep(-1,configs),
    rmseBase=rep(-1,configs),
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
        comparison$medianBaseBayesShr[counter] <- round ( median(subresults$mseBase / subresults$mseBayesShr), digits = 2)
        comparison$medianMintBayesShr[counter] <- round ( median(subresults$mseCombMintShr / subresults$mseBayesShr), digits = 2)
        comparison$medianBaseMint[counter] <- round ( median(subresults$mseBase / subresults$mseCombMintShr), digits = 2)
        comparison$rmseMint[counter] <- median(sqrt(subresults$mseCombMintShr))
        comparison$rmseBayes[counter] <- median(sqrt(subresults$mseBayesShr))
        comparison$rmseBase[counter] <- median(sqrt(subresults$mseBase))    
      }
      counter <- counter + 1
    }
}
filename=paste("results/summaryEachH_",dset,".csv",sep="")
write.table(comparison,file=filename,sep=",",row.names = FALSE)


#creation of the aggregated results and related graphs
# analysis aggregated over h
for (fmethod in fmethods){
    idx = results$fmethod==fmethod & results$h==1
    ggBplot(results[idx,],dset)
  }
 
}
