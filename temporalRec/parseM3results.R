parseTemporalHierResults <- function (type="monthly", fmethod="ets"){
  library(readr)
  source(bayesianSignedRank.R)
  filename <- paste("results/temporalHier","_",type,"_",fmethod,".csv",sep = "")
  results <- read_csv(filename)
  results <- unique(results)
  
  configs<-1
  #correlated method is not used in temporal reconciliation
  favorableProps <- data.frame(
    fmethod=rep(fmethod,configs),
    # propBeatBase=rep(-1,configs),
    # propBeatThief=rep(-1,configs),
    medianBaseBayes=rep(-1,configs),
    medianThiefBayes=rep(-1,configs),
    probThiefBeatBayes=rep(-1,configs),
    probBayesBeatThief=rep(-1,configs),
    stringsAsFactors = FALSE
  )
  counter<-1
  favorableProps$propBeatBase[counter] <- mean (results$mseBase>results$mseBayes)
  favorableProps$propBeatThief[counter] <- mean (results$mseThief>results$mseBayes)
  favorableProps$medianBaseBayes[counter] <- median(results$mseBase / results$mseBayes)
  favorableProps$medianThiefBayes[counter] <- median(results$mseThief / results$mseBayes)
  favorableProps$pValBayesThief[counter] <- wilcox.test(log(results$mseThief / results$mseBayes))$p.value
  bayesThief <- bayesianSignedRank(log(subresults$mseThief/subresults$mseBayes),
                                  rope_min = log(1),

                                    favorableProps$probThiefBeatBayes[counter] <- bayesThief$probSmaller
  
  #generate the bplot with ggplot2
  library(ggplot2)
  pdfname <- paste("results/plotTemporal","_",type,"_",fmethod,".pdf",sep = "")
  denom <- results$mseBase 
  resLenght <- length(results$mseThief/denom)
  relMse <- rbind(matrix(results$mseThief/denom), matrix(results$mseBayes/denom))
  label <-  factor(cbind(matrix(rep("Thief",resLenght)),matrix(rep("Bayes-diag",resLenght))),
                   levels = c("Thief","Bayes-diag"))
  dataPlot <- as.data.frame(relMse)
  dataPlot$label <- label
  currentPlot <- ggplot(dataPlot, aes(x = label, y = log(relMse))) + geom_boxplot()  +
    stat_boxplot(geom = "errorbar", width = 0.5) +  #draw the whiskers
    scale_x_discrete(name = "") +
    scale_y_continuous(name = "Log10 (relative mse) ") 
  ylim1 = boxplot.stats(log(dataPlot$V1))$stats[c(1, 5)]
  currentPlot = currentPlot + coord_cartesian(ylim = ylim1*1.1)  + geom_hline(yintercept = 0, color='darkblue', linetype="dashed")
  print(currentPlot)
  ggsave(pdfname, width = 4, height = 3)
  
  # print(paste(type," ",fmethod))
  # print("summary mseBase / mseBayes")
  # print(summary(results$mseBase/results$mseBayes))
  # print("summary mseThief / mseBayes")
  # print(summary(results$mseThief/results$mseBayes))
  
  return (favorableProps)
}
