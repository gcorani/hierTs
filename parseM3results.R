parseM3Results <- function (type="monthly", fmethod="ets"){
  library(readr)
  filename <- paste("results/temporalHier","_",type,"_",fmethod,".csv",sep = "")
  results <- read_csv(filename)
  if (type=="weekly"){
     warning("the AE experiment is only over 13 time series") 
    freqs <- c("_Weekly","_2-Weekly","_4-Weekly","_Quarterly","_Biannual","_Annual")
  }
  if (type=="monthly"){
    freqs <- c("_Monthly","_2-Monthly","_4-Monthly","_Biannual","_Annual")
  }
  else if (type=="quarterly"){
    freqs <- c("Quarterly","_Biannual","_Annual")
  }
  names <- colnames(results)
  baseMethods <- c("Bu","Thief")
  myMethod <-"Bayes"
  
  comparisons <- length(freqs) * length (baseMethods)
  improvIndicator <- matrix(0, ncol = comparisons, nrow = nrow(results))
  improvNames <- vector(length = comparisons)
  favorableSign <- vector(length = comparisons)
  pValueSign <- vector(length = comparisons)
  meanImprovement <- vector(length = comparisons)
  isMseImproved  <- vector(length = comparisons)
  mseImprovement  <- vector(length = comparisons)
  counter <- 1
  
  for (method in baseMethods) {
    for (freq in freqs){
      improvNames[counter] <- paste(method,freq)
      colBase <- intersect(grep(method,names) , grep(freq,names))
      colMy <- intersect(grep(myMethod,names) , grep(freq,names))
      
      favorableSign[counter] <- mean (results[,colBase] > results[,colMy] )
      
      #computing the improvIndicator
      #unfortunaltely, this yields a single-columns *data frame* which we then have to cast as a vector
      currentImprovement <- ( results[,colBase] - results[,colMy] ) / ((results[,colBase]  + results[,colMy]) /2 )
      #casting as vector
      currentImprovement <- currentImprovement[,1]
      improvIndicator[,counter] <- currentImprovement
      meanImprovement[counter] <- mean (currentImprovement)
      counter <- counter + 1
    }
  }
  
  improvIndicator <- as.data.frame(improvIndicator)  
  colnames(improvIndicator) <- improvNames
  
  meanImprovement <-  as.data.frame(t(meanImprovement))  
  colnames(meanImprovement) <- improvNames
  
  favorableSign <- as.data.frame(t(favorableSign)) 
  colnames(favorableSign) <- improvNames
  
  #code used for AE only and uncommented
  # denom <- sum(results$mseBase) + sum(results$mseBu) + sum(results$mseThief) + sum(results$mseBayes)
  # denom <- denom / 4
  # relMseBase <- sum(results$mseBase)/denom
  # relMseBu <- sum(results$mseBu)/denom
  # relMseThief <- sum(results$mseThief)/denom
  # relMseBayes <- sum(results$mseBayes)/denom

  
  #generate the bplot
  pdfname <- paste("results/temporalHier","_",type,"_",fmethod,".pdf",sep = "")
  pdf(pdfname) 
  denom <- results$mseBase 
  a <-  cbind(results$mseBu/denom, 
              results$mseThief/denom, results$mseBayes/denom)
  boxplot(log10(a),names=c("bu","thief","bayes"), outline=TRUE, ylab="Relative MSE (log10)")
  dev.off()
  
  #generate the bplot with ggplot2
  library(ggplot2)
  pdfname <- paste("results/GGPLOTtemporalHier","_",type,"_",fmethod,".pdf",sep = "")
  denom <- results$mseBase 
  resLenght <- length(results$mseThief/denom)
  relMse <- rbind(matrix(results$mseBu/denom), matrix(results$mseThief/denom), matrix(results$mseBayes/denom))
  label <-  factor(rbind(matrix(rep("Bu",resLenght)),matrix(rep("Thief",resLenght)),matrix(rep("Bayes",resLenght))),
                   levels = c("Bu","Thief","Bayes"))
  dataPlot <- as.data.frame(relMse)
  dataPlot$label <- label
  currentPlot <- ggplot(dataPlot, aes(x = label, y = log10(relMse))) + geom_boxplot()  +
                  stat_boxplot(geom = "errorbar", width = 0.5) +  #draw the whiskers
                           scale_x_discrete(name = "") +
    scale_y_continuous(name = "Log10 (MSE / MSE base) ")
  print(currentPlot)
  ggsave(pdfname, width = 4, height = 3)
  # boxplot(log10(a),names=c("bu","thief","bayes"), outline=TRUE, ylab="Relative MSE (log10)")
  # dev.off()
  
  # Scatter plot, commented out
  # pdfname <- paste("results/scatterBayesThier","_",type,"_",fmethod,".pdf",sep = "")
  # pdf(pdfname)
  # if (type=="monthly"){
  # plot(results$mseThief,results$mseBayes, xlab = "mse (thief)",
  #      ylab="mse(Bayes)", xlim = c(0,1e+9), ylim= c(0,1e+9))
  # }
  # else if (type=="quarterly"){
  #   plot(results$mseThief,results$mseBayes, xlab = "mse (thief)",
  #        ylab="mse(Bayes)", xlim = c(0,1e+8), ylim= c(0,1e+8))
  # }
  # abline(0,1)
  # dev.off()
  # 
  #test comparing thief and bayes
  wilcoxPval <- wilcox.test(a[,2],a[,3],paired = TRUE, alternative = "greater")
  # ttestPval <-  t.test(a[,2],a[,3],paired = TRUE, alternative = "greater")
  
  return (list("favorableSign" = favorableSign, "meanImprovement" = meanImprovement,
               "wilcoxPval" = wilcoxPval, "ttestPval"=ttestPval, "improvIndicator"=improvIndicator ) )
}