parseM3Results <- function (type="monthly"){
  library(readr)
  filename <- paste("temporalHier","_",type,"_",fmethod,".csv",sep = "")
  results <- read_csv(filename)
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
  pValueTtest <- vector(length = comparisons)
  meanImprovement <- vector(length = comparisons)
  counter <- 1
  
  for (method in baseMethods) {
    for (freq in freqs){
      improvNames[counter] <- paste(method,freq)
      colBase <- intersect(grep(method,names) , grep(freq,names))
      colMy <- intersect(grep(myMethod,names) , grep(freq,names))
      
      favorableSign[counter] <- mean (results[,colBase] > results[,colMy] )
      pValueSign[counter] <- binom.test (sum(results[,colBase] > results[,colMy]), nrow(results))$p.value
      
      #computing the improvIndicator
      #unfortunaltely, this yields a single-columns *data frame* which we then have to cast as a vector
      currentImprovement <- ( results[,colBase] - results[,colMy] ) / ((results[,colBase]  + results[,colMy]) /2 )
      #casting as vector
      currentImprovement <- currentImprovement[,1]
      improvIndicator[,counter] <- currentImprovement
      meanImprovement[counter] <- mean (currentImprovement)
      pValueTtest[counter] <- t.test (currentImprovement, alternative = "greater")$p.value
      counter <- counter + 1
    }
  }
  
  improvIndicator <- as.data.frame(improvIndicator)  
  colnames(improvIndicator) <- improvNames
  
  meanImprovement <-  as.data.frame(t(meanImprovement))  
  colnames(meanImprovement) <- improvNames
  
  favorableSign <- as.data.frame(t(favorableSign)) 
  colnames(favorableSign) <- improvNames
  
  pValueSign <- as.data.frame(t(pValueSign)) 
  colnames(pValueSign) <- improvNames
  
  pValueTtest <- as.data.frame(t(pValueTtest)) 
  colnames(pValueTtest) <- improvNames
  
  return (list("favorableSign" = favorableSign, "pValueTtest" = pValueTtest, "pValueSign" = pValueSign, "meanImprovement" = meanImprovement,
               "improvIndicator"=improvIndicator ) )
}