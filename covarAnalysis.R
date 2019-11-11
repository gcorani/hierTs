
covarAnalysis <- function (dset, h=1, fmethod="ets", iTest=1, seed=0)
{
  #The hierTs data set can be ("tourism","infantgts", "synthetic") 
  #fmethod can be "ets" or "arima"
  #iTest allows to parallelize many training/test  with different splits (iTest is comprised between 1 and 50 and controls the separation between training and test) 
  #synth_n and synthCorrel are used only when generating synthetic data (synth_n: number of time points, synthCorrel: correlation between the two bottom time series.)
  #seed is especially important when you run synthetic experiments, to make sure you get different data in each experimetn
  #howManyBottom controls how many bottom synthetic time series (supported: 2 or 4)
  #covariance can be either "sam" of "shr".
  #"sam" implies that the sample covariance is used both by minT and Bayes
  #"shr" implies that the shrunken covariance and the glasso are used respectively by minT and Bayes
  
  library(hts)
  library(huge)#covariance matrix via glasso
  library(SHIP)#shrinkage of covarianca matrix
  source("loadTourism.R")
  set.seed(seed)
  
  if (dset=="tourism"){
    hierTs <- loadTourism()
  }
  else if (dset=="infantgts"){
    hierTs <- infantgts
  }
 
  
  testSize <- 50
  
  
  #if h=1, the possible preds are the whole test size lenght; 
  #if h=2, the possible preds are the (test size lenght -1); etc.
  possiblePreds <- testSize - h + 1
  
  if (iTest>possiblePreds){
    stop("iTest>possiblePreds")
  }
  
  #here the experiment starts
  timeIdx             <- time(hierTs$bts[,1])
  endTrain            <- length(timeIdx) - h - (iTest - 1)
  train               <- window(hierTs, start = timeIdx[1], end = timeIdx[endTrain] )
  test                <- window(hierTs, start =timeIdx[endTrain +1], end=timeIdx[endTrain + h])
  
  
  #recompute predictions to be  accessed by the Bayesian method
  allTsTrain <- allts(train)
  numTs <- ncol(allTsTrain)
  alpha <- 0.2
  sigma <- vector(length = numTs)
  preds <- vector(length = numTs)
  #the residuals for the model fitted on each time series
  residuals <- matrix(nrow=dim(allTsTrain)[1], ncol = numTs)
  fitted <- matrix(nrow=dim(allTsTrain)[1], ncol = numTs)
  actual <- matrix(nrow=dim(allTsTrain)[1], ncol = numTs)
  
  #compute, for each  ts, predictions and sigma (h-steps ahead) 
  for (i in 1:numTs){
    if (fmethod=="ets"){
      model <- ets(allTsTrain[,i], additive.only = TRUE)
    }
    else if (fmethod=="arima"){
      model <- auto.arima(allTsTrain[,i])
    }
    tmp <- forecast(model, h=h, level=1-alpha)
    residuals[,i] <- model$residuals #we could store model$residuals if we allowed  multiplicative errors
    fitted[,i] <- model$fitted 
    actual[,i] <- model$x
    preds[i] <- tmp$mean[h]
    #the reconciliation matrix does always refer to h=1
    sigma[i] <- abs ( (tmp$mean[1] - tmp$upper[1])  / (qnorm(alpha / 2)) )
  }
  mseBase =  mean  ( (allts(test)[h,] - preds)^2 )
  
  
  #shrink to the diagonal
  library(stats)
  library(corrplot)
  
  pdf(file = "CovErrorsShr.pdf")
  shCov <-  shrink.estim(residuals, tar=build.target(residuals,type="D"))[[1]]
  shCor <-  cov2cor(shCov)  
  
  if (dset=="infantgts"){
    colnames(shCor) <- colnames(allTsTrain)
    rownames(shCor) <- colnames(allTsTrain)
  }
  

  
  corrplot(shCor, method = "circle", title = "residual correlation after shrinkage")
  dev.off()
  
  
  S <- smatrix(train)
  bottomIdx <- seq( nrow(S) - ncol(S) +1, nrow(S))
  upperIdx <- setdiff(1:nrow(S),bottomIdx)
  pdf(file = paste0(dset,"CovLowerErrorsShr.pdf"))
  corrplot(shCor[bottomIdx,bottomIdx], method = "circle", title = "botton residuals")
  dev.off()
  
  pdf(file = paste0(dset,"CovUpperErrorsShr.pdf"))
  corrplot(shCor[upperIdx,upperIdx], method = "circle", title = "upper residuals")
  dev.off()
  
  pdf(file = paste0(dset,"CovUpperLowerErrorsShr.pdf"))
  corrplot(shCor[bottomIdx,upperIdx], method = "circle", title = "upper - lower")
  dev.off()
  
  
  pdf(file = paste0(dset,"CovY.pdf"))
  corY <- cor(actual)
  colnames(corY) <- colnames(allTsTrain)
  rownames(corY) <- colnames(allTsTrain)
  corrplot(corY, method = "circle", title = "Y correlations")
  dev.off()
  
  
  
    
 
}

