#TODO: does not run on tourism, breaks on accessing residual

# #temporary for debugging
# dset="infantgts"
# h=1
# fmethod="ets"
# iTest=1
# seed=0

positiveReconc <- function (dset, h=1, fmethod="ets", iTest=1, 
                            seed=0, synth_n=100, synthCorrel=0.5)
{
  #compares truncated and gaussian reconciliation.
  #we do not run minT, as it is equal to our full reconciliation
  #The hierTs data set can be ("tourism","infantgts", "synthetic" (2 bottom time series),"syntheticLarge" (4 bottom time series)) 
  #fmethod can be "ets" or "arima"
  #iTest decised the split between  training/test  (iTest is comprised between 1 and 50 and controls the separation between training and test).
  #synth_n and synthCorrel are used only when generating synthetic data (synth_n: number of time points, synthCorrel: correlation between the two bottom time series.)
  #seed is especially important when you run synthetic experiments, to make sure you get different data in each experimetn
  
  library(hts)
  library(SHIP)
  library(scoringRules)
  library(tmvtnorm)
  library(MASS)
  source("loadTourism.R")
  source("reconcUtils.R")
  set.seed(seed)
  testSize <- 50
  hierTs <- setDset(dset)
  
  #if h=1, the possible preds are the whole test size lenght; 
  #if h=2, the possible preds are the (test size lenght -1); etc.
  possiblePreds <- testSize - h + 1
  
  if (iTest>possiblePreds){
    stop("iTest>possiblePreds")
  }
  
  #here the experiment starts
  
  timeIdx <- time(hierTs$bts[,1])
  startTrain          <- timeIdx[1]
  endTrain            <- length(timeIdx) - h - (iTest - 1)
  train               <- window(hierTs, start = startTrain, end = timeIdx[endTrain] )
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
  y <- allts(test)[h,]
  mseBase <-  mean  ( (y - preds)^2 )
  bayesPost <- bayesRecon(train, preds, covType="shr")
  mseBayesShr <-  mean  ( (y - bayesPost$y_tilde)^2 )
  
  esNormal<-energyScore(y , bayesPost, "normal")
  esTrunc<-energyScore(y , bayesPost, "trunc")
  names(esTrunc)<-paste0("trunc_",names(esTrunc))
  
  if (dset=="synthetic"){
    dataFrame <- data.frame(h, fmethod, synth_n, synthCorrel, corrB2_U, mseBase,
                              mseBayesShr, as.vector(esNormal), as.vector(esTrunc))
    colnames(dataFrame) <- c("h","fmethod","sampleSize","correlB1_U","correlB2_U",
                             "mseBase","mseBayesShr",
                             names(esNormal), names(esTrunc))
    dset <- paste0(dset,"_correl",synthCorrel,"_n",synth_n)
  }
  
  else if (dset=="syntheticLarge"){
    dataFrame <- data.frame(h, fmethod, synth_n, mseBase, mseBayesShr, mseBu, as.vector(esNormal), as.vector(esTrunc))
    colnames(dataFrame) <- c("h","fmethod","sampleSize",
                             "mseBase","mseBayesShr","mseBu",names(esNormal), names(esTrunc))
    dset <- paste0("largeSynthetic_n",synth_n)
  }
  
  else
  {
    dataFrame <- data.frame(h, fmethod, dset, mseBase, mseBayesShr,
                            as.vector(esNormal), as.vector(esTrunc))
  }
  
  
  filename <- paste("results/mse_",dset,".csv",sep="")
  
  writeNames <- TRUE
  if(file.exists(filename)){
    writeNames <- FALSE
  }
  write.table(dataFrame, file=filename, append = TRUE, sep=",", row.names = FALSE, col.names = writeNames)
}

