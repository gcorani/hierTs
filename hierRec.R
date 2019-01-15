hierRec <- function (dset, h=1, fmethod="ets", iTest){
  #The hierTs data set ("tourism","infantgts"), reconciles the h-steps ahead forecast 
  #fmethod can be "ets" or "arima"
  #iTest allows to perform many experiments with rolling origin (iTest is comprised between 1 and 50) 
  
  library(hts)
  library(huge)
  source("loadTourism.R")
  
  if (is.character(dset) == FALSE) {
    stop ("dset should be a string")
  }
  
  bayesRecon <- function (correlation){
    S <- smatrix(train)
    bottomIdx <- seq( nrow(S) - ncol(S) +1, nrow(S))
    upperIdx <- setdiff(1:nrow(S),bottomIdx)
    
    #prior mean and covariance of the bottom time series
    priorMean <- preds[bottomIdx]
    Y_vec <- preds[upperIdx]
    Sigma_y <- matrix(nrow = length(upperIdx), ncol = length(upperIdx))
    
    #prior covariance for the bottom time series
    bottomVar <- sigma[bottomIdx]^2
    priorCov <- diag(bottomVar)
    if (correlation){
      #the covariances are the covariances of the time series
      #the variances are the variances of the forecasts, hence the variances of the residuals
      bottomResiduals <- residuals[,bottomIdx]
      priorCov <- cov(bottomResiduals)
      out.glasso <- huge(bottomResiduals, method = "glasso", cov.output = TRUE)
      out.select <- huge.select(out.glasso, criterion = "stars")
      priorCov <- out.select$opt.cov
    }
    
    
    
    #covariance for the upper time series
    upperVar <- sigma[upperIdx]^2
    Sigma_y <- diag(upperVar)
    if (correlation){
      #get variance and covariance of the residuals
      upperResiduals <- residuals[,upperIdx]
      Sigma_y <- cov(upperResiduals)
      out.glasso <- huge(upperResiduals, method = "glasso", cov.output = TRUE)
      out.select <- huge.select(out.glasso, criterion = "stars")
      Sigma_y <- out.select$opt.cov
    }
    
    
    #==updating
    #A explains how to combin the bottom series in order to obtain the
    # upper series
    A <- t(S[upperIdx,])
    
    M <- ncol ( t(A) %*% priorCov %*% A + Sigma_y )
    correl <- priorCov %*% A %*%
      solve (t(A) %*% priorCov %*% A + Sigma_y + 1e-6*diag(M))
    
    postMean <- priorMean + correl  %*%
      (Y_vec - t(A) %*% priorMean)
    bayesPreds <- buReconcile(postMean, S, predsAllTs = FALSE)
    return(bayesPreds)
  }
  
  
  #The buReconcile function computes the bu prediction given the predictions (1 x tot time series) and the S matrix
  #(tot time series X bottom time series)
  #predsAllTs is a flag: is set to true, the input preds contains predictions for all the hierarchy
  #and the function retrieves the bottom series; if set to false, this is not needed
  #as preds only contains only bottom time series
  buReconcile <- function (preds,S, predsAllTs = FALSE) {
    bottomPreds <- preds
    if (predsAllTs) {
      #retrieves the bottom prediction from all predictions
      upperIdx <- 1 : (nrow(S) - ncol(S))
      bottomIdx <- setdiff (1:nrow(S), upperIdx)
      bottomPreds <- preds [,bottomIdx]
    }
    
    buPreds <- preds
    
    #nrow(S) is the total number of time series
    for (i in 1:nrow(S)){
      buPreds[i] <- S[i,] %*% bottomPreds
    }
    return (buPreds)
  }
  
  #check the calibration of the prediction interval with coverage (1-currentAlpha)
  checkCalibration <- function(preds,sigmas,htsActual,coverage){
    stdQuant <- abs(qnorm((1-coverage)/2))
    included <- vector(length = length(preds))
    actual <- allts(htsActual)[h,]
    for (ii in seq_along(preds)){
      upper <- preds[ii] + stdQuant * sigmas[ii]
      lower <- preds[ii] - stdQuant * sigmas[ii]
      included[ii] <- (actual[ii] > lower) &  (actual[ii] < upper)
    }
    return (mean(included))
  }
  
  
  
  hierMse <- function (htsPred, htsActual, h) {
    #receives two hts objects, containing  forecast and actual value.
    #computes the mse for the whole hierarchy.
    mse <- mean  ( (allts(htsPred)[h,] - allts(htsActual)[h,])^2 )
    return (mse)
  }
  
  #extract the time from the data set to then split into train / test (test set contains 25 or 5 time points)
  set.seed(seed = 0)
  
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
    error("iTest>possiblePreds")
  }
  
  #here the experiment starts
  timeIdx             <- time(hierTs$bts[,1])
  endTrain            <- length(timeIdx) - h - (iTest - 1)
  train               <- window(hierTs, start = timeIdx[1], end = timeIdx[endTrain] )
  test                <- window(hierTs, start =timeIdx[endTrain +1], end=timeIdx[endTrain + h])
  
  fcastCombMint       <- forecast(train, h = h, method = "comb", weights="mint", fmethod=fmethod)
  mseCombMint  <- hierMse(fcastCombMint, test,  h)
  
  #recompute predictions to be  accessed by the Bayesian method
  allTsTrain <- allts(train)
  numTs <- ncol(allTsTrain)
  alpha <- 0.2
  sigma <- vector(length = numTs)
  preds <- vector(length = numTs)
  #the residuals for the model fitted on each time series
  residuals <- matrix(nrow=dim(allTsTrain)[1], ncol = numTs)
  
  #compute, for each  ts, predictions and sigma (h-steps ahead) 
  for (i in 1:numTs){
    if (fmethod=="ets"){
      model <- ets(allTsTrain[,i], additive.only = TRUE)
    }
    else if (fmethod=="arima"){
      model <- auto.arima(allTsTrain[,i])
    }
    tmp <- forecast(model, h=h, level=1-alpha)
    residuals[,i] <- model$x - model$fitted #we cannot store model$residuals due to models with multiplicative errors
    preds[i] <- tmp$mean[h]
    sigma[i] <- abs ( (tmp$mean[h] - tmp$upper[h])  / (qnorm(alpha / 2)) )
  }
  mseBase =  mean  ( (allts(test)[h,] - preds)^2 )
  
  
  calibration50 <- checkCalibration(preds, sigma, test, coverage = 0.5)
  calibration80 <- checkCalibration(preds, sigma, test, coverage = 0.8)
  
  mseBayes =  mean  ( (allts(test)[h,] - bayesRecon(correlation=FALSE))^2 )
  mseBayesCorr =  mean  ( (allts(test)[h,] - bayesRecon(correlation=TRUE))^2 )
  
  #save to file the results, at every iteration
  dataFrame <- data.frame(h, fmethod, dset, calibration50, calibration80, mseBase,mseCombMint,mseBayes,mseBayesCorr)
  filename <- paste("results/mseHierReconc",dset,".csv",sep="")
  writeNames <- TRUE
  if(file.exists(filename)){
    writeNames <- FALSE
  }
  write.table(dataFrame, file=filename, append = TRUE, sep=",", row.names = FALSE, col.names = writeNames)
}
