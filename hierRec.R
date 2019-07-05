hierRec <- function (dset, h=1, fmethod="ets", iTest=1, 
                     seed=0, synth_n=100, synthCorrel=0.5)
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
  library(huge)
  source("loadTourism.R")
  set.seed(seed)
  
  if (is.character(dset) == FALSE) {
    stop ("dset should be a string")
  }
  
  #covariance can be "diagonal", "sam" or "glasso"
  bayesRecon <- function (covariance){
    S <- smatrix(train)
    bottomIdx <- seq( nrow(S) - ncol(S) +1, nrow(S))
    upperIdx <- setdiff(1:nrow(S),bottomIdx)
    
    #prior mean and covariance of the bottom time series
    priorMean <- preds[bottomIdx]
    Y_vec <- preds[upperIdx]
    
    
    #prior covariance for the bottom time series
    bottomVar <- sigma[bottomIdx]^2
    if (covariance=="diagonal"){
      priorCov <- diag(bottomVar)
    }
    else if (covariance=="sam"){
      #the covariances are the covariances of the time series
      #the variances are the variances of the forecasts, hence the variances of the residuals
      bottomResiduals <- residuals[,bottomIdx]
      priorCov <- cov(bottomResiduals)
    }
    else if (covariance=="glasso"){
      #the covariances are the covariances of the time series
      #the variances are the variances of the forecasts, hence the variances of the residuals
      bottomResiduals <- residuals[,bottomIdx]
      out.glasso <- huge(bottomResiduals, method = "glasso", cov.output = TRUE)
      out.select <- huge.select(out.glasso, criterion = "ebic")
      priorCov <- out.select$opt.cov
    }
    
    upperVar <- sigma[upperIdx]^2
    #covariance for the upper time series; we need managing separately the case where only a single time series is present
    #as diag will try to create a matrix of size upperVar instead.
    
    if (length(upperIdx)==1) {
      Sigma_y <- upperVar
    }
    
    else if (covariance=="diagonal"){
      Sigma_y <- diag(upperVar)
    }
    #if we only one upper time series, there is no covariance matrix to be estimated. 
    else if (covariance=="glasso") {
      #get variance and covariance of the residuals
      upperResiduals <- residuals[,upperIdx]
      out.glasso <- huge(upperResiduals, method = "glasso", cov.output = TRUE)
      out.select <- huge.select(out.glasso, criterion = "ebic")
      Sigma_y <- out.select$opt.cov
    }
    
    else if (covariance=="sam") {
      #get variance and covariance of the residuals
      upperResiduals <- residuals[,upperIdx]
      Sigma_y <- cov(upperResiduals)
    }
    
    #==updating
    #A explains how to combin the bottom series in order to obtain the
    # upper series
    
    #if upperIdx contains a single row, R behaves oddily; hence we need to manually manage that case.
    if (length(upperIdx)==1){
      A <- cbind(S[upperIdx,])
    }
    else {
      A <- t(S[upperIdx,])
    }
    
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
  
  
  if (dset=="tourism"){
    hierTs <- loadTourism()
  }
  else if (dset=="infantgts"){
    hierTs <- infantgts
  }
  else if (dset=="synthetic"){
    source("drawSmallHierarchy.R")
    #we generate the hierarchy with two bottom time series
    #training and the test
    listSynth <- artificialTs(n=synth_n + h, correl = synthCorrel)
    synthTs <- listSynth$bottomTs
    corrB2_U <- listSynth$corrB2_U
    y=ts(synthTs, frequency = 1)
    hierTs <- hts(y, bnames = colnames(y), characters=c(1,1))
  }
  else if (dset=="syntheticLarge"){
    source("drawLargeHierarchy.R")
    #we generate the hierarchy with *four* bottom time series
    synthTs <- simulFourBottom(n=synth_n)
    y=ts(synthTs, frequency = 1)
    hierTs <- hts(y, bnames = colnames(y))
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
  
  #sometimes the sample matrix is not positive definite and minT crashes
  #the matrix is computed internally by minT and cannot be controlled from here.
  mseCombMintSample <- NA
  try({
  fcastCombMintSam <- 
    forecast(train, h = h, method = "comb", weights="mint", fmethod=fmethod, 
             covariance="sam")
  mseCombMintSample  <- hierMse(fcastCombMintSam, test,  h)
  })
  fcastCombMintShr <- 
    forecast(train, h = h, method = "comb", weights="mint", fmethod=fmethod, 
             covariance="shr")
  
  
  mseCombMintShr  <- hierMse(fcastCombMintShr, test,  h)
  
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
  
  
  mseBayesDiag =  mean  ( (allts(test)[h,] - bayesRecon(covariance="diagonal"))^2 )
  mseBayesGlasso =  mean  ( (allts(test)[h,] - bayesRecon(covariance="glasso"))^2 )
  
  mseBayesSample <- NA
  try({
  mseBayesSample =  mean  ( (allts(test)[h,] - bayesRecon(covariance="sam"))^2 )
  })
  #save to file the results, at every iteration
  
  if (dset=="synthetic"){
    dataFrame <- data.frame(h, fmethod, synth_n, synthCorrel, corrB2_U, mseBase,mseCombMintSample,
                          mseCombMintShr, mseBayesDiag, mseBayesSample, mseBayesGlasso)
    colnames(dataFrame) <- c("h","fmethod","sampleSize","correlB1_U","correlB2_U",
                             "mseBase","mseMintSample","mseCombMintShr","mseBayesDiag","mseBayesSample",
                             "mseBayesGlasso")
    dset <- paste0(dset,"_correl",synthCorrel,"_n",synth_n)
  }
  
  else if (dset=="syntheticLarge"){
    dataFrame <- data.frame(h, fmethod, synth_n, mseBase,mseCombMintSample,
                            mseCombMintShr, mseBayesDiag, mseBayesSample, mseBayesGlasso)
    colnames(dataFrame) <- c("h","fmethod","sampleSize",
                             "mseBase","mseMintSample","mseCombMintShr","mseBayesDiag","mseBayesSample",
                             "mseBayesGlasso")
    dset <- paste0("largeSynthethic_n",synth_n)
  }
  
  else
  {
    dataFrame <- data.frame(h, fmethod, dset, calibration50, calibration80, 
                            mseBase,mseCombMintSample,mseCombMintShr,mseBayesDiag,mseBayesSample,mseBayesGlasso)
  }
  
  
  filename <- paste("results/",dset,".csv",sep="")
  writeNames <- TRUE
  if(file.exists(filename)){
    writeNames <- FALSE
  }
  write.table(na.omit(dataFrame), file=filename, append = TRUE, sep=",", row.names = FALSE, col.names = writeNames)
}

