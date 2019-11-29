hierRec <- function (dset, h=1, fmethod="ets", iTest=1, 
                     seed=0, synth_n=100, synthCorrel=0.5)
{
  #The hierTs data set can be ("tourism","infantgts", "synthetic" (2 bottom time series),"syntheticLarge" (4 bottom time series)) 
  #fmethod can be "ets" or "arima"
  #iTest decised the split between  training/test  (iTest is comprised between 1 and 50 and controls the separation between training and test).
  #The training set contains the data from the first observation up to the observation in position (length(timeSeries) - h - (iTest - 1)).
  #synth_n and synthCorrel are used only when generating synthetic data (synth_n: number of time points, synthCorrel: correlation between the two bottom time series.)
  #seed is especially important when you run synthetic experiments, to make sure you get different data in each experimetn
  
  library(hts)
  source("loadTourism.R")
  set.seed(seed)
  
  feasibleDset <- c("infantgts", "tourism", "synthetic", "syntheticLarge")
  if (! (dset %in% feasibleDset)){
    print("feasible dset are:")
    print(feasibleDset)
    stop ("wrong dset supplied" )
  }
  
  
  
  #Reconciliation using Bayes' rule
  #covariance can be "diagonal", "sam" (sample estimate), "shr" (shrinkage) 
  #in real applications, the shr is recommended
  bayesRecon <- function (covariance){
    S <- smatrix(train)
    bottomIdx <- seq( nrow(S) - ncol(S) +1, nrow(S))
    upperIdx <- setdiff(1:nrow(S),bottomIdx)
    
    #prior mean and covariance of the bottom time series
    priorMean <- preds[bottomIdx]
    Y_vec <- preds[upperIdx]
    
    
    #prior covariance for the bottom time series
    bottomVar <- sigma[bottomIdx]^2
    bottomResiduals <- residuals[,bottomIdx]
    if (covariance=="diagonal"){
      priorCov <- diag(bottomVar)
    }
    else if (covariance=="sam"){
      #the covariances are the covariances of the time series
      #the variances are the variances of the forecasts, hence the variances of the residuals
      priorCov <- cov(bottomResiduals)
    }
    else if (covariance=="shr"){
      sigmaDiag <- diag(bottomVar)
      priorCov <-  shrink.estim(bottomResiduals, tar=build.target(bottomResiduals,type="D"))[[1]]
    }
    
    upperVar <- sigma[upperIdx]^2
    #covariance for the upper time series; we need managing separately the case where only a single time series is present
    #as diag will try to create a matrix of size upperVar instead.
    upperResiduals <- residuals[,upperIdx]
    if (length(upperIdx)==1) {
      Sigma_y <- upperVar
    }
    
    else if (covariance=="diagonal"){
      Sigma_y <- diag(upperVar)
    }
    #if we only one upper time series, there is no covariance matrix to be estimated. 
    else if (covariance=="sam") {
      #get variance and covariance of the residuals
      Sigma_y <- cov(upperResiduals)
    }
    else if (covariance=="shr") {
      sigma_y_diag <- diag(upperVar)
      Sigma_y <-  shrink.estim(upperResiduals, tar=build.target(upperResiduals,type="D"))[[1]]
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
  
  
  #The buReconcile function computes the bu prediction given the predictions for the bottom time series
  #and the S matrix. We use it after having update the distribution of the bottom time series.
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
    hierTs <- hts(y, bnames = colnames(y))
  }
  else if (dset=="syntheticLarge"){
    source("drawLargeHierarchy.R")
    #we generate the hierarchy with *four* bottom time series
    synthTs <- simulFourBottom(n=synth_n)
    y=ts(synthTs, frequency = 1)
    hierTs <- hts(y, bnames = colnames(y), characters=c(1,1))
  }
  
  testSize <- 50
  
  
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
  
  #sometimes the sample matrix is not positive definite and minT crashes
  #the matrix is computed internally by.libPaths() minT and cannot be controlled from here.
  # mseMintSample <- NA
  mseMintSample <- NA
  try({
    fcastMintSam <-
      forecast(train, h = h, method = "comb", weights="mint", fmethod=fmethod,
               covariance="sam")
    mseMintSample  <- hierMse(fcastMintSam, test,  h)
  })
  fcastMintShr <-
    forecast(train, h = h, method = "comb", weights="mint", fmethod=fmethod, 
             covariance="shr")
  mseMintShr  <- hierMse(fcastMintShr, test,  h)
  
  fcastBu <- NA
  #we run bu only on synthetic data, otherwise becomes too slow
  if (dset=="syntheticLarge"){
    fcastBu <- forecast(train, h = h, method = "bu", fmethod=fmethod)
    mseBu   <- hierMse(fcastBu, test,  h)
  }
  
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
  mseBayesShr =  mean  ( (allts(test)[h,] - bayesRecon(covariance="shr"))^2 )
  
  
  mseBayesSample <- NA
  try({
    mseBayesSample =  mean  ( (allts(test)[h,] - bayesRecon(covariance="sam"))^2 )
  })
  #save to file the results, at every iteration
  
  if (dset=="synthetic"){
    dataFrame <- data.frame(h, fmethod, synth_n, synthCorrel, corrB2_U, mseBase,mseMintSample,
                            mseMintShr, mseBayesDiag, mseBayesSample, mseBayesShr)
    colnames(dataFrame) <- c("h","fmethod","sampleSize","correlB1_U","correlB2_U",
                             "mseBase","mseMintSample","mseMintShr","mseBayesDiag","mseBayesSample", "mseBayesShr")
    dset <- paste0(dset,"_correl",synthCorrel,"_n",synth_n)
  }
  
  else if (dset=="syntheticLarge"){
    dataFrame <- data.frame(h, fmethod, synth_n, mseBase,mseMintSample,
                            mseMintShr, mseBayesSample, mseBayesShr, mseBu)
    colnames(dataFrame) <- c("h","fmethod","sampleSize",
                             "mseBase","mseMintSample","mseMintShr","mseBayesSample",
                             "mseBayesShr","mseBu")
    dset <- paste0("largeSynthetic_n",synth_n)
  }
  
  else
  {
    dataFrame <- data.frame(h, fmethod, dset, calibration50, calibration80, 
                            mseBase,mseMintSample,mseMintShr,mseBayesDiag,mseBayesSample,mseBayesShr)
  }
  
  
  filename <- paste("results/mse_",dset,".csv",sep="")
  
  writeNames <- TRUE
  if(file.exists(filename)){
    writeNames <- FALSE
  }
  write.table(dataFrame, file=filename, append = TRUE, sep=",", row.names = FALSE, col.names = writeNames)
}

