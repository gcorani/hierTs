hier <- function (dset, h=1, fmethod="ets"){
  #given the hierTs data set ("tourism","infantgts"), reconciles the h-steps ahead forecast 
  #fmethod can be "ets" or "arima"
  
  library(hts)
  source("loadTourism.R")
  
  if (is.character(dset) == FALSE) {
    stop ("dset should be a string")
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
  
  hierMae <- function (htsPred, htsActual) {
    #receives two hts objects, containing  forecast and actual value.
    #computes the mae for the relevant forecast horizon only
    #Sets correctly the columns names currently not possible
    maeHts <- abs (allts(htsPred) - allts(htsActual))[h,]
    # colnames(maeHts) <- colnames(allts(htsPred))
    return (maeHts)
  }
  
  hierMse <- function (htsPred, htsActual) {
    #receives two hts objects, containing  forecast and actual value.
    #computes the mse for the whole hierarchy.
    mse <- mean  ( (allts(htsPred) - allts(htsActual)[h,])^2 )
    return (mse)
  }
  
  #extract the time from the data set to then split into train / test (test set contains 25 or 5 time points)
  set.seed(seed = 0)
  
  if (dset=="tourism"){
    hierTs <- loadTourism()
  }
  else if (dset=="htseg1"){
    hierTs <- htseg1
  }
  else if (dset=="htseg2"){
    hierTs <- htseg2
  }
  else if (dset=="infantgts"){
    hierTs <- infantgts
  }
  
  
  
  testSize <- 25
  if (length(hierTs$bts[,1]) < 25){
    testSize <- 5
  }
  
  #if h=1, the possible preds are the whole test size lenght; 
  #if h=2, the possible preds are the (test size lenght -1); etc.
  possiblePreds <- testSize - h + 1
  
  totTs       <- nrow(smatrix(hierTs))
  maeBase     <- matrix(nrow = possiblePreds, ncol = totTs)
  maeBu       <- matrix(nrow = possiblePreds, ncol = totTs)
  maeComb     <- matrix(nrow = possiblePreds, ncol = totTs)
  maeCombWls  <- matrix(nrow = possiblePreds, ncol = totTs)
  maeCombMint <- matrix(nrow = possiblePreds, ncol = totTs)
  maeBayes    <- matrix(nrow = possiblePreds, ncol = totTs)
  
  #These vectors will contain the global mse, summed over all the time series of the hierarchy
  mseBase     <- vector(length   = possiblePreds)
  mseBu       <- vector(length   = possiblePreds)
  mseComb     <- vector(length   = possiblePreds)
  mseCombWls  <- vector(length   = possiblePreds)
  mseCombMint <- vector(length   = possiblePreds)
  mseBayes    <- vector(length   = possiblePreds)
  
  elapsedBase     <- vector(length   = possiblePreds)
  elapsedBu     <- vector(length   = possiblePreds)
  elapsedComb     <- vector(length   = possiblePreds)
  elapsedCombMint  <- vector(length   = possiblePreds)
  elapsedBayes    <- vector(length   = possiblePreds)
  
  
  for (iTest in 1:possiblePreds) {
    timeIdx             <- time(hierTs$bts[,1])
    endTrain            <- length(timeIdx) - h - (iTest - 1)
    train               <- window(hierTs, start = timeIdx[1], end = timeIdx[endTrain] )
    test                <- window(hierTs, start =timeIdx[endTrain +1], end=timeIdx[endTrain + h])
    
    ptm <- proc.time()
    fcastBu             <- forecast(train, h = h, method = "bu", fmethod = fmethod)
    elapsedBu[iTest] <- (proc.time() - ptm)["elapsed"]
    
    ptm <- proc.time()
    fcastComb           <- forecast(train, h = h, method = "comb", weights="ols", fmethod=fmethod)
    elapsedComb[iTest] <- (proc.time() - ptm)["elapsed"]
    
    fcastCombWls        <- forecast(train, h = h, method = "comb", weights="wls", fmethod=fmethod)
    
    ptm <- proc.time()
    fcastCombMint       <- forecast(train, h = h, method = "comb", weights="mint", fmethod=fmethod)
    elapsedCombMint[iTest] <- (proc.time() - ptm)["elapsed"]
    
    mseBu[iTest]        <- hierMse(fcastBu, test )
    mseComb[iTest]      <- hierMse(fcastComb, test )
    mseCombWls[iTest]   <- hierMse(fcastCombWls, test )
    mseCombMint[iTest]  <- hierMse(fcastCombMint, test )
    
    #recompute predictions to be easily accessed by the Bayesian method
    allTsTrain <- allts(train)
    numTs <- ncol(allTsTrain)
    alpha <- 0.2
    sigma <- vector(length = numTs)
    preds <- vector(length = numTs)
    
    #compute, for each  ts, predictions and sigma (h-steps ahead) 
    ptm <- proc.time()
    for (i in 1:numTs){
      if (fmethod=="ets"){
        model <- ets(allTsTrain[,i])
        tmp <- forecast(model, h=h, level=1-alpha)
        # print(paste(as.character(i),"/",as.character(numTs)))
        # print(model$components)
      }
      else if (fmethod=="arima"){
        model <- auto.arima(allTsTrain[,i])
        tmp <- forecast(model, h=h, level=1-alpha)
      }
      preds[i] <- tmp$mean[h]
      sigma[i] <- abs ( (tmp$mean[h] - tmp$upper[h])  / (qnorm(alpha / 2)) )
    }
    elapsedBase[iTest] <- (proc.time() - ptm)["elapsed"]
    mseBase[iTest] =  mean  ( (allts(test)[h,] - preds)^2 )
    
    
    S <- smatrix(train)
    ptm <- proc.time()
    #Bayesian reconciliation
    bottomIdx <- seq( nrow(S) - ncol(S) +1, nrow(S))
    upperIdx <- setdiff(1:nrow(S),bottomIdx)
    #prior mean and covariance of the bottom time series
    priorMean <- preds[bottomIdx]
    #prior mean and covariance of the upper time series
    Y_vec <- preds[upperIdx]
    Sigma_y <- matrix(nrow = length(upperIdx), ncol = length(upperIdx))
    
    #prior covariance for the bottom time series
    bottomVar <- sigma[bottomIdx]^2
    priorCov <- diag(bottomVar)
    
    
    #covariance for the upper time series
    upperVar <- sigma[upperIdx]^2
    Sigma_y <- diag(upperVar)
    
    
    #==update in a single shot
    #A explains how to combin the bottom series in order to obtain the
    # upper series
    A <- t(S[upperIdx,])
    
    correl <- priorCov %*% A %*%
      solve (t(A) %*% priorCov %*% A + Sigma_y)
    
    postMean <- priorMean + correl  %*%
      (Y_vec - t(A) %*% priorMean)
    bayesPreds <- buReconcile(postMean, S, predsAllTs = FALSE)
    elapsedBayes[iTest] <- (proc.time() - ptm)["elapsed"]
    mseBayes[iTest] =  mean  ( (allts(test)[h,] - bayesPreds)^2 )
  }
  
  
  h <- rep(h, possiblePreds)
  fmethod <- rep(fmethod, possiblePreds)
  dset <- rep(dset, possiblePreds)
  dataFrame <- data.frame(h, fmethod, dset, mseBase,mseBu,mseComb,mseCombWls,mseCombMint,mseBayes)
  
  
  filename <- "results/mseHierReconc.csv"
  writeNames <- TRUE
  if(file.exists(filename)){
    writeNames <- FALSE
  }
  
  write.table(dataFrame, file=filename, append = TRUE, sep=",", row.names = FALSE, col.names = writeNames)
  return(dataFrame)
  
}



