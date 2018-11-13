hier <- function (hierTs, h=1, fmethod="ets"){
  #given the hierTs data set, reconciles the h-steps ahead forecast 
  #fmethod can be "ets", "arima" or "rw" (rw still to be implemented)
  
  library(hts)
  
  #computes the bu prediction given the predictions (1 x tot time series) and the S matrix
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
  
  
  #extract the time from the data set to then split into train / test (test set contains 25 or 5 time points)
  testSize <- 25
  if (length(hierTs$bts[,1]) < 25){
    testSize <- 5
  }
  #if h=1, the possible preds are the whole test size lenght; 
  #if h=2, the possible preds are the (test size lenght -1); etc.
  possiblePreds <- testSize - h + 1
  
  totTs <- nrow(smatrix(hierTs))
  maeBu <- matrix(nrow = possiblePreds, ncol = totTs)
  maeComb <- matrix(nrow = possiblePreds, ncol = totTs)
  maeCombWls <- matrix(nrow = possiblePreds, ncol = totTs)
  maeCombMint <- matrix(nrow = possiblePreds, ncol = totTs)
  maeBayes <- matrix(nrow = possiblePreds, ncol = totTs)
  
  for (iTest in 1:possiblePreds) {
    timeIdx <- time(hierTs$bts[,1])
    endTrain <- length(timeIdx) - h - (iTest - 1)
    train <- window(hierTs, start = timeIdx[1], end = timeIdx[endTrain] )
    test <- window(hierTs, start =timeIdx[endTrain +1], end=timeIdx[endTrain + h])
    fcastBu <- forecast(train, h = h, method = "bu", fmethod = fmethod)
    fcastComb <- forecast(train, h = h, method = "comb", weights="ols", fmethod=fmethod)
    fcastCombWls <- forecast(train, h = h, method = "comb", weights="wls", fmethod=fmethod)
    fcastCombMint <- forecast(train, h = h, method = "comb", weights="mint", fmethod=fmethod)
    maeBu[iTest,] <- hierMae(fcastBu, test )
    maeComb[iTest,] <- hierMae(fcastComb, test )
    maeCombWls[iTest,] <- hierMae(fcastCombWls, test )
    maeCombMint[iTest,] <- hierMae(fcastCombMint, test )
    
    #recompute predictions to be easily accessed by the Bayesian method
    allTsTrain <- allts(train)
    numTs <- ncol(allTsTrain)
    alpha <- 0.2
    sigma <- vector(length = numTs)
    preds <- vector(length = numTs)
    #compute, for each  ts, predictions and sigma (h-steps ahead) 
    for (i in 1:numTs){
      if (fmethod=="ets"){
        model <- ets(ts(allTsTrain[,i]))
        tmp <- forecast(model, h=h, level=1-alpha)
      }
      else if (fmethod=="arima"){
        model <- auto.arima(ts(allTsTrain[,i]))
        tmp <- forecast(model, h=h, level=1-alpha)
      }
      else if (fmethod=="rw"){
        tmp <- rwf(ts(allTsTrain[,i]),h=h, level=1-alpha)
       }
      preds[i] <- tmp$mean[h]
      sigma[i] <- abs ( (tmp$mean[h] - tmp$upper[h])  / (qnorm(alpha / 2)) )
    }
    
    S <- smatrix(train)
    
    #Bayesian reconciliation
    bottomIdx <- seq( nrow(S) - ncol(S) +1, nrow(S))
    upperIdx <- setdiff(1:nrow(S),bottomIdx)
    #p is the number of the bottom time series
    p <- length(bottomIdx)
    #prior mean and covariance of the bottom time series
    priorMean <- preds[bottomIdx]
    priorCov <- matrix (nrow = p, ncol = p)
    #prior mean and covariance of the upper time series
    Y_vec <- preds[upperIdx]
    Sigma_y <- matrix(nrow = length(upperIdx), ncol = length(upperIdx))
    
    
    for (i in 1:p) {
      priorCov[i,i] <- sigma[i]^2
    }
    priorCov[is.na(priorCov)] <- 0
    #prior covariance is now instantiated
    
    #next row is necessary to finalize the covariance matrix of y
    for (i in 1:length(upperIdx)) {
      Sigma_y[i,i] <- sigma[upperIdx[i]]^2
    }
    Sigma_y[is.na(Sigma_y)] <- 0  
    #Sigma_y  is now instantiated
    
    
    #==update in a single shot
    #A explains how to combin the bottom series in order to obtain the
    # upper series
    A <- t(S[upperIdx,])
    
    correl <- priorCov %*% A %*%
      solve (t(A) %*% priorCov %*% A + Sigma_y)
    
    postMean <- priorMean + correl  %*%
      (Y_vec - t(A) %*% priorMean)
    bayesPreds <- buReconcile(postMean, S, predsAllTs = FALSE)
    maeBayes[iTest,] = abs (allts(test)[h,] - bayesPreds) 
    
  }
  
  return( list (
    "signBetterBu"=mean(maeBayes<maeBu),
    "signBetterComb"=mean(maeBayes<maeComb),
    "signBetterCombWls"=mean(maeBayes<maeCombWls),
    "signBetterCombMint"=mean(maeBayes<maeCombMint),
    
    "maeImprovementBu"= mean ( (maeBu-maeBayes)/ ((maeBu+maeBayes)/2) ),
    "maeImprovementComb"= mean ( (maeComb-maeBayes)/ ((maeComb+maeBayes)/2) ),
    "maeImprovementWls"= mean ( (maeCombWls-maeBayes)/ ((maeCombWls+maeBayes)/2) ),
    "maeImprovementMint"= mean ( (maeCombMint-maeBayes)/ ((maeCombMint+maeBayes)/2) ),
    
    "maeBu"=maeBu, "maeComb"=maeComb, "maeCombWls"=maeCombWls, "maeCombMint"=maeCombMint, "maeBayes"=maeBayes
  ))
  
  
}

