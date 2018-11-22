thier <- function (ts="usdeaths", fmethod="ets"){
  #given a time series, creates all the possible higher-level aggregations 
  # and then computes the reconciled forecasts
  #fmethod can be "ets" or "arima"
  
  library(fpp2)
  library(hts)
  library(thief)
  source("hier.R")
  
  #computes mae for temporal hierarchies
  #both actual and forecast are temporal hierarchies
  tHierMae <- function (actual, forecast){
    hierMae <- vector(length = length(forecast))
    for (i in seq_along(forecast)) {
      hierMae[i] <- mean( abs (forecast[[i]]$mean - actual[[i]]) )
    }
    return (tHierMae)
  }
  
  #builds the A matrix, which indicates  which bottom time series sum up to each upper time series.
  buildMatrix <- function() {
    A <- matrix(data = 0, nrow = length(upperIdx), ncol = length(bottomIdx))
    maxFreq <- frequency(trainTs[[1]])
    counter <- 1
    for (ii in (2:length(trainTs))){
      currentFreq <- frequency(trainTs[[ii]])
      aggregatedTs <- maxFreq / currentFreq
      offset <- 1
      for (jj in (1:aggregatedTs)) {
        A[counter, offset : (offset + currentFreq - 1)] <- 1
        offset <- offset + currentFreq
      }
    }
  }
  
  
  #coverage of the PI is 0.8
  alpha <- 0.2
  
  #both lines temporary
  ts <- usdeaths
  tsName <- "pippo"
  
  #predict one full season ahead (12 months or 4 quarters etc)
  timeIdx <- time(ts)
  endTrain <- length(timeIdx) - frequency(ts)
  train <- window(ts, start = timeIdx[1], end = timeIdx[endTrain] )
  test <- window(ts, start =timeIdx[endTrain +1])
  
  trainTs <- tsaggregates(train)
  testTs <- tsaggregates(test)
  
  # Compute forecasts one full season ahead
  fc <- list()
  for(i in seq_along(trainTs)){
    if (fmethod == "ets") {
      fc[[i]] <- forecast(ets(trainTs[[i]]), h=frequency(trainTs[[i]]), level = (1-alpha))
    }
    else if (fmethod == "arima") {
      fc[[i]] <- forecast(auto.arima(aggts[[i]]), h=2*frequency(aggts[[i]]))
    }
  }
  
  # Reconcile forecasts using thief
  thiefReconc <- reconcilethief(fc, comb = "struc")
  buReconc <- reconcilethief(fc, comb = "bu")
  
  #Reconcile using the Bayesian approach
  #how many predictions we manage within the hierarchy
  numTs <- 0
  for (i in seq_along(trainTs) ){
    numTs <- numTs + frequency(trainTs[[i]])
  }
  
  #recover sigma and mean of each prediction
  sigma <- vector(length = numTs)
  preds <- vector(length = numTs)
  offset <- 1
  #fill the predictions
  for (i in seq_along(fc)) {
    currentLenght <- length(fc[[i]]$mean)
    currentPreds <- fc[[i]]$mean
    preds[offset: (offset + currentLenght - 1 )] <- currentPreds
    offset <- offset + currentLenght
  }
  #fill the sigma 
  offset <- 0
  for (i in seq_along(fc)) {
    currentLenght <- length(fc[[i]]$mean)
    #the sigma is different for each prediction, even if they are at the same level
    for (j in 1:currentLenght) {
      sigma[j + offset] <- abs ( (fc[[i]]$mean[j] - fc[[i]]$upper[j])  / (qnorm(alpha / 2)) )
    }
    offset <- offset + currentLenght
  }
  
  #setup the S matrix (code from combine.R of thief package)
  freq <- unlist(lapply(trainTs,frequency))
  m <- max(freq)
  nsum <- rev(rep(m/freq, freq))
  unsum <- unique(nsum)
  grps <- matrix(0, nrow=length(unsum)-1, ncol=m)
  for(i in 1:(length(unsum)-1))
  {
    mi <- m/unsum[i]
    grps[i,] <- rep(1:mi, rep(unsum[i],mi))
  }
  ##== do not know if above code helpful
  
  #the time series in the first element of the list are the bottom ones.
  priorMean <- fc[[1]]$mean
  #prior covariance for the bottom time series
  bottomIdx <- 1:length(fc[[1]]$mean)
  bottomVar <- sigma[bottomIdx]^2
  priorCov <- diag(bottomVar)
  
  upperIdx <- setdiff(1:numTs,bottomIdx)
  #prior mean and covariance of the upper time series
  Y_vec <- preds[upperIdx]
  upperVar <- sigma[upperIdx]^2
  Sigma_y <- diag(upperVar)
  A <- buildMatrix()
  
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



