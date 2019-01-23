temporalRec <- function (tsObj, fmethod="ets", periodType="monthly"){
  #given a time series, creates all the possible higher-level aggregations 
  # and then computes the reconciled forecasts
  #fmethod can be "ets" or "arima"
  #tsObj is an object in the M3comp format
  #periodType can be "monthly" or "quarterly"
  
  library(fpp2)
  library(hts)
  library(thief)
  
  # computes mae for temporal hierarchies
  #both actual and forecast are temporal hierarchies
  #it  averages the  different time series having the same frequency
  getHierMae <- function (actual, forecast) {
    hierMae <- vector(length = length(forecast))
    for (i in seq_along(forecast)) {
      hierMae[i] <- mean( abs (forecast[[i]]$mean - actual[[i]]) )
    }
    return (hierMae)
  }
  
  globalMse <- function (actual, forecast) {
    mse <- 0
    for (i in seq_along(forecast)) {
      mse <- mse + sum ( (forecast[[i]]$mean - actual[[i]])^2 )
    }
    return (mse)
  }
  
  calibration <- function (actual, forecast, sigmas) {
    included50 <- as.vector(rep(0,length(actual)))
    included80 <- as.vector(rep(0,length(actual)))
    #implementation to be finalized
    for (ii in seq_along(forecast)) {
      upper50 <- as.numeric(fc[[i]]$upper[,"50%"])
      lower50 <- as.numeric(fc[[i]]$lower[,"50%"])
      included50[ii] <- (actual[ii] > lower) &  (actual[ii] < upper)
      upper80 <- as.numeric(fc[[i]]$upper[,"80%"])
      lower80 <- as.numeric(fc[[i]]$lower[,"80%"])
      included80[ii] <- (actual[ii] > lower) &  (actual[ii] < upper)
    }
    return (mean(included50),mean(included80))
  }
  
  
  
  
  #builds the A matrix, which indicates  which bottom time series sum up to each upper time series.
  buildMatrix <- function() {
    A <- matrix(data = 0, nrow = length(upperIdx), ncol = length(bottomIdx))
    maxFreq <- frequency(trainHier[[1]])
    counter <- 1
    for (ii in (2:length(trainHier))){
      currentFreq <- frequency(trainHier[[ii]])
      aggregatedTs <- currentFreq
      howManyBottomToSum <- maxFreq / currentFreq
      offset <- 1
      for (jj in (1:aggregatedTs)) {
        A[counter, offset : (offset + howManyBottomToSum - 1)] <- 1
        offset <- offset + howManyBottomToSum
        counter <- counter + 1
      }
    }
    return (t(A))
  }
  
  
  #coverage of the PI is 0.8
  alpha <- 0.2
  
  #default test set for the M3 is 18 months for the monthly (tp be adapted) and 8 quarters for the quarterly.
  #for the moment force the test to be as long as exactly one period 
  #we should implement a mechanism for doing repeated predictions
  train <- tsObj$x
  trainHier <- tsaggregates(train)
  timeIdx <- time(tsObj$xx)
  test <- window(tsObj$xx, end=timeIdx[frequency(tsObj$xx)])
  testHier <- tsaggregates(test)
  
  # Compute forecasts one full season ahead
  fc <- list()
  for(i in seq_along(trainHier)){
    #how many test observations are available
    h=length(testHier[[i]])
    if (fmethod == "ets") {
      fc[[i]] <- forecast(ets(trainHier[[i]]), h=h, level = c(0.5,0.8), additive.only = TRUE)
    }
    else if (fmethod == "arima") {
      fc[[i]] <- forecast(auto.arima(trainHier[[i]]), h=h, , level = (1-alpha))
    }
  }
  
  # Reconcile forecasts using thief
  thiefReconc <- reconcilethief(fc, comb = "struc")
  buReconc <- reconcilethief(fc, comb = "bu")
  
  #Reconcile using the Bayesian approach
  #how many predictions we manage within the hierarchy
  numTs <- 0
  for (i in seq_along(trainHier) ){
    numTs <- numTs + frequency(trainHier[[i]])
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
  
  
  #this code copied from hier.R
  correl <- priorCov %*% A %*%
    solve (t(A) %*% priorCov %*% A + Sigma_y)
  
  postMean <- priorMean + correl  %*%
    (Y_vec - t(A) %*% priorMean)
  
  #pay attention: we only overwrite the point forecast for the bottom time series,
  #without managing the covariance
  bayesFc <- fc
  bayesFc[[1]]$mean <- postMean
  bayesReconc <- reconcilethief(bayesFc, comb = "bu")
  
  #prepare the data to be saved
  #Header of the results
  #extract a string describing each frequency
  freqNames <- names(lapply(testHier, frequency))
  methodsNames <- c("Bu","Thief","Bayes")
  #title: mae for each method, for each frequency
  columnNames <- vector(length = length(methodsNames) * length(freqNames))
  
  counter <- 1
  for (mName in methodsNames){
    for (fName in freqNames){
      columnNames[counter] <- paste(mName,fName,sep = "_")
      counter <- counter + 1
    }
  }
  
  #vector containing the mae of each method in each level
  tmp <- as.vector(cbind(getHierMae(testHier, buReconc), getHierMae(testHier, thiefReconc), getHierMae(testHier, bayesReconc)))
  dataFrame <- as.data.frame(t(tmp))
  colnames(dataFrame) <- columnNames
  dataFrame$mseBase <- globalMse(testHier, fc)
  dataFrame$mseBu <- globalMse(testHier, buReconc)
  dataFrame$mseThief <- globalMse(testHier, thiefReconc)
  dataFrame$mseBayes <- globalMse(testHier, bayesReconc)
  #calibration still to be implemented: the easiest seems to ask the forecast function to already produce
  #the relevant intervals
  # dataFrame$calibration50 <- calibration(testHier, preds, sigma, 0.5) 
  # dataFrame$calibration80 <- calibration(testHier, preds, sigma, 0.8) 
  
  
  dataFrame$tsName <- tsObj$sn
  idx <- c(ncol(dataFrame), 1:(ncol(dataFrame)-1))
  dataFrame <- dataFrame[,idx]
  
  filename <- paste("results/temporalHier","_",periodType,"_",fmethod,".csv",sep = "")
  writeNames <- TRUE
  if(file.exists(filename)){
    writeNames <- FALSE
  }
  if(!dir.exists("results/")){
    dir.create("results")
  }
  
  write.table(dataFrame, file=filename, append = TRUE, sep=",", row.names = FALSE, col.names = writeNames)
}



