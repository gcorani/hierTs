example <- function (hierTs){
  library(hts)

  #computes the bu prediction given the predictions (1 x tot time series) and the S matrix
  #(tot time series X bottom time series)
  #predsAllTs is a flag: is set to true, the input preds contains predictions for all the hierarchy
  #and the function retrieves the bottom series; if set to false, this is not needed
  #as preds only contains only bottom time series
  buReconcile <- function (preds,S, predsAllTs = FALSE) {
    bottomPreds <- preds
    if (predsAllTs) {
      #retrive the bottom prediction from all predictions
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
  
  #for simplicity we run only 1-step ahead predictions
  h=1
  
  #extract the time from the data set to then split into train / test (test set contains 25 or 5 time points)
  testSize <- 25
  if (length(hierTs) < 25){
    testSize <- 5
  }
  timeIdx <- time(hierTs$bts[,1])
  train <- window(hierTs, start = timeIdx[1], end = timeIdx[length(timeIdx) - testSize] )
  test <- window(hierTs, start =timeIdx[length(timeIdx) - testSize +1], end=timeIdx[length(timeIdx) - testSize + 1])
  fcastBu <- forecast(train, h = h, method = "bu")
  fcastComb <- forecast(train, h = h, method = "comb", weights="ols")
  fcastCombWls <- forecast(train, h = h, method = "comb", weights="wls")
  fcastCombMint <- forecast(train, h = h, method = "comb", weights="mint")
  maeBu <- accuracy(fcastBu, test)["MAE",]
  maeComb <- accuracy(fcastComb, test)["MAE",]
  maeCombWls <- accuracy(fcastCombWls, test)["MAE",]
  maeCombMint <- accuracy(fcastCombMint, test)["MAE",]
  
  
  allTsTrain <- allts(train)
  numTs <- ncol(allTsTrain)
  alpha <- 0.2
  sigma <- vector(length = numTs)
  preds <- matrix(nrow = h, ncol = numTs)
  #now we compute predictions and their sigma for each trainign ts
  for (i in 1:numTs){
    model <- ets(ts(allTsTrain[,i]))
    tmp <- forecast(model, h=1, level=1-alpha)
    preds[h,i] <- tmp$mean[1]
    #we  need the [1] to access the numerical information within a ts objects
    sigma[i] <- abs ( (tmp$mean[1] - tmp$upper[1])  / (qnorm(alpha / 2)) )
  }
  
  S <- smatrix(train)
  buPreds <- buReconcile(preds, S, predsAllTs=TRUE)
  maeBu_gc = abs (allts(test) - buPreds)
  #unit test: is our bu consistent with hyndmand?
  consistent = matrix(maeBu_gc) == maeBu
  if (mean(consistent) < 1){
    error("bu implementation not consistent")
  }
  
  #let's implement now the bayesian version
  #p is the number of the bottom time series
  bottomIdx <- seq( nrow(S) - ncol(S) +1, nrow(S))
  upperIdx <- setdiff(1:nrow(S),bottomIdx)
  p <- length(bottomIdx)
  priorMean <- preds[bottomIdx]
  priorCov <- matrix (nrow = p, ncol = p)
  
  for (i in 1:p) {
    priorCov[i,i] <- sigma[i]^2
  }
  priorCov[is.na(priorCov)] <- 0
  #prior covariance is now instantiated
  
  #now we start the updating with the forecast of each upper time series
  for (i in 1:length(upperIdx)){
    #get the matrix A of the linear map
    #A states which are the time series to be summed
    #in order to obtain the i-th upper time series
    A <- S[i,]
    mu_y <- preds[i]
    var_y <- sigma[i]^2
    
    correl <- priorCov %*% A %*%
      solve (t(A) %*% priorCov %*% A + var_y)
    
    priorMean <- priorMean + 
      correl  %*%
      (mu_y - t(A) %*% priorMean)
    
    priorCov <- priorCov - correl  %*% t(A) %*% priorCov
  }
  
  bayesPreds <- buReconcile(priorMean, S, predsAllTs = FALSE)
  #now you need to reconstruct the bu predictions
  maeBayes = abs (allts(test) - bayesPreds) 
  
  return( list (
    "percBetterBu"=mean(maeBayes<maeBu),
    "percBetterComb"=mean(maeBayes<maeComb),
    "percBetterCombWls"=mean(maeBayes<maeCombWls),
    "percBetterCombMint"=mean(maeBayes<maeCombMint),
    "maeBu"=maeBu, "maeComb"=maeComb, "maeCombWls"=maeCombWls, "maeCombMint"=maeCombMint, "maeBayes"=as.numeric(maeBayes)
                ))
  
  
}

