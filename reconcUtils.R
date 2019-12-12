library(modeest)  
#computes the energy score for the bottom / upper time 
#series, using the mvn or the truncated distribution
#distribution can be "normal" or "trunc"
#y is the actual vector of hierachical observations
#bayesPost contains posterior mean and var of the bottom ts, S and A.
#it returns a list which contains only the energy scores (normal case)
#and also the mse for mean / median / mode (truncated case)
energyScore <- function(y, bayesPost, distribution){
  n <- 10000
  #mean and variance for the bottom time series
  b_tilde <- bayesPost$b_tilde
  Sigma_b <- bayesPost$Sigma_b
  ySamples <- matrix(ncol = n, nrow = length(y))
  #samples for the bottom time series, each column is a drawn from the mvn
  if (distribution=="normal"){
    bottomSamples <- t(mvrnorm(n=n, mu = b_tilde, Sigma = Sigma_b))
  }
  else if (distribution=="trunc"){
    bottomSamples <- t(rtmvnorm(n=n, mean = as.vector(b_tilde), sigma = Sigma_b,
                               lower = rep(0, length = length(b_tilde))))
  }
  bottomIdx <- rowSums(bayesPost$S) == 1
  upperIdx  <- rowSums(bayesPost$S) != 1
  ySamples[bottomIdx,] <- bottomSamples
  
  #use tapply here, when you refactor
  for (i in 1:ncol(ySamples)){
    ySamples[upperIdx,i] <-
        bayesPost$A %*% ySamples[bottomIdx,i] 
  }
  
  
  esAll <-  es_sample(y, ySamples)
  #the rows of the bottom time series sum up to 1
  
  esUpper   <- es_sample(y[upperIdx], ySamples[upperIdx,])
  esBottom  <- es_sample(y[bottomIdx], ySamples[bottomIdx,])
  
  myList <- list (
    "esUpper"=esUpper,
    "esAll"=esAll,
    "esBottom"=esBottom
  )
  
  #if we are using the truncated, we compute the mse for the mean, median and mode
  if (distribution == "trunc"){
  #in order to guarantee coherence, we get the median of the 
  #bottom and then we multiply it by A
  median_y <- y
  median_y[bottomIdx] <- apply(bottomSamples, 1, median)
  median_y[upperIdx] <- bayesPost$A %*% median_y[bottomIdx]
  # print(paste("error on median", as.character(median_y[1]-sum(b_tilde))))
  
  mean_y <- y
  mean_y[bottomIdx] <- apply(bottomSamples, 1, mean)
  mean_y[upperIdx] <- bayesPost$A %*% mean_y[bottomIdx]
  # print(paste("error on mean", as.character(mean_y[1]-sum(b_tilde))))
  
  #mode computation via density, commented out
  # mode_y <- y
  # for (i in 1:sum(bottomIdx)){
  #   bDensity = density(bottomSamples[i,], n=n)
  #   mode_y[which(bottomIdx==TRUE)[i]] = bDensity$x[which.max(bDensity$y)]
  # }
  # mode_y[upperIdx] <- bayesPost$A %*% mode_y[bottomIdx]
  # print(paste("error on mode, density estimate", as.character(mode_y[1]-sum(b_tilde))))
  
  mode_y <- y
  for (i in 1:sum(bottomIdx)){
    mode_y[which(bottomIdx==TRUE)[i]] = mlv(bottomSamples[i,], method = "hsm")
  }
  mode_y[upperIdx] <- bayesPost$A %*% mode_y[bottomIdx]
  # print(paste("error on mode, mfv estimate", as.character(mode_y[1]-sum(b_tilde))))
  
    myList$mseMean =  mean  (( y - mean_y)^2)
    myList$mseMedian = mean (( y - median_y)^2)
    myList$mseMode = mean   (( y - mode_y)^2)
  }
  return (myList)
}


#Reconciliation using Bayes' rule
#covType can be "diagonal", "sam" (sample estimate), "shr" (shrinkage) 
#in real applications, the shr is recommended
#Function is modified in order to return, besides \tilde{y},  also \tilde{b} and the covariance of B. 
bayesRecon <- function (train, preds, covType){
  S <- smatrix(train)
  bottomIdx <- which(rowSums(S)==1)
  upperIdx  <- which(rowSums(S)!=1)
  
  #prior mean and covariance of the bottom time series
  b_hat <- preds[bottomIdx]
  Y_vec <- preds[upperIdx]
  
  
  #prior covariance for the bottom time series
  bottomVar <- sigma[bottomIdx]^2
  bottomResiduals <- residuals[,bottomIdx]
  if (covType=="diagonal"){
    priorCov <- diag(bottomVar)
  }
  else if (covType=="sam"){
    #the covariances are the covariances of the time series
    #the variances are the variances of the forecasts, hence the variances of the residuals
    priorCov <- cov(bottomResiduals)
  }
  else if (covType=="shr"){
    sigmaDiag <- diag(bottomVar)
    priorCov <-  shrink.estim(bottomResiduals, tar=build.target(bottomResiduals,type="D"))[[1]]
  }
  
  #variance of the upper forecast
  upperVar <- sigma[upperIdx]^2
  #covariance for the upper time series; we need managing separately the case where only a single time series is present
  #as diag will try to create a matrix of size upperVar instead.
  upperResiduals <- residuals[,upperIdx]
  if (length(upperIdx)==1) {
    Sigma_u <- upperVar
  }
  
  else if (covType=="diagonal"){
    Sigma_u <- diag(upperVar)
  }
  #if we only one upper time series, there is no covariance matrix to be estimated. 
  else if (covType=="sam") {
    #get variance and covariance of the residuals
    Sigma_u <- cov(upperResiduals)
  }
  else if (covType=="shr") {
    sigma_y_diag <- diag(upperVar)
    Sigma_u <-  shrink.estim(upperResiduals, tar=build.target(upperResiduals,type="D"))[[1]]
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
  
  M <- ncol ( t(A) %*% priorCov %*% A + Sigma_u )
  #cov_u is the covariance matrix of the upper time series, A Sigma_b A^T + Sigma_u + nugget
  #the nugget is a small diagonal noise which makes it more robust the inversion
  cov_U = t(A) %*% priorCov %*% A + Sigma_u + 1e-6*diag(M)
  G <- priorCov %*% A %*% solve (cov_U)
  
  b_tilde <- b_hat + G  %*% (Y_vec - t(A) %*% b_hat)
  y_tilde <- buReconcile(b_tilde, S, predsAllTs = FALSE)
  sigma_b <- priorCov - G %*% t(A) %*% priorCov
  
  return( list (
    "b_tilde" = b_tilde,
    "y_tilde" = y_tilde,
    "Sigma_b" = sigma_b,
    "A"=t(A),
    "S"=S
  )
  )
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


setDset <- function (dset, synth_n = 100, synth_correl=0.5){
  feasibleDset <- c("infantgts", "tourism", "synthetic", "syntheticLarge")
  if (! (dset %in% feasibleDset)){
    print("feasible dset are:")
    print(feasibleDset)
    stop ("wrong dset supplied" )
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
  return (hierTs)
}