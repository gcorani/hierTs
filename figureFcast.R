figureFcast <- function(seed=0){
  #create figure of forecasts for the paper
  source("hierRec.R")
  library(hts)
  library(SHIP)
  set.seed(seed)
  source("drawLargeHierarchy.R")
  fmethod <- "arima"
  synth_n <- 20
  #we generate the hierarchy with *four* bottom time series
  synthTs <- simulFourBottom(n=synth_n)
  y=ts(synthTs, frequency = 1)
  hierTs <- hts(y, bnames = colnames(y), characters=c(1,1))
  maxH = 4
  
  iTest <- 1
  timeIdx <- time(hierTs$bts[,1])
  startTrain          <- timeIdx[1]
  endTrain            <- length(timeIdx) - maxH - (iTest - 1)
  train               <- window(hierTs, start = startTrain, end = timeIdx[endTrain] )
  test                <- window(hierTs, start =timeIdx[endTrain +1], end=timeIdx[endTrain +   maxH])
  
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
  #we consider 4 steps ahead

  base <- matrix(ncol=numTs, nrow=maxH)
  reconcBayes <- matrix(ncol=numTs, nrow=maxH)
  for (h in (seq(1:maxH))){
    for (i in 1:numTs){
      model <- auto.arima(allTsTrain[,i])
      tmp <- forecast(model, h=h, level=1-alpha)
      residuals[,i] <- model$residuals
      preds[i] <- tmp$mean[h]
    }
    mseBase =  mean  ( (allts(test)[h,] - preds)^2 )
    reconcBayes[h,] <- bayesRecon(covariance="shr")
    mseBayesShr =  mean   (allts(test)[h,] - reconcBayes[h,])^2 
    base[h,] <- preds
    actual <- allts(test)
  }
  
  
  bottomIdx <- 1

  actual <- hierTs$bts[,bottomIdx]
  base   <- ts(base[,3+bottomIdx], start=timeIdx[endTrain +1])
  reconc <- ts(reconcBayes[,3+bottomIdx], start=timeIdx[endTrain +1])
  p <- autoplot(actual) + autolayer(base)   + autolayer(reconc) # + theme(legend.position = c(.95, .95))
  
  #create a 2 by 2 figure
    #the last 4 ts are the useful ones
  # for (idx in (seq(4:8))){
  #   autoplot(actual[,idx])
  # }
  pdf("rplot.pdf")
}