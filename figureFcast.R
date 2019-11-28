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
  
  baseAll <- matrix(ncol=numTs, nrow=maxH+1)
  baseAll[1,] <- tail(allts(train),1)
  fmintAll <- baseAll
  bayesAll <- baseAll
  
  
  for (h in (seq(1:maxH))){
    for (i in 1:numTs){
      model <- auto.arima(allTsTrain[,i])
      tmp <- forecast(model, h=h, level=1-alpha)
      residuals[,i] <- model$residuals
      preds[i] <- tmp$mean[h]
    }
    mseBase =  mean  ( (allts(test)[h,] - preds)^2 )
    bayesAll[2:(h+1),] <- bayesRecon(covariance="shr")
    mseBayesShr =  mean   (allts(test)[h,] - bayesAll[h,])^2 
    baseAll[2:(h+1),] <- preds
    actual <- allts(test)
    fmintAll[2:(h+1),] <-
      allts(forecast(train, h = h, method = "comb", weights="mint", fmethod=fmethod, 
                     covariance="shr"))[h,]
    mseMint <-  mean   (allts(test)[h,] - fmintAll[h,])^2 
  }
  
  for (bottomIdx in seq(1:4)){
    y <- hierTs$bts[,bottomIdx]
    base  <- ts (baseAll[,3+bottomIdx], start=timeIdx[endTrain])
    bayes <- ts (bayesAll[,3+bottomIdx], start=timeIdx[endTrain])
    fmint <- ts (fmintAll[,3+bottomIdx], start=timeIdx[endTrain])
    filename <- paste0("results/plots/plot_seed",seed,"_bottom",bottomIdx,".pdf")
    pdf(filename, width = 4, height = 3)
    p <- autoplot(y, size = 0.75)+ autolayer(base, size = 0.75) + autolayer(bayes, size = 0.75)   + autolayer(fmint, size = 0.5)
    print(p)
    dev.off()
  }
  #create a 2 by 2 figure
  #the last 4 ts are the useful ones
  # for (idx in (seq(4:8))){
  #   autoplot(actual[,idx])
  # }
  
}
