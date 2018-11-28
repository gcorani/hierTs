batchAE <- function(fmethod="arima"){
  #prepares the AEdemand dataset and launches the temporal reconciliation experiment
  #the fmethod is forced to be arima as the expsmoothing does not support such long seasonality
  library(thief)
  library(fpp2)
  source("thier.R")
  
  #finestra di 40 esperimenti su 13-weeks ahead forecast
  
  for (i in (1:dim(AEdemand)[2])) {
    print(paste("AEdemand: ",colnames(AEdemand)[i]))
    freq <- frequency(AEdemand[,i])
    n <- length(AEdemand[,i])
    endTrain <- time(AEdemand[,i])[n - freq]
    startTest <- time(AEdemand[,i])[n - freq +1]
    train <- window(AEdemand[,i], end=endTrain) 
    test <- window(AEdemand[,i], start=startTest) 
    tsObj$x <- train
    tsObj$xx <- test
    tsObj$sn <- colnames(AEdemand)[i]
    thier (tsObj, fmethod=fmethod, periodType="weekly")
  }
  
}
