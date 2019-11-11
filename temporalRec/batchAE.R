batchAE <- function(fmethod="arima"){
  #prepares the AEdemand dataset and launches the temporal reconciliation experiment
  #the fmethod is forced to be arima as the expsmoothing does not support such long seasonality
  library(thief)
  library(fpp2)
  source("temporalRec.R")
  
  #finestra di 40 esperimenti su 13-weeks ahead forecast
  totExp=20
  
  for (i in (1:dim(AEdemand)[2])) {
    for (currentExp in (1:totExp)){
      print(paste("AEdemand: ",colnames(AEdemand)[i]))
      freq <- frequency(AEdemand[,i])
      n <- length(AEdemand[,i])
      endTrain <- time(AEdemand[,i])[n - freq - currentExp ]
      startTest <- time(AEdemand[,i])[n - freq - currentExp + 1]
      train <- window(AEdemand[,i], end=endTrain) 
      test <- window(AEdemand[,i], start=startTest) 
      tsObj<- list()
      tsObj$x <- train
      tsObj$xx <- test
      tsObj$sn <- colnames(AEdemand)[i]
      temporalRec (tsObj, fmethod=fmethod, periodType="weekly")
    }
  }
}
