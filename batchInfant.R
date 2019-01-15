batchInfant <- function(){
  source("hierRec.R")
  for (h in c(2,3,4)){
    for (model in c("ets","arima")){
      for (iTest in (1:(50-h))) {
        print(paste(as.character(h)," ",model,"",as.character(iTest)))
        hierRec("infantgts",h=h, fmethod = model, iTest = iTest)
      }
    }
  }
}
