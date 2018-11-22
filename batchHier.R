batchHier <- function(){
  source("hier.R")
  
  
  
  for (h in c(1,2,3,4)) {
    for (method in c("ets","arima")){
      a <- hier ("htseg2", h=h, fmethod = method)
      print(paste("htseg2 ", as.character(method)))
    }
  }
  
  for (h in c(1,2)) {
    for (method in c("ets","arima")){
      print(paste("htseg1 ", as.character(method)))
      a <- hier ("htseg1", h=h, fmethod = method)
    }
  }
  
  for (h in c(1,2,3,4)) {
    for (method in c("ets","arima")){
      print(paste("infantgts ", as.character(method)))
      a <- hier("infantgts", h=h, fmethod = method)
    }
  }
  
  for (h in c(1,2,3,4)) {
    for (method in c("ets","arima")){
      print(paste("tourism", as.character(method)))
      a <- hier("tourism", h=h, fmethod = method)
    }
  }
  
  
  
}