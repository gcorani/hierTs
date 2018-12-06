checkVar <- function(){
library(fpp2)
library(hts)  
  
  fets <- function(x, h) {
    forecast(ets(x), h = h)
  }
  
  for (i in 1:dim(htseg2[[1]])[2]){
    ts <- htseg2[[1]][,i]
    e <- tsCV(ts, fets, h=1)  
    model <- ets(ts)
  }
  
  
}