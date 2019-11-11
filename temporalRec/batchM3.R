batchM3 <- function(type="monthly",fmethod="ets"){
  library(Mcomp)
  source("thier.R")
  if (type=="monthly"){
    M3.selected <- subset(M3,12)
  }
  else if (type=="quarterly"){
    M3.selected <- subset(M3,4)
  }
  #we do not need any selection: 
  #summary of the lenght quarterly ts:
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 16.00   36.00   44.00   40.95   44.00   64.00 
  #summary of the lenght of monthly ts:
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  # 48.00   78.00  115.00   99.34  116.00  126.00   
  
  for (i in seq_along(M3.selected)) {
    print(M3.selected[[i]]$sn)
    thier(M3.selected[[i]],fmethod=fmethod, periodType=type)
  }
  
}
