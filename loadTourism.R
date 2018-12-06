loadTourism  <- function() {
  # library(tidyverse)
  tourism <- read.csv("TourismData_v3.csv", 
                      na = "empty")
  
  #the first two columns contain the time
  tourism <- tourism[,-1:-2]
  
  #we known that it start in Jan 1998
  y=ts(tourism, frequency = 12, start = c(1998,1))
  hierTourism <- hts(y=y, bnames = colnames(tourism),
                     characters = c(3,3))
  #at this point,   
  #Number of nodes at each level: 1 76 304 
  #Total number of series: 381 
  #Number of observations per series: 228  
  
  return (hierTourism)
}