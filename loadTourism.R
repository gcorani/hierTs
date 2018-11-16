loadTourism  <- function() {
tourism <- read_csv("TourismData_v3.csv", 
                           na = "empty")

#the first two columns contain the time
tourism <- tourism[,-1:-2]
dim(tourism)
hierTourism <- hts(y=tourism, bnames = colnames(tourism),
                   characters = c(3,3))
#at this point,   
#Number of nodes at each level: 1 76 304 
#Total number of series: 381 
#Number of observations per series: 228  

return (hierTourism)
}