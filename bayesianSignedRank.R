bayesianSignedRank <- function(zVector,rope_min,rope_max) {
  #Bayesian signed-rank test as in sec. 4.2.1 of Journal of Machine Learning Research 18 (2017) 1-36
  # "Time for a Change: a Tutorial for Comparing Multiple Classifiers...", Benavoli, Corani et al.
  #Setting: z0=0, s=0.5.
  #It returns the posterior probability  of thetaLeft, thetaRope and thetaRight being the largest.
  #zVector is the vector whose median has to be checked, rope_min and rope_max the upper and lower bound of the rope.
  library(MCMCpack)
  
  samples <- 4000
  #we adopt s=0.5, z=0
  #build the vector 0.5 1 1 ....... 1
  weights <- c(0,rep(1,length(zVector)))
  
  #add the fake first observation in 0
  zVector <- c (0, zVector)
  
  sampledWeights <- rdirichlet(samples,weights)
  
  winLeft <- vector(length = samples)
  winRope <- vector(length = samples)
  winRight <- vector(length = samples)
  
  for (rep in 1:samples){
    currentWeights <- sampledWeights[rep,]
    for (i in 1:length(currentWeights)){
      for (j in i:length(currentWeights)){
        product= currentWeights[i] * currentWeights[j]
        if ((zVector[i]+zVector[j]) > (2*rope_max) ) {
          winRight[rep] <- winRight[rep] + product
        }
        else if ((zVector[i]+zVector[j]) < (2*rope_min) ) {
  
                 winLeft[rep] <- winLeft[rep] + product
        }
        else {
          winRope[rep] <- winRope[rep] + product
        }
      }
    }
    maxWins=max(winRight[rep],winRope[rep],winLeft[rep])
    winners = (winRight[rep]==maxWins)*1 + (winRope[rep]==maxWins)*1 + (winLeft[rep]==maxWins)*1
    winRight[rep] <- (winRight[rep]==maxWins)*1/winners
    winRope[rep] <- (winRope[rep]==maxWins)*1/winners
    winLeft[rep] <- (winLeft[rep]==maxWins)*1/winners
  }
  
  
  results = list ("probSmaller"=mean(winLeft), "probRope"=mean(winRope),
                  "probLarger"=mean(winRight) )
  
  return (results)
  
}