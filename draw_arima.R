

# draw parameters of an invertible and stationary arima model, see https://otexts.com/fpp2/AR.html
# draw_ar <- function (arOrder){
#   #drawing parameters   
#   phi = vector(length = 2)
#   if (arOrder ==1){
#     phi[1]=runif(n=1, min=-1, max=1)
#     phi[2]=0              
#   }
#   
#   if (arOrder == 2){
#     phi[2]=runif(n=1, min=-1, max=1)              
#     upper = 1 - phi[2]
#     lower = phi[2] - 1
#     phi[1]= runif(1, min=lower, max=upper)
#   }
#   return(phi)
# }
# 
# 
# draw_ma <- function (maOrder){
#   #drawing parameters   
#   theta = vector(length = 2)
#   if (maOrder == 1){
#     theta[1]=runif(n=1, min=-1, max=1)
#     theta[2]=0              
#   }
#   
#   if (maOrder == 2){
#     theta[2]=runif(n=1, min=-1, max=1)              
#     upper = 1 + theta[2]
#     lower = -theta[2] - 1
#     theta[1]= runif(1, min=lower, max=upper)
#   }
#   return (theta)
# }

#draw correlated noise for n time steps, return a vector (sampleSize x 2)
#assumes the variance to be 1, the bottomCorrel is set by the user
#in this very specific case, covariance = correlation.
#by default 2 bottom time series are assumed
drawNoise <- function(sampleSize, bottomCovar, bottomTs=2){
  library(MASS)
  #generate an excess of noise observations
  #non-generic code, but might be enough for some limited simulation
  if (bottomTs==2)
    noise <- mvrnorm(n=2*sampleSize, mu = c(0,0), Sigma= rbind(c(1,bottomCovar),c(bottomCovar,1)))
  if (bottomTs==4)
    noise <- mvrnorm(n=2*sampleSize, mu = c(0,0,0,0), Sigma= rbind(c(1,bottomCovar,bottomCovar, bottomCovar),
                                                                   c(bottomCovar,1, bottomCovar, bottomCovar),
                                                                   c(bottomCovar, bottomCovar, 1, bottomCovar),
                                                                   c(bottomCovar, bottomCovar, bottomCovar, 1)
    ))
  
  return (noise)
}


simulVAR1 <- function(n, phi, noise){
  howManyTs <- length(phi)
  simul <- matrix( nrow = (n+1), ncol = howManyTs)
  simul[1,] = rep(0,howManyTs)
  c <- 1
  for (t in 2:(n+1)){
    for (ts in 1:howManyTs ){
      simul[t,ts] <- c + phi[ts] * simul[t-1,ts] + noise[t,ts]
    }
  }
  return (simul[-1,])
}

# simulArma <- function(n, phi, theta, noise){
#   #we hard code the fact that models have order 2.
#   simul <- vector(length = n)
#   
#   #burn-in
#   endIdx <- n+10 #rough trick to have the noise index and the simul index aligned: basically we get the first noise value among the excess noise values generated.
#   simul[1] <- noise[endIdx] * theta[1] + noise[1] #only ma(1)
#   simul[2] <- noise[1] * theta[1] + noise[endIdx] * theta[2] + phi[1] * simul[1] + noise[2]#ma(1,2) ar(1)
#   simul[3] <- noise[2] * theta[1] + noise[1] * theta[2] + phi[1] * simul[2] + phi[2] * simul[1] + noise[3]#full
#   for (i in 4:n){
#     #time index of noise and simulated data are shifted of one
#     simul[i] <- noise[i-1] * theta[1] + noise[i-2] * theta[2] + phi[1] * simul[i-1] + phi[2] * simul[i-2] + noise[i]
#   }
#   return(simul)
# }

#return the two bottom time series and the sum time series
artificialTs <- function(n, correl, howMany=2){
  #we use arma(2,2) only
  # arOrders <- c(2,2)
  # maOrders <- c(2,2) #(runif(2) > 0.5)  + 1
  #the pars referring to the same  time series are in the same column: different columns refer to different models 
  arOrder <- 1
  if (! howMany==2) {
    stop("only two  bottom time series are supported")
  } 
  phi <- runif(howMany)
  noise <- drawNoise(n, correl, howMany)
  
  bottomTs <- simulVAR1(n, phi, noise)
  upperTs <- rowSums(bottomTs)
  hierarchy <- data.frame(cbind(bottomTs, upperTs))
  colnames(hierarchy) <- c("b1","b2","u")
  
  # debug
  print ("phi:")
  print(phi)
  expectedVar_B1 <- 1/(1-phi[1]^2)
  print("expected and empirical var of B1:")
  print (c(expectedVar_B1, var(hierarchy[,1])))
  
  expectedCovarB1B2 <- correl / (1 - prod(phi))
  empiricalCovarB1B2 <- cov(hierarchy)[1,2]
  
  print ("expected and empirical covar of B1 and B2:")
  print (c(expectedCovarB1B2, empiricalCovarB1B2))
  
  print ("expected and empirical correl of B1 and B2:")
  #we assume sigma^2=1, which disappears from below
  expectedCorrelB1B2 = expectedCovarB1B2 * sqrt( (1-phi[1]^2) * (1-phi[2]^2) )  
  empiricalCorrelB1B2 = cor(hierarchy)[1,2]
  print (c(expectedCorrelB1B2, empiricalCorrelB1B2))
  
  print ("expected and empirical variance of U:")
  #we assume sigma^2=1, which disappears from below
  varBottom = 1 / (1 - phi^2)
  expectedVar_U = varBottom[1] + varBottom[2] + 2 * expectedCorrelB1B2 * 
    sqrt(prod (varBottom))
  expectedVar_U_second = varBottom[1] + varBottom[2] + 2 * correl / 
    (1-prod (phi))
  
  empiricalVar_U = cov(hierarchy)[3,3]
  print (c(expectedVar_U, expectedVar_U_second, empiricalVar_U))
  
  print ("expected and empirical covariance (B1,U):")
  expectedCovB1_U = 1/ (1-phi[1]^2) + correl/(1-prod(phi))
  empiricalCovB1_U = cov(hierarchy)[1,3]
  print (c(expectedCovB1_U, empiricalCovB1_U))
  
  print ("expected and empirical correl of B1 and U:")
  expectedCorB1_U =  expectedCovB1_U / sqrt  ( 1/(1-phi[1]^2) * expectedVar_U )
  empiricalCorB1_U = cor(hierarchy)[1,3]
  print (c(expectedCorB1_U, empiricalCorB1_U))
  
  
  
  return (hierarchy)
}

#determines the value of sigma12
get_rho12<- function (phi, corB1_U) {
  
  #we assume everywhere noise with variance 1
  inv_sigmaB1 <- sqrt(1 - phi[1]^2)
  a <- 1/(1-phi[1]^2)
  b <- 1 - prod(phi)
  c <- 1/(1-phi[1]^2)
  d <- 1/(1-phi[2]^2)
  
  #range is +-1 becasue we are talking about a correlation
  inverse = function (f, lower = -1, upper = 1) {
    function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
  }
  
  
  correl_inverse = inverse(function (x) inv_sigmaB1 * (a + x/(1-prod(phi))) * 
                             1/(sqrt(c + d + 2 * x/b)) )
  
  rho12 = correl_inverse(corB1_U)
  return(rho12)
}


