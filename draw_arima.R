
#TODO:
#check if true parameters are eventually found
#which are the experiments? to generate with arma(2,2) and fit arima or ets

# draw parameters of an invertible and stationary arima model, see https://otexts.com/fpp2/AR.html
draw_ar <- function (arOrder){
  #drawing parameters   
  phi = vector(length = 2)
  if (arOrder ==1){
    phi[1]=runif(n=1, min=-1, max=1)
    phi[2]=0              
  }
  
  if (arOrder == 2){
    phi[2]=runif(n=1, min=-1, max=1)              
    upper = 1 - phi[2]
    lower = phi[2] - 1
    phi[1]= runif(1, min=lower, max=upper)
  }
  return(phi)
}


draw_ma <- function (maOrder){
  #drawing parameters   
  theta = vector(length = 2)
  if (maOrder == 1){
    theta[1]=runif(n=1, min=-1, max=1)
    theta[2]=0              
  }
  
  if (maOrder == 2){
    theta[2]=runif(n=1, min=-1, max=1)              
    upper = 1 + theta[2]
    lower = -theta[2] - 1
    theta[1]= runif(1, min=lower, max=upper)
  }
  return (theta)
}

#draw correlated noise for n time steps, return a vector (n x 2)
drawNoise <- function(n, correl){
  library(MASS)
  #generate an excess of noise observations
  noise <- mvrnorm(n=2*n, mu = c(0,0), Sigma= rbind(c(1,correl),c(correl,1)))
  return (noise)
}

simulArma <- function(n, phi, theta, noise){
  #we hard code the fact that models have order 2.
  simul <- vector(length = n)
  
  #burn-in
  endIdx <- n+10 #rough trick to have the noise index and the simul index aligned: basically we get the first noise value among the excess noise values generated.
  simul[1] <- noise[endIdx] * theta[1] + noise[1] #only ma(1)
  simul[2] <- noise[1] * theta[1] + noise[endIdx] * theta[2] + phi[1] * simul[1] + noise[2]#ma(1,2) ar(1)
  simul[3] <- noise[2] * theta[1] + noise[1] * theta[2] + phi[1] * simul[2] + phi[2] * simul[1] + noise[3]#full
  for (i in 4:n){
    #time index of noise and simulated data are shifted of one
    simul[i] <- noise[i-1] * theta[1] + noise[i-2] * theta[2] + phi[1] * simul[i-1] + phi[2] * simul[i-2] + noise[i]
  }
  return(simul)
}

#return the two bottom time series and the sum time series
artificialTs <- function(n, correl){
  #we use arma(2,2) only
  # arOrders <- c(2,2)
  # maOrders <- c(2,2) #(runif(2) > 0.5)  + 1
  #the pars referring to the same  time series are in the same column: different columns refer to different models 
  armaOrder <- 2
  models <- 2
  phi <- matrix(nrow = armaOrder, ncol = models)
  theta <- phi
  phi[,1] <- draw_ar(armaOrder)
  phi[,2] <- draw_ar(armaOrder)
  theta[,1] <- draw_ma(armaOrder)
  theta[,2] <- draw_ma(armaOrder)
  noise <- drawNoise(n, correl)
  
  timeSeries <- matrix(nrow = n, ncol = models)
  for (i in 1:models){
    timeSeries[,i] <- simulArma(n, phi[,i],theta[,i],noise[,i])
  }
  #it is debugged that the true parameters are recovered by autorima if we set n to a large value
  a <- auto.arima(timeSeries[,2])
}



