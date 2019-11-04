# %values to be reconciled
x_0= c(5,-3)
sigma_0 = matrix( data = c(1,0.5,0.5,1), nrow = 2, byrow = TRUE)

library(Matrix)

states <- nrow(sigma_0)
currentSigma <- sigma_0

for (st in 1:states){
  # the SVD: A = U * D * V', where U and V are orthogonal
  decomp <- svd(sigma_0)
  D <-  diag(decomp$d) # d contains the singular values of D
  
  #dcomposition idea:
  # sigma_0 = decomp$u %*% D %*% t(decomp$v) 

  # has been the decomposition succesfull?
  if (norm( (decomp$u %*% t(decomp$u)) - diag(states)) > 1e-7 ){
    stop ("u is not orthogonal")
  }
  
  if (norm( (decomp$u %*%  D %*% t(decomp$v)) - currentSigma) > 1e-7 ){
    stop ("decomposition failed")
  }
  
  }
