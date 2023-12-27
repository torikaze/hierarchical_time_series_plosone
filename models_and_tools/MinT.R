library(R6);library(dplyr); library(corpcor)

MinT = R6Class(
  classname = "MinT",
  
  public = list(
    x = NULL,
    samplecov = NULL,
    samplecor = NULL,
    res_normal = NULL,
    
    initialize = function(pred, train_res){
      self$x = pred
      self$samplecov = cov(train_res)
      self$samplecor = cor(train_res)
      self$res_normal = scale(train_res)
      #print("pred must be samplesize*features matrix")
    },
    
    reconcile = function(smatrix){
      w_inv = solve(self$samplecov)
      P = solve(t(smatrix) %*% w_inv %*% smatrix) %*% t(smatrix) %*% w_inv
      y_recon = smatrix %*% P %*% t(self$x)
      return(y_recon)
    },
    
    reconcile_shrinkage = function(smatrix){
      R = self$samplecor
      lambda_denominator = (sum(R^2) - nrow(R))/2
      
      x = self$res_normal
      W = array(0, dim = c(ncol(x), ncol(x), nrow(x)))
      for (i in 1:nrow(x)) {
        for (j in 1:ncol(x)) {
          for (k in 1:ncol(x)) {
            W[j, k, i] = x[i, j] * x[i, k]
          }
        }
      }
      r_var = matrix(0, nrow = ncol(x), ncol = ncol(x))
      for (i in 1:ncol(x)) {
        for (j in 1:ncol(x)) {
          r_var[i, j] = sum((W[i, j,] - mean(W[i, j,]))^2)
        }
      }
      r_var = r_var * nrow(x) / (nrow(x) - 1)^3
      lambda_numerator = sum(r_var - diag(diag(r_var)))/2
      lambda = lambda_numerator / lambda_denominator
      W_shrinkage = lambda * diag(diag(self$samplecov)) + (1 - lambda) * self$samplecov
      W_inv = solve(W_shrinkage)
      P = solve(t(smatrix) %*% W_inv %*% smatrix) %*% t(smatrix) %*% W_inv
      y_recon_shrink = smatrix %*% P %*% t(self$x)
      return(y_recon_shrink)
    },
    
    OCR = function(smatrix){
      P = solve(t(smatrix) %*% smatrix) %*% t(smatrix)
      y_recon = smatrix %*% P %*% t(self$x)
      return(y_recon)
    }
  )
)
