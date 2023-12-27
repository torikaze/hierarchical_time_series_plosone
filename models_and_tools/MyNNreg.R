library(dplyr)
library(R6)
# feature*sample matrix is required
MyNNreg = R6Class(
  "MyNNreg",
  
  private = list(
    
  ),
  
  public = list(
    # initialization
    x = NULL,
    y = NULL,
    x_hts = NULL,
    hidden = NULL,
    out = NULL,
    hid_w = NULL,
    hid_b = NULL,
    out_w = NULL,
    out_b = NULL,
    unitvec = NULL,
    output = NULL,
    dW = NULL,
    smat = NULL,
    
    # Constructor. x_hts is the data that exclude bottom level series
    # smat is a summing-up matrix excluding bottom level
    initialize = function(x, y, x_hts, hidden_num, out_num, smat){
      self$x = x
      self$y = y
      self$x_hts = x_hts
      self$hidden = hidden_num
      self$out = out_num
      self$unitvec = rep(1, ncol(x))
      self$smat = smat
      # parameter initialization
        # hidden layer
      self$hid_w = matrix(rnorm(nrow(x)*hidden_num, 0, 1), hidden_num, nrow(x))
      self$hid_b = matrix(0, nrow = hidden_num)
        # output layer
      self$out_w = matrix(rnorm(hidden_num*out_num, 0, 1), out_num, hidden_num)
      self$out_b = matrix(0, nrow = out_num)
    },
    
    # sigmoid function and it's derivative
    sigm = function(x){
      return(1/(1 + exp(-x)))
    },
      #derivative
    sigm_deriv = function(x){
      temp = self$sigm(x)*(1-self$sigm(x))
      return(temp)
    },
    
    pred = function(data = self$x){
      # cut the paths to estimate single model
        # hidden layer
      num_each = self$hidden/self$out
      num_each_hid = nrow(data)/self$out
      temp1 = matrix(0, nrow = self$hidden, ncol = nrow(data))
      for(i in 1:self$out){
        n = 1 + num_each*(i-1)
        m = 1 + num_each_hid*(i-1)
        temp1[n:(n+num_each-1),m:(m+1)] = self$hid_w[n:(n+num_each-1),m:(m+1)]
      }
      self$hid_w = temp1
        # output layer
      temp2 = matrix(0, nrow = self$out, ncol = self$hidden)
      for(i in 1:self$out){
        n = 1 + num_each*(i-1)
        temp2[i,n:(n+num_each-1)] = self$out_w[i,n:(n+num_each-1)]
      }
      self$out_w = temp2
      
      # feedforward 
      uvec = rep(1, ncol(data))
      hid_u = self$hid_w %*% data + self$hid_b %*% t(uvec)
      hid_z = self$sigm(hid_u)
      out_u = self$out_w %*% hid_z + self$out_b %*% t(uvec)
      return(out_u)
    },
    
    # prediction and store the value in "output"
    feedforward = function(){
      self$output = self$pred()
      invisible(self)
    },
    
    # backpropagation
    backpropagation = function(eta = 1e-2, lambda_top = 0.1, lambda_mid = 0.1){
      # penalty matrix (diagonal)
      Lambda = diag(c(lambda_top, rep(lambda_mid, nrow(self$smat) - 1)))
    
      # gammma for top layers (differential of penalty terms)
      out_gamma = -1*t(self$smat) %*% Lambda^2 %*% (self$x_hts - self$smat %*% self$output)
    
      # delta (differential of output units + penalty terms)
      out_delta = -1*(self$y - self$output) + out_gamma
      hid_u = self$hid_w %*% self$x + self$hid_b %*% t(self$unitvec)
      hid_delta = self$sigm_deriv(hid_u)*(t(self$out_w) %*% out_delta)
      
      # update parameters
      hid_dw = hid_delta %*% t(self$x)
      hid_db = hid_delta %*% self$unitvec
      out_dw = out_delta %*% t(self$sigm(hid_u))
      out_db = out_delta %*% self$unitvec
      
      self$hid_w = self$hid_w - eta*hid_dw
      self$hid_b = self$hid_b - eta*hid_db
      self$out_w = self$out_w - eta*out_dw
      self$out_b = self$out_b - eta*out_db
      
      self$dW = hid_dw
      invisible(self)
    },
    
    train = function(iterations = 1e4, learn_rate = 1e-2, tolerance = 1e3, lambda_top = 0.1, lambda_mid = 0.1){
      error = matrix(0, nrow = iterations)
      error[1] = 1 # prevent dividing by 0 at first iteration
      error_fin = NULL
      for(i in 2:iterations){
        self$feedforward()$backpropagation(eta = learn_rate, lambda_top = lambda_top, lambda_mid = lambda_mid)
        error[i] = self$error_function(dx = self$hid_w, lambda_top = lambda_top, lambda_mid = lambda_mid)
        diff_percentage = abs((error[i] - error[i-1])/error[i-1])
      
        if(i > 10){
          if(diff_percentage < tolerance){
            error_fin = error[1:i]
            break
          }
        }
      }
      if(is.null(error_fin)){
        #cat("not converged")
        return(error[1:iterations])
      }else{
        #cat("converged")
        return(error_fin)
      }
    },
    
    error_function = function(dx, lambda_top, lambda_mid){
      # prediction error except penalty terms
      hid_u = dx %*% self$x + self$hid_b %*% t(self$unitvec)
      hid_z = self$sigm(hid_u)
      out_u = (self$out_w) %*% hid_z + self$out_b %*% t(self$unitvec)
      res = out_u - self$y
      res_diag = (res %*% t(res)) %>% diag()
      
      # penalty matrix (diagonal)
      Lambda = diag(c(lambda_top, rep(lambda_mid, nrow(self$smat)-1)))
      uvec = rep(1, self$out)
      uvec_penal = rep(1, nrow(Lambda))

      # gammma for top layers (differential of penalty terms)
      res_penal = Lambda %*% (self$x_hts - self$smat %*% out_u)
      
      res_penal_diag = (res_penal %*% t(res_penal)) %>% diag()
      
      return((t(res_diag) %*% uvec)/2 +(t(res_penal_diag) %*% uvec_penal)/2)
    },
    
    # numerical differentiation (to see if gradient is correctly calculated)
    numerical_diff = function(f = self$error_function, lambda_top, lambda_mid){
      h = .Machine$double.eps %>% sqrt()
      d1 = self$hid_w %>% dim()
      nd = matrix(0, nrow = d1[1], ncol = d1[2])
      temp1 = self$hid_w
      temp2 = self$hid_w
      for(i in 1:d1[1]){
        for(j in 1:d1[2]){
          temp1[i,j] = temp1[i,j] + h
          temp2[i,j] = temp2[i,j] - h
          nd[i,j] = (f(temp1, lambda_top = lambda_top, lambda_mid = lambda_mid) - f(temp2, lambda_top = lambda_top, lambda_mid = lambda_mid))/(2*h)
          temp1 = self$hid_w
          temp2 = self$hid_w
        }
      }
      return(nd)
    },
    
    # return accuracy of input data and naive prediction
    accuracy = function(data_pred, naive_pred = NULL, data_test){
      series_num = ncol(data_pred)
      rownames = colnames(data_test)
      rownames_naive = paste(rownames, "_naive", sep = "")
      result = matrix(0, nrow = series_num, ncol = 5)
      result_naive = matrix(0, nrow = series_num, ncol = 5)
      
      if(is.null(naive_pred) != TRUE){
        for(i in 1:series_num){
          result_naive[i,] = accuracy(naive_pred[,i], data_test[,i])
        }
      }
      for(i in 1:series_num){
        result[i,] = accuracy(data_pred[,i], data_test[,i])
      }
      rownames(result) = rownames
      rownames(result_naive) = rownames_naive
      return(rbind(result, result_naive))
    }
  )
)