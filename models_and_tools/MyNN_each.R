library(dplyr)
library(R6)
# feature*sample matrix is required
MyNN_each = R6Class(
  "MyNN_each",
  
  private = list(
    
  ),
  
  public = list(
    # initialization
    x = NULL,
    y = NULL,
    hidden = NULL,
    out = NULL,
    hid_w = NULL,
    hid_b = NULL,
    out_w = NULL,
    out_b = NULL,
    unitvec = NULL,
    output = NULL,
    dW = NULL,

    initialize = function(x, y, hidden, out){
      self$x = x
      self$y = y
      self$hidden = hidden
      self$out = out
      self$unitvec = rep(1, ncol(x))
     
      self$hid_w = matrix(rnorm(nrow(x)*self$hidden, 0, 1), self$hidden, nrow(x))
      self$hid_b = matrix(0, nrow = hidden)
     
      self$out_w = matrix(rnorm(self$hidden*self$out, 0, 1), self$out, self$hidden)
      self$out_b = matrix(0, nrow = self$out)
    },
    
    # sigmoid function and it'S derivative
    sigm = function(x){
      return(1/(1 + exp(-x)))
    },
    #derivative
    sigm_deriv = function(x){
      temp = self$sigm(x)*(1-self$sigm(x))
      return(temp)
    },
    
    # return the predicted value
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
    backpropagation = function(eta = 1e-2){
      out_delta = self$output - self$y
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
    
    train = function(iterations = 1e4, learn_rate = 1e-2, tolerance = 1e3){
      error = matrix(0, nrow = iterations)
      error[1] = 1 # prevent dividing by 0 at first iteration
      error_fin = NULL
      for(i in 2:iterations){
        self$feedforward()$backpropagation(learn_rate)
        error[i] = self$error_function(dx = self$hid_w)
        diff_percentage = abs((error[i] - error[i-1])/error[i-1])
        if(i > 10){
          if(diff_percentage < tolerance){
            error_fin = error[1:i]
            break
          }
        }
      }
      return(error_fin)
    },
    
    error_function = function(dx){
      uvec = rep(1, self$out)
      hid_u = dx %*% self$x + self$hid_b %*% t(self$unitvec)
      hid_z = self$sigm(hid_u)
      out_u = (self$out_w) %*% hid_z + self$out_b %*% t(self$unitvec)
      res = out_u - self$y
      res_diag = (res %*% t(res)) %>% diag()
      return((res_diag %*% uvec)/2)
    },
    
    numerical_diff = function(f = self$error_function){
      h = .Machine$double.eps %>% sqrt()
      d1 = self$hid_w %>% dim()
      nd = matrix(0, nrow = d1[1], ncol = d1[2])
      temp1 = self$hid_w
      temp2 = self$hid_w
      for(i in 1:d1[1]){
        for(j in 1:d1[2]){
          temp1[i,j] = temp1[i,j] + h
          temp2[i,j] = temp2[i,j] - h
          nd[i,j] = (f(temp1) - f(temp2))/(2*h)
          temp1 = self$hid_w
          temp2 = self$hid_w
        }
      }
      return(nd)
    },
    
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