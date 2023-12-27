library(R6);library(dplyr);library(tseries)
TSset = R6Class(
  "TSset",
  private = list(
    
  ),
  
  public = list(
    x = NULL,
    initialize = function(x){
      self$x = as.matrix(x)
    },
    
    arrange = function(nodes){
      library(hts)
      return(as.matrix(aggts(hts(self$x, nodes = nodes))))
    },
    
    lagset = function(lag){
      series_num = ncol(self$x)
      bottom = matrix(0, ncol = series_num*lag, nrow = nrow(self$x))
      outcome_num = rep(0, series_num)
      
      for(i in 1:series_num){
        for(j in 1:lag){
          bottom[,(i-1)*lag+j] = lag(as.matrix(self$x[,i]), j-1)
        }
      }
      name_x = colnames(self$x)
      name = c()
      for(i in 1:series_num){
        for(j in 1:lag){
          name[j+(i-1)*lag] = paste(name_x[i], j, sep = "")
        }
      }
      colnames(bottom) = name
      for(i in 1:series_num){
        outcome_num[i] = 1 + (i-1)*lag
      }
      bottom = na.omit(bottom)
      return(cbind(self$x[-(1:(lag-1)),], bottom[, -outcome_num]))
    },
    
    split_traintest = function(){
      
    },
    
    adf = function(data = self$x){
      result = matrix(0, nrow = ncol(data))
      for(i in 1:ncol(data)){
        result[i] = adf.test(data[,i])$p.value
      }
      rownames(result) = colnames(data)
      return(result)
    }
    
  )
)