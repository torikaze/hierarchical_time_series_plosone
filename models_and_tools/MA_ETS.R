library(dplyr);library(R6)


MovingAverage = function(data, lag, slide = TRUE){
  # this function calculate moving-average prediction
  # the number of fitted vaue is (nr  ow(data) - lag + 1).
  
  row_number = nrow(data)
  column_number = ncol(data)
  endpoint = row_number - lag + 1 # the number of predicted value
  
  # return slide prediction
  if(isTRUE(slide)){
    fitted_value = matrix(0, nrow = endpoint, ncol = column_number)
    if(lag == 1){ 
      # naive prediction
      fitted_value = data[1:endpoint,]
    }else{
      for(i in 1:endpoint){
        fitted_value[i,] = apply(data[i:(i + lag - 1),], 2, mean)
      }
    }
  }
  # return one prediction
  if(isFALSE(slide)){
    fitted_value = apply(data[(row_number - lag + 1):row_number,], 2, mean)
  }
  
  # insert NA to match the original data and fitted value
  insertNA = rep(NA, column_number)
  if(lag != 1){
    for(i in 1:(lag - 1)){
      insertNA = rbind(insertNA, rep(NA, column_number))
    }
  }
  fitted_value_withNA = rbind(insertNA, fitted_value)
  result = data.frame("original_data" = rbind(data, rep(NA, column_number)), "fitted_value" = fitted_value_withNA)
  return(result)
}

ETS = function(data, alpha, initial){
  # this function calculate Simple Exponential Smoothing prediction
  # the number of fitted vaue is (nrow(data) - lag + 1)
  row_number = nrow(data)
  column_number = ncol(data)
  
  # prediction
  fitted_value = rbind(initial, matrix(0, nrow = row_number, ncol = column_number))
  for(i in 1:(row_number)){
    fitted_value[i+1,] = alpha*data[i,] + (1 - alpha)*fitted_value[i,]
  }
  
  # insert NA to match the original data and fitted value
  insertNA = rep(NA, column_number)
  
  result = data.frame("original_data" = rbind(data, rep(NA, column_number)), "fitted_value" = fitted_value)
  return(result)
}
