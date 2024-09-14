library(dplyr);library(xts)
library(tseries);library(hts)
library(forecast); library(tictoc)
library(stats); library(tictoc)
source('MA_ETS.R')
source('MinT.R')
source('TSset.R')
source('MyNNreg.R')
source('MyNN_each.R')

NameOfData = "neg"
wind_part = read.csv(paste("data/hts_sim_", NameOfData, ".csv", sep = ""), fileEncoding = "Shift_JIS")

# seasonal decomposition
### line 16~24 Applies only to unemployment rate data (unemployment_UE_~.csv) ###
# temp = data.frame(rep(0,dim(wind_part)[1]))
# for(i in 1:dim(wind_part)[2]){
#   ts = ts(as.numeric(wind_part[,i]), frequency = 4)
#   a_stl = stl(ts, s.window = "periodic")
#   aaa = a_stl$time.series[,3]
#   temp = cbind(temp, as.data.frame(aaa))
# }
# colnames(temp) = c("", "AA", "AB", "AC", "BA", "BB", "BC", "CA", "CB", "CC")
# wind_part = temp[, -1]


 # see if each series has a unit root
tsset = TSset$new(x = wind_part)
tsset$adf()
smatrix = wind_part %>% hts(nodes = list(3,c(3,3,3))) %>% smatrix()

# compute a lagged series
wind_hts = tsset$arrange(nodes = list(3, c(3,3,3)))
lset = TSset$new(x = wind_hts)
wind_hts_lag = lset$lagset(lag = 3)

# split the data into test and train
split_point = floor(nrow(wind_hts)*2/3)
data_train_mint = wind_hts[1:(nrow(wind_hts)*2/3),]
data_test_mint = wind_hts[(nrow(wind_hts)*2/3-1):nrow(wind_hts),]
data_train = wind_part[1:(nrow(wind_hts)*2/3),]
data_test  = wind_part[(nrow(wind_hts)*2/3-1):nrow(wind_part),]

naive_test  = wind_hts_lag[(nrow(wind_hts)*2/3-1):nrow(wind_hts_lag),c(14,16,18,20,22,24,26,28,30,32,34,36,38)]
test = data_test_mint[-c(1:2),]

# standardization
# standardize after aggregation for MinT
mean_mint = apply(data_train_mint, 2, mean)
sd_mint = apply(data_train_mint, 2, sd)
train_mean_mint = (mean_mint %*% t(rep(1, nrow(data_test_mint)))) %>% t()
train_sd_mint = (sd_mint %*% t(rep(1, nrow(data_test_mint)))) %>% t()

train_mean_mint_temp = (mean_mint %*% t(rep(1, nrow(data_train_mint)))) %>% t()
train_sd_mint_temp = (sd_mint %*% t(rep(1, nrow(data_train_mint)))) %>% t()

train_mint_scale = scale(data_train_mint)
test_mint_scale = (data_test_mint - train_mean_mint)/train_sd_mint

set = TSset$new(x = train_mint_scale)
hts_scale_train_lag_mint = set$lagset(lag = 3)
set = TSset$new(x = test_mint_scale)
hts_scale_test_lag_mint = set$lagset(lag = 3)

# standardize before aggregation
data_mean = apply(data_train, 2, mean)
data_sd = apply(data_train, 2, sd)  
train_mean =(data_mean %*% t(rep(1, nrow(data_test)))) %>% t() 
train_sd = (data_sd %*% t(rep(1, nrow(data_test)))) %>% t()

data_train_scale = apply(data_train, 2, scale)
data_test_scale =  (data_test - train_mean)/train_sd

scale_set = TSset$new(x = data_train_scale)
hts_scale_train = scale_set$arrange(nodes = list(3, c(3,3,3)))
scale_set = TSset$new(x = data_test_scale)
hts_scale_test = scale_set$arrange(nodes = list(3, c(3,3,3)))

set = TSset$new(x = hts_scale_train)
hts_scale_train_lag = set$lagset(lag = 3)
set = TSset$new(x = hts_scale_test)
hts_scale_test_lag = set$lagset(lag = 3)

# hts_scale_train_lag_mint:  hts -> standardization -> lag (for MinT, which requires all nodes' value)
# hts_scale_train_lag  bottom: -> standardization -> sum-up -> lag (for reguralization, bottom must be added up to upper series)
# train_mean;train_sd:  for standardization
# train_mean_mint;train_sd_mint:  for standardization of Mint
# naive_test:  naive prediction
# test: test data

# predictions
data_train_t = hts_scale_train_lag %>% t()
data_test_t = hts_scale_test_lag %>% t()
d = dim(data_train_t)

mint_train_t = hts_scale_train_lag_mint %>% t()
mint_test_t = hts_scale_test_lag_mint %>% t()
d2 = dim(mint_train_t)


learn_rate <- 1e-5
tolerance <- 5e-5
# 2-folds cross-validation to decide penalty parameters
lambda_top_range = seq(0, 0.4, by = 0.04)
lambda_mid_range = seq(0, 3, by = 0.3)
rmse_cv = c()
cv_length = (dim(data_train_t)[2]*0.3) %>% round(digits = 0)
cv_score_1val = data.frame()
k = cv_length
for(i in lambda_top_range){
  for(j in lambda_mid_range){
    # fit the model using first N-k data
    set.seed(100)
    model_Mnnreg_cv = MyNNreg$new(x = data_train_t[22:39, 1:(dim(data_train_t)[2] - k)], y = data_train_t[5:13, 1:(dim(data_train_t)[2] - k)], x_hts =  data_train_t[1:4, 1:(dim(data_train_t)[2] - k)], hidden = 36, out = 9, smat = smatrix[1:4,])
    suppressMessages(model_Mnnreg_cv$train(iterations = 4e3, learn_rate = 1e-5, tolerance = tolerance, lambda_top = i, lambda_mid = j))
    # 1 step ahead forecast
    pred_Mnnreg_cv = model_Mnnreg_cv$pred(data = data_train_t[22:39,]) # predict all
    scale_set = TSset$new(x = t(pred_Mnnreg_cv[, (dim(data_train_t)[2] - k + 1):dim(data_train_t)[2]])) # pick up 1 step ahead forecast
    cv_hts = suppressMessages(scale_set$arrange(nodes = list(3, c(3,3,3)))) # sum up
    # considering the size of values, we might have to scale the top and the middle forecasts
    # but calculate rmse without scaling this time
    rmse_cv = ((t(data_train_t[1:13, (dim(data_train_t)[2] - k + 1):dim(data_train_t)[2]]) - cv_hts)^2 %*% rep(1, 13))/13
    cv_score_1val = rbind(cv_score_1val, cbind(mean(rmse_cv), i, j))
  }
  cat("\n ###", i/lambda_top_range[length(lambda_top_range)]*100, "% ### \n")
}
# write.csv(cv_score_1val, paste0('tables/cv_score_', NameOfData, '.csv'), row.names = FALSE)
# cv_score_1val = read.csv(paste("tables/cv_score_", NameOfData, ".csv", sep = ""))
(lt = (cv_score_1val %>% filter(V1 == min(V1)))[,2])
(lm = (cv_score_1val %>% filter(V1 == min(V1)))[,3])

# calculate RMSE with confidence intervals
rep_num = 30
result_NN = c()
result_NNreg = c()
result_Mint = c()
result_Mint_shrinkage = c()
for(i in 1:rep_num){
  # base forecasts
  set.seed(i)
  model_Mnn_bottom = MyNNreg$new(x = data_train_t[22:39,], y = data_train_t[5:13,], x_hts = data_train_t[1:4,], hidden = 36, out = 9, smat = smatrix[1:4,])
  model_Mnn_bottom$train(learn_rate = learn_rate, iterations = 4000, tolerance = tolerance, lambda_top = 0, lambda_mid = 0) %>% length()
  pred_Mnn_test  = model_Mnn_bottom$pred(data = data_test_t[22:39,])
  # base forecasts with penalty
  set.seed(i)
  model_Mnnreg_bottom = MyNNreg$new(x = data_train_t[22:39,], y = data_train_t[5:13,], x_hts =  data_train_t[1:4,], hidden = 36, out = 9, smat = smatrix[1:4,])
  model_Mnnreg_bottom$train(iterations = 4000, learn_rate = learn_rate, tolerance = tolerance, lambda_top = lt, lambda_mid = lm) %>% length()
  pred_Mnnreg_test  = model_Mnnreg_bottom$pred(data = data_test_t[22:39,])
  # MinT
  set.seed(i)
  model_Mnn_mint = MyNN_each$new(x = data_train_t[14:39,], y = data_train_t[1:13,], hidden = 52, out = 13)
  model_Mnn_mint$train(iterations = 40000, learn_rate = learn_rate, tolerance = tolerance) 
  pred_mint_test  = model_Mnn_mint$pred(data = data_test_t[14:39,])
  
  # MyNN
  output = t(pred_Mnn_test)*train_sd[-c(1:2),] + train_mean[-c(1:2),]
  output_set = TSset$new(x = output)
  output_hts = suppressMessages(output_set$arrange(nodes = list(3, c(3,3,3))))
  # NN with penalty
  output_p = t(pred_Mnnreg_test)*train_sd[-c(1:2),] + train_mean[-c(1:2),]
  output_set_p = TSset$new(x = output_p)
  output_hts_p = suppressMessages(output_set_p$arrange(nodes = list(3, c(3,3,3))))
  # MinT
  mint_base = t(pred_mint_test)*train_sd_mint[-c(1:2),] + train_mean_mint[-c(1:2),]
  train_res = t(model_Mnn_mint$pred())*train_sd_mint_temp[-(1:2),] + train_mean_mint_temp[-(1:2),] - data_train_mint[-(1:2),]
  mint = MinT$new(pred = mint_base, train_res = train_res)
  output_mint_shrinkage = mint$reconcile_shrinkage(smatrix = smatrix) %>% t()
    # accuracy
  result = model_Mnn_bottom$accuracy(data_pred = output_hts, naive_pred = naive_test, data_test = test)
  result_p = model_Mnn_bottom$accuracy(data_pred = output_hts_p, naive_pred = naive_test, data_test = test)
  result_mint_shrinkage = model_Mnn_bottom$accuracy(data_pred = output_mint_shrinkage, naive_pred = naive_test, data_test = test)
  
  result_NN = cbind(result_NN, result[1:13,2])
  result_NNreg = cbind(result_NNreg, result_p[1:13,2])
  result_Mint_shrinkage = cbind(result_Mint_shrinkage, result_mint_shrinkage[1:13, 2])
  cat("process:", i/rep_num*100, "%\n")
}
CI <- qt(p = 0.975, df = rep_num -1)
NN_mid_sd = (apply(result_NN[2:4,], 2, mean) %>% sd())/sqrt(rep_num) * CI
NN_bottom_sd = (apply(result_NN[5:13,], 2, mean) %>% sd())/sqrt(rep_num) * CI
NN_all_sd = (apply(result_NN, 2, mean) %>% sd())/sqrt(rep_num) * CI
NNreg_mid_sd = (apply(result_NNreg[2:4,], 2, mean) %>% sd())/sqrt(rep_num) * CI
NNreg_bottom_sd = (apply(result_NNreg[5:13,], 2, mean) %>% sd())/sqrt(rep_num) * CI
NNreg_all_sd = (apply(result_NNreg, 2, mean) %>% sd())/sqrt(rep_num) * CI
Mint_shrinkage_mid_sd = (apply(result_Mint_shrinkage[2:4,], 2, mean) %>% sd())/sqrt(rep_num) * CI
Mint_shrinkage_bottom_sd = (apply(result_Mint_shrinkage[5:13,], 2, mean) %>% sd())/sqrt(rep_num) * CI
Mint_shrinkage_all_sd = (apply(result_Mint_shrinkage, 2, mean) %>% sd())/sqrt(rep_num) * CI

result_CI_mean = data.frame(NN_mid_CI = NN_mid_sd, NN_bottom_CI = NN_bottom_sd, NN_all_CI = NN_all_sd, 
                            NNreg_mid_CI = NNreg_mid_sd, NNreg_bottom_CI = NNreg_bottom_sd, NNreg_all_CI = NNreg_all_sd, 
                            Mint_shrinkage_mid_CI = Mint_shrinkage_mid_sd, Mint_shrinkage_bottom_CI = Mint_shrinkage_bottom_sd, Mint_shrinkage_all_CI = Mint_shrinkage_all_sd) %>% t()
result_data = round(data.frame(NN = apply(result_NN, 1, mean), NN_sd = apply(result_NN, 1, sd), 
                               NNreg = apply(result_NNreg, 1, mean), NNreg_sd = apply(result_NNreg, 1, sd), 
                               result_MinT_shrinkage = apply(result_Mint_shrinkage, 1, mean), result_MinT_shrinkage_sd = apply(result_Mint_shrinkage, 1, sd)), digits = 2)


# MovingAverage
MALagDecide = function(data, slide = TRUE, lag_max, lag_min){
  # @return lag number and rmse by cross-validation
  cv_result = data.frame()
  cv_range = seq(lag_min, lag_max, by = 1)
  for(lag in cv_range){
    pred = MovingAverage(data = data, lag = lag, slide = slide) %>% as.matrix()
    range_compare = (lag_max + 1):nrow(data - 1)
    mean_squared_error = mean((pred[range_compare, 1:ncol(data)] - pred[range_compare, (ncol(data) + 1):ncol(pred)])^2)
    rmse = sqrt(mean_squared_error)
    cv_result = rbind(cv_result, cbind(rmse, lag))
  }
  return(cv_result)
}
MA_lag_cv = MALagDecide(data = data_train_mint, slide = TRUE, lag_max = 20, lag_min = 1)
MA_lag = (MA_lag_cv %>% filter(rmse == min(rmse)))$lag
data_for_MA = rbind(data_train_mint, data_test_mint[-c(1:2),])
pred_MA = MovingAverage(data = data_for_MA[(split_point - MA_lag + 1):(nrow(data_for_MA)),], lag = MA_lag, slide = TRUE)
pred_MA_test = pred_MA[-c(1:MA_lag, nrow(pred_MA)), -c(1:ncol(data_test_mint))]
# Exponential Smoothing
ETSAlphaDecide = function(data, by = 0.01, initial){
  # cv for ETS. 0<=alpha<=1
  alpha_range = seq(0, 1, by = by)
  cv_result = data.frame()
  for(alpha in alpha_range){
    pred = ETS(data = data, alpha = alpha, initial = initial) %>% as.matrix()
    range_compare = 2:(nrow(pred) - 1)
    mean_squared_error = mean((pred[range_compare, 1:ncol(data)] - pred[range_compare, (ncol(data) + 1):ncol(pred)])^2)
    rmse = sqrt(mean_squared_error)
    cv_result = rbind(cv_result, cbind(rmse, alpha))
  }
  return(cv_result)
}
ETS_alpha_cv = ETSAlphaDecide(data = data_train_mint, by = 0.01, initial = data_train_mint[1,])
ETS_alpha = (ETS_alpha_cv %>% filter(rmse == min(rmse)))$alpha
# due to it's prediction algorithm, use all data
# because of index-related problem (line 30~35), use wind_hts[-nrow(wind_hts),]
pred_ETS = ETS(data = wind_hts[-nrow(wind_hts),], alpha = ETS_alpha, initial = wind_hts[1,])
pred_ETS_test = pred_ETS[(split_point + 1):(nrow(pred_ETS) - 1), (ncol(data_test_mint) + 1):ncol(pred_ETS)]

# accuracy
result_MA = model_Mnn_bottom$accuracy(data_pred = pred_MA_test, naive_pred = naive_test, data_test = test)
result_ETS = model_Mnn_bottom$accuracy(data_pred = pred_ETS_test, naive_pred = naive_test, data_test = test)
ma_round = round(result_MA, digits = 2)
ets_round = round(result_ETS, digits = 2)
MA_ETS_result = cbind(ma_round[1:13, 2], ets_round[1:13, 2])
ETS_mid_average = apply(MA_ETS_result[2:4,], 2, mean)
ETS_bottom_average = apply(MA_ETS_result[5:13,], 2, mean)
ETS_average = apply(MA_ETS_result, 2, mean)
MA_ETS_result_withAverage = rbind(MA_ETS_result[1:4,], ETS_mid_average, MA_ETS_result[5:13,], ETS_bottom_average, ETS_average) %>% round(digits = 2)
colnames(MA_ETS_result_withAverage) = c(paste("MA(", MA_lag, ")", sep = ""), paste("ETS(", ETS_alpha, ")", sep = ""))


# accuracy
result_data
MA_ETS_result_withAverage

