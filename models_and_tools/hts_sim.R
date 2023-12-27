### for better simulation data ###
library(dplyr);library(tseries)

# each bottom series is generated from common factor "sim_commonN"
# control correlation by the coefficients of "sim_commonN"


# common factor 1, 2, 3 and all
# common factor 1, 2, 3 affects 3 nodes respectively
# common_all affects all series
set.seed(10)
sim_common1 = rep(0, 100)
for(i in 1:99){
  sim_common1[i+1] = 0.3*sim_common1[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(11)
sim_common2 = rep(0, 100)
for(i in 1:99){
  sim_common2[i+1] = 0.3*sim_common2[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(12)
sim_common3 = rep(0, 100)
for(i in 1:99){
  sim_common3[i+1] = 0.3*sim_common3[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(13)
sim_common_all = rep(0, 100)
for(i in 1:99){
  sim_common_all[i+1] = 0.3*sim_common_all[i] + rnorm(1, mean = 0, sd = 0.3)
}

### bottom series ###
  ## negative correlation ##
set.seed(1)
sim1 = rep(0,100)
for(i in 1:99){
  sim1[i+1] = 0.1*sim_common_all[i+1] + 1*sim_common1[i+1] + 0.3*sim1[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(2)
sim2 = rep(0,100)
for(i in 1:99){
  sim2[i+1] = -0.1*sim_common_all[i+1] + -1*sim_common1[i+1] + 0.3*sim2[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(3)
sim3 = rep(0,100)
for(i in 1:99){
  sim3[i+1] = 1*sim_common_all[i+1] + 0.1*sim_common1[i+1] + 0.3*sim3[i] + rnorm(1, mean = 0, sd = 0.3)
}


set.seed(4)
sim4 = rep(0,100)
for(i in 1:99){
  sim4[i+1] = 0.1*sim_common_all[i+1] + 1*sim_common2[i+1] + 0.3*sim4[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(5)
sim5 = rep(0,100)
for(i in 1:99){
  sim5[i+1] = -0.1*sim_common_all[i+1] + -1*sim_common2[i+1] + 0.3*sim5[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(6)
sim6 = rep(0,100)
for(i in 1:99){
  sim6[i+1] = -1*sim_common_all[i+1] + 0.1*sim_common2[i+1] + 0.3*sim6[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(7)
sim7 = rep(0,100)
for(i in 1:99){
  sim7[i+1] = 0.1*sim_common_all[i+1] + 1*sim_common3[i+1] + 0.3*sim7[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(8)
sim8 = rep(0,100)
for(i in 1:99){
  sim8[i+1] = -0.1*sim_common_all[i+1] + -1*sim_common3[i+1] + 0.3*sim8[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(9)
sim9 = rep(0,100)
for(i in 1:99){
  sim9[i+1] = 1*sim_common_all[i+1] + 0.1*sim_common3[i+1] + 0.3*sim9[i] + rnorm(1, mean = 0, sd = 0.3)
}

data_sim_neg = cbind(AA = sim1, AB = sim2, AC = sim3, BA = sim4, BB = sim5, BC = sim6, CA = sim7, CB = sim8, CC = sim9)
write.csv(data_sim_neg, "C:/Users/Owner/OneDrive/ドキュメント/研究室/追加実験資料_表/hts_sim_neg.csv", row.names = FALSE)
 
  ## uncorrelated ##
set.seed(1)
sim1 = rep(0,100)
for(i in 1:99){
  sim1[i+1] = 0.1*sim_common_all[i+1] + 0.1*sim_common1[i+1] + 0.3*sim1[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(2)
sim2 = rep(0,100)
for(i in 1:99){
  sim2[i+1] = 0.1*sim_common_all[i+1] + 0.1*sim_common1[i+1] + 0.3*sim2[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(3)
sim3 = rep(0,100)
for(i in 1:99){
  sim3[i+1] = 0.1*sim_common_all[i+1] + 0.1*sim_common1[i+1] + 0.3*sim3[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(4)
sim4 = rep(0,100)
for(i in 1:99){
  sim4[i+1] = 0.1*sim_common_all[i+1] + 0.1*sim_common2[i+1] + 0.3*sim4[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(5)
sim5 = rep(0,100)
for(i in 1:99){
  sim5[i+1] = 0.1*sim_common_all[i+1] + 0.1*sim_common2[i+1] + 0.3*sim5[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(6)
sim6 = rep(0,100)
for(i in 1:99){
  sim6[i+1] = 0.1*sim_common_all[i+1] + 0.1*sim_common2[i+1] + 0.3*sim6[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(7)
sim7 = rep(0,100)
for(i in 1:99){
  sim7[i+1] = 0.1*sim_common_all[i+1] + 0.1*sim_common3[i+1] + 0.3*sim7[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(8)
sim8 = rep(0,100)
for(i in 1:99){
  sim8[i+1] = 0.1*sim_common_all[i+1] + 0.1*sim_common3[i+1] + 0.3*sim8[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(9)
sim9 = rep(0,100)
for(i in 1:99){
  sim9[i+1] = 0.1*sim_common_all[i+1] + 0.1*sim_common3[i+1] + 0.3*sim9[i] + rnorm(1, mean = 0, sd = 0.3)
}

data_sim_ind = cbind(AA = sim1, AB = sim2, AC = sim3, BA = sim4, BB = sim5, BC = sim6, CA = sim7, CB = sim8, CC = sim9)
write.csv(data_sim_ind, "C:/Users/Owner/OneDrive/ドキュメント/研究室/追加実験資料_表/hts_sim_ind.csv", row.names = FALSE)


### positive correlation ###
set.seed(1)
sim1 = rep(0,100)
for(i in 1:99){
  sim1[i+1] = sim_common_all[i+1] + sim_common1[i+1] + 0.3*sim1[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(2)
sim2 = rep(0,100)
for(i in 1:99){
  sim2[i+1] = sim_common_all[i+1] + sim_common1[i+1] + 0.3*sim2[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(3)
sim3 = rep(0,100)
for(i in 1:99){
  sim3[i+1] = 1*sim_common_all[i+1] + 1*sim_common1[i+1] + 0.3*sim3[i] + rnorm(1, mean = 0, sd = 0.3)
}


set.seed(4)
sim4 = rep(0,100)
for(i in 1:99){
  sim4[i+1] = 1*sim_common_all[i+1] + 1*sim_common2[i+1] + 0.3*sim4[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(5)
sim5 = rep(0,100)
for(i in 1:99){
  sim5[i+1] = 1*sim_common_all[i+1] + 1*sim_common2[i+1] + 0.3*sim5[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(6)
sim6 = rep(0,100)
for(i in 1:99){
  sim6[i+1] = 1*sim_common_all[i+1] + 1*sim_common2[i+1] + 0.3*sim6[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(7)
sim7 = rep(0,100)
for(i in 1:99){
  sim7[i+1] = 1*sim_common_all[i+1] + 1*sim_common3[i+1] + 0.3*sim7[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(8)
sim8 = rep(0,100)
for(i in 1:99){
  sim8[i+1] = 1*sim_common_all[i+1] + 1*sim_common3[i+1] + 0.3*sim8[i] + rnorm(1, mean = 0, sd = 0.3)
}

set.seed(9)
sim9 = rep(0,100)
for(i in 1:99){
  sim9[i+1] = 1*sim_common_all[i+1] + 1*sim_common3[i+1] + 0.3*sim9[i] + rnorm(1, mean = 0, sd = 0.3)
}

data_sim_pos = cbind(AA = sim1, AB = sim2, AC = sim3, BA = sim4, BB = sim5, BC = sim6, CA = sim7, CB = sim8, CC = sim9)
write.csv(data_sim_pos, "C:/Users/Owner/OneDrive/ドキュメント/研究室/追加実験資料_表/hts_sim_pos.csv", row.names = FALSE)
