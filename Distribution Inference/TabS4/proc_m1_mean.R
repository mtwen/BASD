## MengtaoWEN
## 2024-07-17
## Executed file for different methods

rm(list = ls()); gc()
set.seed(20240722)
source('basis.R')
library(doSNOW)

n = 1000
m = 10*n
num_simu = 200
ticks = 100 - 1
taus = seq(1/(ticks + 1), 1 - 1/(ticks + 1), by = 1/(ticks + 1))
alpha = 0.1
y = qnorm(taus, 0, sd = sqrt(1+1/4)) 
mean.gold = 0

################################################################################
#### p = 100
################################################################################
p = 100
p1 = ceiling(0.1*p)
beta0 = c(rep(1, p1), rep(0, p-p1))
filename.true <- paste0('RES_dist_M1_trueCCDF_n_', n, '_p_', p, '_m_', m, '_p1_', p1, 
                        '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

filename.gamlss <- paste0('RES_dist_M1_gamlss_n_', n, '_p_', p, '_m_', m, '_p1_', p1, 
                          '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

filename.drf <- paste0('RES_dist_M1_drf_n_', n, '_p_', p, '_m_', m, '_p1_', p1, 
                       '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

filename.eng <- paste0('RES_dist_M1_engression_n_', n, '_p_', p, '_m_', m, '_p1_', p1, 
                       '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

filename.mean <- paste0('RES_mean_M1_n_', n, '_p_', p, '_m_', m, '_p1_', p1,
                   '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

#### Processing BASD*

load(file = filename.true)
RES.true <- RES
rm(RES)

MSE.true = rowMeans(sapply(1:num_simu, function(num){
  (apply(RES.true[[1, num]], 2, function(X){
    sum(c(X[1], diff(X, lag = 1, differences = 1))*y)
  }) - mean.gold)^2
}))
MSEr.true <- MSE.true/MSE.true[2]

#### Processing GAMLSS

load(file = filename.gamlss)
RES.GAMLSS <- RES
rm(RES)

MSE.gamlss = rowMeans(sapply(1:num_simu, function(num){
  (apply(RES.GAMLSS[[1, num]], 2, function(X){
    sum(c(X[1], diff(X, lag = 1, differences = 1))*y)
  }) - mean.gold)^2
}))
MSEr.gamlss <- MSE.gamlss/MSE.gamlss[5]

#### Processing DRF

load(file = filename.drf)
RES.DRF <- RES
rm(RES)

MSE.drf = rowMeans(sapply(1:num_simu, function(num){
  (apply(RES.DRF[[1, num]], 2, function(X){
    sum(c(X[1], diff(X, lag = 1, differences = 1))*y)
  }) - mean.gold)^2
}))
MSEr.drf <- MSE.drf/MSE.drf[5]

#### Processing ENG

load(file = filename.eng)
RES.ENG <- RES
rm(RES)

MSE.eng = rowMeans(sapply(1:num_simu, function(num){
  (apply(RES.ENG[[1, num]], 2, function(X){
    sum(c(X[1], diff(X, lag = 1, differences = 1))*y)
  }) - mean.gold)^2
}))
MSEr.eng <- MSE.eng/MSE.eng[4]

#### Processing ZhangY and ZhangA

load(file = filename.mean)
RES.Mean <- RES.F
rm(RES.F)

MSE.zhang <- rowMeans((RES.Mean - mean.gold)^2)
MSEr.zhang <- MSE.zhang/MSE.zhang[6]

################################################################################

out = rbind(rep(MSEr.true[1], 4), MSEr.gamlss[c(2, 1, 3, 4)], c(MSEr.eng[c(2, 1, 3)], -1), 
      MSEr.drf[c(2, 1, 3, 4)], rep(MSEr.zhang[1], 4), MSEr.zhang[c(1, 2, 3, 4)+1])
rownames(out) <- c('BASD*', 'BASD-GAMLSS', 'BASD-ENG', 'BASD-DRF', 'ZhangA', 'ZhangY')
colnames(out) <- paste('K = ', c(2, 5, 10, 20), sep = '')
out100 = out[, 2]

################################################################################
#### p = 500
################################################################################
p = 500
p1 = ceiling(0.1*p)
beta0 = c(rep(1, p1), rep(0, p-p1))
filename.true <- paste0('RES_dist_M1_trueCCDF_n_', n, '_p_', p, '_m_', m, '_p1_', p1, 
                        '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

filename.gamlss <- paste0('RES_dist_M1_gamlss_n_', n, '_p_', p, '_m_', m, '_p1_', p1, 
                          '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

filename.drf <- paste0('RES_dist_M1_drf_n_', n, '_p_', p, '_m_', m, '_p1_', p1, 
                       '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

filename.eng <- paste0('RES_dist_M1_engression_n_', n, '_p_', p, '_m_', m, '_p1_', p1, 
                       '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

filename.mean <- paste0('RES_mean_M1_n_', n, '_p_', p, '_m_', m, '_p1_', p1,
                        '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

#### Processing BASD*

load(file = filename.true)
RES.true <- RES
rm(RES)

MSE.true = rowMeans(sapply(1:num_simu, function(num){
  (apply(RES.true[[1, num]], 2, function(X){
    sum(c(X[1], diff(X, lag = 1, differences = 1))*y)
  }) - mean.gold)^2
}))
MSEr.true <- MSE.true/MSE.true[2]

#### Processing GAMLSS

load(file = filename.gamlss)
RES.GAMLSS <- RES
rm(RES)

MSE.gamlss = rowMeans(sapply(1:num_simu, function(num){
  (apply(RES.GAMLSS[[1, num]], 2, function(X){
    sum(c(X[1], diff(X, lag = 1, differences = 1))*y)
  }) - mean.gold)^2
}))
MSEr.gamlss <- MSE.gamlss/MSE.gamlss[5]

#### Processing DRF

load(file = filename.drf)
RES.DRF <- RES
rm(RES)

MSE.drf = rowMeans(sapply(1:num_simu, function(num){
  (apply(RES.DRF[[1, num]], 2, function(X){
    sum(c(X[1], diff(X, lag = 1, differences = 1))*y)
  }) - mean.gold)^2
}))
MSEr.drf <- MSE.drf/MSE.drf[5]

#### Processing ENG

load(file = filename.eng)
RES.ENG <- RES
rm(RES)

MSE.eng = rowMeans(sapply(1:num_simu, function(num){
  (apply(RES.ENG[[1, num]], 2, function(X){
    sum(c(X[1], diff(X, lag = 1, differences = 1))*y)
  }) - mean.gold)^2
}))
MSEr.eng <- MSE.eng/MSE.eng[4]

#### Processing ZhangY and ZhangA

load(file = filename.mean)
RES.Mean <- RES.F
rm(RES.F)

MSE.zhang <- rowMeans((RES.Mean - mean.gold)^2)
MSEr.zhang <- MSE.zhang/MSE.zhang[6]

################################################################################

out = rbind(rep(MSEr.true[1], 4), MSEr.gamlss[c(2, 1, 3, 4)], c(MSEr.eng[c(2, 1, 3)], -1), 
            MSEr.drf[c(2, 1, 3, 4)], rep(MSEr.zhang[1], 4), MSEr.zhang[c(1, 2, 3, 4)+1])
rownames(out) <- c('BASD*', 'BASD-GAMLSS', 'BASD-ENG', 'BASD-DRF', 'ZhangA', 'ZhangY')
colnames(out) <- paste('K = ', c(2, 5, 10, 20), sep = '')
out500 = out[, 2]

################################################################################
#### p = 10
################################################################################
p = 10
p1 = ceiling(0.1*p)
beta0 = c(rep(1, p1), rep(0, p-p1))
filename.basd <- paste0('RES_dist_M1_method_n_', n, '_p_', p, '_m_', m, '_p1_', p1, 
                                      '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

filename.mean <- paste0('RES_mean_M1_n_', n, '_p_', p, '_m_', m, '_p1_', p1,
                   '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

load(filename.mean)
RES.Mean <- RES.F
rm(RES.F)

load(filename.basd)
RES.basd <- RES
rm(RES)

MSE.basd = rowMeans(sapply(1:num_simu, function(num){
  (apply(RES.basd[[1, num]], 2, function(X){
    sum(c(X[1], diff(X, lag = 1, differences = 1))*y)
  }) - mean.gold)^2
}))
MSE.zhang <- rowMeans((RES.Mean - mean.gold)^2)
out10 <- c((MSE.basd/MSE.basd[length(MSE.basd)])[c(1, 2, 4, 3)], (MSE.zhang/MSE.zhang[length(MSE.zhang)])[c(1, 2)])

################################################################################
#### p = 1000
################################################################################
p = 1000
p1 = ceiling(0.1*p)
beta0 = c(rep(1, p1), rep(0, p-p1))
filename.basd <- paste0('RES_dist_M1_method_n_', n, '_p_', p, '_m_', m, '_p1_', p1, 
                        '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

filename.mean <- paste0('RES_mean_M1_n_', n, '_p_', p, '_m_', m, '_p1_', p1,
                        '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

load(filename.mean)
RES.Mean <- RES.F
rm(RES.F)

load(filename.basd)
RES.basd <- RES
rm(RES)

MSE.basd = rowMeans(sapply(1:num_simu, function(num){
  (apply(RES.basd[[1, num]], 2, function(X){
    sum(c(X[1], diff(X, lag = 1, differences = 1))*y)
  }) - mean.gold)^2
}))
MSE.zhang <- rowMeans((RES.Mean - mean.gold)^2)
out1000 <- c((MSE.basd/MSE.basd[length(MSE.basd)])[c(1, 2, 4, 3)], -1, (MSE.zhang/MSE.zhang[length(MSE.zhang)])[1])


################################################################################
#### Pooling
################################################################################

OUT = rbind(out10, out100, out500, out1000)
paste('$p =', c(10, 100, 500, 1000), '$ && ', apply(OUT, 1, function(X){
  paste(format(round(X, 2)), collapse = ' & ')
}), '\\ ', collapse = '', sep = '')
