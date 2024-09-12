rm(list = ls()); gc()
set.seed(20240722)
library(doSNOW)
source('basis.R', local = TRUE)

################################################################################
#### Functions
################################################################################

sampleQuant <- function(Y, tau){ # Y: a vector; tau: the tau-th quantile (0, 1]
  n = length(Y)
  sy = sort(Y, decreasing = FALSE)
  return(sy[ceiling(n*tau)])
}

zhangYsemi <- function(Y, X, X.un, K = 10){
  n = nrow(X)
  m = nrow(X.un)
  
  folds_n = lapply(1:K, function(k){ seq(k, n, by = K) })
  folds_m = lapply(1:K, function(k){ seq(k, m, by = K) })
  P = vector(length = n+m)
  for (k in 1:K){
    idx_n = folds_n[[k]]
    idx_m = folds_m[[k]]
    
    fit.rf = ranger::ranger(tY~., data = data.frame(tY = Y[-idx_n], tX = X[-idx_n, ]), write.forest = TRUE)
    P[idx_n] = predict(fit.rf, data.frame(tX = X[idx_n, ]))$predictions
    P[n+idx_m] = predict(fit.rf, data.frame(tX = X.un[idx_m, ]))$predictions
  }
  P
}


################################################################################
#### Data generation
################################################################################
n = 1000
p = 1000
m = 10*n
p1 = ceiling(0.1*p)
num_simu = 200
ticks = 100 - 1
beta0 = c(rep(1, p1), rep(0, p-p1))
taus = seq(1/(ticks + 1), 1 - 1/(ticks + 1), by = 1/(ticks + 1))
alpha = 0.1
y = qnorm(taus, 0, sd = sqrt(1+1/4)) 
filename <- paste0('RES_mean_M1_n_', n, '_p_', p, '_m_', m, '_p1_', p1,
                   '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

################################################################################
#### Methods Evaluation
################################################################################

progress <- function(nfin){
  cat(sprintf('%s: tasks completed: %d.\n', Sys.time(), nfin))
}
opts <- list(progress = progress)

cat(sprintf('%s: tasks begin...\n', Sys.time()))
cl <- makeSOCKcluster(5)
registerDoSNOW(cl)
RES.F = foreach (num=1:num_simu, .packages = c('mgcv'), .combine = 'cbind', 
                 .multicombine = TRUE, .options.snow = opts) %dopar% 
  {
    set.seed(20240723 + num*1000)
    X = matrix(rnorm(n*p), n, p)
    Y = X %*% beta0 / sqrt(p1) + rnorm(n)/2
    Xa = matrix(rnorm(m*p), m, p)
    
    mean.sample = mean(Y)
    
    # ZhangA 2019 for p < n
    beta.A <- lm(Y~X+1)$coefficients
    mean.zhA = mean(cbind(1, rbind(X, Xa)) %*% beta.A)
    
    # ZhangY 2021
    # pred.2 = zhangYsemi(Y, X, Xa, K = 2)
    pred.5 = zhangYsemi(Y, X, Xa, K = 5)
    # pred.10 = zhangYsemi(Y, X, Xa, K = 10)
    # pred.20 = zhangYsemi(Y, X, Xa, K = 20)
    # mean.zhY.2 = mean(pred.2) + mean(Y) - mean(pred.2[1:n])
    mean.zhY.5 = mean(pred.5) + mean(Y) - mean(pred.5[1:n])
    # mean.zhY.10 = mean(pred.10) + mean(Y) - mean(pred.10[1:n])
    # mean.zhY.20 = mean(pred.20) + mean(Y) - mean(pred.20[1:n])
    
    # MSet = c(mean.zhA, mean.zhY.2, mean.zhY.5, mean.zhY.10, mean.zhY.20, mean.sample) # p = 100, 500
    MSet = c(mean.zhA, mean.zhY.5, mean.sample) # p = 10
    # MSet = c(mean.zhY.5, mean.sample) # p = 1000
    MSet
  }
stopCluster(cl)
save(RES.F, file = filename)
cat(sprintf('%s: tasks end and results saved.\n', Sys.time()))

