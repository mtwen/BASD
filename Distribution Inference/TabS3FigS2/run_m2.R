## MengtaoWEN
## 2024-07-28
## Executed file for different methods

rm(list = ls()); gc()
set.seed(20240730)
source('basis.R')
library(doSNOW)

etaFun <- function(X){
  3*cos(sum(X[1:3]))
}

sigFun <- function(X){
  exp(sum(X[4:7]/2))/3
}

n = 1000
p = 50
# mSet = c(1*n, 5*n)
mSet = c(0.5*n, 1*n, 5*n, 10*n)
mMax = max(mSet)
num_simu = 200
ticks = 100 - 1
taus = seq(1/(ticks + 1), 1 - 1/(ticks + 1), by = 1/(ticks + 1))
alpha = 0.1
y = qunif(taus, min = -3.5, max = 4) 
filename <- paste0('RES_dist_M2_n_', n, '_p_', p, '_mMax_', mMax, 
                   '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

progress <- function(nfin){
  cat(sprintf('%s: tasks completed: %d.\n', Sys.time(), nfin))
}
opts <- list(progress = progress)

cat(sprintf('%s: tasks begin...\n', Sys.time()))
cl <- makeSOCKcluster(5)
registerDoSNOW(cl)
RES = foreach (num=1:num_simu, .packages = c('mgcv'), .combine = 'cbind', 
               .multicombine = TRUE, .options.snow = opts) %dopar% 
  {
    set.seed(20240730 + num*1000)
    X = matrix(rnorm(n*p), n, p)
    Y = apply(X, 1, etaFun) + apply(X, 1, sigFun) * rnorm(n)
    XaM = matrix(rnorm(mMax*p), mMax, p)
    
    ecdf.BASD.true = vector('list', length = length(mSet))    
    ecdf.BASD.LSS.5 = vector('list', length = length(mSet))
    ecdf.BASD.DRF.5 = vector('list', length = length(mSet))
    ecdf.BASD.ENG.5 = vector('list', length = length(mSet))
    for (mm in 1:length(mSet)){
      m = mSet[mm]
      Xa = XaM[1:m, ]
      
      P.true = matrix(pnorm( (rep(y, each = n+m) - rep(apply(rbind(X, Xa), 1, etaFun), length(y)))/rep(apply(rbind(X, Xa), 1, sigFun), length(y)) ), n+m, length(y))
      P.LSS.5 = CondCDF(y, Y, X, Xa, K = 5, method = 'gamlss', control = mboost::boost_control(mstop = 170))
      P.DRF.5 = CondCDF(y, Y, X, Xa, K = 5, method = 'drf', control = list(num.trees = 10000, num.features = 7))
      P.ENG.5 = CondCDF(y, Y, X, Xa, K = 5, method = 'engression', control = list(hidden_dim = 50, num_layer = 2, num_epochs = 500))
      ecdf.BASD.true[[mm]] = BASD(y, Y, P.true, inference = TRUE, alpha = alpha)
      ecdf.BASD.LSS.5[[mm]] = BASD(y, Y, P.LSS.5, inference = TRUE, alpha = alpha)
      ecdf.BASD.DRF.5[[mm]] = BASD(y, Y, P.DRF.5, inference = TRUE, alpha = alpha)
      ecdf.BASD.ENG.5[[mm]] = BASD(y, Y, P.ENG.5, inference = TRUE, alpha = alpha)
    }
    ecdf.classic = empiricalCDF(y, Y, inference = TRUE, alpha = alpha)
    
    out.ecdf = sapply(c(ecdf.BASD.true, ecdf.BASD.LSS.5, ecdf.BASD.DRF.5, ecdf.BASD.ENG.5,
                        list(ecdf.classic)), function(ecdf.BASD){
                          ecdf.BASD$estCDF
                        })
    out.L = sapply(c(ecdf.BASD.true, ecdf.BASD.LSS.5, ecdf.BASD.DRF.5, ecdf.BASD.ENG.5,
                     list(ecdf.classic)), function(ecdf.BASD){
                       ecdf.BASD$L
                     })
    list(ECDF = out.ecdf, L = out.L)
  }
stopCluster(cl)
save(RES, file = filename)
cat(sprintf('%s: tasks end and results saved.\n', Sys.time()))

