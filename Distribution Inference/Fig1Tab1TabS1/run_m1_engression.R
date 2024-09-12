## MengtaoWEN
## 2024-07-17
## Executed file for different methods

rm(list = ls()); gc()
set.seed(20240722)
source('basis.R')
library(doSNOW)

n = 1000
p = 500
m = 10*n
p1 = ceiling(0.1*p)
num_simu = 200
ticks = 100 - 1
beta0 = c(rep(1, p1), rep(0, p-p1))
taus = seq(1/(ticks + 1), 1 - 1/(ticks + 1), by = 1/(ticks + 1))
alpha = 0.1
y = qnorm(taus, 0, sd = sqrt(1+1/4)) 
filename <- paste0('RES_dist_M1_engression_n_', n, '_p_', p, '_m_', m, '_p1_', p1, 
                   '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')

# hidden_dim = 50, num_layer = 2, num_epochs = 5000; p = 100
# hidden_dim = 50, num_layer = 3, num_epochs = 5000; p = 500
# hidden_dim = 50, num_layer = 2, num_epochs = 5000; p = 1000

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
    set.seed(20240723 + num*1000)
    X = matrix(rnorm(n*p), n, p)
    Y = X %*% beta0 / sqrt(p1) + rnorm(n)/2
    Xa = matrix(rnorm(m*p), m, p)
    
    # quantile(Y, probs = taus, names = FALSE)
    P.ENG.5 = CondCDF(y, Y, X, Xa, K = 5, method = 'engression', control = list(hidden_dim = 50, num_layer = 3, num_epochs = 5000))
    P.ENG.2 = CondCDF(y, Y, X, Xa, K = 2, method = 'engression', control = list(hidden_dim = 50, num_layer = 3, num_epochs = 5000))
    P.ENG.10 = CondCDF(y, Y, X, Xa, K = 10, method = 'engression', control = list(hidden_dim = 50, num_layer = 3, num_epochs = 5000))
    # P.ENG.20 = CondCDF(y, Y, X, Xa, K = 20, method = 'engression', control = list(hidden_dim = 50, num_layer = 2, num_epochs = 5000))
    ecdf.BASD.ENG.5 = BASD(y, Y, P.ENG.5, inference = TRUE, alpha = alpha)
    ecdf.BASD.ENG.2 = BASD(y, Y, P.ENG.2, inference = TRUE, alpha = alpha)
    ecdf.BASD.ENG.10 = BASD(y, Y, P.ENG.10, inference = TRUE, alpha = alpha)
    # ecdf.BASD.ENG.20 = BASD(y, Y, P.ENG.20, inference = TRUE, alpha = alpha)
    ecdf.classic = empiricalCDF(y, Y, inference = TRUE, alpha = alpha)
    
    out.ecdf = sapply(list(ecdf.BASD.ENG.5, ecdf.BASD.ENG.2, ecdf.BASD.ENG.10, ecdf.classic), 
                      function(ecdf.BASD){
                        ecdf.BASD$estCDF
                      })
    out.L = sapply(list(ecdf.BASD.ENG.5, ecdf.BASD.ENG.2, ecdf.BASD.ENG.10, ecdf.classic), 
                   function(ecdf.BASD){
                     ecdf.BASD$L
                   })
    list(ECDF = out.ecdf, L = out.L)
  }
stopCluster(cl)
save(RES, file = filename)
cat(sprintf('%s: tasks end and results saved.\n', Sys.time()))
