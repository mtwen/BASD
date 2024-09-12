rm(list = ls()); gc()
set.seed(1234)
library(doSNOW)
source("funcShare.R")

## Settings
n = 100
mSet = c(100, 500, 1000, 5000)
M = max(mSet)
p = 2
sigma = 2
maxX = 5
minX = 0
num_simu = 1000
cSet = c(0, 5, 10, 15, 20, 30)
mthds = c('JC', 'BS.ora', 'BS.dim')

X.teO = matrix(c(4, 4), 1, p)
# X.teO = matrix(c(3, 1), 1, p)
filename = paste('res_cpv_split_x_', X.teO[1, 1], '_y_', X.teO[1, 2], '_num_', num_simu, '.RData', sep = '')

X.tr = matrix(runif(n*p, minX, maxX), n, p)
Y.tr = muFun(X.tr)  + sigma*sigmaFun(X.tr)*rnorm(n, 0, 1)
# fit = nnet::nnet(tY~., data=data.frame(tX = X.tr, tY=Y.tr), size = 100, linout=TRUE, maxit=2000, trace=FALSE)
# fit = lm(tY~., data = data.frame(tX = X.tr, tY = Y.tr))
fit = randomForest::randomForest(tY~., data=data.frame(tX = X.tr, tY=Y.tr))

mu.hatO = as.vector(predict(fit, data.frame(tX = X.teO)))  
Vs.teO = cSet - mu.hatO

N = 2000000
X.ad =  matrix(runif(N*p, minX, maxX), N, p)
Y.ad = muFun(X.ad) + sigma*sigmaFun(X.ad)*rnorm(N, 0, 1)
mu.hatN = as.vector(predict(fit, data.frame(tX = X.ad)))
V.ad = Y.ad - mu.hatN
pval.gd = (1 + n * colMeans(matrix(rep(V.ad, length(cSet)) <= rep(Vs.teO, each = N), N, length(cSet))))/(n+1)

progress <- function(nfin){
  cat(sprintf('%s: tasks completed: %d.\n', Sys.time(), nfin))
}
opts <- list(progress = progress)

cat(sprintf('%s: tasks begin...\n', Sys.time()))
cl <- makeSOCKcluster(5)
registerDoSNOW(cl)
res = foreach (num=1:num_simu, .packages = c('nnet', 'randomForest'), .combine = 'cbind', .multicombine = TRUE, .options.snow = opts) %dopar% {
  
  set.seed(20231124 + num*1000)
  
  X.ca = matrix(runif(n*p, minX, maxX), n, p)
  Y.ca = muFun(X.ca)  + sigma*sigmaFun(X.ca)*rnorm(n, 0, 1)
  mu.hatn = as.vector(predict(fit, data.frame(tX = X.ca)))
  V.ca = Y.ca - mu.hatn
  
  pval.JC = (colSums(matrix(rep(V.ca, length(cSet)) <= rep(Vs.teO, each = n), n, length(cSet))) + 1)/(n+1)
  
  X.unM = matrix(runif(M*p, minX, maxX), M, p)
  mu.hatM = as.vector(predict(fit, data.frame(tX = X.unM)))
  pval.BS.ora = matrix(NA, length(cSet), length(mSet))
  pval.BS.est = matrix(NA, length(cSet), length(mSet))
  for (mm in 1:length(mSet)){
    m = mSet[mm]
    X.un = X.unM[1:m, ]
    mu.hatm = mu.hatM[1:m]
    P.ora = matrix(pnorm((rep(Vs.teO, each = n+m) + rep(c(mu.hatn, mu.hatm) - muFun(rbind(X.ca, X.un)), length(cSet)))/sigma/sigmaFun(rbind(X.ca, X.un)), 0, 1), n+m, length(cSet))
    pval.BS.ora[, mm] = (1 + n * BAS(Vs.teO, V.ca, P.ora))/(n+1)
    P.est = ConCDF(Vs.teO, V.ca, X.ca, X.un, K = 10, method = 'DRF')
    pval.BS.est[, mm] = (1 + n * BAS(Vs.teO, V.ca, P.est))/(n+1)  
  }
  list(JC = pval.JC, BS.ora = pval.BS.ora, BS.est = pval.BS.est)
}
stopCluster(cl)
# save(res, pval.gd, file = filename)
cat(sprintf('%s: tasks end and results saved.\n', Sys.time()))

load(file = filename)
JC = rowMeans((sapply(1:num_simu, function(num){
  res[[1, num]]
}) - pval.gd)^2)

ORA = matrix(rowMeans((sapply(1:num_simu, function(num){
  as.vector(res[[2, num]])
}) - pval.gd)^2), length(cSet), length(mSet))

EST = matrix(rowMeans((sapply(1:num_simu, function(num){
  as.vector(res[[3, num]])
}) - pval.gd)^2), length(cSet), length(mSet))

pval.gd
ORA/JC
EST/JC
