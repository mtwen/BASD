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

X.teO = matrix(c(3, 1), 1, p)
filename = paste('res_cpv_split_x_', X.teO[1, 1], '_y_', X.teO[1, 2], '_num_', num_simu, '.RData', sep = '')
load(file = filename)
pval.gd = (pval.gd * (n+1) - 1)/n
pnorm((cSet - muFun(X.teO))/sigmaFun(X.teO)/sigma)

JC = rowMeans((sapply(1:num_simu, function(num){
  res[[1, num]]
}) - pval.gd)^2)

# ORA = matrix(rowMeans((sapply(1:num_simu, function(num){
#   as.vector(res[[2, num]])
# }) - pval.gd)^2), length(cSet), length(mSet))

EST = matrix(rowMeans((sapply(1:num_simu, function(num){
  as.vector(res[[3, num]])
}) - pval.gd)^2), length(cSet), length(mSet))

pval1 = pval.gd
out1 = EST/JC
rm(pval.gd, res)

X.teO = matrix(c(4, 4), 1, p)
filename = paste('res_cpv_split_x_', X.teO[1, 1], '_y_', X.teO[1, 2], '_num_', num_simu, '.RData', sep = '')
load(file = filename)
pval.gd = (pval.gd * (n+1) - 1)/n
pnorm((cSet - muFun(X.teO))/sigmaFun(X.teO)/sigma)

JC = rowMeans((sapply(1:num_simu, function(num){
  res[[1, num]]
}) - pval.gd)^2)

# ORA = matrix(rowMeans((sapply(1:num_simu, function(num){
#   as.vector(res[[2, num]])
# }) - pval.gd)^2), length(cSet), length(mSet))

EST = matrix(rowMeans((sapply(1:num_simu, function(num){
  as.vector(res[[3, num]])
}) - pval.gd)^2), length(cSet), length(mSet))

pval2 = pval.gd
out2 = EST/JC
rm(pval.gd, res)

scSet = c(3:6)
smSet = c(2:4)
paste(sapply(1:length(scSet), function(c){
  paste0(cSet[scSet[c]], ' && ', format(round(pval1[scSet[c]], 3)), ' & ', paste(format(round(out1[scSet[c], smSet], 3)), collapse = ' & '), 
         ' && ', format(round(pval2[scSet[c]], 3)), ' & ', paste(format(round(out2[scSet[c], smSet], 3)), collapse = ' & '), '\\')
}), collapse = ' ')

paste0(paste(format(round(pval1, 3)), collapse = ' & '), ' && ', paste(format(round(pval2, 3)), collapse = ' & '))
