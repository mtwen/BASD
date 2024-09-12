#### Mengtao WEN
#### 2023-09-21

rm(list = ls()); gc()
library(doSNOW)

n = 1000
p = 5
sigma = 1
rho = 0.5
deltaSet = seq(0, 0.6, 0.2)
# delta = 0
B = 2000
alpha = 0.1

K = 10 # number of folds for cross-fitting
num_simu = 1000
filename = paste('SSF_LDTE_n_', n, '_p_', p, '_rho_', rho, '_sigma_', sigma, 
                 '_Sigma_I_K_', K, '_B_', B, '_alpha_', alpha, '_sim_', num_simu, '.RData', sep = '')

load(filename)

mthdSet = c('Ab', 'S3')
res = sapply(1:num_simu, function(num){
  out = array(NA, c(2, length(deltaSet), length(mthdSet)))
  for(dd in 1:length(deltaSet)){
    out[1, dd, ] = RES.LDTE[[dd, num]]$rj.eq
    out[2, dd, ] = RES.LDTE[[dd, num]]$rj.fsd
  }
  as.vector(out)
})
ave = array(rowMeans(res), c(2, length(deltaSet), length(mthdSet)))
ave

paste(sapply(1:length(mthdSet), function(mm){
  mthd = mthdSet[mm]
  paste(mthd, '&&', paste(sapply(1:2, function(ii){paste(format(ave[ii, , mm]), collapse = ' & ')}), collapse = ' && '), '\\')
}), collapse = ' ')

rowMeans(sapply(1:num_simu, function(num){
  RES.LDTE[[3, num]]$rj.fsd
}))
