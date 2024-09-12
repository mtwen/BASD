## MengtaoWEN
## 2024-07-17
## Executed file for different methods

rm(list = ls()); gc()
set.seed(20240722)
source('basis.R')
library(doSNOW)

n = 1000
pSet = c(100, 500)
m = 10*n
num_simu = 200
ticks = 100 - 1
taus = seq(1/(ticks + 1), 1 - 1/(ticks + 1), by = 1/(ticks + 1))
alpha = 0.1
y = qnorm(taus, 0, sd = sqrt(1+1/4))

for (pp in 1:length(pSet)){
  p = pSet[pp]
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
  
  
  total_tick = 1000000
  Z = qnorm(seq(1/total_tick, 1 - 1/total_tick, by = 1/total_tick), 0, 1)
  gamman = m/(n+m)
  var.cCDF <- sapply(y, function(loc){
    mean((pnorm(2*loc - 2*Z) - mean(pnorm(2*loc - 2*Z)))^2)
  })
  var.y <- taus*(1 - taus)
  LB = 1 - gamman * var.cCDF/var.y

  ################################################################################
  #### Processing BASD*
  ################################################################################
  
  load(file = filename.true)
  RES.true <- RES
  rm(RES)
  
  MSE.true = matrix(rowMeans(sapply(1:num_simu, function(num){
    as.vector(RES.true[[1, num]] - taus)^2
  })), ticks, 2)
  MSEr.true <- MSE.true/MSE.true[, 2]
  
  ################################################################################
  #### Processing GAMLSS
  ################################################################################
  
  load(file = filename.gamlss)
  RES.GAMLSS <- RES
  rm(RES)
  
  MSE.gamlss = matrix(rowMeans(sapply(1:num_simu, function(num){
    as.vector(RES.GAMLSS[[1, num]] - taus)^2
  })), ticks, 5)
  MSEr.gamlss <- MSE.gamlss/MSE.gamlss[, 5]
  
  
  ################################################################################
  #### Processing DRF
  ################################################################################
  
  load(file = filename.drf)
  RES.DRF <- RES
  rm(RES)
  
  MSE.drf = matrix(rowMeans(sapply(1:num_simu, function(num){
    as.vector(RES.DRF[[1, num]] - taus)^2
  })), ticks, 5)
  MSEr.drf <- MSE.drf/MSE.drf[, 5]
  
  ################################################################################
  #### Processing ENG
  ################################################################################
  
  load(file = filename.eng)
  RES.ENG <- RES
  rm(RES)
  
  MSE.eng = matrix(rowMeans(sapply(1:num_simu, function(num){
    as.vector(RES.ENG[[1, num]] - taus)^2
  })), ticks, 4)
  MSEr.eng <- MSE.eng/MSE.eng[, 4]
  
  ################################################################################
  #### Pooling all results
  ################################################################################
  cv.idx = 3
  nmethod = 5
  MSEr = cbind(MSEr.true[, 1], MSEr.gamlss[, cv.idx], MSEr.drf[, cv.idx], MSEr.eng[, cv.idx], 1)
  
  eval(parse(text = paste0('pMSE', pp, ' <- data.frame(location = rep(y, nmethod+1), 
                     value = c(as.vector(MSEr), LB),
                     method = rep(c(\'BASD*\', \'BASD-GAMLSS\', \'BASD-DRF\', \'BASD-Engression\', \'ECDF\', \'LB\'), each = length(y)))', collapse = '')))
}

pMSE = rbind(pMSE1, pMSE2)
pMSE$method <- factor(pMSE$method, levels = c('BASD*', 'BASD-GAMLSS', 'BASD-DRF', 'BASD-Engression', 'ECDF', 'LB'))
pMSE$p <- factor(rep(paste('p =', pSet), each = nrow(pMSE1)), levels = paste('p =', pSet))

require(ggplot2)
require(extrafont)
man_linetype = c(1, 2, 4, 6, 5, 5)
man_color = c(RColorBrewer::brewer.pal(8, 'Set1')[c(3, 1, 2, 4)], "#878787", "#4D4D4D")
plt <- ggplot(data = pMSE, aes(x = location, y = value, color = method, linetype = method)) +
  geom_line(linewidth = 0.9) + 
  facet_wrap(vars(p)) + 
  ylim(0.3, 1.1) + ylab(expression('pMSE'(y)/'pMSE'[0](y))) + xlab(expression(y)) +
  theme_bw() +
  theme(legend.position = 'none',
        plot.title.position = 'plot') + 
  scale_color_manual(values = man_color) + 
  scale_linetype_manual(values = man_linetype) 
plt

setEPS()
postscript(file = paste0('FIG_dist_M1.eps'), horizontal = FALSE, width = 8, height = 3)
plt
dev.off()


