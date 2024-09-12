## MengtaoWEN
## 2024-07-17
## Executed file for different methods

rm(list = ls()); gc()
set.seed(20240722)
source('basis.R')
library(doSNOW)

n = 1000
pSet = c(10, 1000)
m = 10*n
num_simu = 200
ticks = 100 - 1
taus = seq(1/(ticks + 1), 1 - 1/(ticks + 1), by = 1/(ticks + 1))
alpha = 0.1
y = qnorm(taus, 0, sd = sqrt(1+1/4)) 

total_tick = 1000000
Z = qnorm(seq(1/total_tick, 1 - 1/total_tick, by = 1/total_tick), 0, 1)
gamman = m/(n+m)
var.cCDF <- sapply(y, function(loc){
  mean((pnorm(2*loc - 2*Z) - mean(pnorm(2*loc - 2*Z)))^2)
})
var.y <- taus*(1 - taus)
LB = 1 - gamman * var.cCDF/var.y

for (pp in 1:length(pSet)){
  p = pSet[pp]
  p1 = ceiling(0.1*p)
  beta0 = c(rep(1, p1), rep(0, p-p1))
  filename <- paste0('RES_dist_M1_method_n_', n, '_p_', p, '_m_', m, '_p1_', p1, 
                     '_ticks_', ticks, '_alpha_', alpha, '_numSimu_', num_simu, '.RData')
  
  load(file = filename)
  nmethod = 5
  MSE = matrix(rowMeans(sapply(1:num_simu, function(num){
    as.vector(RES[[1, num]] - taus)^2
  })), ticks, nmethod)
  MSEr = MSE/MSE[, nmethod]
  oMSEr = colMeans(MSE)/mean(MSE[, nmethod])
  cover = 1 - rowMeans(ceiling(sapply(1:num_simu, function(num){
    colMeans(matrix(abs(as.vector(RES[[1, num]]) - taus) > rep(RES[[2, num]]/sqrt(n), each = length(y)), length(y), nmethod))
  })))
  len = rowMeans(sapply(1:num_simu, function(num){
    RES[[2, num]]/sqrt(n)*2
  }))
  
  eval(parse(text = paste0('pMSE', pp, ' <- data.frame(location = rep(y, nmethod+1), 
                     value = c(as.vector(MSEr), LB),
                     method = rep(c(\'BASD*\', \'BASD-GAMLSS\', \'BASD-DRF\', \'BASD-Engression\', \'ECDF\', \'LB\'), each = length(y)))', collapse = '')))
  eval(parse(text = paste0('oMSEr', pp, ' <- oMSEr', collapse = '')))
  eval(parse(text = paste0('cover', pp, ' <- cover', collapse = '')))
  eval(parse(text = paste0('len', pp, ' <- len', collapse = '')))
}

################################################################################
#### Figures
################################################################################

pMSE = rbind(pMSE1, pMSE2)
pMSE$method <- factor(pMSE$method, levels = c('BASD*', 'BASD-GAMLSS', 'BASD-DRF', 'BASD-Engression', 'ECDF', 'LB'))
pMSE$p <- rep(paste('p =', pSet), each = nrow(pMSE1))


require(ggplot2)
require(extrafont)
man_linetype = c(1, 2, 4, 6, 5, 5)
man_color = c(RColorBrewer::brewer.pal(8, 'Set1')[c(3, 1, 2, 4)], "#878787", "#4D4D4D")
plt <- ggplot(data = pMSE, aes(x = location, y = value, color = method, linetype = method)) +
  geom_line(linewidth = 0.8) + 
  facet_wrap(vars(p)) + 
  ylim(0.3, 1.2) + ylab(expression('pMSE'(y)/'pMSE'[0](y))) + xlab(expression(y)) +
  theme_bw() +
  theme(legend.position = 'none',
        plot.title.position = 'plot') + 
  scale_color_manual(values = man_color) + 
  scale_linetype_manual(values = man_linetype) 
plt


setEPS()
postscript(file = paste0('FIG_dist_M1_supp.eps'), horizontal = FALSE, width = 8, height = 3)
plt
dev.off()

################################################################################
#### Numbers
################################################################################

out = cbind(oMSEr1, cover1, len1*100, oMSEr2, cover2, len2*100)
methodnames <- c('BASD*', 'BASD-GAMLSS', 'BASD-DRF', 'BASD-Engression', 'ECDF')
paste(methodnames, ' && ', apply(out, 1, function(X){
  paste(apply(matrix(format(round(X, 2)), nrow = 3), 2, function(Y){
    paste(Y, collapse = ' & ')
  }), collapse = ' && ')
}), '\\', collapse = ' ')

