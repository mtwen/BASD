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

pMSE <- data.frame(location = rep(y, nmethod+1), 
                   value = c(as.vector(MSEr), LB),
                   method = rep(c('BASD*', 'BASD-GAMLSS', 'BASD-DRF', 'BASD-Engression', 'ECDF', 'LB'), each = length(y)))
pMSE$method <- factor(pMSE$method, levels = c('BASD*', 'BASD-GAMLSS', 'BASD-DRF', 'BASD-Engression', 'ECDF', 'LB'))

require(ggplot2)
require(extrafont)
man_linetype = c(1, 6, 2, 4, 5, 5)
man_color = c(RColorBrewer::brewer.pal(8, 'Set1')[c(3, 1, 2, 4)], "#878787", "#4D4D4D")
plt <- ggplot(data = pMSE, aes(x = location, y = value, color = method, linetype = method)) +
  geom_line(linewidth = 0.8) + 
  ylim(0.3, 1.1) + ylab(expression('pMSE'(y)/'pMSE'[0](y))) + xlab(expression(y)) +
  theme_bw() +
  ggtitle('(b)') +
  theme(legend.position = 'none',
        plot.title.position = 'plot') + 
  scale_color_manual(values = man_color) + 
  scale_linetype_manual(values = man_linetype) 
plt


setEPS()
postscript(file = paste0('FIG_dist_p', p, '.eps'), horizontal = FALSE, width = 5, height = 3.5)
plt
dev.off()

################################################################################
#### Numbers
################################################################################

mser.true = colMeans(MSE.true)/mean(MSE.true[, 2])
mser.gamlss = colMeans(MSE.gamlss)/mean(MSE.gamlss[, 5])
mser.drf = colMeans(MSE.drf)/mean(MSE.drf[, 5])
mser.eng = colMeans(MSE.eng)/mean(MSE.eng[, 4])

cover.true = 1 - rowMeans(ceiling(sapply(1:num_simu, function(num){
  colMeans(matrix(abs(as.vector(RES.true[[1, num]]) - taus) > rep(RES.true[[2, num]]/sqrt(n), each = length(y)), length(y), 2))
})))
len.true = rowMeans(sapply(1:num_simu, function(num){
  RES.true[[2, num]]/sqrt(n)*2
}))

cover.gamlss = 1 - rowMeans(ceiling(sapply(1:num_simu, function(num){
  colMeans(matrix(abs(as.vector(RES.GAMLSS[[1, num]]) - taus) > rep(RES.GAMLSS[[2, num]]/sqrt(n), each = length(y)), length(y), 5))
})))
len.gamlss = rowMeans(sapply(1:num_simu, function(num){
  RES.GAMLSS[[2, num]]/sqrt(n)*2
}))

cover.drf = 1 - rowMeans(ceiling(sapply(1:num_simu, function(num){
  colMeans(matrix(abs(as.vector(RES.DRF[[1, num]]) - taus) > rep(RES.DRF[[2, num]]/sqrt(n), each = length(y)), length(y), 5))
})))
len.drf = rowMeans(sapply(1:num_simu, function(num){
  RES.DRF[[2, num]]/sqrt(n)*2
}))

cover.eng = 1 - rowMeans(ceiling(sapply(1:num_simu, function(num){
  colMeans(matrix(abs(as.vector(RES.ENG[[1, num]]) - taus) > rep(RES.ENG[[2, num]]/sqrt(n), each = length(y)), length(y), 4))
})))
len.eng = rowMeans(sapply(1:num_simu, function(num){
  RES.ENG[[2, num]]/sqrt(n)*2
}))


out = data.frame(MSEr = c(mser.true[1], mser.gamlss[cv.idx], mser.eng[cv.idx], mser.drf[cv.idx], 1),
           Cover = c(cover.true[1], cover.gamlss[cv.idx], cover.eng[cv.idx], cover.drf[cv.idx], mean(c(cover.true[2], cover.gamlss[5], cover.eng[4], cover.drf[5]))), 
           Lens = c(len.true[1], len.gamlss[cv.idx], len.eng[cv.idx], len.drf[cv.idx], mean(c(len.true[2], len.gamlss[5], len.eng[4], len.drf[5])))*100)

apply(out, 1, function(X){
  paste(format(round(X, 2)), collapse = ' & ')
})


## Number of cross-fitting folds

MSEr.cv = rbind(rep(mser.true[1], 4), mser.gamlss[c(2, 1, 3, 4)], mser.eng[c(2, 1, 3, 4)], mser.drf[c(2, 1, 3, 4)], rep(1, 4))
rownames(MSEr.cv) <- c('BASD*', 'BASD-GAMLSS', 'BASD-Engression', 'BASD-DRF', 'ECDF')
colnames(MSEr.cv) <- paste('K = ', c(2, 5, 10, 20), sep = '')
MSEr.cv['BASD-Engression', 4] = -1
MSEr.cv

Cov.cv = rbind(rep(cover.true[1], 4), cover.gamlss[c(2, 1, 3, 4)], c(cover.eng[c(2, 1, 3)], -1), cover.drf[c(2, 1, 3, 4)],
               rep(mean(c(cover.true[2], cover.gamlss[5], cover.eng[4], cover.drf[5])), 4))
rownames(Cov.cv) <- c('BASD*', 'BASD-GAMLSS', 'BASD-Engression', 'BASD-DRF', 'ECDF')
colnames(Cov.cv) <- paste('K = ', c(2, 5, 10, 20), sep = '')
Cov.cv

Len.cv = rbind(rep(len.true[1], 4), len.gamlss[c(2, 1, 3, 4)], c(len.eng[c(2, 1, 3)], -0.01), len.drf[c(2, 1, 3, 4)],
               rep(mean(c(len.true[2], len.gamlss[5], len.eng[4], len.drf[5])), 4))*100
rownames(Len.cv) <- c('BASD*', 'BASD-GAMLSS', 'BASD-Engression', 'BASD-DRF', 'ECDF')
colnames(Len.cv) <- paste('K = ', c(2, 5, 10, 20), sep = '')
Len.cv

