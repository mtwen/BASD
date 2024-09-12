
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

load(file = filename)

Z = matrix(rnorm(1000000*7), 1000000, 7)
etaZ = apply(Z, 1, etaFun)
sigZ = apply(Z, 1, sigFun)
rm(Z)
Fy.true = sapply(y, function(yy){ mean(pnorm((yy - etaZ)/sigZ)) })

var.y = Fy.true * (1 - Fy.true)
var.cCDF = sapply(y, function(yy){ var(pnorm( (yy - etaZ)/sigZ )) })
LB = matrix(NA, ticks, length(mSet))
for (mm in 1:length(mSet)){
  gamman = mSet[mm]/(n+mSet[mm])
  LB[, mm] = 1 - gamman * var.cCDF/var.y
}


nmethod = 4
nrecord = length(mSet) * nmethod + 1
MSE = matrix(rowMeans(sapply(1:num_simu, function(num){
  as.vector((RES[[1, num]]  - Fy.true)^2)
})), ticks, nrecord)
MSEr = MSE/MSE[, nrecord]

################################################################################
#### FIGURE plot
################################################################################

pMSE1 <- data.frame(location = rep(y, nmethod+2), 
                   value = c(as.vector(MSEr[, c(0:(nmethod-1))*length(mSet) + 1]), MSEr[, nrecord], LB[, 1]),
                   method = rep(c('BASD*', 'BASD-GAMLSS', 'BASD-DRF', 'BASD-Engression',
                                  'ECDF', 'LB'), each = length(y)))
pMSE1$method <- factor(pMSE1$method, levels = c('BASD*', 'BASD-GAMLSS', 'BASD-DRF', 'BASD-Engression',
                                                'ECDF', 'LB'))


pMSE2 <- data.frame(location = rep(y, nmethod+2), 
                    value = c(as.vector(MSEr[, c(0:(nmethod-1))*length(mSet) + 2]), MSEr[, nrecord], LB[, 2]),
                    method = rep(c('BASD*', 'BASD-GAMLSS', 'BASD-DRF', 'BASD-Engression',
                                   'ECDF', 'LB'), each = length(y)))
pMSE2$method <- factor(pMSE2$method, levels = c('BASD*', 'BASD-GAMLSS', 'BASD-DRF', 'BASD-Engression',
                                                'ECDF', 'LB'))

pMSE3 <- data.frame(location = rep(y, nmethod+2), 
                    value = c(as.vector(MSEr[, c(0:(nmethod-1))*length(mSet) + 3]), MSEr[, nrecord], LB[, 3]),
                    method = rep(c('BASD*', 'BASD-GAMLSS', 'BASD-DRF', 'BASD-Engression',
                                   'ECDF', 'LB'), each = length(y)))
pMSE3$method <- factor(pMSE3$method, levels = c('BASD*', 'BASD-GAMLSS', 'BASD-DRF', 'BASD-Engression',
                                                'ECDF', 'LB'))

pMSE4 <- data.frame(location = rep(y, nmethod+2), 
                    value = c(as.vector(MSEr[, c(0:(nmethod-1))*length(mSet) + 4]), MSEr[, nrecord], LB[, 4]),
                    method = rep(c('BASD*', 'BASD-GAMLSS', 'BASD-DRF', 'BASD-Engression',
                                   'ECDF', 'LB'), each = length(y)))
pMSE4$method <- factor(pMSE4$method, levels = c('BASD*', 'BASD-GAMLSS', 'BASD-DRF', 'BASD-Engression',
                                                'ECDF', 'LB'))

require(ggplot2)
require(extrafont)
man_linetype = c(1, 2, 4, 6, 5, 5)
man_color = c(RColorBrewer::brewer.pal(8, 'Set1')[c(3, 1, 2, 4)], "#878787", "#4D4D4D")
plt1 <- ggplot(data = pMSE1, aes(x = location, y = value, color = method, linetype = method)) +
  geom_line(linewidth = 0.9) + 
  ylim(0.1, 1.3) + ylab(expression('MSE'(y)/'MSE'[0](y))) + xlab(expression(y)) +
  theme_bw() +
  ggtitle(paste('(a) m =', mSet[1])) +
  theme(legend.position = 'none',
        plot.title.position = 'plot') +
  scale_color_manual(values = man_color) + 
  scale_linetype_manual(values = man_linetype) 
plt2 <- ggplot(data = pMSE2, aes(x = location, y = value, color = method, linetype = method)) +
  geom_line(linewidth = 0.9) + 
  ylim(0.1, 1.3) + ylab(expression('MSE'(y)/'MSE'[0](y))) + xlab(expression(y)) +
  theme_bw() +
  ggtitle(paste('(b) m =', mSet[2])) +
  theme(legend.position = 'none',
        plot.title.position = 'plot') +
  scale_color_manual(values = man_color) + 
  scale_linetype_manual(values = man_linetype) 
plt3 <- ggplot(data = pMSE3, aes(x = location, y = value, color = method, linetype = method)) +
  geom_line(linewidth = 0.9) + 
  ylim(0.1, 1.3) + ylab(expression('MSE'(y)/'MSE'[0](y))) + xlab(expression(y)) +
  theme_bw() +
  ggtitle(paste('(c) m =', mSet[3])) +
  theme(legend.position = 'none',
        plot.title.position = 'plot') +
  scale_color_manual(values = man_color) + 
  scale_linetype_manual(values = man_linetype) 
plt4 <- ggplot(data = pMSE4, aes(x = location, y = value, color = method, linetype = method)) +
  geom_line(linewidth = 0.9) + 
  ylim(0.1, 1.3) + ylab(expression('MSE'(y)/'MSE'[0](y))) + xlab(expression(y)) +
  theme_bw() +
  ggtitle(paste('(d) m =', mSet[4])) +
  theme(legend.position = 'none',
        plot.title.position = 'plot') +
  scale_color_manual(values = man_color) + 
  scale_linetype_manual(values = man_linetype) 
plt1
plt2
plt3
plt4

pMSE <- rbind(pMSE1, pMSE2, pMSE3, pMSE4)
pMSE$m <- c(rep(paste('m = ', mSet, sep = ''), each = length(y)*(nmethod + 2)))
pMSE$m <- factor(pMSE$m, levels = paste('m = ', mSet, sep = ''))
require(ggplot2)
require(extrafont)
man_linetype = c(1, 2, 4, 6, 5, 5)
man_color = c(RColorBrewer::brewer.pal(8, 'Set1')[c(3, 1, 2, 4)], "#878787", "#4D4D4D")
plt <- ggplot(data = pMSE, aes(x = location, y = value, color = method, linetype = method)) +
  geom_line(linewidth = 0.9) + 
  facet_wrap(vars(m), dir = 'h') +
  ylim(0.1, 1.3) + ylab(expression('pMSE'(y)/'pMSE'[0](y))) + xlab(expression(y)) +
  theme_bw() +
  theme(legend.position = 'none',
        plot.title.position = 'plot') +
  scale_color_manual(values = man_color) + 
  scale_linetype_manual(values = man_linetype) 
plt

setEPS()
postscript(file = 'FIG_dist_M2.eps', horizontal = FALSE, width = 7, height = 4.9)
plt
dev.off()

################################################################################
#### TABLE plot
################################################################################

out.mse <- cbind(matrix(colMeans(MSE)[-nrecord]/mean(MSE[, nrecord]), 4), 1)
colnames(out.mse) <- c('BASD*', 'BASD-GAMLSS', 'BASD-DRF', 'BASD-Engression', 'ECDF')
out.mse

out.len.v <- rowMeans(sapply(1:num_simu, function(num){
  RES[[2, num]]/sqrt(n)*2*100
}))
out.len <- cbind(matrix(out.len.v[-nrecord], 4), out.len.v[nrecord])

out.cov.v <- rowMeans(sapply(1:num_simu, function(num){
  ceiling(rowMeans(t(abs(RES[[1, num]] - Fy.true)) > RES[[2, num]]/sqrt(n)))
}))
out.cov <- 1 - cbind(matrix(out.cov.v[-nrecord], 4), out.cov.v[nrecord])

out = t(rbind(out.mse, out.cov, out.len))

methodnames <- c('BASD*', 'BASD-GAMLSS', 'BASD-DRF', 'BASD-Engression', 'ECDF')
paste(sapply(1:nrow(out), function(rr){
  tmp = matrix(format(round(out[rr, ], 2)), nrow = 4)
  paste0(methodnames[rr], ' && ', paste(apply(tmp, 2, function(Y){ paste(Y, collapse = ' & ') }), collapse = ' && '), '\\')
}), collapse = ' ')


