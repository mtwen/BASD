### Mengtao WEN
### 2023-09-11
### LDTE with covariates

rm(list = ls()); gc()
library(doSNOW)

################################################################################
#### Functions
################################################################################

sampleQuant <- function(Y, tau){ # Y: a vector; tau: the tau-th quantile (0, 1]
  n = length(Y)
  sy = sort(Y, decreasing = FALSE)
  return(sy[ceiling(n*tau)])
}

genFun <- function(X){
  rowSums(sin(X[, 1:2]) + X[, 1:2]) + rowSums((X[, 1:5])^2)
}

### Bias-correct framework
BAS <- function(y, Y, P){ 
  # P should be a (n+m) x q matrix where q = length(y)
  q = length(y)
  n = length(Y)
  
  Q = matrix(rep(Y, q) <= rep(y, each = n), n)
  colMeans(P) + colMeans(Q) - colMeans(P[1:n, ])
}

### Estimators for P matrix
ConCDF <- function(y, Y, X, Xa, K, method = c('QR', 'DR', 'DIM', 'DRF'), taus = NULL){
  q = length(y)
  n = length(Y)
  m = nrow(Xa)
  method = match.arg(method)
  if (method == 'QR') { 
    if (is.null(taus)){ stop('Please specify taus when using method=QR') }
    s = length(taus)
    tauSet = sort(taus)
  }
  
  folds_n = lapply(1:K, function(k){ seq(k, n, by = K) })
  folds_m = lapply(1:K, function(k){ seq(k, m, by = K) })
  P = matrix(nrow = n+m, ncol = q)
  for (k in 1:K){
    idx_n = folds_n[[k]]
    idx_m = folds_m[[k]]
    
    if (method == 'QR'){
      beta_qr = sapply(tauSet, function(t){
        quantreg::rq.fit.br(cbind(1, X[-idx_n, ]), Y[-idx_n], tau = t)$coefficients
      })
      Qh_n = cbind(1, as.matrix(X[idx_n, ])) %*% beta_qr
      Qh_m = cbind(1, as.matrix(Xa[idx_m, ])) %*% beta_qr
      P[idx_n, ] = t(apply(Qh_n, 1, function(Qn){ colMeans(matrix(rep(Qn, q) <= rep(y, each = s), s)) }))
      P[n + idx_m, ] = t(apply(Qh_m, 1, function(Qm){ colMeans(matrix(rep(Qm, q) <= rep(y, each = s), s)) }))
    } else if (method == 'DR'){
      beta_dr = sapply(1:q, function(j){
        glm.fit(cbind(1, X[-idx_n, ]), (Y[-idx_n]<=y[j]), family=binomial(link="logit"), 
                control = list(maxit = 50))$coefficients
      })
      Fh_n = plogis(cbind(1, as.matrix(X[idx_n, ])) %*% beta_dr)
      Fh_m = plogis(cbind(1, as.matrix(Xa[idx_m, ])) %*% beta_dr)
      
      P[idx_n, ] = t(apply(Fh_n, 1, sort))
      P[n + idx_m, ] = t(apply(Fh_m, 1, sort))
    } else if (method == 'DIM'){
      p = ncol(X)
      xnam = paste0("s(X.", 1:p, ", bs = \'cr\')")
      fmla = as.formula(paste('Y ~ ', paste(xnam, collapse = "+")))
      distrIndexMod = isodistrreg::dindexm(data = data.frame(Y = Y[-idx_n], X = X[-idx_n, ]),
                                           indexfit = mgcv::gam, response = 'Y',
                                           formula = fmla)
      preds_n <- predict(distrIndexMod, data.frame(X = X[idx_n, ]))
      preds_m <- predict(distrIndexMod, data.frame(X = Xa[idx_m, ]))
      P[idx_n, ] = t(sapply(preds_n, function(pred_i){ 
        sfun = stepfun(pred_i[, 1], c(0, pred_i[, 3])) 
        sfun(y)
      }))
      P[n + idx_m, ] = t(sapply(preds_m, function(pred_i){ 
        sfun = stepfun(pred_i[, 1], c(0, pred_i[, 3])) 
        sfun(y)
      }))
    } else if (method == 'DRF'){
      drf.forest = drf::drf(X[-idx_n, ], Y[-idx_n], num.trees = 5000)
      W_n = predict(drf.forest, newdata = X[idx_n, ])$weights
      W_m = predict(drf.forest, newdata = Xa[idx_m, ])$weights
      P[idx_n, ] = as.matrix(W_n %*% matrix(rep(Y[-idx_n], q) <= rep(y, each = n - length(idx_n)), n - length(idx_n)))
      P[idx_m + n, ] = as.matrix(W_m %*% matrix(rep(Y[-idx_n], q) <= rep(y, each = n - length(idx_n)), n - length(idx_n)))
    }
  }
  return(P)
}

BAS_LDTE <- function(y, Y1, X1, Y0, X0, method=c('QR', 'DR', 'DIM', 'DRF'), B, K = 10, taus=NULL){
  gtick = length(y)
  n1 = length(Y1)
  n0 = length(Y0)
  n = n1+n0
  mthd = match.arg(method)
  
  P1 = ConCDF(y, Y1, X1, X0, K, mthd, taus)
  P0 = ConCDF(y, Y0, X0, X1, K, mthd, taus)
  P0.inv = rbind(P0[(1+n0):n, ], P0[1:n0, ])
  P = P1 * n0/n + P0.inv * n1/n
  Q = matrix(rep(c(Y1, Y0), gtick) <= rep(y, each = n), n)
  ks.eq = sqrt(n1*n0/n) * max(abs( colMeans((Q - P)[1:n1, ]) - colMeans((Q - P)[(n1+1):n, ]) ))
  ks.fsd = sqrt(n1*n0/n) * max( colMeans((Q - P)[(n1+1):n, ]) - colMeans((Q - P)[1:n1, ]) )
  Sig = cov(Q - P)
  M = MASS::mvrnorm(B, rep(0, gtick), Sig)
  Fn.eq = apply(M, 1, function(X) max(abs(X)))
  Fn.fsd = apply(M, 1, function(X) max(X))
  return(list(ks.eq = ks.eq, ks.fsd = ks.fsd, Fn.eq = Fn.eq, Fn.fsd = Fn.fsd))
}

Abe_LDTE <- function(y, Y1, Y0, B){
  gtick = length(y)
  n1 = length(Y1)
  n0 = length(Y0)
  n = n1 + n0
  
  Q1 = matrix(rep(Y1, gtick) <= rep(y, each = n1), n1)
  Q0 = matrix(rep(Y0, gtick) <= rep(y, each = n0), n0)
  ks.eq = sqrt(n1*n0/n) * max(abs(colMeans(Q1) - colMeans(Q0)))
  ks.fsd = sqrt(n1*n0/n) * max( colMeans(Q0) - colMeans(Q1) )
  Sig = cov(rbind(Q1, Q0))
  M = MASS::mvrnorm(B, rep(0, gtick), Sig)
  Fn.eq = apply(M, 1, function(X){ max(abs(X)) })
  Fn.fsd = apply(M, 1, function(X){ max(X) })
  return(list(ks.eq = ks.eq, ks.fsd = ks.fsd, Fn.eq = Fn.eq, Fn.fsd = Fn.fsd))
}


################################################################################
#### Data generation
################################################################################
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

################################################################################
#### Methods Evaluation
################################################################################

progress <- function(nfin){
  cat(sprintf('%s: tasks completed: %d.\n', Sys.time(), nfin))
}
opts <- list(progress = progress)

cat(sprintf('%s: tasks begin...\n', Sys.time()))
cl <- makeSOCKcluster(48)
registerDoSNOW(cl)
RES.LDTE = foreach (num=1:num_simu, .packages = c('mgcv'), .combine = 'cbind', .multicombine = TRUE, .options.snow = opts) %dopar% {
  set.seed(20230920 + num*10)
  
  Z = rbinom(n, 1, rho)
  n1 = sum(Z)
  n0 = n - n1
  tmp = rmultinom(n, 1, c(1/5, 1/5, 3/5))
  D0 = rep(0, n)
  D1 = rep(1, n)
  D0[tmp[2, ]] = 1
  D1[tmp[1, ]] = 0
  D = D0 + Z * (D1 - D0)
  X = matrix(rnorm(n*p), n, p)
  Y00 = genFun(X) + sigma*rnorm(n) # Y0 ~ N(0, 6.25)
  
  grid.tick = 2000 - 1
  taus = seq(0.0005, 0.9995, length = grid.tick)
  y0 = quantile(Y00, taus, names = FALSE)
  
  res = vector('list', length = length(deltaSet))
  for (dd in 1:length(deltaSet)){
    delta = deltaSet[dd]
    
    Y = Y00 + D*delta
    y1 = y0 + delta
    y = sort(unique(c(y0, y1)), decreasing = FALSE)
    
    Y1 = Y[which(Z == 1)]
    Y0 = Y[which(Z == 0)]
    X1 = X[which(Z == 1), ]
    X0 = X[which(Z == 0), ]
    
    KS.bas = BAS_LDTE(y, Y1, X1, Y0, X0, method = 'DIM', B, K = K)
    cn.bas.eq = quantile(KS.bas$Fn.eq, 1-alpha, names = FALSE)
    cn.bas.fsd = quantile(KS.bas$Fn.fsd, 1-alpha, names = FALSE)
    
    # Abadie 2002
    KS.Ab = Abe_LDTE(y, Y1, Y0, B)
    cn.Ab.eq = quantile(KS.Ab$Fn.eq, 1-alpha, names = FALSE)
    cn.Ab.fsd = quantile(KS.Ab$Fn.fsd, 1-alpha, names = FALSE)
    
    # res[[dd]] = list(ks = c(KS.Ab, KS.s0, KS.s3), cn = c(cn.Ab, cn.s0, cn.s3), 
    #      rj = c(KS.Ab > cn.Ab, KS.s0 > cn.s0, KS.s3 > cn.s3))
    cat(sprintf('\n %s: tasks with delta = %.1f completed: %d.\n', Sys.time(), delta, num))
    
    res[[dd]] = list(ks.eq = c(KS.Ab$ks.eq, KS.bas$ks.eq), 
                     ks.fsd = c(KS.Ab$ks.fsd, KS.bas$ks.fsd), 
                     cn.eq = c(cn.Ab.eq, cn.bas.eq), 
                     cn.fsd = c(cn.Ab.fsd, cn.bas.fsd),
                     rj.eq = c(KS.Ab$ks.eq > cn.Ab.eq, KS.bas$ks.eq > cn.bas.eq),
                     rj.fsd = c(KS.Ab$ks.fsd > cn.Ab.fsd, KS.bas$ks.fsd > cn.bas.fsd))
  }
  res
}
stopCluster(cl)
save(RES.LDTE, file = filename)
cat(sprintf('%s: tasks end and results saved.\n', Sys.time()))

