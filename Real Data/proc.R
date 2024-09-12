#### Mengtao Wen 
#### 2021-08-19
rm(list = ls()); gc()
set.seed(20240812)

################################################################################
#### Bias-correct framework
################################################################################

sampleQuant <- function(Y, tau){ # Y: a vector; tau: the tau-th quantile (0, 1]
  n = length(Y)
  sy = sort(Y, decreasing = FALSE)
  return(sy[ceiling(n*tau)])
}

BASD <- function(y, Y, P, inference = FALSE, alpha = 0.05, B = NULL){ 
  # y should be an increasing vector
  # P should be a (n+m) x q matrix where q = length(y)
  q = length(y)
  n = length(Y)
  m = nrow(P) - n
  if (inference == TRUE & is.null(B)){
    B = 1000
  }
  
  Q = matrix(rep(Y, q) <= rep(y, each = n), nrow = n, byrow = FALSE)
  funValue <- colMeans(P) + colMeans(Q) - colMeans(P[1:n, ])
  if(inference == TRUE){
    Sig.y = cov(Q)
    Sig.eps = cov(Q - P[1:n, ])
    gamma = m/(n+m)
    Sig = (1-gamma)*Sig.y + gamma*Sig.eps
    Fb = MASS::mvrnorm(B, rep(0, q), Sig)
    gb = apply(Fb, 1, function(bb){ max(abs(bb)) })
    L = quantile(gb, 1-alpha, names = FALSE)
  } else {
    L = NULL
    alpha = NULL
  }
  
  return(list(estCDF = funValue, L = L, alpha = alpha))
}

CondCDF.gamlss <- function(y, Y, X, Xa, K = 5){
  q = length(y)
  n = nrow(X)
  m = nrow(Xa)
  p = ncol(X)
  
  folds_n = lapply(1:K, function(k){ seq(k, n, by = K) })
  folds_m = lapply(1:K, function(k){ seq(k, m, by = K) })
  P = matrix(nrow = n+m, ncol = q)
  for (k in 1:K){
    idx_n = folds_n[[k]]
    idx_m = folds_m[[k]]
    
    require(gamboostLSS)
    dat.fit = data.frame(Y = Y[-idx_n], X = X[-idx_n, ])
    fml <- as.formula(paste0('Y~', paste('bbs(X.', 1:p, ')', sep = '', collapse = '+')))
    gamModel <- gamboostLSS::gamboostLSS(fml, data = dat.fit, method = 'noncyclic', 
                                         families = gamboostLSS::GaussianLSS(stabilization = 'MAD'),
                                         control = boost_control(mstop = 40))
    mu_n = predict(gamModel, parameter = 'mu', newdata = data.frame(X = X[idx_n, ]), type = 'response')
    sig_n = predict(gamModel, parameter = 'sigma', newdata = data.frame(X = X[idx_n, ]), type = 'response')
    P[idx_n, ] = t(sapply(1:length(idx_n), function(ii){ pnorm(y, mu_n[ii], sig_n[ii]) }))
    mu_m = predict(gamModel, parameter = 'mu', newdata = data.frame(X = Xa[idx_m, ]), type = 'response')
    sig_m = predict(gamModel, parameter = 'sigma', newdata = data.frame(X = Xa[idx_m, ]), type = 'response')
    P[n + idx_m, ] = t(sapply(1:length(idx_m), function(ii){ pnorm(y, mu_m[ii], sig_m[ii]) }))
  }
  return(P)
}

CondCDF.drf <- function(y, Y, X, Xa, K = 5){
  q = length(y)
  n = nrow(X)
  m = nrow(Xa)
  p = ncol(X)

  folds_n = lapply(1:K, function(k){ seq(k, n, by = K) })
  folds_m = lapply(1:K, function(k){ seq(k, m, by = K) })
  P = matrix(nrow = n+m, ncol = q)
  for (k in 1:K){
    idx_n = folds_n[[k]]
    idx_m = folds_m[[k]]
    
    drf.forest = drf::drf(X[-idx_n, ], Y[-idx_n], num.trees = 2000, num.features = 1)
    W_n = predict(drf.forest, newdata = X[idx_n, ])$weights
    W_m = predict(drf.forest, newdata = Xa[idx_m, ])$weights
    P[idx_n, ] = as.matrix(W_n) %*% matrix(rep(Y[-idx_n], q) <= rep(y, each = n - length(idx_n)), nrow = n - length(idx_n))
    P[n + idx_m, ] = as.matrix(W_m) %*% matrix(rep(Y[-idx_n], q) <= rep(y, each = n - length(idx_n)), nrow = n - length(idx_n))
  }
  return(P)
}

CondCDF.engression <- function(y, Y, X, Xa, K = 5){
  q = length(y)
  n = nrow(X)
  m = nrow(Xa)
  p = ncol(X)

  folds_n = lapply(1:K, function(k){ seq(k, n, by = K) })
  folds_m = lapply(1:K, function(k){ seq(k, m, by = K) })
  P = matrix(nrow = n+m, ncol = q)
  for (k in 1:K){
    idx_n = folds_n[[k]]
    idx_m = folds_m[[k]]
    
    ns = 10
    engModel = engression::engression(X[-idx_n, ], Y[-idx_n], standardize = FALSE,
                                      hidden_dim = 50, num_layer = 2,
                                      num_epochs = 3000, silent = FALSE)
    Fs_n = predict(engModel, X[idx_n, ], type = 'sample', nsample = ns*q)
    P[idx_n, ] = t(sapply(1:nrow(Fs_n), function(rr){ colMeans(matrix(rep(as.vector(Fs_n[rr, ]), q) <= rep(y, each = ns*q), ns*q, q)) }))
    Fs_m = predict(engModel, Xa[idx_m, ], type = 'sample', nsample = ns*q)
    P[n + idx_m, ] = t(apply(Fs_m, 1, function(samp){ colMeans(matrix(rep(as.vector(samp), q) <= rep(y, each = ns*q), ns*q, q)) }))
  }
  return(P)
}

empiricalCDF <- function(y, Y, inference = FALSE, alpha = 0.05, B = NULL){
  q = length(y)
  n = length(Y)
  if (inference == TRUE & is.null(B)){
    B = 1000
  }
  
  Q = matrix(rep(Y, q) <= rep(y, each = n), nrow = n, byrow = FALSE)
  funValue <- colMeans(Q)
  if(inference == TRUE){
    Sig.y = cov(Q)
    Fb = MASS::mvrnorm(B, rep(0, q), Sig.y)
    gb = apply(Fb, 1, function(bb){ max(abs(bb)) })
    L = quantile(gb, 1-alpha, names = FALSE)
  } else {
    L = NULL
    alpha = NULL
  }
  
  return(list(estCDF = funValue, L = L, alpha = alpha))
}

zhangYsemi <- function(Y, X, X.un, K = 10){
  n = nrow(X)
  m = nrow(X.un)
  
  folds_n = lapply(1:K, function(k){ seq(k, n, by = K) })
  folds_m = lapply(1:K, function(k){ seq(k, m, by = K) })
  P = vector(length = n+m)
  for (k in 1:K){
    idx_n = folds_n[[k]]
    idx_m = folds_m[[k]]
    
    fit.rf = ranger::ranger(tY~., data = data.frame(tY = Y[-idx_n], tX = X[-idx_n, ]), write.forest = TRUE)
    P[idx_n] = predict(fit.rf, data.frame(tX = X[idx_n, ]))$predictions
    P[n+idx_m] = predict(fit.rf, data.frame(tX = X.un[idx_m, ]))$predictions
  }
  P
}

Trunct <- function(x, n, ymax, ymin){
  theta = 0
  if (x > (n+1)*ymax - n*ymin){
    theta = (n+1)*ymax - n*ymin
  } else if (x < (n+1)*ymin - n*ymax){
    theta = (n+1)*ymin - n*ymax
  } else {
    theta = x
  }
  return(theta)
}

data0 = read.csv("HomelessData.csv")
preidx = which(data0$selection. == 2)
ranidx = which(data0$selection. == 1)
ulidx = which(data0$selection. == 0)
data = data0[, c(1, 2, 3, 4, 6, 7, 10, 11)]
data$MedianHouseholdIncome = data$MedianHouseholdIncome/1000
summary(data)

pre.data = data[preidx, ]
ran.data = data[ranidx, ]
ul.data = data[ulidx, -1]
n = nrow(ran.data)
p = ncol(ran.data) - 1
m = nrow(ul.data)
X = as.matrix(ran.data[, -1])
colnames(X) <- NULL
Y = ran.data[, 1]
Xa = as.matrix(ul.data)
colnames(Xa) <- NULL 
# ticks = 99
# taus = seq(1/(ticks+1), 1 - 1/(ticks + 1), by = 1/(1+ticks))
# y = qgamma(taus, shape = mean(Y)^2/var(Y), scale = var(Y)/mean(Y))
alpha = 0.05

################################################################################
#### Estimation and Inference
################################################################################
y = 0:max(Y)
P.LSS = CondCDF.gamlss(y, Y, X, Xa, K = 10)
ecdf.BASD.LSS = BASD(y, Y, P.LSS, inference = TRUE, alpha = alpha)
P.DRF = CondCDF.drf(y, Y, X, Xa, K = 10)
ecdf.BASD.DRF = BASD(y, Y, P.DRF, inference = TRUE, alpha = alpha)
P.ENG = CondCDF.engression(y, Y, X, Xa, K = 10)
ecdf.BASD.ENG = BASD(y, Y, P.ENG, inference = TRUE, alpha = alpha)
ecdf.classic = empiricalCDF(y, Y, inference = TRUE, alpha = alpha)

mean(Y)
mean.basd = c(sum(c(ecdf.BASD.LSS$estCDF[1], diff(ecdf.BASD.LSS$estCDF, lag = 1, differences = 1))*y),
  sum(c(ecdf.BASD.DRF$estCDF[1], diff(ecdf.BASD.DRF$estCDF, lag = 1, differences = 1))*y),
  sum(c(ecdf.BASD.ENG$estCDF[1], diff(ecdf.BASD.ENG$estCDF, lag = 1, differences = 1))*y))
var.basd = c(m/(n+m) * var(Y - apply(P.LSS[1:n, ], 1, function(prob){ sum( c(prob[1], diff(prob, lag = 1, differences = 1))*y ) })) + n/(n+m) * var(Y),
             m/(n+m) * var(Y - apply(P.DRF[1:n, ], 1, function(prob){ sum( c(prob[1], diff(prob, lag = 1, differences = 1))*y ) })) + n/(n+m) * var(Y),
             m/(n+m) * var(Y - apply(P.ENG[1:n, ], 1, function(prob){ sum( c(prob[1], diff(prob, lag = 1, differences = 1))*y ) })) + n/(n+m) * var(Y))

# ZhangA 2019 for p < n
beta.A <- lm(Y~X+1)$coefficients
mean.zhA = mean(cbind(1, rbind(X, Xa)) %*% beta.A)
var.zhA <- m/(n+m) * sum((Y - cbind(1, X) %*% beta.A)^2)/(n - p - 1) + n/(n+m) * var(Y)

# ZhangY 2021
pred = zhangYsemi(Y, X, Xa, K = 10)
mean.zhY = mean(pred) + mean(Y) - mean(pred[1:n])
var.zhY = m/(n+m) * var(Y - pred[1:n]) + n/(n+m) * var(Y)

c(mean.basd, mean.zhA, mean.zhY, mean(Y))
c(sqrt(var(Y)), sqrt(var.basd), sqrt(var.zhA), sqrt(var.zhY))*qnorm(1 - alpha/2, 0, 1)/sqrt(n)*2

paste(format(round(c(mean.basd, mean.zhA, mean.zhY, mean(Y)), 2)), collapse = ' & ')
paste(format(round(c(sqrt(var.basd), sqrt(var.zhA), sqrt(var.zhY), sqrt(var(Y)))*qnorm(1 - alpha/2, 0, 1)/sqrt(n)*2, 2)), collapse = ' & ')

taus = c(0.1, 0.3, 0.5, 0.7, 0.9)
quants = cbind(sapply(taus, function(tau){
        y[min(which(ecdf.BASD.LSS$estCDF >= tau))]
      }), 
      sapply(taus, function(tau){
        y[min(which(ecdf.BASD.DRF$estCDF >= tau))]
      }), 
      sapply(taus, function(tau){
        y[min(which(ecdf.BASD.ENG$estCDF >= tau))]
      }), 
      sampleQuant(Y, taus)
)
paste(apply(quants, 1, function(X){
  paste('&&', paste(X[1:3], collapse = ' & '), '&&', X[4])
}), collapse = '\\ ')
