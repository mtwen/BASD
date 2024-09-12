### Mengtao WEN
### 2023-09-26
### Functions of Conformal Prediction

muFun <- function(X){
  rowSums(2*(X^2 - 3))
}

sigmaFun <- function(X){
  abs((5.5 - abs(muFun(X)))/2)
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
  if (is.null(dim(X))){ X = matrix(X, ncol = 1); Xa = matrix(Xa, ncol = 1) }
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

## Conditional CDF estimation with variation for P(eps < v/sigma(X)|X)
ConCDF_var <- function(y, V, X, Xa, method = c('DIM', 'DRF'), K = 10){
  gtick = length(y)
  n = length(V)
  p = ncol(X)
  m = nrow(Xa)
  method = match.arg(method)
  
  folds_n = lapply(1:K, function(k){ seq(k, n, by = K) })
  folds_m = lapply(1:K, function(k){ seq(k, m, by = K) })
  P = matrix(NA, nrow = n+m, ncol = gtick)
  for (k in 1:K){
    idx.n = folds_n[[k]]
    idx.m = folds_m[[k]]
    
    V.t = V[-idx.n]
    X.t = X[-idx.n, ]
    U.t = sign(V.t)
    W.t = log(abs(V.t))
    
    if (method == 'DIM'){
      xnam = paste0("s(X.", 1:p, ", bs = \'cr\')")
      fmla = as.formula(paste('W ~ ', paste(xnam, collapse = "+"), '+ U'))
      DistIM = isodistrreg::dindexm(data = data.frame(W = W.t, X = X.t, U = U.t), indexfit = mgcv::gam, response = 'W', formula = fmla)
      GeneLM = glm(B~.+1, family = binomial, data = data.frame(B = (U.t+1)/2, X = X.t))
      
      pred_n1 = predict(DistIM, data.frame(X = X[idx.n, ], U = 1))
      pred_n0 = predict(DistIM, data.frame(X = X[idx.n, ], U = -1))
      pred_m1 = predict(DistIM, data.frame(X = Xa[idx.m, ], U = 1))
      pred_m0 = predict(DistIM, data.frame(X = Xa[idx.m, ], U = -1))
      pred_ns = predict(GeneLM, data.frame(X = X[idx.n, ]), type = 'response')
      pred_ms = predict(GeneLM, data.frame(X = Xa[idx.m, ]), type = 'response')
      
      idx.y = which(y >= 0)
      P[idx.n, ] = t(sapply(1:length(idx.n), function(ii){
        out = rep(NA, gtick)
        sfun1 = stepfun(pred_n1[[ii]][, 1], c(0, pred_n1[[ii]][, 3]))
        sfun0 = stepfun(pred_n0[[ii]][, 1], c(0, pred_n0[[ii]][, 3]))
        out[idx.y] = 1 - pred_ns[ii] + pred_ns[ii]*sfun1(log(y[idx.y]))
        out[-idx.y] = (1 - pred_ns[ii])*(1-sfun0(log(-y[-idx.y])))
        out
      }))
      P[n+idx.m, ] = t(sapply(1:length(idx.m), function(ii){
        out = rep(NA, gtick)
        sfun1 = stepfun(pred_m1[[ii]][, 1], c(0, pred_m1[[ii]][, 3]))
        sfun0 = stepfun(pred_m0[[ii]][, 1], c(0, pred_m0[[ii]][, 3]))
        out[idx.y] = 1 - pred_ms[ii] + pred_ms[ii]*sfun1(log(y[idx.y]))
        out[-idx.y] = (1 - pred_ms[ii])*(1-sfun0(log(-y[-idx.y])))
        out
      }))
    } else if (method == 'DRF'){
      DistRF = drf::drf(cbind(X.t, U.t), W.t, num.trees = 1000, num.features = 4, seed = 0, honesty.prune.leaves = FALSE) 
      GeneLM = glm(B~.+1, family = binomial, data = data.frame(B = (U.t+1)/2, X = X.t))
      
      W.n1 = predict(DistRF, newdata = cbind(X[idx.n, ], 1))$weights
      W.n0 = predict(DistRF, newdata = cbind(X[idx.n, ], -1))$weights
      W.m1 = predict(DistRF, newdata = cbind(Xa[idx.m, ], 1))$weights
      W.m0 = predict(DistRF, newdata = cbind(Xa[idx.m, ], -1))$weights
      pred.ns = as.vector(predict(GeneLM, data.frame(X = X[idx.n, ]), type = 'response'))
      pred.ms = as.vector(predict(GeneLM, data.frame(X = Xa[idx.m, ]), type = 'response'))
      
      idx.y = which(y >= 0)
      P[idx.n, idx.y] = as.matrix(W.n1 %*% matrix(rep(W.t, length(idx.y)) <= rep(log(y[idx.y]), each = length(W.t)), length(W.t), length(idx.y))) * pred.ns + (1 - pred.ns)
      P[idx.n, -idx.y] = (1 - as.matrix(W.n0 %*% matrix(rep(W.t, gtick - length(idx.y)) <= rep(log(-y[-idx.y]), each = length(W.t)), length(W.t), gtick - length(idx.y))))*(1-pred.ns)
      P[n+idx.m, idx.y] = as.matrix(W.m1 %*% matrix(rep(W.t, length(idx.y)) <= rep(log(y[idx.y]), each = length(W.t)), length(W.t), length(idx.y))) * pred.ms + (1 - pred.ms)
      P[n+idx.m, -idx.y] = (1 - as.matrix(W.m0 %*% matrix(rep(W.t, gtick - length(idx.y)) <= rep(log(-y[-idx.y]), each = length(W.t)), length(W.t), gtick - length(idx.y))))*(1-pred.ms)
    }
  }
  P
}


BH <- function(pval, alpha){
  n = length(pval)
  sp = sort(pval)
  thrsh = ifelse(sum(sp <= (1:n)/n*alpha), sp[max(which(sp <= (1:n)/n*alpha))], 0)
  which(pval <= thrsh)
}
