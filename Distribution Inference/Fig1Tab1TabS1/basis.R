## MengtaoWEN
## 2024-07-17
## Functions related to BASD

### Bias-correct framework
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

### Estimators for P matrix
CondCDF <- function(y, Y, X, Xa, K, method = c('dr', 'dim', 'drf', 'qrf', 'engression', 'gamlss'), control = NULL){
  # y should be an increasing vector
  q = length(y)
  n = nrow(X)
  m = nrow(Xa)
  p = ncol(X)
  method = match.arg(method)
  
  folds_n = lapply(1:K, function(k){ seq(k, n, by = K) })
  folds_m = lapply(1:K, function(k){ seq(k, m, by = K) })
  P = matrix(nrow = n+m, ncol = q)
  for (k in 1:K){
    idx_n = folds_n[[k]]
    idx_m = folds_m[[k]]
    
    if (method == 'dr'){
      if (is.null(control)) { control = list(maxit = 100) }
      idxl = which(sapply(1:q, function(jj){ sum(Y[-idx_n] <= y[jj]) }) <= 1)
      idxu = which(sapply(1:q, function(jj){ sum(Y[-idx_n] > y[jj]) }) <= 1)
      beta_dr = sapply(1:q, function(j){
        if (j %in% union(idxl, idxu)){ return(rep(NA, p+1)) }
        glm.fit(cbind(1, X[-idx_n, ]), (Y[-idx_n]<=y[j]), family=binomial(link="logit"), 
                intercept = FALSE, control = control)$coefficients
      })
      Fh_n = plogis(cbind(1, as.matrix(X[idx_n, ])) %*% beta_dr)
      Fh_m = plogis(cbind(1, as.matrix(Xa[idx_m, ])) %*% beta_dr)
      for (jj in idxl){ Fh_n[, jj] = 0; Fh_m[, jj] = 0 }
      for (jj in idxu){ Fh_n[, jj] = 1; Fh_m[, jj] = 1 }
      P[idx_n, ] = t(apply(Fh_n, 1, sort)) # rearrangement step
      P[n + idx_m, ] = t(apply(Fh_m, 1, sort)) # rearrangement step
    } else if (method == 'gamlss'){
      require(gamboostLSS)
      if (is.null(control)) { control = mboost::boost_control(mstop = 200) }
      dat.fit = data.frame(Y = Y[-idx_n], X = X[-idx_n, ])
      fml <- as.formula(paste0('Y~', paste('bbs(X.', 1:p, ')', sep = '', collapse = '+')))
      # formu <- as.formula(paste0('Y ~ ', paste('X.', 1:p, sep = '', collapse = ' + ')))
      gamModel <- gamboostLSS::gamboostLSS(fml, data = dat.fit, method = 'noncyclic', 
                                           families = GaussianLSS(stabilization = 'MAD'),
                                           control = control)
      mu_n = predict(gamModel$mu, newdata = data.frame(X = X[idx_n, ]), type = 'response')
      sig_n = predict(gamModel$sigma, newdata = data.frame(X = X[idx_n, ]), type = 'response')
      P[idx_n, ] = t(sapply(1:length(idx_n), function(ii){ pnorm(y, mu_n[ii], sig_n[ii]) }))
      mu_m = predict(gamModel$mu, newdata = data.frame(X = Xa[idx_m, ]), type = 'response')
      sig_m = predict(gamModel$sigma, newdata = data.frame(X = Xa[idx_m, ]), type = 'response')
      P[n + idx_m, ] = t(sapply(1:length(idx_m), function(ii){ pnorm(y, mu_m[ii], sig_m[ii]) }))
    } else if (method == 'engression'){
      if (is.null(control)){ control = list(hidden_dim = 100, num_layer = 4, num_epochs = 2000) }
      ns = 10
      engModel = engression::engression(X[-idx_n, ], Y[-idx_n], standardize = FALSE,
                                        hidden_dim = control$hidden_dim, num_layer = control$num_layer,
                                        num_epochs = control$num_epochs, silent = FALSE)
      Fs_n = predict(engModel, X[idx_n, ], type = 'sample', nsample = ns*q)
      P[idx_n, ] = t(sapply(1:nrow(Fs_n), function(rr){ colMeans(matrix(rep(as.vector(Fs_n[rr, ]), q) <= rep(y, each = ns*q), ns*q, q)) }))
      Fs_m = predict(engModel, Xa[idx_m, ], type = 'sample', nsample = ns*q)
      P[n + idx_m, ] = t(apply(Fs_m, 1, function(samp){ colMeans(matrix(rep(as.vector(samp), q) <= rep(y, each = ns*q), ns*q, q)) }))
    } else if (method == 'dim'){
      p = ncol(X)
      xnam = paste0("s(X.", 1:p, ", bs = \'cr\', sp = ", floor(n*(K-1)/K/p - 1), ')')
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
    } else if (method == 'drf'){
      if (is.null(control)){ control = list(num.trees = 2000, num.features = 5) }
      drf.forest = drf::drf(X[-idx_n, ], Y[-idx_n], num.trees = control$num.trees, num.features = control$num.features)
      W_n = predict(drf.forest, newdata = X[idx_n, ])$weights
      W_m = predict(drf.forest, newdata = Xa[idx_m, ])$weights
      P[idx_n, ] = as.matrix(W_n) %*% matrix(rep(Y[-idx_n], q) <= rep(y, each = n - length(idx_n)), nrow = n - length(idx_n))
      P[n + idx_m, ] = as.matrix(W_m) %*% matrix(rep(Y[-idx_n], q) <= rep(y, each = n - length(idx_n)), nrow = n - length(idx_n))
    } else if (method == 'qrf'){
      qtick = 2000
      probs = seq(1/qtick, 1- 1/qtick, by = 1/qtick)
      grfModel = grf::quantile_forest(X[-idx_n, ], Y[-idx_n], quantiles = probs)
      Qn = predict(grfModel, newdata = X[idx_n, ])$predictions
      P[idx_n, ] = t(apply(Qn, 1, function(margin){ sf <- stepfun(margin, c(0, probs)); sf(y) }))
      Qm = predict(grfModel, newdata = Xa[idx_m, ])$predictions
      P[n + idx_m, ] = t(apply(Qm, 1, function(margin){ sf<-stepfun(margin, c(0, probs)); sf(y) }))
    }
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