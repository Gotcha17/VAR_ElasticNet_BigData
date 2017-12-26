setClass("enetVAR.enetVAR", representation = "list", S3methods = T)
setClass("enetVAR.GroupEnetVAR", representation = "list", S3methods = T)
setClass("enetVAR.VARZ", representation = "list", S3methods = T)


# VAR elastic net ---------------------------------------------------------
require(glmnet)


.enetVAR = function(j, y, p, Z, intercept = TRUE, alpha.vec, lambda.vec, parallel) {
  colIndex = j[1]
  if (length(alpha.vec)==1){
    alpha <- alpha.vec
  } else {
    alpha <- alpha.vec[colIndex]
  }
  if (!missing(lambda.vec)){
    if (length(lambda.vec)==1){
      lambda <- lambda.vec
    } else if (length(lambda.vec)>1){
      lambda <- unlist(lambda.vec)[colIndex]
    }
    return (glmnet(Z, j[-1], family = "gaussian", intercept = intercept,
                   alpha=alpha, lambda = seq((lambda*2), (lambda/2), length.out = 10),
                   standardize = TRUE))
  } else {
    # same fold structure is created for comparability of the results
    k.fold <- 10
    folds <- floor(NROW(Z)/k.fold)           # 10-fold
    last.fold <- NROW(Z) - folds*k.fold
    foldid <- as.vector(sapply(seq(1,(folds+1)), function(x){rep(x, k.fold)}))
    foldid <- foldid[1:(length(foldid)-(k.fold-last.fold))]
    return (cv.glmnet(Z, j[-1], family = "gaussian", intercept = intercept,
                      alpha=alpha, standardize = TRUE, foldid = foldid,
                      parallel = parallel))
  }
}

#' Vector Autoregression Elastic Net
#'
#' Fit a vector autoregressive model with elastic net penalty.
#' The VAR model is estimated using a multiresponse linear regression.
#' The VAR Elastic Net fits multiple uniresponse linear regressions with elastic net penalty.
#' @param y A matrix where each column represents an individual time series
#' @param p Number of lags to include in the design matrix
#' @param alpha numeric, vector. Elastic net mixing parameter, as defined in 'glmnet'
#'  (vector for each equation or single value that is used in all equations)
#' @param intercept logical.  If true, include a 1 vector for the intercept term
#' @param lambda numeric, vector. elastic net penalty as defined in 'glmnet' (same as in alpha)
#' @param parallel logical. If true, parallel backend must be registered prior
#' @export
enetVAR = function(y, p, alpha = rep(0.4,ncol(y)), intercept = F, lambda,
                  parallel = F, ...) {
  if (!length(alpha)==0){
    alpha.rtn <- alpha[1]
  } else {
    alpha.rtn <- alpha
  }
  if (p < 1) stop("p must be a positive integer")
  if (!is.matrix(y)) {
    stop("y must be a matrix (not a data frame).  Consider using as.matrix(y)")
  }
  var.z = VAR.Z(y,p,intercept = intercept)
  Z = var.z$Z
  y.augmented = rbind(1:ncol(y),var.z$y.p)
  
  var.lasso = apply(y.augmented, 2, .enetVAR, y=y,p=p,
                    Z=Z, alpha=alpha, lambda.vec = lambda,
                    intercept = intercept, parallel = parallel)

  return(structure(list(
    model = var.lasso,
    var.z = var.z,
    alpha = alpha.rtn),
    class="enetVAR.enetVAR"))
}

#' Coefficients of a enetVAR model
#'
#' The underlying library, glmnet, computes the full path to the lasso.
#' This means it is computationally easy to compute the lasso solution
#' for any penalty term.  This function allows you to pass in the desired
#' l1 penalty and will return the coefficients
#' @param enetVAR an object of class enetVAR.enetVAR
#' @param l1penalty The l1 penalty to be applied to the enetVAR
#'   model.
#' @method coef enetVAR.enetVAR
#' @S3method coef enetVAR.enetVAR
coef.enetVAR.enetVAR = function(enetVAR, l1penalty) {
  if (missing(l1penalty)) {
    B = data.frame(lapply(enetVAR$model, function(model) {
      as.vector(coef(model, model$lambda.min))
    }))
  }
  else if (length(l1penalty) == 1) {
    B = data.frame(lapply(enetVAR$model, function(model) {
      as.vector(coef(model, as.numeric(l1penalty)))
    }))
  } else {
    B = matrix(0, nrow=(ncol(enetVAR$var.z$Z)+1), ncol=ncol(enetVAR$var.z$y.p))
    for (i in 1:length(l1penalty)) {
      B[,i] = as.vector(coef(enetVAR$model[[i]], as.numeric(l1penalty[i])))
    }
  }
  if (enetVAR$var.z$intercept){
    B <- B[-2,]
  } else {
    B <- B[-1,]
  }
  colnames(B) = colnames(enetVAR$var.z$y.orig)
  rownames(B) = colnames(enetVAR$var.z$Z)
  
  return (as.matrix(B))
}

#' enetVAR Predict
#'
#' Predict n steps ahead from a enetVAR.enetVAR object
#' @param enetVAR an object of class enetVAR.enetVAR returned from enetVAR
#' @param n.ahead number of steps to predict
#' @param ... extra parameters to pass into the coefficients method
#'   for objects of type enetVAR.enetVAR
#' @examples
#'   data(Canada)
#'   predict(enetVAR(Canada, p = 3), 1)
#' @method predict enetVAR.enetVAR
#' @S3method predict enetVAR.enetVAR
predict.enetVAR.enetVAR = function(enetVAR, n.ahead=1, lambda.vec, ...) {
  y.pred = matrix(nrow=n.ahead, ncol=ncol(enetVAR$var.z$y.orig))
  colnames(y.pred) = colnames(enetVAR$var.z$y.orig)
  y.orig = enetVAR$var.z$y.orig
  for (i in 1:n.ahead) {
    if (enetVAR$var.z$intercept) {
      Z.ahead = c(1,as.vector(t(y.orig[
        ((nrow(y.orig)):
           (nrow(y.orig)-enetVAR$var.z$p+1))
        ,])))
    } else {
      Z.ahead = as.vector(t(y.orig[
        ((nrow(y.orig)):
           (nrow(y.orig)-enetVAR$var.z$p+1))
        ,]))
    }
    if (missing(lambda.vec)){
      y.ahead = Z.ahead %*% coef.enetVAR.enetVAR(enetVAR, ...)
    } else {
      y.ahead = Z.ahead %*% coef.enetVAR.enetVAR(enetVAR, lambda.vec, ...)
    }
    y.pred[i,] = y.ahead
    if (i == n.ahead) break
    y.orig = rbind(y.orig, y.ahead)
  }
 return (y.pred)
}


# Exctracting residuals from enetVAR.enetVAR
#' ElnetVAR Residuals
#'
#' Exctracting residuals from enetVAR.enetVAR
#' @param enetVAR an object of class enetVAR.GroupEnetVAR
#' 
#' @method residuals enetVAR.enetVAR
#' @S3method residuals enetVAR.enetVAR
residuals.enetVAR <- function(enetVAR, lambda){
  co <- coef.enetVAR.enetVAR(enetVAR, lambda)
  Y <- enetVAR$var.z$y.p
  Z <- enetVAR$var.z$Z
  resid <- Y - Z %*% co
  colnames(resid) <- colnames(enetVAR$var.z$y.orig)
  rownames(resid) <- rownames(Y)
  
  return(resid)
}

### Calculation of information criteria (FPE,AIC,HQ,SC) for Elnet Penalty with adopted dof
infCrit <- function(enetVAR){
  T <- NROW(enetVAR$var.z$y.p)
  K <- NCOL(enetVAR$var.z$y.orig)
  m <- enetVAR$var.z$p
  alpha <- enetVAR$alpha
  lambda <- mean(sapply(enetVAR$model, function(x){x$lambda.min}))
  U_hat <- t(residuals.enetVAR(enetVAR))
  Sigma_hat <- (U_hat %*% t(U_hat))/T
  det_Sigma_hat <- det(Sigma_hat)
  ## Negative determinant of the covariance matrix indicates an ill-conditioned model,
  #  thus this model should not be considered and gets a large penalty
  if (det_Sigma_hat<0) det_Sigma_hat <- 1000
  # y = (Z' kron I_k)*beta+u
  X <- kronecker(enetVAR$var.z$Z, diag(K))                         # (Z' kron I_k)
  B <- t(coef.enetVAR.enetVAR(enetVAR))                             # B=(nu, A1,...,Ap) coeffs
  beta <- t(as.vector(B))                                         # VAR(p) for LS estimation form
  A <- (beta>0 | beta<0)                                          # Support of beta
  X_A <- X[,A]                                                    # Formed from X and non-zero coeffs
  # Degrees of freedom for Elastic Net in accordance to TIBSHIRANI AND TAYLOR (2011)
  dof <- sum(diag(X_A%*%solve(t(X_A)%*%X_A+lambda*0.5*(1-alpha)*diag(sum(A)))%*%t(X_A)))
  FPE <- (1 + dof/T)/(1 - dof/T) * det_Sigma_hat
  AIC <- log(det_Sigma_hat) + 2/T * dof
  HQ <- log(det_Sigma_hat) + 2*log(log(T))/T * dof
  SC <- log(det_Sigma_hat) + (log(T))/T * dof
  return(list(FPE = FPE, AIC = AIC, HQ = HQ, SC = SC, dof = dof))
}

enetVARselect <- function(data, max.lag = 30, alpha = 0.25, intercept = F){
  iteration <- 0
  result <- c()
  for (p in (1:max.lag)){
    fit.enetVAR <- enetVAR(data, p = p, alpha = alpha, intercept = intercept)
    ic <- infCrit(fit.enetVAR)
    result <- cbind(result, ic)
    iteration <- iteration + 1
    tic <- matrix(unlist(result), ncol = iteration, nrow = 4)
    if (iteration > 3){
      inf.sum <- sum(tic[,iteration]==(-Inf))
      if (inf.sum>2) break
      ch <- 0
      for (i in 1:4){
        min <- min(tic[i,])
        ch <- ch + (sum(result[i,(iteration-3):iteration]>min)>3)
      }
      if (ch==4) break
    }
  }
  FPE_p <- which(tic[1,]==min(tic[1,]))
  AIC_p <- which(tic[2,]==min(tic[2,]))
  HQ_p <- which(tic[3,]==min(tic[3,]))
  SC_p <- which(tic[4,]==min(tic[4,]))
  colnames(result) <- seq(1,iteration)
  ic_matrix <- matrix(c(FPE_p, AIC_p, HQ_p, SC_p), ncol = 4, nrow = 1)
  colnames(ic_matrix) <- c("FPE", "AIC", "HQ", "SC")
  return(list(IC_lag = ic_matrix, IC_value = result))
}

######### Enet Variables Preselection
enetVARpreselection <- function(data, lag = 1, alpha=0.25, maxnrvar){
  # target variable is supposed to be in the first column
  targetvars <- colnames(data)[1]
  varlist <- colnames(data)[-1]
  while (length(targetvars)<maxnrvar){
    ic.tmp <- c()
    for (i in 1:length(varlist)){
      selcvar <- varlist[i]
      test.matrix <- data[,c(targetvars, selcvar)]
      fit.enetVAR <- enetVAR(test.matrix, p = lag, alpha = alpha)
      ic.sc <- infCrit(fit.enetVAR)[[4]]
      ic.tmp <- c(ic.tmp, ic.sc)
    }
    min.ic <- which(min(ic.tmp)==ic.tmp)
    min.var <- varlist[min.ic]
    varlist <- varlist[-min.ic]
    targetvars <- c(targetvars, min.var)
  }
  return(targetvars)
}
# VAR design matrix -------------------------------------------------------


#' VAR Design Matrix
#'
#' Compute the Design Matrix for a standard VAR model
#'
#' @param y a multivariate time series matrix where rows represent
#'   time and columns represent different time series
#' @param p the lag order of the VAR model
#' @param intercept logical.  If true, include a 1 vector for the intercept term
#' @return
#'  \item{n}{Number of endogenous time series that are being measured}
#'  \item{T }{The number of time points in the reduced response matrix}
#'  \item{k }{The total number of predictor variables used to model each endogenous time series}
#'  \item{dof }{The degrees of freedom of the residuals}
#'  \item{y.p }{The reduced response matrix}
#'  \item{Z }{The design matrix}
#'  \item{y.orig}{The original input matrix}
#'  \item{p}{The lag order of the VAR model}
#'  \item{intercept}{logical.  If true, include a 1 vector for the intercept term}
#' @export
VAR.Z = function(y, p, intercept=F) {
  if(p < 1) stop("p must be a positive integer")
  
  if(is.null(colnames(y))) {
    colnames(y) = sapply(1:ncol(y), function(j) {
      paste('y',j,sep='')
    })
  }
  
  n = ncol(y) #numcols in response
  T = nrow(y) #numrows in response
  k = n*p     #numcols in design matrix
  if(intercept) k = k+1
  dof = T - p - k
  
  y.p = y[(1+p):T,]
  Z = do.call('cbind', (sapply(0:(p-1), function(i) {
    y[(p-i):(T-1-i),]
  }, simplify=F)))
  Z.names = colnames(Z)
  colnames(Z) = as.vector(sapply(1:p, function(i) {
    startIndex = (i-1)*n + 1
    endIndex = i*n
    paste(Z.names[startIndex:endIndex], '.l', i, sep='')
  }))
  
  if(intercept) {
    Z = cbind(1, Z)
    colnames(Z)[1] = "intercept"
  }
  
  return ( structure ( list (
    n=n,
    T=nrow(y.p),
    k=k,
    dof=dof,
    y.p=as.matrix(y.p),
    Z=as.matrix(Z),
    y.orig=y,
    p = p,
    intercept = intercept
  ), class="enetVAR.VARZ"))
}

# Group VAR elastic net -------------------------------------------------------


#' Group Vector Autoregression via Group Lasso
#'
#' Fit a VAR model by creating the lagged design matrix
#' and fitting a multivariate response matrix to it.  Note that
#' the VAR response matrix omits the first p responses from the input
#' matrix Y.  This is because observations in Y are modeled by the
#' p previous values, so the first p observations cannot be modeled.
#'
#' While multivariate response regressions can be solved as multiple
#' univariate response regressions, this multivariate response problem
#' can better be solved by using Group Lasso.  Instead of seeking sparsity
#' in the coefficients for each univariate response, Group Lasso attempts to
#' find sparsity in the coefficient matrix as a whole.
#' @param y A matrix where each column represents an individual time series
#' @param freq only used if the time series are periodic.  freq is a vector of
#'   frequencies for each of the time series, as in 'ts(y, freq = ...)'.
#'   If the time series are not periodic, then this vector can be a vector of NA
#' @param p the number of lags to include in the design matrix
#' @param alpha the mixing parameter between group lasso and quadratic, as in 'glmnet'
#' @export
GroupEnetVAR = function(y, p=1, alpha = 0.4, intercept = F,
                        parallel = F, ...) {
  if (!is.matrix(y)){
    stop("y must be a matrix (not a data frame).  Consider using as.matrix(y)")
  }
  if (p < 1) {
    stop("p must be a positive integer")
  }
  var.z = VAR.Z(y, p, intercept = intercept)
  k.fold <- 10
  folds <- floor(NROW(var.z$Z)/k.fold)           # 10-fold
  last.fold <- NROW(var.z$Z) - folds*k.fold
  foldid <- as.vector(sapply(seq(1,(folds+1)), function(x){rep(x, k.fold)}))
  foldid <- foldid[1:(length(foldid)-(k.fold-last.fold))]
  model = cv.glmnet(var.z$Z, var.z$y.p, family = "mgaussian", alpha = alpha,
                    intercept = intercept, standardize = TRUE,
                    standardize.response = TRUE, parallel = parallel,
                    foldid = foldid, ...)
  result = structure(list(model = model,
                          var.z = var.z),
                     class = "enetVAR.GroupEnetVAR")
  return (result)
}

#' GroupEnetVAR Coefficients
#'
#' @param GroupEnetVAR an object of class enetVAR.GroupEnetVAR
#' @param ... extra parameters to pass into the coefficients method
#'   for objects of type enetVAR.GroupEnetVAR
#' @return The coefficients for the VAR model
#' @method coef enetVAR.GroupEnetVAR
#' @S3method coef enetVAR.GroupEnetVAR
coef.enetVAR.GroupEnetVAR = function(GroupEnetVAR, ...) {
  coefs = as.matrix(do.call('cbind', lapply(coef(GroupEnetVAR$model$glmnet.fit,
                                                 GroupEnetVAR$model$lambda.min,
                                                 ...), as.matrix)))
  if (GroupEnetVAR$var.z$intercept){
    coefs <- coefs[-2,]
  } else {
    coefs <- coefs[-1,]
  }
  colnames(coefs) = colnames(GroupEnetVAR$var.z$y.orig)
  return (coefs)
}

#' GroupEnetVAR Predict
#'
#' Predict n steps ahead from a enetVAR.GroupEnetVAR object
#' @param GroupEnetVAR an object of class enetVAR.GroupEnetVAR returned from GroupEnetVAR
#' @param n.ahead number of steps to predict
#' @param ... extra parameters to pass into the coefficients method
#'   for objects of type enetVAR.GroupEnetVAR
#' @method predict enetVAR.GroupEnetVAR
#' @S3method predict enetVAR.GroupEnetVAR
predict.enetVAR.GroupEnetVAR = function(GroupEnetVAR, n.ahead, ...) {
  y.pred = matrix(nrow=n.ahead, ncol=ncol(GroupEnetVAR$var.z$y.orig))
  colnames(y.pred) = colnames(GroupEnetVAR$var.z$y.orig)
  y.orig = GroupEnetVAR$var.z$y.orig
  for (i in 1:n.ahead) {
    if (GroupEnetVAR$var.z$intercept) {
      Z.ahead = c(1,as.vector(t(y.orig[
        ((nrow(y.orig)):
           (nrow(y.orig)-GroupEnetVAR$var.z$p+1))
        ,])))
      y.ahead = Z.ahead %*% coef.enetVAR.GroupEnetVAR(GroupEnetVAR, ...)
    } else {
      Z.ahead = as.vector(t(y.orig[
        ((nrow(y.orig)):
           (nrow(y.orig)-GroupEnetVAR$var.z$p+1))
        ,]))
      y.ahead = Z.ahead %*% coef.enetVAR.GroupEnetVAR(GroupEnetVAR, ...)
    }
    y.pred[i,] = y.ahead
    if (i == n.ahead) break
    y.orig = rbind(y.orig, y.ahead)
  }
  return (y.pred)
}



# Model Training ----------------------------------------------------------

modeltrain <- function(data, start.pred="2000 Q1", step=1, h=8,
                       method="enet", alpha=0.4, lambda,
                       lag=1, const=FALSE, ...){
  if (!is.zoo(data)){
    data <- zoo(data)
  }
  dates <- time(data)
  startdate <- dates[1]
  window.size <- which(dates==start.pred)-h
  sequence <- seq(from = window.size, to = (length(dates)-1), by = step)
  pred.ind <- c(1, 2, 2*seq(from = 2, to = h/2,  by = 2))
  var.err <- c()
  forecasts <- c()
  u_1 <- 0
  u_2 <- 0
  for (i in sequence){
    traindata <- as.matrix(window(data, start = startdate, end = dates[i]))
    if (method == "enet"){
      train.fit <- enetVAR(traindata, p = lag, alpha = alpha, lambda = lambda,
                          intercept = const)
      train.pred <- predict.enetVAR.enetVAR(train.fit, n.ahead = h, lambda = lambda)
    } else if (method == "genet"){
      train.fit <- GroupEnetVAR(traindata, p = lag, alpha = alpha,
                                intercept = const, ...)
      train.pred <- predict.enetVAR.GroupEnetVAR(train.fit, n.ahead = h)
    }
    var.pred.values <- train.pred[pred.ind,1]
    true.values <- as.vector((window(data, index = dates[(i+pred.ind)]))[,1])
    while(length(true.values)<4){
      true.values <- c(true.values, 0)
    }
    forecasts <- cbind(forecasts, var.pred.values)
    var.err <- cbind(var.err, (var.pred.values-true.values))
    y_t <- as.vector((window(data, index = dates[(i+pred.ind-1)]))[,1])
    while(length(y_t)<4){
      y_t <- c(y_t, 0)
    }
    u_2 <- cbind(u_2, (true.values-y_t))
  }
  h1.ind <- h:ncol(var.err)
  h2.ind <- (h-1):(ncol(var.err)-1)
  h4.ind <- (h-3):(ncol(var.err)-3)
  h8.ind <- 1:(ncol(var.err)-h+1)
  h.ind.len <- length(h1.ind)
  for.err <- list(h1 = var.err[1,h1.ind], h2 = var.err[2,h2.ind],
                  h4 =  var.err[3,h4.ind], h8 = var.err[4,h8.ind])
  forecast <- list(h1 = forecasts[1,h1.ind], h2 = forecasts[2,h2.ind],
                  h4 =  forecasts[3,h4.ind], h8 = forecasts[4,h8.ind])
  u_1 <- list(h1 = sum((for.err$h1)^2), h2 = sum((for.err$h2)^2),
              h4 = sum((for.err$h4)^2), h8 = sum((for.err$h8)^2))
  u_2 <- list(h1 = sum(u_2[1,h1.ind]^2), h2 = sum(u_2[2,h2.ind]^2),
              h4 = sum(u_2[3,h4.ind]^2), h8 = sum(u_2[4,h8.ind]^2))
  theils_u <- list(h1 = sqrt(u_1$h1/u_2$h1), h2 = sqrt(u_1$h2/u_2$h2),
                   h4 = sqrt(u_1$h4/u_2$h4), h8 = sqrt(u_1$h8/u_2$h8))
  var.msfe <- list(h1 = (u_1$h1)/h.ind.len, h2 = (u_1$h2)/h.ind.len,
                   h4 =  (u_1$h4)/h.ind.len, h8 =  (u_1$h8)/h.ind.len)
  theils_u_ar1 <- list(h1 = theils_u_ar1(data[,1], mse_pred = var.msfe$h1, h="h1"), 
                       h2 = theils_u_ar1(data[,1], mse_pred = var.msfe$h2, h="h2"),
                       h4 = theils_u_ar1(data[,1], mse_pred = var.msfe$h4, h="h4"),
                       h8 = theils_u_ar1(data[,1], mse_pred = var.msfe$h8, h="h8"))
  residuals <- residuals.enetVAR(train.fit, lambda)
  return(structure(list(
    forecast = forecast,
    for.err = for.err,
    var.msfe = var.msfe,
    theils_u.r.w. = theils_u,
    theils_u_ar1 = theils_u_ar1,
    residuals = residuals,
    model = train.fit)))
}

modeltrain.slim <- function(data, start.pred="2000 Q1", step=1, h=8,
                            alpha=0.4, lambda, lag=1, const=FALSE){
  if (!is.zoo(data)){
    data <- zoo(data)
  }
  dates <- time(data)
  startdate <- dates[1]
  window.size <- which(dates==start.pred)-h
  sequence <- seq(from = window.size, to = (length(dates)-1), by = step)
  pred.ind <- c(1, 2, 2*seq(from = 2, to = h/2,  by = 2))
  var.err <- c()
  for (i in sequence){
    traindata <- as.matrix(window(data, start = startdate, end = dates[i]))
    train.fit <- enetVAR(traindata, p = lag, alpha = alpha, lambda = lambda,
                        intercept = const)
    train.pred <- predict.enetVAR.enetVAR(train.fit, n.ahead = h, lambda = lambda)
    var.pred.values <- train.pred[pred.ind,1]
    true.values <- as.vector((window(data, index = dates[(i+pred.ind)]))[,1])
    var.err <- cbind(var.err, (var.pred.values-true.values))
  }
  h1.ind <- h:ncol(var.err)
  h2.ind <- (h-1):(ncol(var.err)-1)
  h4.ind <- (h-3):(ncol(var.err)-3)
  h8.ind <- 1:(ncol(var.err)-h+1)
  h.ind.len <- length(h1.ind)
  for.err <- list(h1 = var.err[1,h1.ind], h2 = var.err[2,h2.ind],
                  h4 =  var.err[3,h4.ind], h8 = var.err[4,h8.ind])
  var.msfe <- list(h1 = sum((for.err$h1)^2)/h.ind.len, h2 = sum((for.err$h2)^2)/h.ind.len,
                   h4 =  sum((for.err$h4)^2)/h.ind.len, h8 =  sum((for.err$h8)^2)/h.ind.len)
  return(structure(list(
    for.err = for.err,
    var.msfe = var.msfe)))
}



require(caret);
#library(doParallel);
#clcaret <- makeCluster(7) # create a cluster with x=? cores
#registerDoParallel(cores=clcaret)
enetVARtune <- function(data, lag, init.window, horizon, parallel = F){
  oldw <- getOption("warn")
  options(warn = -1)
  
  caretdata <- VAR.Z(data, p = lag, intercept = FALSE)
  caret.y <- as.matrix(caretdata$y.p)
  caret.x <- as.matrix(caretdata$Z)
  # Control of the training procedure settings
  myTimeControl <- trainControl(method = "timeSlice",
                                initialWindow = init.window,
                                horizon = horizon,
                                fixedWindow = FALSE,
                                allowParallel = parallel)

  # Train sequence with settings for the method
  var.tune = apply(caret.y, 2, train, x=caret.x, method = "glmnet",
                   family = "gaussian", trControl = myTimeControl,
                   standardize = TRUE, maxit = 100000, intercept = FALSE,
                   tuneGrid = expand.grid(.alpha = seq(0.05, 0.95, by = 0.05),
                                          .lambda = 10^seq(1,-4,length = 200)))

  # Extraction of best tune values
  best.tune = data.frame(lapply(var.tune, function(model) {
    t(model$bestTune)}))
  # returning alpha and lambda tune vectors
  options(warn = oldw)
  return(best.tune)
}

  
ar1_train <- function(data, start.pred="2000 Q1", step=1, h=8,
                       lag=1, const=FALSE){
  if (!is.zoo(data)){
    data <- zoo(data)
  }
  mse.sum <- 0
  dates <- time(data)
  startdate <- dates[1]
  window.size <- which(dates==start.pred)-h
  sequence <- seq(from = window.size, to = (length(dates)-1), by = step)
  pred.ind <- c(1, 2, 2*seq(from = 2, to = h/2,  by = 2))
  var.err <- c()
  forecast <- c()
  for (i in sequence){
    traindata <- as.matrix(window(data, start = startdate, end = dates[i]))
    train.fit <- arima(traindata, order = c(1,0,0), 
                       include.mean = const, transform.pars = FALSE)
    train.pred <- predict(train.fit, n.ahead = h)
    var.pred.values <- as.vector(train.pred$pred)[pred.ind]
    true.values <- as.vector(window(data, index = dates[(i+pred.ind)]))
    while(length(true.values)<4){
      true.values <- c(true.values, 0)
    }
    var.err <- cbind(var.err, (var.pred.values-true.values))
    forecast <- cbind(forecast, var.pred.values)
  }
  h1.ind <- h:ncol(var.err)
  h2.ind <- (h-1):(ncol(var.err)-1)
  h4.ind <- (h-3):(ncol(var.err)-3)
  h8.ind <- 1:(ncol(var.err)-h+1)
  h.ind.len <- length(h1.ind)
  forecast <- list(h1 = forecast[1,h1.ind], h2 = forecast[2,h2.ind],
                  h4 =  forecast[3,h4.ind], h8 = forecast[4,h8.ind])
  for.err <- list(h1 = var.err[1,h1.ind], h2 = var.err[2,h2.ind],
                  h4 =  var.err[3,h4.ind], h8 = var.err[4,h8.ind])
  var.msfe <- list(h1 = sum((for.err$h1)^2)/h.ind.len, h2 = sum((for.err$h2)^2)/h.ind.len,
                   h4 =  sum((for.err$h4)^2)/h.ind.len, h8 =  sum((for.err$h8)^2)/h.ind.len)
  return(structure(list(
    forecast = forecast,
    for.err = for.err,
    var.msfe = var.msfe)))
}


# variable selection ------------------------------------------------------

require(caret)
require(glmnet)
## Returns variables from lasso variable selection, use alpha=0 for ridge
ezlasso <- function(y,x,alpha=0, maxnrvar=10 ){
  # alpha: don't choose to high alpha, otherwise there might stay to few non-zero
  # coefficient values for the individual regressors (alpha<0.5)
  
  myTimeControl <- trainControl(method = "timeSlice",
                                initialWindow = 159,
                                horizon = 1,
                                fixedWindow = FALSE)
  
  glmFitTime <- train(x, y,
                      method = "glmnet",
                      family = "gaussian",
                      trControl = myTimeControl,
                      intercept = FALSE,
                      standardize = TRUE, maxit = 10000,
                      tuneGrid = expand.grid(.alpha = alpha,
                                             .lambda = 10^seq(2,-2,length = 100)))
  
  co <- coef(glmFitTime$finalModel, glmFitTime$bestTune$lambda)
  inds <- (order(co, decreasing = T))[1:maxnrvar]
  variables <- row.names(co)[inds]
  variables <- variables[!(variables %in% '(Intercept)')];
  
  return(c("GDP", variables));
}


## Selects variables with highest mean squared autocorrelation for given lag
# First naive version, without accounting for correlation of included variables
acf.var.selection <- function(y, lag=5, maxnrvar=10){
  y <- na.omit(y)
  set.seed(777)
  acf.fit <- acf(y, lag.max = (lag+1), plot = FALSE)
  # it is assumed that the first variable in the data is of interest
  ac.values.sq <- acf.fit$acf[2:(lag+1),1,]^2
  ms.ac <- apply(ac.values.sq, 2, mean)
  var.selection <- acf.fit$snames[(order(ms.ac, decreasing = T))[1:maxnrvar]]
  if (!any(var.selection=="GDP")){
    var.selection <- c("GDP", var.selection)
  } else {
    var.selection <- var.selection[-(var.selection=="GDP")]
    var.selection <- c("GDP", var.selection)
  }
  return(var.selection)
}

## Selects variables with highest mean squared autocorrelation for given lag
# and accounts for similar patters in the correlation of the lags
acf.var.selection2 <- function(y, lag=10, maxnrvar=10){
  # @ lag: works for lag>1 and variables selection should perform better for
  # hiher lags.
  y <- na.omit(y)
  set.seed(777)
  acf.fit <- acf(y, lag.max = (lag+1), plot = FALSE)
  # it is assumed that the first variable in the data is of interest
  ac.values.sq <- acf.fit$acf[2:(lag+1),1,]^2
  ms.ac <- apply(ac.values.sq, 2, mean)
  var.selection <- acf.fit$snames[(order(ms.ac, decreasing = T))[1]]
  if (any(var.selection=="GDP")){
    var.selection <- acf.fit$snames[(order(ms.ac, decreasing = T))[2]]
  }
  select <- var.selection
  select.ind.total <- c()
  while (length(select.ind.total)<NCOL(y)){
    select.ind <- which(acf.fit$snames==select)
    select.ind.total <- c(select.ind.total, select.ind)
    # Squared distance to last selected variable is calculated
    ac.values.sq.dist <- (acf.fit$acf[2:(lag+1),1,]-acf.fit$acf[2:(lag+1),1,select.ind])^2
    # Values of already selected variables are set to zero
    ac.values.sq.dist[,select.ind.total] <- 0
    ms.ac <- apply(ac.values.sq.dist, 2, mean)
    select <- acf.fit$snames[(order(ms.ac, decreasing = T))[1]]
    if (any(select=="GDP")){
      select <- acf.fit$snames[(order(ms.ac, decreasing = T))[2]]
    }
    # "total variables"/"maxnrvar" of total variables with least scores are removed
    select.ind.total <- c(select.ind.total, (order(ms.ac[-select.ind.total],
                                                   decreasing = F))[1:floor((NCOL(y)/maxnrvar))])
    var.selection <- c(var.selection, select)
    if (length(var.selection)==(maxnrvar-1)){
      break
    }
  }
  return(c("GDP", var.selection))
}


pacf.var.selection <- function(y, lag=8, maxnrvar=10){
  # @ lag: works for 9>lag>1 and variables selection should perform better for
  # higher lags.
  y <- na.omit(y)
  
  set.seed(777)
  # pacf function is limited in the simultaneous calculation...
  pacf.fit <- pacf(y[,1:4], lag.max = (lag+1), plot = FALSE)
  pacf.fit$acf <- pacf.fit$acf[-1,1,]
  for (i in 1:(ceiling(NCOL(y)/4)-1)){
    if (NCOL(y)<(i+1)*4){
      ind <- (i+1)*4-1
    } else {
      ind <- ((1+i)*4)
    }
    y.m <- merge(y[,1], y[,(1+i*4):ind])
    set.seed(777)
    pacf.fit.new <- pacf(y.m, lag.max = (lag+1), plot = FALSE)
    pacf.fit$acf <- cbind(pacf.fit$acf, pacf.fit.new$acf[-1,1,-1])
    pacf.fit$snames <- c(pacf.fit$snames, pacf.fit.new$snames[-1])
  }
  
  # it is assumed that the first variable in the data is of interest
  ac.values.sq <- pacf.fit$acf^2
  ms.ac <- apply(ac.values.sq, 2, mean)
  var.selection <- pacf.fit$snames[(order(ms.ac, decreasing = T))[1]]
  if (any(var.selection=="GDP")){
    var.selection <- pacf.fit$snames[(order(ms.ac, decreasing = T))[2]]
  }
  select <- var.selection
  select.ind.total <- c()
  while (length(select.ind.total)<NCOL(y)){
    select.ind <- which(pacf.fit$snames==select)
    select.ind.total <- c(select.ind.total, select.ind)
    # Squared distance to last selected variable is calculated
    ac.values.sq.dist <- (pacf.fit$acf-pacf.fit$acf[,select.ind])^2
    # Values of already selected variables are set to zero
    ac.values.sq.dist[,select.ind.total] <- 0
    ms.ac <- apply(ac.values.sq.dist, 2, mean)
    select <- pacf.fit$snames[(order(ms.ac, decreasing = T))[1]]
    if (any(select=="GDP")){
      select <- pacf.fit$snames[(order(ms.ac, decreasing = T))[2]]
    }
    # "total variables"/"maxnrvar" of total variables with least scores are removed
    select.ind.total <- c(select.ind.total, (order(ms.ac[-select.ind.total],
                                                   decreasing = F))[1:floor((NCOL(y)/maxnrvar))])
    var.selection <- c(var.selection, select)
    if (length(var.selection)==(maxnrvar-1)){
      break
    }
  }
  return(c("GDP", var.selection))
}

# Tests -------------------------------------------------------------------

### Function for stationarity for p-th lag with Augmented Dickey-Fuller test
aug_dick_fuller <- function(x, crit=0.01){
  oldw <- getOption("warn")
  options(warn = -1)
  var_number <- NCOL(x)
  adf_result <- c()
  for (i in 1:var_number){
    adf_result[i] <- adf.test(na.omit(x[,i]))$p.value                             # From 23rd, because of NAs
  }
  non_stat <- names(variables_diff_quart[0, which(adf_result>crit)])
  options(warn = oldw)
  return(non_stat)                                                                        # List of non-stationary variables       
}

### "MSPE-adjusted" test of Clark and West
CW_test <- function(e1,e2,yf1,yf2,nwlag){
  # Clark-West test
  # Alternative: larger model has smaller MSPE 
  # e1, e2: forecast errors
  # yf1 yf2: forecasts
  # nwlag: Newey-West lag window (forecast step: h)
  # output: CWstat statistic is t distributed with lag degrees of freedom
  # to obtain a p-value of the test simply use: pt(CWstat, df = lag, lower.tail = FALSE)
  
  P <- length(e1)
  froll_adj <- e1^2-( e2^2-(yf1-yf2)^2 )
  varfroll_adj <- nw(froll_adj,nwlag)
  CWstat <- sqrt(P)*(mean(froll_adj))/sqrt(varfroll_adj)
  p_val <- pt(abs(CWstat), df = nwlag, lower.tail = FALSE)
  return(structure(
    list(CWStat = CWstat,
         p.value = p_val)))
} 
  
nw <- function(y,qn){
#input: y is a T*k vector and qn is the truncation lag (forecast step: h)
#output: the newey west HAC covariance estimator
#Formulas are from Hayashi
  t_caps <- NROW(y); k <- NCOL(y); ybar <- matrix(1,t_caps,1)*((sum(y))/t_caps)
  dy <- y-ybar
  G0 <- t(dy)%*%dy/t_caps
  for (j in 1:(qn-1)){
    gamma <- (t(dy[(j+1):t_caps,])%*%dy[1:(t_caps-j),])/(t_caps-1)
    G0 <- G0+(gamma+t(gamma))*(1-abs(j/qn))
  }
  return(G0)
}



##### Diebold-Mariano test statistic
DMtest <- function(d,l){
  
  # Diebold-Mariano test of predictive accuracy, JBES 1995, 13(3), p.253-263
  #
  # H0: equality of forecast accuracy of two forecasts
  # H1: difference of forecast accuracy of two forecasts
  #
  # d:		(T x 1) loss-differential series, usually e1^2-e2^2, where e1 and e2
  #				are the errors of two competing forecasts
  # l:        lag window (forecast step: h)
  # DMstat	Diebold-Mariano test statistic ~ N(0,1) distributed!!!
  # use: pnorm(DMstat, lower.tail = F) to get the p-value of the test!!!
  
  
  t_caps <- length(d)
  m <- mean(d)
  
  # Newey-West variance estimator
  e <- d-m
  lag <- seq(-l,l,1)
  gamma <- vector(length = (2*l+1))
  for (j in 1:(2*l+1)){
    gamma[j] <- t(e[(abs(lag[j])+1):t_caps])%*%e[1:(t_caps-abs(lag[j]))]/t_caps
  }
  weights <- 1-abs(t(lag)/(l+1))
  s2 <- sum(gamma*weights)/t_caps
  # test statistic
  DMstat <- m / sqrt(s2)
  p_val <- pnorm(abs(DMstat), lower.tail = FALSE)
  return(structure(
    list(DMStat = DMstat,
         p.value = p_val)))
}


# Theils U with AR(1) as Benchmark
theils_u_ar1 <- function(y_data, mse_pred, h="h1"){
  # y_data: whole Tx1 data set
  # mse_pred: mspe of the model that AR(1) should be compared to
  fit.ar1 <- ar1_train(na.omit(y_data))
  ar1.rmse <- sqrt(fit.ar1$var.msfe[[h]])
  model.rmse <- sqrt(mse_pred)
  U <- model.rmse/ar1.rmse
  return(U)
}

# Other functions ------------------------------------------------------


# Unscaling of scaled objects function
unscale <- function(x, scaled, columns=0){
  attr <- attributes(scaled)
  scale <- attr$`scaled:scale`
  center <- attr$`scaled:center`
  unscaled.x <- x
  if (columns==0){
    columns = length(scale)
  }
  for (i in 1:columns){
    unscaled.x[,i] <- x[,i] * scale[i] + center[i]
  }
  return(unscaled.x)
}


# Maximal lag to consider because of k>>t restriction for regular VAR(p) models
max.lag <- function(y){
  t <- NROW(y)
  k <- NCOL(y)
  lag <- floor(t/(k+1))-1 
  return(lag)
}


# Transform log differences to initial values
diff_log2norm <- function(logdiffs, init.val, accumulate = T){
  return(Reduce(f = function(x,y){x*exp(y)}, x = logdiffs, init = init.val,
                accumulate = accumulate))
}

