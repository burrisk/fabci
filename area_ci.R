#' Log-likelhood evaluation of sampling variances
#'
#' Evaluates the negative log-likelihood, assuming \eqn{(\frac{n_j - 1)/s_j^2}{\sigma_j^2} \sim 
#' \Chi^2_{n_j - 1}} and \eqn{1/\sigma_j^2 \sim 
#' \text{Gamma}(a, b)}
#'
#' @param theta A vector of length 2, corresponding to \code{c(a, b)}
#' @param samp_size A vector of length m, corresponding to the sample sizes for each of the groups
#' @param S A vector of length m, corresponding to plug-in estimates of the variance
#'
#' @return The negative log-likelihood
#'
ab_f <- function(theta, samp_size, S){
  if (length(S) != length(samp_size)){
    stop("You must specify the sample size for each group.")
  }
  m <- length(S); a <- theta[1]; b <- theta[2]
  omega <- (samp_size - 1) / 2
  -sum(lgamma(omega + a) + a * log(b) -  (omega + a) * log(omega * S + b) - lgamma(a))
}

#' Log-likelihood gradient evaluation of sampling variances
#'
#' Evaluates the negative log-likelihood, assuming \eqn{(\frac{n_j - 1)/s_j^2}{\sigma_j^2} \sim 
#' \Chi^2_{n_j - 1}} and \eqn{1/\sigma_j^2 \sim 
#' \text{Gamma}(a, b)}
#'
#' @param theta A vector of length 2, corresponding to \code{c(a, b)}
#' @param samp_size A vector of length m, corresponding to the sample sizes for each of the groups
#' @param S A vector of length m, corresponding to plug-in estimates of the variance
#'
#' @return The gradient of the negative log-likelihood
#'
ab_g <- function(theta, samp_size, S){
  if (length(S) != length(samp_size)){
    stop("You must specify the sample size for each group.")
  }
  m <- length(S); a <- theta[1]; b <- theta[2]
  omega <- (samp_size - 1) / 2
  gradient <- rep(0, 2)
  gradient[1] <- -sum(digamma(omega + a) + log(b)  - log(omega * S + b) - digamma(a))
  gradient[2] <- -sum(a / b - (omega + a) / (omega * S + b))
  gradient
}

#' Maximum likelihood estimates for sampling variance hyperparameters
#'
#'  Finds the MLEs \eqn{\hat{a}} and \eqn{\hat{b}}, assuming  \eqn{(\frac{n_j - 1)/s_j^2}{\sigma_j^2} \sim 
#' \Chi^2_{n_j - 1}} and \eqn{1/\sigma_j^2 \sim 
#' \text{Gamma}(a, b)}.  Used in construction of FAB t-intervals
#'
#' @param S Plug-in estimates of the sampling variances of each area mean
#' @param samp_size The sample size for each group
#'
#' @return The maximum likelihood estimates of a and b
#'
optim_ab <- function(S, samp_size){
  S <- S * samp_size
  # Method of moments initialization
  mg <- mean(1 / S); vg <- var(1 / S) 
  a <- mg ^ 2 / vg;  b <- mg / vg
  return(optim(par = c(a, b), S = S, samp_size = samp_size, fn = ab_f, gr = ab_g,
               method = "L-BFGS-B", lower = c(1e-6, 1e-6),
               upper = c(1e4, 1e4))$par)
}

a <- 3
b <- 5
N <- 100000
samp_size <- sample(5:100, N, replace = T)
phi <- rgamma(N, a, b)
x <- rchisq(N, samp_size - 1)
s <- x / (phi * (samp_size - 1))

optim_ab(s, samp_size)


#' An iterative procedure to find MLEs of hyperparameters for the spatial Fay-Herriot model
#'
#'  Finds the MLEs \eqn{\rho} and \eqn{\tau^2}, which are hyperparameters of a SAR Fay-Herriot model
#'
#' @param y Direct estimates of small area means
#' @param X The design matrix of covariates; should include a column of ones
#' @param W A row-standardized proximity matrix to represent spatial dependence
#' @param direct_var A vector of known sampling variances for each area
#' @param maxiter Maximum number of iterations until stopping, defaults to 100
#' @param precision Precision needed for convergence to be assessed  
#'
#' @return A list consisting of the MLEs of \eqn{\rho} and \eqn{\tau^2}, as well as the number of iterations
#' until convergence
#'
fisherScoringSpatial <- function(y, X, W, direct_var, maxiter = 100, precision = 1e-4){
  m <- length(y); yt <- t(y); W_t <- t(W); I <- diag(m); Xt = t(X)
  par_old <- matrix(0, 2, 1); par_new <- matrix(0, 2, 1)
  gradient <- matrix(0, 2, 1); hessian <- matrix(0, 2, 2)
  tau2_path <- rep(0, maxiter); rho_path <- rep(0, maxiter); log_lik <- rep(-Inf, maxiter)
  tau2_path[1] <- median(direct_var); rho_path[1] <- 0
  k <- 0; diff.S <- precision + 1
  while ((diff.S > precision) & (k < maxiter)) {
    k <- k + 1
    dtau2 <- solve((I - rho_path[k] * W_t) %*% (I - rho_path[k] * W), 
                   tol = 1e-17)
    dRho <- 2 * rho_path[k] * W_t %*% W -  W - W_t
    dVRho <- (-1) * tau2_path[k] * (dtau2 %*% dRho %*% dtau2)
    V <- tau2_path[k] * dtau2 + diag(direct_var); Vi <- solve(V)
    XtVi <- Xt %*% Vi
    Q <- solve(XtVi %*% X)
    P <- Vi - t(XtVi) %*% Q %*% XtVi
    beta <- Q %*% XtVi %*% y
    log_lik[k] <- dmvnorm(y, X %*% beta, V, log = T)
    PD <- P %*% dtau2; PR <- P %*% dVRho; Pdir <- P %*% y
    ViD <- Vi %*% dtau2
    ViR <- Vi %*% dVRho
    gradient[1, 1] <- (-0.5) * sum(diag(ViD)) + (0.5) * (yt %*% 
                                                           PD %*% Pdir)
    gradient[2, 1] <- (-0.5) * sum(diag(ViR)) + (0.5) * (yt %*% 
                                                           PR %*% Pdir)
    hessian[1, 1] <- (0.5) * sum(diag(ViD %*% ViD))
    hessian[1, 2] <- (0.5) * sum(diag(ViD %*% ViR))
    hessian[2, 1] <- hessian[1, 2]
    hessian[2, 2] <- (0.5) * sum(diag(ViR %*% ViR))
    par_old[1, 1] <- tau2_path[k]
    par_old[2, 1] <- rho_path[k]
    step_size <- 1; feasibleStep = F;
    while (!(feasibleStep)){
      par_new <- par_old + step_size * solve(hessian) %*% gradient
      if (par_new[1, 1] <= 0.001 | par_new[2, 1] <= -0.999 | par_new[2, 1] >= 0.999){
        step_size <- step_size / 2
      } else{
        feasibleStep = T
      }
    }
    tau2_path[k + 1] <- par_new[1, 1]
    rho_path[k + 1] <- par_new[2, 1]
    diff.S <- max(abs(par_new - par_old) / par_old)
  }
  k = which.max(log_lik)
  list(rho = rho_path[k+1], tau2 = tau2_path[k+1], iterations = k)
}

#' 
#'
#' Approximates the mean-squared-error of the EBLUP estimate in the spatial Fay-Herriot model, using the formula
#' given in Rao and Molina (2015).
#'
#' @param y Direct estimates of small area means
#' @param X The design matrix of covariates; should include a column of ones
#' @param W A row-standardized proximity matrix to represent spatial dependence
#' @param direct_var A vector of known sampling variances for each area
#' @param rho The maximum likelihood estimate of \eqn{\rho}, the spatial autocorrelation parameter
#' @param tau2 The maximum likelihood estimate of \eqn{\tau^2}, the dispersion parameter 
#'
#' @return A list consisting of EBLUP estimates, MSE of the EBLUP estimates, and the estimated regression coefficients
#' \eqn{\beta}
#'
mseSpatial <- function(y, X, W, direct_var, rho, tau2){
  m <- dim(X)[1]; p <- dim(X)[2]
  g1d <- rep(0, m); g2d <- rep(0, m); g3d <- rep(0, m); g4d <- rep(0, m)
  mse2d <- rep(0, m); eblup <- rep(0, m); I <- diag(m); Xt <- t(X); Wt <- t(W)
  Ci <- solve((I - rho * Wt) %*% (I - rho * W)); G <- tau2 * Ci
  V <- G + diag(direct_var); Vi <- solve(V)
  XtVi <- Xt %*% Vi
  Q <- solve(XtVi %*% X)
  beta <- Q %*% XtVi %*% y; Xbeta <- X %*% beta
  res <- y - Xbeta
  eblup <- Xbeta + G %*% Vi %*% res
  Ga <- G - G %*% Vi %*% G; Gb <- G %*% Vi %*% X
  Xa <- matrix(0, 1, p)
  for (i in 1:m) {
    g1d[i] <- Ga[i, i]
    Xa[1, ] <- X[i, ] - Gb[i, ]
    g2d[i] <- Xa %*% Q %*% t(Xa)
  }
  dRho <- 2 * rho * Wt %*% W - W - Wt
  tau2mat <- (-1) * tau2 * (Ci %*% dRho %*% Ci)
  P <- Vi - t(XtVi) %*% Q %*% XtVi
  PCi <- P %*% Ci
  Ptau2mat <- P %*% tau2mat
  hessian <- matrix(0, 2, 2)
  hessian[1, 1] <- (0.5) * sum(diag((PCi %*% PCi)))
  hessian[1, 2] <- (0.5) * sum(diag((PCi %*% Ptau2mat)))
  hessian[2, 1] <- hessian[1, 2]
  hessian[2, 2] <- (0.5) * sum(diag((Ptau2mat %*% Ptau2mat)))
  hessian_inv <- solve(hessian)
  ViCi <- Vi %*% Ci
  Vitau2mat <- Vi %*% tau2mat
  l1 <- ViCi - tau2 * ViCi %*% ViCi
  l1t <- t(l1)
  l2 <- Vitau2mat - tau2 * Vitau2mat %*% ViCi
  l2t <- t(l2)
  L <- matrix(0, 2, m)
  for (i in 1:m) {
    L[1, ] <- l1t[i, ]
    L[2, ] <- l2t[i, ]
    g3d[i] <- sum(diag(L %*% V %*% t(L) %*% hessian_inv))
  }
  QXtVi <- Q %*% XtVi
  ViX <- Vi %*% X
  h1 <- (-1) * sum(diag(QXtVi %*% Ci %*% ViX))
  h2 <- (-1) * sum(diag(QXtVi %*% tau2mat %*% ViX))
  h <- matrix(c(h1, h2), nrow = 2, ncol = 1)
  bML <- (hessian_inv %*% h)/2
  tbML <- t(bML)
  GVi <- G %*% Vi
  GViCi <- GVi %*% Ci
  GVitau2mat <- GVi %*% tau2mat
  ViCi <- Vi %*% Ci
  dg1_dtau2 <- Ci - 2 * GViCi + tau2 * GViCi %*% ViCi
  dg1_dp <- tau2mat - 2 * GVitau2mat + tau2 * GVitau2mat %*% ViCi
  gradg1d <- matrix(0, nrow = 2, ncol = 1)
  bMLgradg1 <- rep(0, m)
  for (i in 1:m) {
    gradg1d[1, 1] <- dg1_dtau2[i, i]
    gradg1d[2, 1] <- dg1_dp[i, i]
    bMLgradg1[i] <- tbML %*% gradg1d
  }
  mse2d <- g1d + g2d + 2 * g3d - bMLgradg1
  list(est = eblup, mse = mse2d, beta = beta)
}

#' An iterative procedure to find MLEs of hyperparameters for the vanilla Fay-Herriot model
#'
#'  Finds the MLE of \eqn{\tau^2}, the dispersion hyperparameter in the Fay-Herriot model
#'
#' @param y Direct estimates of small area means
#' @param X The design matrix of covariates; should include a column of ones
#' @param direct_var A vector of known sampling variances for each area
#' @param maxiter Maximum number of iterations until stopping, defaults to 100
#' @param precision Precision needed for convergence to be assessed  
#'
#' @return A list consisting of the MLE of \eqn{\tau^2}
#'
fisherScoringNonSpatial <- function(y, X, direct_var, maxiter = 100, precision = 1e-4){
  m <- length(y); Xt <- t(X); tau2_path <- rep(0, maxiter)
    tau2 <- 0; tau2_path[1] <- median(direct_var)
    k <- 0; diff <- precision + 1
    while ((diff > precision) & (k < maxiter)) {
      k <- k + 1
      Vi <- 1/(tau2_path[k] + direct_var)
      XtVi <- t(Vi * X)
      Q <- solve(XtVi %*% X)
      P <- diag(Vi) - t(XtVi) %*% Q %*% XtVi
      Py <- P %*% y
      gradient <- (-0.5) * sum(Vi) + 0.5 * (t(Py) %*% Py)
      hessian <- 0.5 * sum(Vi^2)
      step_size <- 1; feasibleStep = F;
      while (!(feasibleStep)){
        tau2_path[k + 1] <- tau2_path[k] + step_size * gradient/hessian
        if (tau2_path[k + 1] <= 0.001){
          step_size <- step_size / 2
        } else{
          feasibleStep = T
        }
      }
      diff <- abs((tau2_path[k + 1] - tau2_path[k])/tau2_path[k])
    }
    list(tau2 = tau2_path[k+1])
}

#' 
#'
#' Approximates the mean-squared-error of the EBLUP estimate in the vanilla Fay-Herriot model, using the formula
#' given in Rao and Molina (2015).
#'
#' @param y Direct estimates of small area means
#' @param X The design matrix of covariates; should include a column of ones
#' @param direct_var A vector of known sampling variances for each area
#' @param tau2 The maximum likelihood estimate of \eqn{\tau^2}, the dispersion parameter 
#'
#' @return A list consisting of EBLUP estimates, MSE of the EBLUP estimates, and the estimated regression coefficients
#' \eqn{\beta}
#'
mseNonSpatial <- function(y, X, direct_var, tau2){
  m <- dim(X)[1]; p <- dim(X)[2]
  g1d <- rep(0, m); g2d <- rep(0, m); g3d <- rep(0, m); 
  eblup <- rep(0, m); mse2d <- rep(0, m)
  Vi <- 1/(tau2 + direct_var)
  Bd <- direct_var/(tau2 + direct_var)
  trV <- sum(Vi^2)
  XtVi <- t(Vi * X)
  Q <- solve(XtVi %*% X)
  beta <- Q %*% XtVi %*% y; Xbeta <- X %*% beta
  res <- y - Xbeta
  eblup <- Xbeta + tau2 * Vi * res
  Vartau2 <- 2/trV
  b <- (-1) * sum(diag(Q %*% (t((Vi^2) * X) %*% X)))/trV
  for (d in 1:m) {
    g1d[d] <- direct_var[d] * (1 - Bd[d])
    xd <- matrix(X[d, ], nrow = 1, ncol = p)
    g2d[d] <- (Bd[d]^2) * xd %*% Q %*% t(xd)
    g3d[d] <- (Bd[d]^2) * Vartau2/(tau2 + direct_var[d])
    mse2d[d] <- g1d[d] + g2d[d] + 2 * g3d[d] - b * (Bd[d]^2)
  }
  list(est = eblup, mse = mse2d, beta = beta)
}

#' Calculates a spatial Fay-Herriot FAB interval for a given area mean 
#'
#' Uses data from other groups to adaptively estimate the hyperparameters of a SAR Fay-Herriot linking model,
#' leading to a FAB interval
#'
#' @param y Direct estimates of small area means
#' @param X The design matrix of covariates; should include a column of ones
#' @param W A row-standardized proximity matrix to represent spatial dependence
#' @param direct_var A vector of sampling variances (either known or plug-in estimates) for each area
#' @param index The index of the area of interest
#' @param alpha \eqn{1 - \alpha} is the coverage probability of the interval
#' @param varknown A logical indicating whether the sampling variances are known 
#' @param sample_size The sample size for each area
#' @param maxiter The maximum number of iterations for the Fisher scoring algorithm
#' @param precision The precision needed for algorithm to converge
#'
#' @return A vector of length 2, corresponding to the FAB z or t interval for the area mean of the specified index
#'
fabSFH <- function(y, X, W, direct_var, index, alpha = 0.05, varknown = T,
                   sample_size = NULL, maxiter = 100, 
                   precision = 1e-4){
  y_ind <- y[-index]; X_ind <- as.matrix(X[-index, ]); v_ind <- direct_var[-index]
  W_ind <- W[-index, -index]; I <- diag(length(y))
  if (abs(sum(rowSums(W_ind) - 1)) > 1e-6){
    W_ind <- sweep(W_ind, 1, rowSums(W_ind), "/")
  }
  maxlik_hyper <- fisherScoringSpatial(y_ind, X_ind, W_ind, v_ind, maxiter,
                                       precision)
  rho <- maxlik_hyper$rho; tau2 <- maxlik_hyper$tau2
  eb_est <- mseSpatial(y_ind, X_ind, W_ind, v_ind, rho, tau2)
  eblup <- eb_est$est; mse <- eb_est$mse; beta <- eb_est$beta
  marg_mean <- X %*% beta
  marg_var <- tau2 * solve((I - rho * t(W)) %*% (I - rho * W)) + diag(direct_var)
  Sigma12 <- matrix(marg_var[index, -index], nrow = 1)
  Sigma22_inv <- solve(marg_var[-index, -index])
  cond_mean <- marg_mean[index] + Sigma12 %*% Sigma22_inv %*%
    matrix(eblup - marg_mean[-index], ncol = 1)
  cond_var <- marg_var[index, index] - Sigma12 %*% Sigma22_inv %*% t(Sigma12) 
  if (varknown){
    interval <- fabzCI(y[index], cond_mean[1, 1], cond_var[1, 1],
                       direct_var[index], alpha = alpha)
  } else{
    ab <- optim_ab(direct_var[-index], sample_size[-index])
    nu0 <- ab[1] * 2
    s20 <- (ab[2] / ab[1])
    interval <- fabt_ci(y[index], sqrt(sample_size[index] * direct_var[index]), sample_size[index], 
            cond_mean[1, 1], cond_var[1, 1], nu0, s20, alpha = alpha)
  }
  interval
}

#' Calculates a Fay-Herriot FAB interval for a given area mean 
#'
#' Uses data from other groups to adaptively estimate the hyperparameters of a Fay-Herriot linking model,
#' leading to a FAB interval
#'
#' @param y Direct estimates of small area means
#' @param X The design matrix of covariates; should include a column of ones
#' @param direct_var A vector of sampling variances (either known or plug-in estimates) for each area
#' @param index The index of the area of interest
#' @param alpha \eqn{1 - \alpha} is the coverage probability of the interval
#' @param varknown A logical indicating whether the sampling variances are known 
#' @param sample_size The sample size for each area
#' @param maxiter The maximum number of iterations for the Fisher scoring algorithm
#' @param precision The precision needed for algorithm to converge
#'
#' @return A vector of length 2, corresponding to the FAB z or t interval for the area mean of the specified index
#'
fabFH <- function(y, X, direct_var, index, alpha = 0.05, varknown = T,
                  sample_size = NULL, maxiter = 100, precision = 1e-4){
  y_ind <- y[-index]; X_ind <- as.matrix(X[-index, ]); v_ind <- direct_var[-index]
  maxlik_hyper <- fisherScoringNonSpatial(y_ind, X_ind, v_ind, maxiter,
                                       precision)
  tau2 <- maxlik_hyper$tau2
  eb_est <- mseNonSpatial(y_ind, X_ind, v_ind, tau2)
  eblup <- eb_est$est; mse <- eb_est$mse; beta <- eb_est$beta
  marg_mean <- X %*% beta
  marg_var <- tau2 + direct_var
  cond_mean <- marg_mean[index]
  cond_var <- marg_var[index]
  if (varknown){
    interval <- fabzCI(y[index], cond_mean, cond_var,
                       direct_var[index], alpha = alpha)
  } else{
    ab <- optim_ab(direct_var[-index], sample_size[-index])
    nu0 <- ab[1] * 2
    s20 <- (ab[2] / ab[1])
    interval <- fabt_ci(y[index], sqrt(sample_size[index] * direct_var[index]), sample_size[index], 
                        cond_mean, cond_var, nu0, s20, alpha = alpha)
  }
  interval
}

#' Model-based approach to calculate confidence intervals for each area
#' 
#' Calculates FAB, Empirical Bayes, and Direct confidence intervals for all area means, enabling the use of
#'  either a Fay-Herriot or a spatial Fay-Herriot linking model 
#'
#'
#' @param y A object of class formula; a symbolic description of the model to be fitted
#' @param data The dataset used for the model
#' @param direct_var A vector of sampling variances (either known or plug-in estimates) for each area
#' @param alpha \eqn{1 - \alpha} is the coverage probability of the interval
#' @param method One of "FAB", "EB", or "Direct", the type of confidence interval desired
#' @param type One of "FH" or "SAR", corresponding to the linking model used for interval construction
#' @param varknown A logical indicating whether the sampling variances are known 
#' @param sample_size The sample size for each area
#' @param maxiter The maximum number of iterations for the Fisher scoring algorithm
#' @param precision The precision needed for algorithm to converge
#'
#' @return A \eqn{n \times 2} matrix of confidence intervals, where n is the number of areas
#' 
#' TODO: Error handling
#' 
#' @export
areaCI <- function(formula, data, direct_var, alpha = 0.05, 
                   method = c("FAB", "EB", "Direct"), type = c("FH", "SAR"),
                   W = NULL, varknown = T, 
                   sample_size = NULL, maxiter = 100, precision = 1e-4){
  formuladata <- model.frame(formula, na.action = na.omit, data)
  y <- formuladata[, 1]
  X <- model.matrix(formula, data)
  m <- nrow(data)
  
  if (varknown){
    q_star <- qnorm(alpha / 2)
  } else{
    q_star <- qt(alpha / 2, df = sample_size - 1)
  }
  
  if (method == "FAB"){
    if (type == "FH"){
      result <- lapply(1:m, function(i){
        fabFH(y, X, direct_var, i, alpha, varknown, sample_size,
              maxiter, precision)
      })
      result <- do.call(rbind, result)
    } else if (type == "SAR"){
      result <- lapply(1:m, function(i){
        fabSFH(y, X, W, direct_var, i, alpha, varknown, sample_size,
              maxiter, precision)
      })
      result <- do.call(rbind, result)
    } else{
      stop("Type must be FH or SAR")
    }
  } else if (method == "EB"){
    result <- matrix(0, nrow = m, ncol = 2)
    if (type == "FH"){
      maxlik_hyper <- fisherScoringNonSpatial(y, X, direct_var, maxiter,
                                           precision)
      tau2 <- maxlik_hyper$tau2
      eb_est <- mseNonSpatial(y, X, direct_var, tau2)
      eblup <- eb_est$est; mse <- eb_est$mse
      result[,1] <- eblup + q_star * sqrt(mse)
      result[,2] <- eblup - q_star * sqrt(mse)
    } else if (type == "SAR"){
      maxlik_hyper <- fisherScoringSpatial(y, X, W, direct_var, maxiter,
                                              precision)
      tau2 <- maxlik_hyper$tau2; rho <- maxlik_hyper$rho
      eb_est <- mseSpatial(y, X, W, direct_var, rho, tau2)
      eblup <- eb_est$est; mse <- eb_est$mse
      result[,1] <- eblup + q_star * sqrt(mse)
      result[,2] <- eblup - q_star * sqrt(mse)
    } else{
      stop("Type must be FH or SAR")
    }
  } else if (method == "Direct"){
      result <- matrix(0, nrow = m, ncol = 2)
      result[,1] <- y + q_star * sqrt(direct_var)
      result[,2] <- y - q_star * sqrt(direct_var)
  } else{
      stop("Method must be FAB, EB, or Direct")
  }
  colnames(result) <- c("lower", "upper")
  return(result)
}


