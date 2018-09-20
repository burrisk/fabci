#' Create a neighborhood matrix
#'
#' Creates a binary contiguity neighborhood matrix on a \eqn{n \times m} lattice
#'
#' @param n The number of rows in the lattice
#' @param m The number of columns in the lattice
#'
#' @return A \eqn{nm \times nm} matrix H, such that \eqn{H_{ii}} corresponds to the number of neighbors for
#'  vertex \eqn{i}, and \eqn{-H_{ij}} is a binary indicator for whether vertex i and j are neighbors.
#'
#' @examples
#' H <- createNeighborHoodMat(5, 4)
#' dim(H)
createNeighborhoodMat <- function(n, m){
  nMat <- matrix(0, nrow = n * m, ncol = n * m)
  for (i in 1:(n * m)){
    if (i %% m != 0){
      nMat[i, i + 1] = -1
      nMat[i + 1, i] = -1
      nMat[i, i] = nMat[i, i] + 1
      nMat[i + 1, i + 1] = nMat[i + 1, i + 1] + 1
    }
    if (ceiling(i / m) != n){
      nMat[i, i + m ] = -1
      nMat[i + m, i] = -1
      nMat[i, i] = nMat[i, i] + 1
      nMat[i + m, i + m] = nMat[i + m, i + m] + 1
    }
  }
  nMat
}

#' Simulate data from the Spatial Fay-Herriot model
#'
#' Simulates a dataset of observed area means via the spatial Fay-Herriot model, where the underlying linking model is a 
#' Simultaneous Autoregressive (SAR) Model, the sampling variances are known, and auxillary covariates are uniformly
#'  generated and then standardized
#'
#' @param p The number of areas in the dataset
#' @param W A row standardized proximity matrix with diagonal elements equal to zero 
#' @param beta A vector of regression coefficients 
#' @param intercept The intercept for use in the Fay-Herriot Model
#' @param rho A value in \eqn(-1, 1) that specifies the level of spatial autocorrelation
#' @param tau2 The dispersion of the random effects
#' @param var The known (length-p) vector of sampling variances
#'
#' @return A \eqn{p \times (d + 4)} dataset consisting of auxillary covariates, direct estimates of area means, 
#' the sampling variances, the sample sizes, and the true area means
#'
#' @examples
#' rows <- 7
#' cols <- 7
#' d <- 1
#' p <- rows*cols
#' R <- createNeighborhoodMat(rows, cols)
#' W <- R
#' diag(W) <- 0
#' W <- sweep(W, 1, rowSums(W), "/")
#' genDataSAR(p = p, W = W, beta = c(1, -1), intercept = 0, rho = 0.8, tau2 = 10, var = rep(1, p))
#' @export
genDataSAR <- function(p, W, beta, intercept, rho = 0, tau2 = 1, var = rep(1, p)){
  d <- length(beta)
  X <- matrix(scale(runif(p * d)), ncol = d)
  mean <- X %*% matrix(beta, ncol = 1) + intercept
  iminusw <- diag(p) - rho * W
  theta <- rmvnorm(1, mean, tau2 * solve(iminusw %*% t(iminusw)))
  theta_mat <- t(theta)
  Y <- rnorm(p, as.vector(theta_mat), sqrt(var))
  data <- data.frame(X, Y, var, as.numeric(theta))
  colnames(data) <- c(paste0("X", 1:d), "dir_est", "var_est", "theta")
  return(data)
}
