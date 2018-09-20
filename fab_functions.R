# Main function for FAB z-interval ----

#' Numerically calculates the FAB z-interval
#'
#' Uses numerical methods to approximate the Bayes-optimal spending function and endpoints of the FAB z-interval with
#' a normal prior on the mean
#'
#' @param y An observation drawn from \eqn{\text{Normal}(\theta, s2)}
#' @param mu The prior mean for \eqn{\theta = E[y]}
#' @param t2 The prior variance of \eqn{\theta}
#' @param s2 The known variance of y
#' @param alpha One minus the coverage rate
#'
#' @return A vector of length 2, corresponding to the lower and upper bounds of the FAB z-interval
#'

fabzCI <- function(y, mu, t2, s2, alpha = 0.05) {
  if (!is.numeric(y) | length(y) > 1) {
    stop("y must be a numeric scalar")
  }
  if (!is.numeric(mu) | length(mu) > 1) {
    stop("mu must be a numeric scalar")
  }
  if (!is.numeric(t2) | length(t2) > 1 | t2 < 0) {
    stop("t2 must be a positive numeric scalar")
  }
  if (!is.numeric(s2) | length(s2) > 1 | s2 < 0) {
    stop("s2 must be a positive numeric scalar")
  }
  if (!is.numeric(alpha) | length(alpha) > 1) {
    stop("alpha must be a numeric scalar")
  }
  if (alpha <= 0 | alpha >= 1) {
    stop("alpha must be between 0 and 1")
  }
  s <- sqrt(s2)
  ubroot <- function(theta) {
    (y + s * qnorm(min(1, 1 - alpha + pnorm((y - theta)/s))) + 
       mu * 2 * s2/t2)/(1 + 2 * s2/t2) - theta
  }
  a <- b <- y + s * qnorm(1 - alpha)
  while (ubroot(a) < 0) {
    a <- a - 1e-12
  }
  while (ubroot(b) > 0) {
    b <- b + s * qnorm(1 - alpha) * 0.25
  }
  thetaU <- uniroot(ubroot, c(a, b))$root
  lbroot <- function(theta) {
    (y + s * qnorm(max(0, alpha - pnorm((theta - y)/s))) + 
       mu * 2 * s2/t2)/(1 + 2 * s2/t2) - theta
  }
  a <- b <- y + s * qnorm(alpha)
  while (lbroot(a) < 0) {
    a <- a + s * qnorm(alpha) * 0.25
  }
  while (lbroot(b) > 0) {
    b <- b + 1e-12
  }
  thetaL <- uniroot(lbroot, c(a, b))$root
  c(thetaL, thetaU)
}

# Helper functions for FAB t-interval ----
paccept <- function(w, theta, mu, t2, nu0, s20, n, alpha){
  a <- nu0 / 2; b <- nu0 * s20 / 2

  f<-function(s2){
    c <- sqrt(s2 / n) / sqrt(s2 / n + t2)
    ncp<- c * (mu - theta) / (sqrt(s2 / n))
    suppressWarnings(
    pacc <- pt(c * qt(1 - alpha * (1 - w), n - 1), n - 1, ncp) -
           pt(c * qt(alpha * w, n - 1), n - 1, ncp)
    )

    a <- nu0 / 2 ; b <- nu0 * s20 / 2
    ds2 <- exp(a * log(b) - lgamma(a) - (a + 1) * log(s2) - b / s2 )

    pacc*ds2
  }

  int<-try(integrate(f,lower=0,upper=Inf)$val,silent=TRUE)
  if(is.character(int)){int<-1}
  int
}

wfabt<-function(theta, mu, t2, nu0, s20, n, alpha){
  optimize(paccept, lower=0, upper=1,
            theta = theta, mu = mu, t2 = t2, nu0 = nu0, s20 = s20, n = n, alpha = alpha)$min
}

# Main t-interval function ----
#' Approximates the FAB t-interval
#'
#' Uses numerical methods to approximate the Bayes-optimal spending function and endpoints of the FAB t-interval with
#' a normal prior on the mean and an inverse-gamma prior on the variance
#'
#' @param ybar The sample mean of \eqn{\{y_i, \cdots y_n\}}
#' @param s A consistent estimate of \eqn{\sigma} the standard deviation of y, not assumed known
#' @param n The sample size
#' @param mu The prior mean for \eqn{\theta = E[y]}
#' @param t2 The prior variance of \eqn{\theta}
#' @param nu0 parameter of inverse-gamma prior for \eqn{\sigma^2}, where \eqn{1/\sigma2 \sim \text{IG}(nu_0/2, s_0^2/2)}
#' @param s20 parameter of inverse-gamma prior for \eqn{\sigma^2}
#' @param alpha One minus the coverage rate
#'
#' @return A vector of length 2, corresponding to the lower and upper bounds of the FAB t-interval
#'
fabt_ci <- function(ybar, s, n, mu, t2, nu0, s20, alpha = 0.05){
  
  ubroot<-function(theta){
    w <- wfabt(theta, mu, t2, nu0, s20, n, alpha)
    ybar + s * qt(1 - alpha * w, n - 1) / sqrt(n) - theta
  }
  a <- b <- ybar + 0.99 * (s / sqrt(n)) * qt(1 - alpha, n - 1)
  while (ubroot(b) > 0) {
    b <- b + (s / sqrt(n)) * qnorm(1 - alpha) * n / (n + 4) 
  }
  thetaU <- uniroot(ubroot,c(a, b))$root
  
  lbroot<-function(theta){
    w<-wfabt(theta, mu, t2, nu0, s20, n, alpha)
    ybar + s * qt(alpha * (1 - w), n - 1)/ sqrt(n) - theta
  }
  a <- b <- ybar + 0.99 * (s / sqrt(n)) * qt(alpha, n - 1)
  while (lbroot(a) < 0) {
    a <- a + (s / sqrt(n)) * qnorm(alpha) * n / (n + 4) 
    }
  thetaL <- uniroot(lbroot, c(a, b))$root
  
  c(thetaL, thetaU)
}

# Other Miscellaneous Functions ----
wfabz<-function(theta,mu,t2,alpha=.05)
{
  igfun<-function(x,alpha)
  {
    gwmx <-function(w){ qnorm(alpha*w) - qnorm(alpha*(1-w)) - x }
    uniroot(gwmx,interval=c(0,1),maxiter=2000,tol=.Machine$double.eps^0.5)$root
  }
  igfun( 2*(theta-mu)/t2,alpha)
}

# Expected width of FAB z-interval
fabz_ew<-function(theta, sigma2, mu, t2,alpha=.05){

  hx<-function(xs, alpha = 0.05){
    hxs<-NULL
    for (x in xs) {
    ih_mx <- function(w){
      qnorm(alpha*w) - qnorm(alpha*(1-w)) - x
      }
    hxs <- c(hxs, uniroot(ih_mx, interval = c(0, 1), maxiter=2000,
               tol = .Machine$double.eps ^ 0.5)$root)
    }
    hxs
  }

  H1<-function(x, theta, t2,alpha = 0.05) {
    h1 <- NULL
    for (xi in x) {
      h1<-c(h1, .5 * t2 * pnorm(theta - .5 * t2 * xi + qnorm(1 - alpha * hx(xi, alpha))))
    }
    h1
  }

  H2<-function(x, theta, t2, alpha = 0.05){
    h2 <- NULL
    for (xi in x) {
      h2 <- c(h2, .5 * t2 * pnorm(theta - .5 * t2 * xi - qnorm(1 - alpha * hx(-xi, alpha))))
    }
    h2
  }

  H<-function(x){ H1(x, theta - mu, t2, alpha) - H2(x, theta - mu, t2, alpha) }
  sqrt(sigma2) * integrate(H, lower= -Inf, upper=Inf,
            subdivisions = 800L)$value
}

# Expected reduction in interval width when $\sigma^2 = 1$
width_reduction_z <- function(mu, t2, alpha = .05){
  direct_width <- 2 * qnorm(1 - alpha/2)
  prior_dens <- function(theta, mu, t2){
    dnorm(theta, mu, sqrt(t2))
  }
  f <- function(theta){
    prior_dens(theta, mu, t2) * sapply(theta, function(x){
      fabz_ew(x, 1, mu, t2, alpha)})
  }
  (integrate(f, lower= -Inf, upper=Inf, subdivisions = 100L)$value) / direct_width
}











