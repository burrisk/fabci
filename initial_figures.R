library(tidyverse)

coverage <- function(theta, mu, sigma2, tau2, alpha = 0.05){
  z_star <- qnorm(1 - alpha/2)
  gamma <- tau2 / (sigma2 + tau2)
  lower <- (((theta - (1 - gamma) * mu - z_star * sqrt(sigma2 * gamma))/ gamma) - theta) / sqrt(sigma2)
  upper <- (((theta - (1 - gamma) * mu + z_star * sqrt(sigma2 * gamma))/ gamma) - theta) / sqrt(sigma2)
  pnorm(upper) - pnorm(lower)
}

coverage2 <- function(theta, mu, sigma2, tau2, alpha = 0.05){
  sigma <- sqrt(sigma2)
  tau <- sqrt(tau2)
  z_star <- qnorm(1 - alpha / 2)
  a <- sigma / tau2 * (theta - mu)
  b <- sqrt(1 + sigma2 / tau2)
  pnorm(a + z_star * b) - pnorm(a - z_star * b)
}

theta_seq <- seq(-5, 5, length = 200)
y_theta <- sapply(theta_seq, function(theta){
  coverage(theta, mu = 0, sigma2 = 1, tau2 = 1, alpha = 0.05)
})
y2_theta <- sapply(theta_seq, function(theta){
  coverage2(theta,  mu = 0, sigma2 = 1, tau2 = 1, alpha = 0.05)
})

# Figure 1
coverage_plot <- tibble(theta = theta_seq, coverage = y_theta)
ggplot(coverage_plot, aes(x = theta, y = coverage)) +
  geom_line(aes(lty = "Bayes")) +
  labs(x = expression(theta[j] - x[j]^T * beta), y = "Coverage Probability") +
  geom_hline(aes(yintercept = 0.95, lty = "Direct")) +
  scale_linetype_manual(values = c("dotted", "dashed"), name = "Interval Type")

ggsave("figures/bayes_cov.pdf", device = "pdf", width = 6, height = 5)

# Figure 2
source("code/fab_functions.R")
ltau_seq <- seq(-2, 5, length = 20)
exp_length <- sapply(ltau_seq, function(ltau){
  width_reduction_z(0, exp(ltau))
})
length_plot <- tibble(tau = ltau_seq, length = exp_length)
ggplot(length_plot, aes(x = tau, y = length)) +
  geom_line(aes(lty = "FAB")) +
  labs(x = expression(log(tau[j]^2/sigma[j]^2)), y = "Relative Expected Length") +
  geom_hline(aes(yintercept = 1, lty = "Direct")) +
  scale_linetype_manual(values = c("dotted", "solid"), name = "Interval Type")

ggsave("figures/tau_width.pdf", device = "pdf", width = 6, height = 5)
