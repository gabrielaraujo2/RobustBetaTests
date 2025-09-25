#***Score vector for MLE****
score_beta_regression <- function(y, X, Z, theta_hat, linkmu = "logit", linkphi = "log") {
  n <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  beta <- theta_hat[1:p]
  gamma <- theta_hat[(p + 1):(p + q)]
  mean_info <- mean_link(beta = beta, y = y, X = X, linkmu = linkmu)
  mu <- mean_info$mu
  T_1_diag <- diag(mean_info$T_1)
  precision_info <- precision_link(gama = gamma, Z = Z, linkphi = linkphi)
  phi <- precision_info$phi
  T_2_diag <- diag(precision_info$T_2)
  a <- mu * phi
  b <- (1 - mu) * phi
  log_y <- log(y)
  log_1_y <- log(1 - y)
  dig_mu_phi <- digamma(a)
  dig_1mu_phi <- digamma(b)
  dig_phi <- digamma(phi)
  grad_beta <- rep(0, p)
  grad_gamma <- rep(0, q)
  
  for (i in 1:n) {
    score_mu <- phi[i] * (log_y[i] - log_1_y[i] - dig_mu_phi[i] + dig_1mu_phi[i])
    grad_beta <- grad_beta + score_mu * T_1_diag[i] * X[i, ]
    score_phi <- mu[i] * (log_y[i] - dig_mu_phi[i]) +
      (1 - mu[i]) * (log_1_y[i] - dig_1mu_phi[i]) +
      dig_phi[i]
    grad_gamma <- grad_gamma + score_phi * T_2_diag[i] * Z[i, ]
  }
  return(c(grad_beta, grad_gamma))
}