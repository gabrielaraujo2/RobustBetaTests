#****Fisher information matrix for MLE****
fisher_info_beta_regression <- function(y, X, Z, theta_hat, linkmu = "logit", linkphi = "log") {
  n <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  beta <- theta_hat[1:p]
  gamma <- theta_hat[(p+1):(p+q)]
  mu <- mean_link(beta, y, X, linkmu)$mu
  phi <- precision_link(gamma, Z, linkphi)$phi
  dmu_deta <- diag(mean_link(beta, y, X, linkmu)$T_1)
  dphi_deta <- diag(precision_link(gamma, Z, linkphi)$T_2)
  psi1_mu_phi <- trigamma(mu * phi)
  psi1_1mu_phi <- trigamma((1 - mu) * phi)
  psi1_phi <- trigamma(phi)
  I_bb <- matrix(0, nrow = p, ncol = p)
  I_bg <- matrix(0, nrow = p, ncol = q)
  I_gg <- matrix(0, nrow = q, ncol = q)
  for (i in 1:n) {
    xi <- matrix(X[i, ], ncol = 1)
    zi <- matrix(Z[i, ], ncol = 1)
    wi_bb <- phi[i]^2 * (psi1_mu_phi[i] + psi1_1mu_phi[i]) * dmu_deta[i]^2
    I_bb <- I_bb + wi_bb * (xi %*% t(xi))
    wi_bg <- phi[i] * (mu[i] * psi1_mu_phi[i] - (1 - mu[i]) * psi1_1mu_phi[i]) * dphi_deta[i] * dmu_deta[i]
    I_bg <- I_bg + wi_bg * (xi %*% t(zi))
    wi_gg <- (mu[i]^2 * psi1_mu_phi[i] + (1 - mu[i])^2 * psi1_1mu_phi[i] - psi1_phi[i]) * dphi_deta[i]^2
    I_gg <- I_gg + wi_gg * (zi %*% t(zi))
  }
  I <- rbind(cbind(I_bb, I_bg), cbind(t(I_bg), I_gg))
  return(I)
}