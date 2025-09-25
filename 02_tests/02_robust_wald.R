robust_wald_statistic <- function(theta_hat,
                                  Hn,
                                  Jn,
                                  theta0,
                                  indices_H0
){
  Hn_inv <- tryCatch(solve(Hn), error = function(e) MASS::ginv(Hn))
  Vq <- Hn_inv %*% Jn %*% Hn_inv
  R <- diag(length(theta_hat))[indices_H0, , drop = FALSE]
  Sigmaq <- R %*% Vq %*% t(R)
  delta <- R %*% matrix(theta_hat, ncol = 1) - matrix(theta0, ncol = 1)
  Sigmaq_inv <- tryCatch(solve(Sigmaq), error = function(e) MASS::ginv(Sigmaq))
  W <- t(delta) %*% Sigmaq_inv %*% delta
  list(statistic = as.numeric(W), p_value = 1 - pchisq(W, df = length(indices_H0)))
}