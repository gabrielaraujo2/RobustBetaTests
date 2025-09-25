score_statistic <- function(Un,
                            Kn,
                            indices_H0) {
  p <- length(Un)
  R  <- diag(p)[indices_H0, , drop = FALSE]
  Kn_inv <- tryCatch(solve(Kn), error = function(e) MASS::ginv(Kn))
  K_psi_psi <- R %*% Kn_inv %*% t(R)
  U_psi <- matrix(Un[indices_H0], ncol = 1)
  S  <- t(U_psi) %*% K_psi_psi %*% U_psi
  p_value <- 1 - pchisq(S, df = length(indices_H0))
  list(statistic = as.numeric(S), p_value = as.numeric(p_value))
}