robust_gradient_statistic <- function(U,
                                      theta_hat,
                                      theta0,
                                      Hn,
                                      Jn,
                                      indices_H0
) {
  U_psi <- matrix(U[indices_H0], ncol = 1)
  delta <- theta_hat[indices_H0] - theta0
  delta <- as.matrix(delta, ncol = 1)
  Jn_inv <- tryCatch(solve(Jn), error = function(e) MASS::ginv(Jn))
  M <- Jn_inv %*% Hn
  M_psi_psi <- M[indices_H0, indices_H0, drop = FALSE]
  T_stat <- -t(U_psi) %*% M_psi_psi %*% delta
  p_value <- 1 - pchisq(T_stat, df = length(indices_H0))
  list(statistic = as.numeric(T_stat), p_value = as.numeric(p_value))
}