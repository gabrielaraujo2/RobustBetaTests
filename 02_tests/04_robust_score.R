robust_score_statistic <- function(U, Hn, Jn, idx) {
  H_inv <- tryCatch(solve(Hn), error = function(e) MASS::ginv(Hn))
  H_psi_psi <- H_inv[idx, idx, drop = FALSE]
  Jn_inv <- tryCatch(solve(Jn), error = function(e) MASS::ginv(Jn))
  M <- tryCatch(solve(Hn %*% Jn_inv %*% Hn), error = function(e) MASS::ginv(Hn %*% Jn_inv %*% Hn))
  M_psi_psi <- M[idx, idx, drop = FALSE]
  M_psi_inv <- tryCatch(solve(M_psi_psi), error = function(e) MASS::ginv(M_psi_psi))
  U_psi <- matrix(U[idx], ncol = 1)
  S <- as.numeric(t(U_psi) %*% H_psi_psi %*% M_psi_inv %*% H_psi_psi %*% U_psi)
  p_val <- 1 - pchisq(S, df = length(idx))
  list(statistic = S, p_value = p_val)
}