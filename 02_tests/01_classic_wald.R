wald_statistic <- function(theta_hat, # Estimated parameters
                           Kn, # Fisher information matrix
                           theta0, # Hypothesized values of parameters under H0
                           indices_H0 # Indices of parameters under H0
                           ){
  # Number of parameters
  p <- length(theta_hat)
  
  # Matrix for the hypothesis test
  R <- diag(p)[indices_H0, , drop = FALSE]
  
  # Difference between estimated and hypothesized parameters
  delta <- R %*% theta_hat - matrix(theta0, ncol = 1)
  
  # Inverse of the Fisher information matrix
  Kn_inv <- tryCatch(solve(Kn), error = function(e) {
    MASS::ginv(Kn)
  })
  
  # Variance of the linear combination of parameters
  K_psi_psi <- R %*% Kn_inv %*% t(R)
  K_psi_inv <- tryCatch(
    solve(K_psi_psi),
    error = function(e) MASS::ginv(K_psi_psi)
  )
  
  # Wald statistic
  W <- t(delta) %*% K_psi_inv %*% delta
  
  # p-value for the Wald statistic
  p_value <- 1 - pchisq(W, df = length(indices_H0))
  
  # Return the Wald statistic and p-value
  list(statistic = as.numeric(W), p_value = as.numeric(p_value))
  }