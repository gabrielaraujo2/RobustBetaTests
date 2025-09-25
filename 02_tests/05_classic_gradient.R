gradient_statistic <- function(Un, # Score statistic for testing hypotheses on parameters
                               theta_hat, # Estimated parameters
                               theta0, # Hypothesized values of parameters under H0
                               indices_H0 # Indices of parameters under H0
                               ) {
  
  # Number of parameters
  p <- length(Un)
  
  # Matrix for the hypothesis test
  R <- diag(p)[indices_H0, , drop = FALSE]
  
  # Difference between estimated and hypothesized parameters
  delta <- R %*% theta_hat - matrix(theta0, ncol = 1)
  
  # U_psi is the score vector for the parameters under H0
  U_psi <- matrix(Un[indices_H0], ncol = 1)
  
  # Gradient Statistics
  G <- t(U_psi) %*% delta
  
  # p-value for the gradient statistic
  p_value <- 1 - pchisq(G, df = length(indices_H0))
  
  # Return the gradient statistic and p-value
  list(statistic = as.numeric(G), p_value = as.numeric(p_value))
}