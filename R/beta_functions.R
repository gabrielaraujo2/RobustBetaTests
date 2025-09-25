#' Beta Regression Log-Likelihood Function
#'
#' Computes the log-likelihood function for beta regression models.
#'
#' @param params Vector of parameters (beta coefficients and precision parameter)
#' @param y Response vector (values between 0 and 1, exclusive)
#' @param X Design matrix
#' @param offset Optional offset vector
#' @return Log-likelihood value
#' @export
beta_loglik <- function(params, y, X, offset = NULL) {
  k <- ncol(X)
  beta <- params[1:k]
  phi <- exp(params[k + 1])  # precision parameter (positive)
  
  eta <- as.vector(X %*% beta)
  if (!is.null(offset)) eta <- eta + offset
  
  mu <- plogis(eta)  # mean parameter
  
  # Beta distribution parameters
  alpha <- mu * phi
  beta_param <- (1 - mu) * phi
  
  # Log-likelihood
  sum(dbeta(y, alpha, beta_param, log = TRUE))
}

#' Beta Regression Score Function
#'
#' Computes the score function (gradient of log-likelihood) for beta regression.
#'
#' @param params Vector of parameters
#' @param y Response vector
#' @param X Design matrix
#' @param offset Optional offset vector
#' @return Score vector
#' @export
beta_score <- function(params, y, X, offset = NULL) {
  k <- ncol(X)
  beta <- params[1:k]
  phi <- exp(params[k + 1])
  
  eta <- as.vector(X %*% beta)
  if (!is.null(offset)) eta <- eta + offset
  
  mu <- plogis(eta)
  dmu_deta <- mu * (1 - mu)  # derivative of logistic function
  
  # Beta distribution parameters
  alpha <- mu * phi
  beta_param <- (1 - mu) * phi
  
  # Score components
  # Using digamma function for score calculation
  d_loglik_dmu <- phi * (digamma(alpha) - digamma(beta_param) + log(y) - log(1 - y))
  
  # Score for beta parameters (using chain rule)
  score_beta <- as.vector(t(X) %*% (d_loglik_dmu * dmu_deta))
  
  # Score for log(phi) parameter
  score_logphi <- phi * sum(mu * (digamma(alpha) - digamma(phi)) + 
                            (1 - mu) * (digamma(beta_param) - digamma(phi)) + 
                            log(y) - log(1 - y))
  
  c(score_beta, score_logphi)
}

#' Beta Regression Hessian Matrix
#'
#' Computes the Hessian matrix (second derivatives of log-likelihood).
#'
#' @param params Vector of parameters
#' @param y Response vector
#' @param X Design matrix
#' @param offset Optional offset vector
#' @return Hessian matrix
#' @export
beta_hessian <- function(params, y, X, offset = NULL) {
  k <- ncol(X)
  beta <- params[1:k]
  phi <- exp(params[k + 1])
  
  eta <- as.vector(X %*% beta)
  if (!is.null(offset)) eta <- eta + offset
  
  mu <- plogis(eta)
  dmu_deta <- mu * (1 - mu)
  
  n <- length(y)
  H <- matrix(0, k + 1, k + 1)
  
  # Approximate Hessian using Fisher Information
  # H_beta_beta
  W <- diag(phi * dmu_deta^2)
  H[1:k, 1:k] <- -t(X) %*% W %*% X
  
  # H_beta_logphi and H_logphi_beta (symmetric)
  H[1:k, k + 1] <- H[k + 1, 1:k] <- rep(0, k)
  
  # H_logphi_logphi
  H[k + 1, k + 1] <- -n * phi^2 / 4  # approximate
  
  H
}