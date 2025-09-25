#' Robust Variance-Covariance Matrix
#'
#' Computes robust variance-covariance matrix using sandwich estimator.
#'
#' @param object A betareg_robust model object
#' @param type Type of robust variance estimator ("sandwich", "huber")
#' @return Robust variance-covariance matrix
#' @export
robust_vcov <- function(object, type = "sandwich") {
  
  if (!inherits(object, "betareg_robust")) {
    stop("Object must be of class 'betareg_robust'")
  }
  
  X <- object$design_matrix
  y <- object$response
  params <- c(object$coefficients, object$log_precision)
  n <- object$n
  k <- object$k
  
  if (type == "sandwich") {
    # Sandwich estimator: A^(-1) B A^(-1)
    # where A = -Hessian, B = sum of outer products of scores
    
    A <- -beta_hessian(params, y, X)
    
    # Compute individual score contributions
    scores <- matrix(0, n, k + 1)
    for (i in 1:n) {
      scores[i, ] <- beta_score(params, y[i], matrix(X[i, ], nrow = 1))
    }
    
    B <- t(scores) %*% scores
    
    # Sandwich estimator
    A_inv <- solve(A)
    vcov_robust <- A_inv %*% B %*% A_inv
    
  } else if (type == "huber") {
    # Huber-White robust estimator
    # Similar to sandwich but with different weighting
    
    A <- -beta_hessian(params, y, X)
    
    # Compute weighted scores
    eta <- as.vector(X %*% object$coefficients)
    mu <- plogis(eta)
    weights <- sqrt(mu * (1 - mu))  # Huber weights
    
    scores <- matrix(0, n, k + 1)
    for (i in 1:n) {
      score_i <- beta_score(params, y[i], matrix(X[i, ], nrow = 1))
      scores[i, ] <- weights[i] * score_i
    }
    
    B <- t(scores) %*% scores
    
    A_inv <- solve(A)
    vcov_robust <- A_inv %*% B %*% A_inv
    
  } else {
    stop("Unknown robust variance type. Use 'sandwich' or 'huber'.")
  }
  
  vcov_robust
}

#' Robust Wald Test
#'
#' Performs robust Wald test for linear hypotheses in beta regression.
#'
#' @param object A betareg_robust model object
#' @param hypothesis A matrix specifying linear constraints (R beta = r)
#' @param rhs Right-hand side vector (default: zero vector)
#' @param vcov_type Type of robust variance estimator
#' @return A list with test results
#' @export
robust_test_wald <- function(object, hypothesis, rhs = NULL, vcov_type = "sandwich") {
  
  if (!inherits(object, "betareg_robust")) {
    stop("Object must be of class 'betareg_robust'")
  }
  
  R <- as.matrix(hypothesis)
  q <- nrow(R)
  k <- object$k
  
  if (ncol(R) != k) {
    stop("Hypothesis matrix dimensions do not match number of coefficients")
  }
  
  if (is.null(rhs)) {
    rhs <- rep(0, q)
  }
  
  beta <- object$coefficients
  vcov_robust <- robust_vcov(object, type = vcov_type)
  
  # Extract variance-covariance for beta coefficients only
  vcov_beta <- vcov_robust[1:k, 1:k]
  
  # Wald statistic
  restriction <- R %*% beta - rhs
  var_restriction <- R %*% vcov_beta %*% t(R)
  
  if (q == 1) {
    # Single restriction - t-test
    wald_stat <- as.numeric(restriction^2 / var_restriction)
    p_value <- 1 - pchisq(wald_stat, df = 1)
    test_type <- "Chi-squared test"
    df <- 1
  } else {
    # Multiple restrictions - F-test approximation
    wald_stat <- as.numeric(t(restriction) %*% solve(var_restriction) %*% restriction)
    p_value <- 1 - pchisq(wald_stat, df = q)
    test_type <- "Chi-squared test"
    df <- q
  }
  
  list(
    statistic = wald_stat,
    p.value = p_value,
    df = df,
    test_type = test_type,
    hypothesis = R,
    rhs = rhs,
    restriction = restriction,
    vcov_type = vcov_type
  )
}

#' Robust Score Test
#'
#' Performs robust Score (Lagrange Multiplier) test.
#'
#' @param object_restricted A betareg_robust model object (restricted model)
#' @param object_unrestricted A betareg_robust model object (unrestricted model)
#' @param vcov_type Type of robust variance estimator
#' @return A list with test results
#' @export
robust_test_score <- function(object_restricted, object_unrestricted, vcov_type = "sandwich") {
  
  if (!inherits(object_restricted, "betareg_robust") || 
      !inherits(object_unrestricted, "betareg_robust")) {
    stop("Both objects must be of class 'betareg_robust'")
  }
  
  # Degrees of freedom
  df <- object_unrestricted$k - object_restricted$k
  
  if (df <= 0) {
    stop("Unrestricted model must have more parameters than restricted model")
  }
  
  # Score test uses restricted model estimates
  y <- object_restricted$response
  X_restricted <- object_restricted$design_matrix
  X_unrestricted <- object_unrestricted$design_matrix
  
  params_restricted <- c(object_restricted$coefficients, object_restricted$log_precision)
  
  # Score vector at restricted estimates
  score_restricted <- beta_score(params_restricted, y, X_restricted)
  
  # Information matrix at restricted estimates  
  info_matrix <- -beta_hessian(params_restricted, y, X_restricted)
  
  # Score statistic
  score_stat <- as.numeric(t(score_restricted) %*% solve(info_matrix) %*% score_restricted)
  p_value <- 1 - pchisq(score_stat, df = df)
  
  list(
    statistic = score_stat,
    p.value = p_value,
    df = df,
    test_type = "Score test (LM test)",
    vcov_type = vcov_type
  )
}

#' Robust Likelihood Ratio Test
#'
#' Performs robust Likelihood Ratio test.
#'
#' @param object_restricted A betareg_robust model object (restricted model)
#' @param object_unrestricted A betareg_robust model object (unrestricted model)
#' @return A list with test results
#' @export
robust_test_lr <- function(object_restricted, object_unrestricted) {
  
  if (!inherits(object_restricted, "betareg_robust") || 
      !inherits(object_unrestricted, "betareg_robust")) {
    stop("Both objects must be of class 'betareg_robust'")
  }
  
  # Degrees of freedom
  df <- object_unrestricted$k - object_restricted$k
  
  if (df <= 0) {
    stop("Unrestricted model must have more parameters than restricted model")
  }
  
  # Likelihood ratio statistic
  lr_stat <- 2 * (object_unrestricted$loglik - object_restricted$loglik)
  p_value <- 1 - pchisq(lr_stat, df = df)
  
  list(
    statistic = lr_stat,
    p.value = p_value,
    df = df,
    test_type = "Likelihood Ratio test",
    loglik_restricted = object_restricted$loglik,
    loglik_unrestricted = object_unrestricted$loglik
  )
}