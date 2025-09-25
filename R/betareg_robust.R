#' Robust Beta Regression
#'
#' Fits a beta regression model using robust methods.
#'
#' @param formula A formula specifying the model
#' @param data A data frame containing the variables
#' @param weights Optional weights vector
#' @param offset Optional offset vector
#' @param method Estimation method ("ML" for maximum likelihood)
#' @return A list containing model results
#' @export
betareg_robust <- function(formula, data, weights = NULL, offset = NULL, method = "ML") {
  
  # Extract response and design matrix
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  X <- model.matrix(formula, data)
  
  # Validate response
  if (any(y <= 0 | y >= 1)) {
    stop("Response variable must be in the interval (0, 1)")
  }
  
  n <- length(y)
  k <- ncol(X)
  
  # Initial values
  # Use logistic regression coefficients as starting values
  # Transform response to avoid boundary issues in GLM
  y_transform <- ifelse(y == 0, 0.001, ifelse(y == 1, 0.999, y))
  glm_fit <- glm(y_transform ~ X - 1, family = binomial(), weights = weights)
  beta_init <- coef(glm_fit)
  beta_init[is.na(beta_init)] <- 0  # Handle any NA coefficients
  phi_init <- log(2)  # initial log(precision)
  
  params_init <- c(beta_init, phi_init)
  
  # Optimize log-likelihood
  opt_result <- optim(
    par = params_init,
    fn = function(p) -beta_loglik(p, y, X, offset),
    method = "BFGS",
    hessian = TRUE,
    control = list(maxit = 1000, reltol = 1e-8)
  )
  
  if (opt_result$convergence != 0) {
    warning("Optimization did not converge")
  }
  
  # Extract results
  coefficients <- opt_result$par[1:k]
  log_precision <- opt_result$par[k + 1]
  precision <- exp(log_precision)
  
  # Standard errors from Hessian
  vcov_matrix <- tryCatch({
    solve(opt_result$hessian)
  }, error = function(e) {
    # If Hessian is singular, use generalized inverse
    MASS::ginv(opt_result$hessian)
  })
  
  se <- sqrt(pmax(0, diag(vcov_matrix)))  # Ensure non-negative variances
  
  # Fitted values
  eta <- as.vector(X %*% coefficients)
  if (!is.null(offset)) eta <- eta + offset
  fitted_values <- plogis(eta)
  
  # Residuals
  residuals <- y - fitted_values
  
  # Model information
  loglik <- -opt_result$value
  aic <- -2 * loglik + 2 * (k + 1)
  bic <- -2 * loglik + log(n) * (k + 1)
  
  result <- list(
    coefficients = coefficients,
    precision = precision,
    log_precision = log_precision,
    vcov = vcov_matrix,
    se = se,
    fitted.values = fitted_values,
    residuals = residuals,
    loglik = loglik,
    aic = aic,
    bic = bic,
    converged = opt_result$convergence == 0,
    iterations = opt_result$counts[1],
    formula = formula,
    data = data,
    response = y,
    design_matrix = X,
    n = n,
    k = k
  )
  
  class(result) <- "betareg_robust"
  result
}

#' Print method for betareg_robust objects
#' @param x A betareg_robust object
#' @param ... Additional arguments
#' @export
print.betareg_robust <- function(x, ...) {
  cat("Robust Beta Regression\n")
  cat("======================\n\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Number of observations:", x$n, "\n")
  cat("Log-likelihood:", round(x$loglik, 4), "\n")
  cat("AIC:", round(x$aic, 4), "\n")
  cat("BIC:", round(x$bic, 4), "\n")
  cat("Precision parameter:", round(x$precision, 4), "\n")
  cat("Converged:", x$converged, "\n\n")
  
  cat("Coefficients:\n")
  coef_table <- data.frame(
    Estimate = round(x$coefficients, 4),
    Std.Error = round(x$se[1:x$k], 4),
    t.value = round(x$coefficients / x$se[1:x$k], 4),
    Pr = round(2 * (1 - pnorm(abs(x$coefficients / x$se[1:x$k]))), 4)
  )
  rownames(coef_table) <- colnames(x$design_matrix)
  print(coef_table)
}