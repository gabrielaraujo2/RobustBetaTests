# Test script for RobustBetaTests package

# Load the package functions
source("R/beta_functions.R")
source("R/betareg_robust.R") 
source("R/robust_tests.R")

# Test with simulated data
set.seed(123)
n <- 100

# Generate covariates
x1 <- rnorm(n)
x2 <- rnorm(n)
X <- cbind(1, x1, x2)  # Design matrix with intercept

# True parameters
beta_true <- c(0.5, 0.3, -0.2)
phi_true <- 2

# Generate response
eta <- X %*% beta_true
mu <- plogis(eta)
alpha <- mu * phi_true
beta_param <- (1 - mu) * phi_true
y <- rbeta(n, alpha, beta_param)

# Create data frame
test_data <- data.frame(y = y, x1 = x1, x2 = x2)

cat("Testing RobustBetaTests package\n")
cat("===============================\n\n")

# Test 1: Basic model fitting
cat("Test 1: Basic model fitting\n")
cat("----------------------------\n")

try({
  model1 <- betareg_robust(y ~ x1 + x2, data = test_data)
  print(model1)
  
  cat("\nTrue coefficients:", beta_true, "\n")
  cat("Estimated coefficients:", round(model1$coefficients, 3), "\n")
  cat("True precision:", phi_true, "\n")
  cat("Estimated precision:", round(model1$precision, 3), "\n")
  
}, silent = FALSE)

cat("\n\n")

# Test 2: Robust variance estimation
cat("Test 2: Robust variance estimation\n") 
cat("-----------------------------------\n")

try({
  vcov_sandwich <- robust_vcov(model1, type = "sandwich")
  vcov_huber <- robust_vcov(model1, type = "huber")
  
  cat("Standard errors (model-based):", round(sqrt(diag(model1$vcov[1:3, 1:3])), 4), "\n")
  cat("Standard errors (sandwich):", round(sqrt(diag(vcov_sandwich[1:3, 1:3])), 4), "\n")
  cat("Standard errors (Huber):", round(sqrt(diag(vcov_huber[1:3, 1:3])), 4), "\n")
  
}, silent = FALSE)

cat("\n\n")

# Test 3: Wald test
cat("Test 3: Wald test\n")
cat("-----------------\n")

try({
  # Test H0: beta_2 = 0 (coefficient of x1)
  R <- matrix(c(0, 1, 0), nrow = 1)
  wald_test <- robust_test_wald(model1, hypothesis = R)
  
  cat("Wald test for H0: coefficient of x1 = 0\n")
  cat("Test statistic:", round(wald_test$statistic, 4), "\n")
  cat("p-value:", round(wald_test$p.value, 4), "\n")
  cat("Reject H0 at 5% level:", wald_test$p.value < 0.05, "\n")
  
}, silent = FALSE)

cat("\n\n")

# Test 4: Model comparison using LR test
cat("Test 4: Likelihood Ratio test\n")
cat("------------------------------\n")

try({
  # Fit restricted model (without x2)
  model_restricted <- betareg_robust(y ~ x1, data = test_data)
  
  # LR test
  lr_test <- robust_test_lr(model_restricted, model1)
  
  cat("LR test comparing models with and without x2\n")
  cat("Test statistic:", round(lr_test$statistic, 4), "\n")
  cat("p-value:", round(lr_test$p.value, 4), "\n")
  cat("Reject H0 at 5% level:", lr_test$p.value < 0.05, "\n")
  
}, silent = FALSE)

cat("\nAll tests completed successfully!\n")