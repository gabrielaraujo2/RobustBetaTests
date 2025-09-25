# Example: Beta Regression Analysis with Robust Methods
# ====================================================

# Load the RobustBetaTests package functions
source("R/beta_functions.R")
source("R/betareg_robust.R") 
source("R/robust_tests.R")

# Example 1: Basic Beta Regression Analysis
# ------------------------------------------

cat("Example 1: Basic Beta Regression Analysis\n")
cat("==========================================\n\n")

# Generate synthetic dataset
set.seed(42)
n <- 200

# Covariates
age <- rnorm(n, mean = 50, sd = 10)
education <- sample(c(0, 1), n, replace = TRUE, prob = c(0.4, 0.6))
income <- rnorm(n, mean = 50000, sd = 15000)

# Standardize continuous variables
age_std <- scale(age)[,1]
income_std <- scale(income)[,1]

# True relationship for proportion (e.g., savings rate)
eta_true <- -0.5 + 0.4*age_std + 0.6*education + 0.3*income_std
mu_true <- plogis(eta_true)
phi_true <- 3

# Generate beta-distributed response
savings_rate <- rbeta(n, mu_true * phi_true, (1 - mu_true) * phi_true)

# Create data frame
data <- data.frame(
  savings_rate = savings_rate,
  age_std = age_std,
  education = education,
  income_std = income_std
)

# Fit beta regression model
cat("Fitting beta regression model...\n")
model <- betareg_robust(savings_rate ~ age_std + education + income_std, data = data)
print(model)

cat("\n")

# Example 2: Robust Variance Estimation Comparison
# ------------------------------------------------

cat("Example 2: Robust Variance Estimation\n")
cat("=====================================\n\n")

# Compare different variance estimators
vcov_model <- model$vcov[1:model$k, 1:model$k]
vcov_sandwich <- robust_vcov(model, type = "sandwich")[1:model$k, 1:model$k]
vcov_huber <- robust_vcov(model, type = "huber")[1:model$k, 1:model$k]

# Extract standard errors
se_model <- sqrt(diag(vcov_model))
se_sandwich <- sqrt(diag(vcov_sandwich))
se_huber <- sqrt(diag(vcov_huber))

# Create comparison table
comparison <- data.frame(
  Coefficient = rownames(vcov_model),
  Model_based = round(se_model, 4),
  Sandwich = round(se_sandwich, 4),
  Huber = round(se_huber, 4)
)

cat("Standard Error Comparison:\n")
print(comparison, row.names = FALSE)

cat("\n")

# Example 3: Hypothesis Testing
# -----------------------------

cat("Example 3: Hypothesis Testing\n")
cat("=============================\n\n")

# Test 1: Individual coefficient test (Wald test)
cat("Test 1: Wald test for education effect\n")
cat("---------------------------------------\n")
R_education <- matrix(c(0, 0, 1, 0), nrow = 1)
wald_test <- robust_test_wald(model, hypothesis = R_education)

cat("H0: Education coefficient = 0\n")
cat("Test statistic:", round(wald_test$statistic, 4), "\n")
cat("p-value:", round(wald_test$p.value, 4), "\n")
cat("Significant at 5% level:", wald_test$p.value < 0.05, "\n\n")

# Test 2: Joint test for multiple coefficients
cat("Test 2: Joint test for age and income effects\n")
cat("---------------------------------------------\n")
R_joint <- matrix(c(0, 1, 0, 0,
                    0, 0, 0, 1), nrow = 2, byrow = TRUE)
wald_joint <- robust_test_wald(model, hypothesis = R_joint)

cat("H0: Age and Income coefficients both = 0\n")
cat("Test statistic:", round(wald_joint$statistic, 4), "\n")
cat("p-value:", round(wald_joint$p.value, 4), "\n")
cat("Significant at 5% level:", wald_joint$p.value < 0.05, "\n\n")

# Test 3: Model comparison using Likelihood Ratio test
cat("Test 3: Model comparison (LR test)\n")
cat("----------------------------------\n")

# Fit reduced model without income
model_reduced <- betareg_robust(savings_rate ~ age_std + education, data = data)
lr_test <- robust_test_lr(model_reduced, model)

cat("Comparing full model vs. model without income\n")
cat("Test statistic:", round(lr_test$statistic, 4), "\n") 
cat("p-value:", round(lr_test$p.value, 4), "\n")
cat("Include income in model:", lr_test$p.value < 0.05, "\n\n")

# Example 4: Model Diagnostics
# ----------------------------

cat("Example 4: Model Diagnostics\n")
cat("============================\n\n")

# Fitted values vs residuals
fitted_vals <- model$fitted.values
residuals <- model$residuals

# Summary statistics
cat("Model Summary:\n")
cat("Log-likelihood:", round(model$loglik, 4), "\n")
cat("AIC:", round(model$aic, 4), "\n")
cat("BIC:", round(model$bic, 4), "\n")
cat("Precision parameter:", round(model$precision, 4), "\n")
cat("Number of observations:", model$n, "\n")
cat("Convergence:", model$converged, "\n\n")

# Residual summary
cat("Residual Summary:\n")
cat("Min:", round(min(residuals), 4), "\n")
cat("Q1:", round(quantile(residuals, 0.25), 4), "\n")
cat("Median:", round(median(residuals), 4), "\n")
cat("Q3:", round(quantile(residuals, 0.75), 4), "\n")
cat("Max:", round(max(residuals), 4), "\n")
cat("Mean absolute residual:", round(mean(abs(residuals)), 4), "\n")

cat("\nExample analysis completed!\n")