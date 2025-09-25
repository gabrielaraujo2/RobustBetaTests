# RobustBetaTests

Computational routines for robust hypothesis testing in beta regression models.

## Description

This R package provides functions for fitting beta regression models with robust methods and performing robust hypothesis tests including Wald, Score, and Likelihood Ratio tests with robust variance estimation.

Beta regression models are useful for analyzing continuous responses that are bounded between 0 and 1, such as rates, proportions, and percentages. This package implements robust estimation methods that provide reliable inference even when model assumptions are violated.

## Installation

```r
# Install from source (development version)
# First, make sure you have the required dependencies
install.packages(c("stats", "MASS"))

# Then install the package from source
install.packages("path/to/RobustBetaTests", repos = NULL, type = "source")
```

## Features

- **Beta Regression Model Fitting**: Robust estimation of beta regression parameters
- **Robust Variance Estimation**: Sandwich and Huber-White robust variance estimators
- **Hypothesis Testing**: 
  - Robust Wald tests for linear hypotheses
  - Score (Lagrange Multiplier) tests
  - Likelihood Ratio tests
- **Model Diagnostics**: Fitted values, residuals, and model statistics

## Quick Start

```r
library(RobustBetaTests)

# Simulate some data
set.seed(123)
n <- 100
x1 <- rnorm(n)
x2 <- rnorm(n)

# True parameters
eta <- 0.5 + 0.3*x1 - 0.2*x2
mu <- plogis(eta)
phi <- 2

# Generate response
y <- rbeta(n, mu*phi, (1-mu)*phi)
data <- data.frame(y = y, x1 = x1, x2 = x2)

# Fit beta regression model
model <- betareg_robust(y ~ x1 + x2, data = data)
print(model)

# Robust variance estimation
vcov_sandwich <- robust_vcov(model, type = "sandwich")
vcov_huber <- robust_vcov(model, type = "huber")

# Wald test for H0: coefficient of x1 = 0
R <- matrix(c(0, 1, 0), nrow = 1)
wald_test <- robust_test_wald(model, hypothesis = R)
print(wald_test)

# Model comparison using likelihood ratio test
model_restricted <- betareg_robust(y ~ x1, data = data)
lr_test <- robust_test_lr(model_restricted, model)
print(lr_test)
```

## Functions

### Core Functions

- `betareg_robust()`: Fit beta regression model with robust methods
- `robust_vcov()`: Compute robust variance-covariance matrix
- `robust_test_wald()`: Perform robust Wald test
- `robust_test_score()`: Perform robust Score (LM) test  
- `robust_test_lr()`: Perform Likelihood Ratio test

### Utility Functions

- `beta_loglik()`: Beta regression log-likelihood function
- `beta_score()`: Beta regression score function
- `beta_hessian()`: Beta regression Hessian matrix

## Theory

Beta regression models assume the response variable Y follows a beta distribution:

Y ~ Beta(μφ, (1-μ)φ)

where μ is the mean parameter and φ is the precision parameter. The mean is related to covariates through a logit link:

logit(μ) = X'β

This package implements robust methods for:

1. **Parameter Estimation**: Maximum likelihood with robust optimization
2. **Variance Estimation**: Sandwich estimators that are robust to model misspecification
3. **Hypothesis Testing**: Robust versions of classical tests

## License

GPL-3

## Authors

Gabriel Araujo

## References

- Ferrari, S. L. P., & Cribari-Neto, F. (2004). Beta regression for modelling rates and proportions. Journal of Applied Statistics, 31(7), 799-815.
- White, H. (1980). A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity. Econometrica, 817-838.