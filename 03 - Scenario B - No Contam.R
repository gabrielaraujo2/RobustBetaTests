#*** We will be storing the following quantities***
#* 1. alpha_lsmle for irrestricted space
#* 2. alpha_lmdpde for irrestricted space
#* 3. alpha_lsmle_h1 for restricted space H1 (beta2 = -5.5)
#* 4. alpha_lmdpde_h1 for restricted space H1 (beta2 = -5.5)
#* 5. beta1_mle for irestricted space
#* 6. beta1_lsmle for irestricted space
#* 7. beta1_lmdpde for irestricted space
#* 8. beta2_mle for irestricted space
#* 9. beta2_lsmle for irrestricted space
#* 10. beta2_lmdpde for irrestricted space
#* 11. gamma1_mle for irestricted space
#* 12. gamma1_lsmle for irestricted space
#* 13. gamma1_lmdpde for irestricted space
#* 14. beta1_mle_h1 for restricted space H1 (beta2 = -5.5)
#* 15. beta1_lsmle_h1 for restricted space H1 (beta2 = -5.5)
#* 16. beta1_lmdpde_h1 for restricted space H1 (beta2 = -5.5)
#* 17. gamma1_mle_h1 for restricted space H1 (beta2 = -5.5)
#* 18. gamma1_lsmle_h1 for restricted space H1 (beta2 = -5.5)
#* 19. gamma1_lmdpde_h1 for restricted space H1 (beta2 = -5.5)
#* 20. beta1_mle_h2 for restricted space H2 (gamma1 = 5)
#* 21. beta1_lsmle_h2 for restricted space H2 (gamma1 = 5)
#* 22. beta1_lmdpde_h2 for restricted space H2 (gamma1 = 5)
#* 23. beta2_mle_h2 for restricted space H2 (gamma1 = 5)
#* 24. beta2_lsmle_h2 for restricted space H2 (gamma1 = 5)
#* 25. beta2_lmdpde_h2 for restricted space H2 (gamma1 = 5)
#* 26. gamma1_mle_h3 for restricted space H3 (beta = (-1, -5.5))
#* 27. gamma1_lsmle_h3 for restricted space H3 (beta = (-1, -5.5))
#* 28. gamma1_lmdpde_h3 for restricted space H3 (beta = (-1, -5.5))
#* 29. All statistics and p-values

# ---- Cleaning environment and loading functions ----
rm(list = ls())

setwd("G:/Meu Drive/03 - Sim/03_case_study_1/")

# Functions files
files_functions <- arquivos <- list.files(path = "../01_functions", full.names = TRUE, recursive = TRUE)

# Testing files
files_test <- list.files(path = "../02_tests", full.names = TRUE, recursive = TRUE)

# Source all functions and tests
lapply(files_functions, function(file) try(source(file), silent = TRUE))
lapply(files_test, function(file) try(source(file), silent = TRUE))

# Loading the X matrix
X_matrix <- read.table("../X_matrix_40.txt", header = TRUE)

# ---- Simulation setup ----
N <- 10000
# sample_sizes <- c(40, 80, 160, 320)
sample_sizes <- c(1200)
beta <- c(-1, -5.5)
gamma <- c(5)
theta <- c(beta, gamma)

# ---- Start loop over sample sizes ----
for (n in sample_sizes) {
  
  cat("\n>>> Simulation for n =", n, " <<<\n")
  
  # -----------> Tuning parameters <-----------
  alpha_lsmle <- numeric(N)
  alpha_lmdpde <- numeric(N)
  alpha_lsmle_h1 <- numeric(N)
  alpha_lmdpde_h1 <- numeric(N)
  
  # -----------> Estimated Parameters <-----------
  beta1_mle         <- numeric(N)
  beta1_lsmle       <- numeric(N)
  beta1_lmdpde      <- numeric(N)
  beta2_mle         <- numeric(N)
  beta2_lsmle      <- numeric(N)
  beta2_lmdpde     <- numeric(N)
  gamma1_mle       <- numeric(N)
  gamma1_lsmle     <- numeric(N)
  gamma1_lmdpde    <- numeric(N)
  beta1_mle_h1     <- numeric(N)
  beta1_lsmle_h1   <- numeric(N)
  beta1_lmdpde_h1  <- numeric(N)
  gamma1_mle_h1    <- numeric(N)
  gamma1_lsmle_h1  <- numeric(N)
  gamma1_lmdpde_h1 <- numeric(N)
  beta1_mle_h2     <- numeric(N)
  beta1_lsmle_h2   <- numeric(N)
  beta1_lmdpde_h2  <- numeric(N)
  beta2_mle_h2     <- numeric(N)
  beta2_lsmle_h2   <- numeric(N)
  beta2_lmdpde_h2  <- numeric(N)
  gamma1_mle_h3    <- numeric(N)
  gamma1_lsmle_h3  <- numeric(N)
  gamma1_lmdpde_h3 <- numeric(N)
  
  # -----------> # Statistics and p-values <-----------
  statistics <- matrix(NA, nrow = N, ncol = 36)
  pvalues <- matrix(NA, nrow = N, ncol = 36)

  # Names of tests
  test_names <- c("wald_h1", "score_h1", "grad_h1", "wald_h2", "score_h2", "grad_h2", "wald_h3", "score_h3", "grad_h3", "wald_h4", "score_h4", "grad_h4")
  colnames(statistics) <- paste0("stat_", rep(test_names, each = 3), rep(c("_mle", "_lsmle", "_lmdpde"), times = 12))
  colnames(pvalues)   <- paste0("pval_",  rep(test_names, each = 3), rep(c("_mle", "_lsmle", "_lmdpde"), times = 12))

  if(n == 40){X = do.call(rbind, replicate(1, X_matrix, simplify = FALSE))}
  if(n == 80){X = do.call(rbind, replicate(2, X_matrix, simplify = FALSE))}
  if(n == 160){X = do.call(rbind, replicate(4, X_matrix, simplify = FALSE))}
  if(n == 320){X = do.call(rbind, replicate(8, X_matrix, simplify = FALSE))}
  if(n == 1200){X = do.call(rbind, replicate(30, X_matrix, simplify = FALSE))}
  
  X <- as.matrix(X)
  Z <- as.matrix(rep(1, n))
  
  # -----------> Initialization the Looping <-----------
  for (i in 1:N) {
    
    cat("Simulation", i, "of", N, "| n =", n, "\n")
    set.seed(123 + i)
    
    repeat {

      # -----------> Generating data <-----------
      mu <- exp(X %*% beta) / (1 + exp(X %*% beta))
      phi <- exp(Z %*% gamma)
      y <- rbeta(n, mu * phi, (1 - mu) * phi)
      
      # -----------> Irrestricted Model <-----------
      fit_mle      <- try(MLE_BETA(y, X, Z), silent = TRUE)
      fit_lsmle    <- try(LSMLE_BETA(y, X, Z), silent = TRUE)
      fit_lmdpde   <- try(LMDPDE_BETA(y, X, Z), silent = TRUE)

      # -----------> Model under H1 (beta2 = -5.5) <-----------
      fit_mle_h1      <- try(MLE_BETA(y, X, Z, beta_fix = beta[2], idx_beta = c(2)), silent = TRUE)
      fit_lsmle_h1    <- try(LSMLE_BETA(y, X, Z, beta_fix = beta[2], idx_beta = c(2), qoptimal = FALSE, q0 = 1 - fit_lmdpde$alpha_const), silent = TRUE)
      fit_lmdpde_h1   <- try(LMDPDE_BETA(y, X, Z, beta_fix = beta[2], idx_beta = c(2), qoptimal = FALSE, q0 = 1 - fit_lmdpde$alpha_const), silent = TRUE)
      fit_lsmle_h1_grad  <- try(LSMLE_BETA(y, X, Z, beta_fix = beta[2], idx_beta = c(2), qoptimal = FALSE, q0 = 1 - fit_lsmle$alpha_const), silent = TRUE)
      fit_lmdpde_h1_grad <- try(LMDPDE_BETA(y, X, Z, beta_fix = beta[2], idx_beta = c(2), qoptimal = FALSE, q0 = 1 - fit_lmdpde$alpha_const), silent = TRUE)

      # -----------> Model under H2 (gamma1 = 5) <-----------
      fit_mle_h2      <- try(MLE_BETA(y, X, Z, gamma_fix = 5, idx_gamma = c(1)), silent = TRUE)
      fit_lsmle_h2    <- try(LSMLE_BETA(y, X, Z, gamma_fix = 5, idx_gamma = c(1), qoptimal = FALSE, q0 = 1 - fit_lsmle$alpha_const), silent = TRUE)
      fit_lmdpde_h2   <- try(LMDPDE_BETA(y, X, Z, gamma_fix = 5, idx_gamma = c(1), qoptimal = FALSE, q0 = 1 - fit_lmdpde$alpha_const), silent = TRUE)

      # -----------> Model under H3 (beta = (-1, -5.5)) <-----------
      fit_mle_h3      <- try(MLE_BETA(y, X, Z, beta_fix = beta, idx_beta = c(1,2)), silent = TRUE)
      fit_lsmle_h3    <- try(LSMLE_BETA(y, X, Z, beta_fix = beta, idx_beta = c(1,2), qoptimal = FALSE, q0 = 1 - fit_lsmle$alpha_const), silent = TRUE)
      fit_lmdpde_h3   <- try(LMDPDE_BETA(y, X, Z, beta_fix = beta, idx_beta = c(1,2), qoptimal = FALSE, q0 = 1 - fit_lmdpde$alpha_const), silent = TRUE)

      # H4: beta and gamma fixed
      fit_mle_h4      <- try(MLE_BETA(y, X, Z, beta_fix = beta, idx_beta = c(1,2), gamma_fix = gamma, idx_gamma = c(1)), silent = TRUE)
      fit_lsmle_h4    <- try(LSMLE_BETA(y, X, Z, beta_fix = beta, idx_beta = c(1,2), gamma_fix = gamma, idx_gamma = c(1), q0 = 1 - fit_lsmle$alpha_const), silent = TRUE)
      fit_lmdpde_h4   <- try(LMDPDE_BETA(y, X, Z, beta_fix = beta, idx_beta = c(1,2), gamma_fix = gamma, idx_gamma = c(1), q0 = 1 - fit_lmdpde$alpha_const), silent = TRUE)

      # -----------> List of all models <-----------
      all_tests <- list(fit_mle, fit_lsmle, fit_lmdpde, # Irrestricted space
                        fit_mle_h1, fit_lsmle_h1, fit_lmdpde_h1, # \beta2 = -5.5
                        fit_lsmle_h1_grad, fit_lmdpde_h1_grad, # \beta2 = -5.5
                        fit_mle_h2, fit_lsmle_h2, fit_lmdpde_h2, # \gamma1 = 5
                        fit_mle_h3, fit_lsmle_h3, fit_lmdpde_h3, # \beta = (-1, -5.5)
                        fit_mle_h4, fit_lsmle_h4, fit_lmdpde_h4 # \beta = (-1, -5.5), \gamma1 = 5
                        )
      
      if (any(sapply(all_tests, inherits, "try-error"))) next else break
    }

    # -----------> Saving the alpha <-----------
    alpha_lsmle[i]  <- fit_lsmle$alpha_const
    alpha_lmdpde[i] <- fit_lmdpde$alpha_const
    alpha_lsmle_h1[i]  <- fit_lsmle_h1$alpha_const
    alpha_lmdpde_h1[i] <- fit_lmdpde_h1$alpha_const

    # -----------> Saving the estimated parameters <-----------
    beta1_mle[i] <- fit_mle$theta_hat_full[1]
    beta1_lsmle[i] <- fit_lsmle$theta_hat[1]
    beta1_lmdpde[i] <- fit_lmdpde$theta_hat[1]
    beta2_mle[i] <- fit_mle$theta_hat_full[2]
    beta2_lsmle[i] <- fit_lsmle$theta_hat[2]
    beta2_lmdpde[i] <- fit_lmdpde$theta_hat[2]
    gamma1_mle[i] <- fit_mle$theta_hat_full[3]
    gamma1_lsmle[i] <- fit_lsmle$theta_hat[3]
    gamma1_lmdpde[i] <- fit_lmdpde$theta_hat[3]

    # Model under H1 (beta2 = -5.5)
    beta1_mle_h1[i] <- fit_mle_h1$theta_hat_full[1]
    beta1_lsmle_h1[i] <- fit_lsmle_h1$theta_hat[1]
    beta1_lmdpde_h1[i] <- fit_lmdpde_h1$theta_hat[1]
    gamma1_mle_h1[i] <- fit_mle_h1$theta_hat_full[3]
    gamma1_lsmle_h1[i] <- fit_lsmle_h1$theta_hat[3]
    gamma1_lmdpde_h1[i] <- fit_lmdpde_h1$theta_hat[3]

    # Model under H2 (gamma1 = 5)
    beta1_mle_h2[i] <- fit_mle_h2$theta_hat_full[1]
    beta1_lsmle_h2[i] <- fit_lsmle_h2$theta_hat[1]
    beta1_lmdpde_h2[i] <- fit_lmdpde_h2$theta_hat[1]
    beta2_mle_h2[i] <- fit_mle_h2$theta_hat_full[2]
    beta2_lsmle_h2[i] <- fit_lsmle_h2$theta_hat[2]
    beta2_lmdpde_h2[i] <- fit_lmdpde_h2$theta_hat[2]

    # Model under H3 (\beta = (-1, -5.5))
    gamma1_mle_h3[i] <- fit_mle_h3$theta_hat_full[3]
    gamma1_lsmle_h3[i] <- fit_lsmle_h3$theta_hat[3]
    gamma1_lmdpde_h3[i] <- fit_lmdpde_h3$theta_hat[3]

    # -----------> Saving the results of hypothesis tests <-----------
    tests <- list(
      # Hypothesis H1 \beta2 = -5.5 x \beta2 =!= -5.5
      wald_statistic(fit_mle$theta_hat_full, fit_mle$fisher_info, beta[2], c(2)),
      robust_wald_statistic(fit_lsmle$theta_hat, fit_lsmle$Vq$H, fit_lsmle$Vq$J, beta[2], c(2)),
      robust_wald_statistic(fit_lmdpde$theta_hat, fit_lmdpde$Vq$H, fit_lmdpde$Vq$J, beta[2], c(2)),
      
      score_statistic(fit_mle_h1$score, fit_mle_h1$fisher_info, c(2)),
      robust_score_statistic(fit_lsmle_h1$score, fit_lsmle_h1$Vq$H, fit_lsmle_h1$Vq$J, c(2)),
      robust_score_statistic(fit_lmdpde_h1$score, fit_lmdpde_h1$Vq$H, fit_lmdpde_h1$Vq$J, c(2)),
      
      gradient_statistic(fit_mle_h1$score, fit_mle$theta_hat_full, beta[2], c(2)),
      robust_gradient_statistic(fit_lsmle_h1_grad$score, fit_lsmle$theta_hat, beta[2], fit_lsmle_h1_grad$Vq$H, fit_lsmle_h1_grad$Vq$J, c(2)),
      robust_gradient_statistic(fit_lmdpde_h1_grad$score, fit_lmdpde$theta_hat, beta[2], fit_lmdpde_h1_grad$Vq$H, fit_lmdpde_h1_grad$Vq$J, c(2)),
      
      # Hypothesis H2 \gamma1 = 5 x gamma1 =!= 5
      wald_statistic(fit_mle$theta_hat_full, fit_mle$fisher_info, gamma[1], c(3)),
      robust_wald_statistic(fit_lsmle$theta_hat, fit_lsmle$Vq$H, fit_lsmle$Vq$J, gamma[1], c(3)),
      robust_wald_statistic(fit_lmdpde$theta_hat, fit_lmdpde$Vq$H, fit_lmdpde$Vq$J, gamma[1], c(3)),

      score_statistic(fit_mle_h2$score, fit_mle_h2$fisher_info, c(3)),
      robust_score_statistic(fit_lsmle_h2$score, fit_lsmle_h2$Vq$H, fit_lsmle_h2$Vq$J, c(3)),
      robust_score_statistic(fit_lmdpde_h2$score, fit_lmdpde_h2$Vq$H, fit_lmdpde_h2$Vq$J, c(3)),

      gradient_statistic(fit_mle_h2$score, fit_mle$theta_hat_full, gamma[1], c(3)),
      robust_gradient_statistic(fit_lsmle_h2$score, fit_lsmle$theta_hat, gamma[1], fit_lsmle_h2$Vq$H, fit_lsmle_h2$Vq$J, c(3)),
      robust_gradient_statistic(fit_lmdpde_h2$score, fit_lmdpde$theta_hat, gamma[1], fit_lmdpde_h2$Vq$H, fit_lmdpde_h2$Vq$J, c(3)),

      # Hypothesis H3 \beta = (-1, -5.5) x \beta =!= (-1, -5.5)
      wald_statistic(fit_mle$theta_hat_full, fit_mle$fisher_info, beta, c(1,2)),
      robust_wald_statistic(fit_lsmle$theta_hat, fit_lsmle$Vq$H, fit_lsmle$Vq$J, beta, c(1,2)),
      robust_wald_statistic(fit_lmdpde$theta_hat, fit_lmdpde$Vq$H, fit_lmdpde$Vq$J, beta, c(1,2)),
      
      score_statistic(fit_mle_h3$score, fit_mle_h3$fisher_info, c(1,2)),
      robust_score_statistic(fit_lsmle_h3$score, fit_lsmle_h3$Vq$H, fit_lsmle_h3$Vq$J, c(1,2)),
      robust_score_statistic(fit_lmdpde_h3$score, fit_lmdpde_h3$Vq$H, fit_lmdpde_h3$Vq$J, c(1,2)),

      gradient_statistic(fit_mle_h3$score, fit_mle$theta_hat_full, beta, c(1,2)),
      robust_gradient_statistic(fit_lsmle_h3$score, fit_lsmle$theta_hat, beta, fit_lsmle_h3$Vq$H, fit_lsmle_h3$Vq$J, c(1,2)),
      robust_gradient_statistic(fit_lmdpde_h3$score, fit_lmdpde$theta_hat, beta, fit_lmdpde_h3$Vq$H, fit_lmdpde_h3$Vq$J, c(1,2)),
      
      # Hyphothesis 4: beta and gamma fixed
      wald_statistic(fit_mle$theta_hat_full, fit_mle$fisher_info, theta, c(1,2,3)),
      robust_wald_statistic(fit_lsmle$theta_hat, fit_lsmle$Vq$H, fit_lsmle$Vq$J, theta, c(1,2,3)),
      robust_wald_statistic(fit_lmdpde$theta_hat, fit_lmdpde$Vq$H, fit_lmdpde$Vq$J, theta, c(1,2,3)),
      
      score_statistic(fit_mle_h4$score, fit_mle_h4$fisher_info, c(1,2,3)),
      robust_score_statistic(fit_lsmle_h4$score, fit_lsmle_h4$Vq$H, fit_lsmle_h4$Vq$J, c(1,2,3)),
      robust_score_statistic(fit_lmdpde_h4$score, fit_lmdpde_h4$Vq$H, fit_lmdpde_h4$Vq$J, c(1,2,3)),

      gradient_statistic(fit_mle_h4$score, fit_mle$theta_hat_full, theta, c(1,2,3)),
      robust_gradient_statistic(fit_lsmle_h4$score, fit_lsmle$theta_hat, theta, fit_lsmle_h4$Vq$H, fit_lsmle_h4$Vq$J, c(1,2,3)),
      robust_gradient_statistic(fit_lmdpde_h4$score, fit_lmdpde$theta_hat, theta, fit_lmdpde_h4$Vq$H, fit_lmdpde_h4$Vq$J, c(1,2,3))
    )
    
    statistics[i, ] <- sapply(tests, `[[`, "statistic")
    pvalues[i, ]    <- sapply(tests, `[[`, "p_value")
  }
  
  final_data <- data.frame(
    n = n,
    sim = 1:N,
    alpha_lsmle,
    alpha_lmdpde,
    alpha_lsmle_h1,
    alpha_lmdpde_h1,
    beta1_mle,
    beta1_lsmle,
    beta1_lmdpde,
    beta2_mle,
    beta2_lsmle,
    beta2_lmdpde,
    gamma1_mle,
    gamma1_lsmle,
    gamma1_lmdpde,
    beta1_mle_h1,
    beta1_lsmle_h1,
    beta1_lmdpde_h1,
    gamma1_mle_h1,
    gamma1_lsmle_h1,
    gamma1_lmdpde_h1,
    beta1_mle_h2,
    beta1_lsmle_h2,
    beta1_lmdpde_h2,
    beta2_mle_h2,
    beta2_lsmle_h2,
    beta2_lmdpde_h2,
    gamma1_mle_h3,
    gamma1_lsmle_h3,
    gamma1_lmdpde_h3,
    statistics,
    pvalues)

  # Save results as txt
  file_name <- paste0("G:/Meu Drive/03 - Sim/03_case_study_1/Results/ScenarioB_NoCont_n", n, ".txt")
  write.table(final_data, file = file_name, row.names = FALSE, col.names = TRUE)
}

# # Plotting the Graph and saving as pdf
# pdf("../ScenarioB_NoCont_n40.pdf", width = 8, height = 6)
#plot(X[, 2], y, xlab = "X", ylab = "y", ylim = c(0, 0.25))
# dev.off()
