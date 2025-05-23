library(mvtnorm)
library(cluster)
library(future.apply)

func_5cov_3g <- function(seed, samples, R, R_burnin, eta_0, eta_1, 
                         sigma_0, sigma_1, p_T) {
  
  # Set the seed for the random number generation at the beginning of each function call
  set.seed(seed)
  
  # Default values for L_0 and L_1
  n = 500
  L_0 = 12
  L_1 = 12
  
  # Generate treatment assignment and exposure
  erand <- runif(n)  
  T <- (erand > (1 - p_T)) * 1
  exposure <- qnorm(erand, mean = 15, sd = 3)
  
  # Simulate covariates
  X1 <- rnorm(n, mean = 20, sd = 5)
  X2 <- rnorm(n, mean = 50, sd = 3)
  X3 <- rnorm(n, mean = 5, sd = 0.3)
  X4 <- rbinom(n, size = 1, prob = 0.6)
  X5 <- rbinom(n, size = 1, prob = 0.4)
  X <- cbind(X1, X2, X3, X4, X5)
  
  X_df <- data.frame(X)
  X_df[, 4:5] <- lapply(X_df[, 4:5], factor)  # Binomial as factor
  
  # Gower distance for clustering
  gower_dist <- daisy(X_df, metric = "gower")
  pam_result <- pam(gower_dist, k = 3)
  S_cluster <- pam_result$clustering
  print(table(S_cluster))
  
  # Simulate potential outcomes
  Y <- sapply(1:n, function(i) rmvnorm(1, 
                                       mean = c(eta_0[S_cluster[i]], eta_1[S_cluster[i]]),
                                       sigma = diag(c(sigma_0[S_cluster[i]], sigma_1[S_cluster[i]]))))
  T_level <- T
  Y_obs <- unlist(sapply(1:n, function(i) Y[(T_level[i] + 1), i]))
  
  # Generate new covariates
  X1NEW <- generate_X1new(X1 = X1, Y_obs = Y_obs, S_cluster = S_cluster)
  X2NEW <- generate_X2new(X2 = X2, Y_obs = Y_obs, S_cluster = S_cluster)
  
  exposureNEW <- generate_exposure_new(X1_new_full = X1NEW$X1_new_full, 
                                       exposure = exposure, S_cluster = S_cluster)
  
  X3NEW <- generate_X3new(X3 = X3, exposure_new_full = exposureNEW$exposure_new_full,
                          S_cluster = S_cluster)
  
  X_NEW <- cbind(X1NEW$X1_new_full, X2NEW$X2_new_full, X3NEW$X3_new_full, 
                 X4, X5, exposureNEW$exposure_new_full) 
  
  # Set the future plan with future.seed = TRUE (ensure future.apply handles random number generation)
  plan(multisession, workers = detectCores() - 1) 
  
  # Run the model with parallel processing, now future.seed is not passed here
  CDBMM_scenario_result <- future_lapply(seq_len(samples), function(c) {
   
    print(paste("Running task", c))
    
    CDBMM_Gibbs(c=c, T_level = T_level, X = X_NEW, Y_obs = Y_obs, R = R, R_burnin = R_burnin, 
                L_0 = L_0, L_1 = L_1, n=n)  # Explicitly passing L_0, L_1, and n
  },
  future.seed = TRUE
  )
  
  # Reset to sequential after parallel processing
  plan(sequential)
  
  return(list(
    CDBMM_results = CDBMM_scenario_result,
    S_cluster = S_cluster,
    
    input_data = list(
      X_NEW = X_NEW,
      T_level = T_level,
      Y = Y
      ),
    
    plot_data = list(
      X1NEW = X1NEW,
      X2NEW = X2NEW,
      exposureNEW = exposureNEW,
      X3NEW = X3NEW)
  ))
}




