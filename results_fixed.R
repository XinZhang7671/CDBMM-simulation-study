source("advanced_gencor.R")
source("corr_plot.R")
source("model_CDBMM.R")
source("func_5cov_3g.R")

# 1. Create a save directory
dir.create("scenario_results", showWarnings = FALSE)

# 2. Set all scenario parameter lists
params_list <- list(
  list(seed = 1, samples = 100, R = 3000, R_burnin = 2000,
       eta_0 = c(2, 2, 2), eta_1 = c(3, 3, 3), 
       sigma_0 = rep(0.3, 3), sigma_1 = rep(0.3, 3), p_T = 0.5),
  
  list(seed = 2, samples = 100, R = 3000, R_burnin = 2000, 
       eta_0 = c(2, 2, 2), eta_1 = c(3, 3, 3), 
       sigma_0 = rep(0.3, 3), sigma_1 = rep(0.3, 3), p_T = 0.3),
  
  list(seed = 3, samples = 100, R = 3000, R_burnin = 2000,
       eta_0 = c(0, 1, 2), eta_1 = c(1, 3, 5), 
       sigma_0 = rep(0.3, 3), sigma_1 = rep(0.3, 3), p_T = 0.5),
  
  list(seed = 4, samples = 100, R = 3000, R_burnin = 2000, 
       eta_0 = c(0, 1, 2), eta_1 = c(1, 3, 5), 
       sigma_0 = rep(0.3, 3), sigma_1 = rep(0.3, 3), p_T = 0.3),
  
  list(seed = 5, samples = 100, R = 3000, R_burnin = 2000, 
       eta_0 = c(1, 3, 5), eta_1 = c(0, 1, 3), 
       sigma_0 = rep(0.3, 3), sigma_1 = rep(0.3, 3), p_T = 0.5),
  
  list(seed = 6, samples = 100, R = 3000, R_burnin = 2000, 
       eta_0 = c(1, 3, 5), eta_1 = c(0, 1, 3), 
       sigma_0 = rep(0.3, 3), sigma_1 = rep(0.3, 3), p_T = 0.3),
  
  list(seed = 7, samples = 100, R = 3000, R_burnin = 2000,
       eta_0 = c(1,1.75,2.5), eta_1 = c(1.5,2,2.5), 
       sigma_0 = rep(0.3, 3), sigma_1 = rep(0.3, 3), p_T = 0.5),
  
  list(seed = 8, samples = 100, R = 3000, R_burnin = 2000,
       eta_0 = c(1,1.75,2.5), eta_1 = c(1.5,2,2.5), 
       sigma_0 = rep(0.3, 3), sigma_1 = rep(0.3, 3), p_T = 0.3)
)

# 3. Initialize result storage
results_5cov_3g <- list()
error_seeds <- c()

# 4. Traverse each scenario and support breakpoint continuation
for (i in seq_along(params_list)) {
  seed_val <- params_list[[i]]$seed
  file_path <- paste0("scenario_results/Scenario_", seed_val, ".RDS")
  
  # Skip if the file already exists
  if (file.exists(file_path)) {
    cat("Skipping seed =", seed_val, "(already exists)\n")
    result <- readRDS(file_path)  
  } else {
    cat("Running task for seed =", seed_val, "\n")
    
    result <- tryCatch({
      res <- do.call(func_5cov_3g, params_list[[i]])
      cat(" -> Success for seed", seed_val, "\n")
      saveRDS(res, file = file_path)
      res
    }, error = function(e) {
      cat(" -> Error for seed", seed_val, ":", conditionMessage(e), "\n")
      error_seeds <- c(error_seeds, seed_val)
      NULL
    })
  }
  
  results_5cov_3g[[paste0("Scenario_", seed_val)]] <- result
}


save(results_5cov_3g, file = "results_5cov_3g.RData")

# 6. Display operation status
cat("✅  all Scenario  finish\n")
if (length(error_seeds) > 0) {
  cat("⚠️ wrong seed label：", paste(error_seeds, collapse = ", "), "\n")
}

plot_X1NEW(results_5cov_3g$Scenario_1$plot_data$X1NEW)
plot_X2NEW(results_5cov_3g$Scenario_1$plot_data$X2NEW)
plot_exposureNEW(results_5cov_3g$Scenario_1$plot_data$exposureNEW)
plot_X3NEW(results_5cov_3g$Scenario_1$plot_data$X3NEW)

 

