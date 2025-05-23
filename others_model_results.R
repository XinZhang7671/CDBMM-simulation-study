source("BART_BCF_CART_model.R")

library(future)
library(future.apply)


plan(multisession, workers = 6)


data_sample <- lapply(results_5cov_3g, function(scen) {
  X <- scen$input_data$X_NEW
  
  colnames(X) <- c("X1", "X2", "X3", "X4", "X5", "exposure")
  
  list(
    data = list(
      X = X,
      T = scen$input_data$T_level,
      Y = scen$input_data$Y
    )
  )
})

samples <- 100 

#########################################################################
#        ---  BART    ----
#########################################################################

BART_results <- lapply(1:8, function(i) {
  scenario_i <- data_sample[[i]]
  future_lapply(1:samples, function(s) {
    bart_sample(s, data_sample = scenario_i)
  }, future.seed = TRUE)
})

names(BART_results) <- paste0("BART_scenario_", 1:8)
save(BART_results, file = "BART_results.RData")

#########################################################################
#        ---  BCF    ----
#########################################################################

plan(multisession, workers = availableCores() - 1)
BCF_results <- vector("list", 8)

for (i in 1:8) {
  save_file <- paste0("BCF_result_scenario", i, ".RData")
  

  if (file.exists(save_file)) {
    load(save_file)  
  } else {
    result_i <- vector("list", samples)
  }
  
  scenario_i <- data_sample[[i]]
  
  not_done <- which(sapply(result_i, is.null))
  
  if (length(not_done) > 0) {
    new_results <- future_lapply(not_done, function(s) {
      tryCatch({
        BCF_sample(s, data_sample = scenario_i)
      }, error = function(e) NULL)
    }, future.seed = TRUE)
    
    for (j in seq_along(not_done)) {
      result_i[[not_done[j]]] <- new_results[[j]]
    }
    
    save(result_i, file = save_file)
  }
  
  BCF_results[[i]] <- result_i
}

save(BCF_results, file = "BCF_results_all.RData")


file.remove(list.files(pattern = "mod_trees.*\\.txt"))
file.remove(list.files(pattern = "con_trees.*\\.txt"))

#########################################################################
#        ---  CART    ----
#########################################################################

CART_results <- lapply(1:8, function(i) {
  scenario_i <- data_sample[[i]]
  estimated_Y_i <- BART_results[[i]]  # or BART_results[[i]]
  
  future_lapply(1:samples, function(s) {
    CART(s, data_sample = scenario_i, estimated_Y = estimated_Y_i)
  }, future.seed = TRUE)
})

names(CART_results) <- paste0("CART_scenario_", 1:8)

save(CART_results, file = "CART_results.RData")

