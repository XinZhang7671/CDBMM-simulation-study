generate_X1new <- function(X1, Y_obs, S_cluster) {
  
  g <- length(unique(S_cluster))
  result_list <- list()
  
  X1_new_full <- rep(NA, length(S_cluster))
  
  for (group in 1:g) {
    idx <- which(S_cluster == group)
    
    y_group <- Y_obs[idx]   
    x1_obs <- X1[idx]
    
    target_cor <- runif(1, 0, 1)
    
    y_corr <- sort(rnorm(length(x1_obs)) + rnorm(length(x1_obs), 0, 1e-8))
    x_corr <- target_cor * y_corr + rnorm(length(x1_obs)) * sqrt(1 - target_cor^2)
    
    sy <- sort.int(y_group, index.return = TRUE)
    xtemp <- x1_obs[sy$ix]
    
    rank_xcorr <- rank(x_corr, ties.method = "random")
    rank_xtemp <- rank(xtemp, ties.method = "random")
    
    x1_new_group <- numeric(length(x_corr))
    for (i in 1:length(x_corr)) {
      nearest_idx <- which.min(abs(rank_xtemp - rank_xcorr[i]))
      x1_new_group[i] <- xtemp[nearest_idx]
    }
    
    achieved_cor <- cor(x1_new_group, sort(y_group))
    
    result_list[[paste0("Group_", group)]] <- list(
      x_new = x1_new_group,
      y_new = sort(y_group),
      target_cor = target_cor,
      achieved_cor = achieved_cor
    )
    
    X1_new_full[idx] <- x1_new_group
  }
  
  result_list[["X1_new_full"]] <- X1_new_full
  return(result_list)
}

generate_X2new <- function(X2, Y_obs, S_cluster) {
  
  g <- length(unique(S_cluster))
  result_list <- list()
  
  X2_new_full <- rep(NA, length(S_cluster))
  
  for (group in 1:g) {
    idx <- which(S_cluster == group)
    
    y_group <- Y_obs[idx]  
    x2_obs <- X2[idx]
    
    target_cor <- runif(1, 0, 1)
    
    y_corr <- sort(rnorm(length(x2_obs)) + rnorm(length(x2_obs), 0, 1e-8))
    x_corr <- target_cor * y_corr + rnorm(length(x2_obs)) * sqrt(1 - target_cor^2)
    
    sy <- sort.int(y_group, index.return = TRUE)
    xtemp <- x2_obs[sy$ix]
    
    rank_xcorr <- rank(x_corr, ties.method = "random")
    rank_xtemp <- rank(xtemp, ties.method = "random")
    
    x2_new_group <- numeric(length(x_corr))
    for (i in 1:length(x_corr)) {
      nearest_idx <- which.min(abs(rank_xtemp - rank_xcorr[i]))
      x2_new_group[i] <- xtemp[nearest_idx]
    }
    
    achieved_cor <- cor(x2_new_group, sort(y_group))
    
    result_list[[paste0("Group_", group)]] <- list(
      x_new = x2_new_group,
      y_new = sort(y_group),
      target_cor = target_cor,
      achieved_cor = achieved_cor
    )
    
    X2_new_full[idx] <- x2_new_group
  }
  
  result_list[["X2_new_full"]] <- X2_new_full
  return(result_list)
}


generate_exposure_new <- function(X1_new_full, exposure, S_cluster) {
  
  g <- length(unique(S_cluster))
  result_list <- list()
  
  exposure_new_full <- rep(NA, length(S_cluster))
  
  for (group in 1:g) {
    idx <- which(S_cluster == group)
    
    x1_group <- X1_new_full[idx]  # Real X1_new
    exposure_group <- exposure[idx]  # Real exposure
    
    target_cor <- runif(1, 0, 1)  # 随机指定每组的目标相关性
    
    # Copula processing: matching with exposure_group and x1_group sorting
    y_corr <- sort(rnorm(length(x1_group)) + rnorm(length(x1_group), 0, 1e-8))
    x_corr <- target_cor * y_corr + rnorm(length(x1_group)) * sqrt(1 - target_cor^2)
    
    sy <- sort.int(x1_group, index.return = TRUE)
    xtemp <- exposure_group[sy$ix]
    
    rank_xcorr <- rank(x_corr, ties.method = "random")
    rank_xtemp <- rank(xtemp, ties.method = "random")
    
    exposure_new_group <- numeric(length(x_corr))
    for (i in 1:length(x_corr)) {
      nearest_idx <- which.min(abs(rank_xtemp - rank_xcorr[i]))
      exposure_new_group[i] <- xtemp[nearest_idx]
    }
    
    achieved_cor <- cor(exposure_new_group, sort(x1_group))
    
    result_list[[paste0("Group_", group)]] <- list(
      exposure_new = exposure_new_group,
      x1_new = sort(x1_group),
      target_cor = target_cor,
      achieved_cor = achieved_cor
    )
    
    exposure_new_full[idx] <- exposure_new_group
  }
  
  result_list[["exposure_new_full"]] <- exposure_new_full
  return(result_list)
}



generate_X3new <- function(exposure_new_full, X3, S_cluster) {
  
  g <- length(unique(S_cluster))
  result_list <- list()
  
  X3_new_full <- rep(NA, length(S_cluster))
  
  for (group in 1:g) {
    idx <- which(S_cluster == group)
    
    exposure_group <- exposure_new_full[idx]  # exposure_new
    x3_group <- X3[idx]  # real X3
    
    target_cor <- runif(1, 0, 1)
    
    y_corr <- sort(rnorm(length(exposure_group)) + rnorm(length(exposure_group), 0, 1e-8))
    x_corr <- target_cor * y_corr + rnorm(length(exposure_group)) * sqrt(1 - target_cor^2)
    
    sy <- sort.int(exposure_group, index.return = TRUE)
    xtemp <- x3_group[sy$ix]
    
    rank_xcorr <- rank(x_corr, ties.method = "random")
    rank_xtemp <- rank(xtemp, ties.method = "random")
    
    x3_new_group <- numeric(length(x_corr))
    for (i in 1:length(x_corr)) {
      nearest_idx <- which.min(abs(rank_xtemp - rank_xcorr[i]))
      x3_new_group[i] <- xtemp[nearest_idx]
    }
    
    achieved_cor <- cor(x3_new_group, sort(exposure_group))
    
    result_list[[paste0("Group_", group)]] <- list(
      x_new = x3_new_group,
      exposure_new = sort(exposure_group),
      target_cor = target_cor,
      achieved_cor = achieved_cor
    )
    
    X3_new_full[idx] <- x3_new_group
  }
  
  result_list[["X3_new_full"]] <- X3_new_full
  return(result_list)
}
