#' @title
#' Estimation of BART, BCF, and BCF+CART combo
#'
#' @description
#' 3 function for the estimation of BART, BCF, and BCF+CART combo, respectively
#'
#' @param c : seed 
#' @param data_sample : list of data
#' @param estimated_Y (only for CART function): vector of estimated outcomes
#'
#' @return
#' For BART and BCF functions:
#' - *tau* : vector of individual treatment effect (ITE)
#' For CART function:
#' - *partition* : point estimation of group partition
#'
#' @import bartCause
#' @import bcf
#' @import rpart.plot
#' @import rattle

#########################################################################
# libraries
library(bartCause)
library(bcf)
library(rpart.plot)
library(rattle)

#########################################################################
#    ---     BART  ----
#########################################################################

# function to estimate the BART and to extrapolate the quantities of interest
bart_sample <- function(c, data_sample) {
  set.seed(c)
  
  # 正确访问数据
  T_level <- data_sample$data$T
  X <- data_sample$data$X
  Y_mat <- data_sample$data$Y
  n <- length(T_level)
  
  Y_obs <- unlist(sapply(1:n, function(i) {
    Y_mat[(T_level[i] + 1), i]
  }))
  
  # BART 估计
  bart_fit <- bartCause::bartc(
    response = as.matrix(Y_obs),
    treatment = as.matrix(T_level),
    confounders = as.matrix(X),
    n.samples = 1000,
    n.burn = 1000
  )
  
  # potential outcome
  E_Y_obs <- apply(bart_fit$mu.hat.obs[1,,], 2, mean)
  E_Y_cf <- apply(bart_fit$mu.hat.cf[1,,], 2, mean)
  
  tau <- rep(NA, n)
  tau[T_level == 0] <- E_Y_cf[T_level == 0] - E_Y_obs[T_level == 0]
  tau[T_level == 1] <- E_Y_obs[T_level == 1] - E_Y_cf[T_level == 1]
  
  return(list(tau = tau))
}


#########################################################################
#    ---    BCF    ----
#########################################################################

BCF_sample <- function(c, data_sample) {
  set.seed(c)
  
  T_level <- data_sample$data$T
  X <- as.matrix(data_sample$data$X)
  Y_mat <- data_sample$data$Y
  n <- length(T_level)
  
  Y_obs <- unlist(sapply(1:n, function(i) {
    Y_mat[(T_level[i] + 1), i]
  }))
  
  # 倾向性得分
  p.score <- glm(T_level ~ ., data = as.data.frame(X), family = binomial)
  pihat <- predict(p.score, newdata = as.data.frame(X), type = "response")
  
  result <- invisible(capture.output({
    bcf_fit <- bcf(Y_obs, T_level, X, X, pihat,
                   nburn = 500, nsim = 500)  # ✅ 改小采样数
  }))
  
  return(list(tau = apply(bcf_fit$tau, 2, mean)))
}



#########################################################################
#    ---  Clustering with CART    ----
#########################################################################

CART <- function(c, data_sample, estimated_Y) {
  set.seed(c)
  
  # 正确提取数据
  T_level <- data_sample$data$T
  X <- data_sample$data$X
  tau <- estimated_Y[[c]]$tau  # BCF 或 BART 的第 c 次样本的 tau
  
  # 组合成数据框
  variables <- as.data.frame(cbind(tau, X))
  
  # CART 决策树拟合
  fit.tree <- rpart::rpart(tau ~ ., data = variables, cp = 0.01, minsplit = 10)
  
  # 分组分配
  GroupsAllocation_cart <- c(fit.tree$where)
  
  return(list(partition = GroupsAllocation_cart))
}

