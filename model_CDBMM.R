#' @title
#' Model Confounder-Dependent Bayesian Mixture Model (CDBMM)
#'
#' @description
#' Gibbs sampler for the CDBMM estimation
#'
#' @param c : seed 
#' @param data_sample : list of data
#' @param n : sample size
#'
#' @return
#' A list composed by:
#' - *partition* : point estimation of group partition
#' - *atoms* : mean of the cluster for the treatment and control level
#' - *tau* : vector of individual treatment effect (ITE)
#'
#' @import mvtnorm
#' @import CholWishart
#' @import parallel
#' @import truncnorm
#' @import invgamma
#' @import BNPmix


#########################################################################
#        ---  ESTIMATING CDBMM    ----
#########################################################################

# libraries
library(mvtnorm)
library(CholWishart)
library(parallel)
library(truncnorm)
library(invgamma)
library(BNPmix)

#########################################################################
#    ---  Confonder-Dependent Bayesian Mixture Model    ----
#########################################################################

CDBMM_Gibbs <- function(c, n, T_level, X, Y_obs, R, R_burnin, L_0, L_1) {

  
  # ------ prearing variables ------
  T_level = T_level              
  X = X                   
  Y_obs = Y_obs
  
  # number covariates (X)
  n_X=dim(X)[2]+1
  
  # dividing observation in the treatment levels
  # level t=1
  T1=which(T_level==1)  
  n1=length(T1)
  Y_obs_1=Y_obs[T1]
  # level t=0
  T0=which(T_level==0)
  n0=length(T0)
  Y_obs_0=Y_obs[T0]
  
  # ------   hyperparameters   -----
  p_beta=c(0,20)                 # mean and variance of normal distribution for beta param.
  p_sigma=c(5,1)                # shape parameters of inv-gamma for sigma param.
  p_eta=c(0,20)                  # mean and variance of normal distribution for eta param.
  
  # ------   initialitation   -----
  # parameters
  sigma_0=rep(1,L_0)
  sigma_1=rep(1,L_1)
  beta_0=rep(0,n_X*(L_0-1))
  beta_1=rep(0,n_X*(L_1-1))
  eta_0=seq(-3,-3+0.5*(L_0-1),0.5)
  eta_1=seq(-3,-3+0.5*(L_1-1),0.5)
  # cluster allocation variables
  xi_0=sample(1:L_0,n0,replace=TRUE)
  xi_1=sample(1:L_1,n1,replace=TRUE)
  
  # ------   useful quantities   ------   
  # number of units in each cluster
  units_for_clusters_0=rep(0,L_0)
  # Count the occurrences of each clusters in `xi_0`
  units_for_clusters_0[as.numeric(names(table(xi_0)))]=table(xi_0)
  units_for_clusters_1=rep(0,L_1)
  # Count the occurrences of each clusters in `xi_1`
  units_for_clusters_1[as.numeric(names(table(xi_1)))]=table(xi_1) 
  
  # eta corresponding to each unit
  # Take each index value from `xi_0` and extract the corresponding value from `eta_0`
  eta_sogg_0=eta_0[xi_0]
  # Take each index value from `xi_1` and extract the corresponding value from `eta_1`.
  eta_sogg_1=eta_1[xi_1]
  
  # latent variable for data augmentation in probit regression 
  Z_0=rep(list(rep(NA,L_0)),n0)
  Z_1=rep(list(rep(NA,L_1)),n1)
  
  # ------   useful functions   -----
  
  # probit stick-breaking weights 
  omega = function(X, beta) {
    L = length(beta) / n_X
    alpha_vals = sapply(1:L, function(l) {
      val = pnorm(beta[(l*n_X - (n_X - 1)):(l*n_X)] %*% c(1, X))
      val = min(max(val, 1e-6), 1 - 1e-6)  # 防止极端值
      return(val)
    })
    alpha = c(alpha_vals, 1)
    
    omega_raw = c(alpha[1],
                  sapply(2:(L + 1), function(l) alpha[l] * prod(1 - alpha[1:(l - 1)])))
    
    omega_safe = omega_raw / sum(omega_raw)
    return(omega_safe)
  }
  
  # point estimate parition
  estimation_partition<-function(xi_0,xi_1){
    
    # using Variation of Information as loss function
    clusters_0=partition.BNPdens(list(clust=t(xi_0)),dist = "VI")$partitions[1,]
    clusters_1=partition.BNPdens(list(clust=t(xi_1)),dist = "VI")$partitions[1,]
    
    # alternativly: Binder's loss
    #clusters_0=partition.BNPdens(list(clust=t(xi_0)),dist = "Binder")$partitions
    #clusters_1=partition.BNPdens(list(clust=t(xi_1)),dist = "Binder")$partitions
    
    clusters=cbind(clusters_0,clusters_1)
    dimnames(clusters)=NULL
    return(clusters)
  }

  # funtion to impute outcomes
  posterior_atoms <- function(partition_iterations, R, R_burnin) {
    
    p = list(
      p_0 = rep(NA, max(partition_iterations[,1])),
      p_1 = rep(NA, max(partition_iterations[,2]))
    )
    
    p_MC = list(
      p_0_MC = vector("list", max(partition_iterations[,1])),
      p_1_MC = vector("list", max(partition_iterations[,2]))
    )
    
    # T0 eta
    for (g in 1:max(partition_iterations[,1])) {
      units = intersect(T0, which(partition_iterations[,1] == g))
      # mean T0 eta
      p$p_0[g] = mean(sapply(1:length(units), function(x) 
        mean(sapply(1:(R - R_burnin), function(j) 
          post_eta[cluster_allocation_0[units[x], j], R_burnin + j]))))
      
      p_MC$p_0_MC[[g]] = sapply(1:length(units), function(x) 
        sapply(1:(R - R_burnin), function(j) 
          post_eta[cluster_allocation_0[units[x], j], R_burnin + j]))
    }
    
    # T1 eta
    for (g in 1:max(partition_iterations[,2])) {
      units = intersect(T1, which(partition_iterations[,2] == g))
      # mean T1 eta
      p$p_1[g] = mean(sapply(1:length(units), function(x) 
        mean(sapply(1:(R - R_burnin), function(j) 
          post_eta[L_0 + cluster_allocation_1[units[x], j], R_burnin + j]))))
      
      p_MC$p_1_MC[[g]] = sapply(1:length(units), function(x) 
        sapply(1:(R - R_burnin), function(j) 
          post_eta[L_0 + cluster_allocation_1[units[x], j], R_burnin + j]))
    }
    
    return(list(mean_p = p, full_p = p_MC))
  }
  
  
  # ------   saving informations   -----
  
  # empty matrix where save all the informations for each iteration
  post_eta=matrix(NA,ncol=R,nrow=L_0+L_1)
  post_var=matrix(NA,ncol=R,nrow=L_0+L_1)
  cluster_allocation_0=matrix(NA,ncol=R-R_burnin,nrow=n)
  cluster_allocation_1=matrix(NA,ncol=R-R_burnin,nrow=n)
  post_beta=matrix(NA,ncol=R, nrow=length(beta_0)+length(beta_1))
  Y0_imp=matrix(NA,ncol=R-R_burnin,nrow=n)
  Y1_imp=matrix(NA,ncol=R-R_burnin,nrow=n)
  
  # -----   updating parameters and variables at each itearation   -----
  tau_samples <- numeric(R - R_burnin)
  
  for (r in 1:R){
    
    #######################################################
    # ----------- Cluster Specific Parameters  ------------
    #######################################################
    
    # -----   ETA: mean of normal distributions (atoms)   -----
    # eta_0 = eta for treatment level 0
    eta_0=sapply(1:L_0,function(s) 
      rnorm(1,mean=1/(units_for_clusters_0[s]/sigma_0[s]+1/p_eta[2])*(sum(Y_obs_0[xi_0==s])/sigma_0[s]+p_eta[1]/p_eta[2]),
            sd=sqrt(1/(units_for_clusters_0[s]/sigma_0[s]+1/p_eta[2]))))
    
    # eta_1 = eta for treatment level 1
    eta_1=sapply(1:L_1,function(s) 
      rnorm(1,mean=1/(units_for_clusters_1[s]/sigma_1[s]+1/p_eta[2])*(sum(Y_obs_1[xi_1==s])/sigma_1[s]+p_eta[1]/p_eta[2]),
            sd=sqrt(1/(units_for_clusters_1[s]/sigma_1[s]+1/p_eta[2]))))
    
    eta_sogg_0=eta_0[xi_0]
    eta_sogg_1=eta_1[xi_1]
    
    # -----   SIGMA: variance of normal distributions (atoms)   -----
    # sigma_0 = sigma for treatment level 0
    sigma_0=sapply(1:L_0, function(l) rinvgamma(1,p_sigma[1]+sum(xi_0==l)/2,
                                                p_sigma[2]+sum((Y_obs_0[xi_0==l]-eta_0[l])^2)/2))
    # sigma_1 = sigma for treatment level 1
    sigma_1=sapply(1:L_0, function(l) rinvgamma(1,p_sigma[1]+sum(xi_1==l)/2,
                                                p_sigma[2]+sum((Y_obs_1[xi_1==l]-eta_1[l])^2)/2))
    
    ##############################################
    # ----------- Cluster Allocation  ------------
    ##############################################
    
    # -----   OMEGA: weights treatment-specific   -----
    # omega_0 = omega for treatment level 0
    omega_0=sapply(1:n0, function(i) omega(X=X[T0[i],],beta=beta_0))
    # omega_1 = omega for treatment level 1
    omega_1=sapply(1:n1, function(i) omega(X=X[T1[i],],beta=beta_1))
    
    # -----   LATENT VARIABLE for cluster allocation   -----
    # xi_0 = latent variable for treatment level 0
    dmn_0=sapply(1:n0, function(i) sapply(1:L_0, function(l) dnorm(Y_obs_0[i], eta_0[l], sqrt(sigma_0[l]), log=TRUE)))+log(omega_0)
    dmn_0[which(is.nan(dmn_0))]=-100
    xi=sapply(1:n0, function(i) rmultinom(1,1,exp(dmn_0[,i])))
    xi_0=sapply(1:n0, function(i) xi[,i]%*%(1:L_0))
    units_for_clusters_0=apply(xi, 1, sum)
    
    # xi_1 = latent variable for treatment level 1
    dmn_1=sapply(1:n1, function(i) sapply(1:L_1, function(l) dnorm(Y_obs_1[i], eta_1[l], sqrt(sigma_1[l]), log=TRUE)))+log(omega_1)
    dmn_1[which(is.nan(dmn_1))]=-100
    xi=sapply(1:n1, function(i) rmultinom(1,1,exp(dmn_1[,i])))
    xi_1=sapply(1:n1, function(i) xi[,i]%*%(1:L_1))
    units_for_clusters_1=apply(xi, 1, sum)
    
    
    ###############################################
    # ----------- Augmentation Scheme  ------------
    ###############################################
    
    # -----   LATENT VARIABLE: Z for probit regression   -----
    # Z_0 = latent variable for treatment level 0
    # building intermediate values
    pesi_0 = t(omega_0)
    mu_z = cbind((pesi_0[, 1]), (pesi_0[, 2] / (1 - pesi_0[, 1])))
    
    if (L_0 > 3) {
      mu_z = cbind(mu_z, sapply(3:(L_0 - 1), function(l) (pesi_0[, l] / (1 - apply(pesi_0[, 1:(l - 1)], 1, sum)))))
    }
    
    # 确保mu_z的值在合理范围
    mu_z[is.nan(mu_z)] = 1
    mu_z[mu_z > 1] = 1
    mu_z = mu_z - 9.9e-15 * (mu_z > (1 - 1e-16))
    
    # 加入平滑处理，确保没有为零的概率
    mu_z = pmax(pmin(mu_z, 1 - 1e-10), 1e-10)
    
    # updating Z_0
    for (i in 1:n0) {
      for (l in 1:(min(xi_0[i], L_0 - 1))) {
        if (l > 1) {
          if (l < xi_0[i]) {
            Z_0[[i]][l] = rtruncnorm(1, b = 0, mean = qnorm(mu_z[i, l]))
          } else {
            Z_0[[i]][l] = rtruncnorm(1, a = 0, mean = qnorm(mu_z[i, l]))
          }
        } else {
          if (l < xi_0[i]) {
            Z_0[[i]][l] = rtruncnorm(1, b = 0, mean = qnorm(mu_z[i, l]))
          } else {
            Z_0[[i]][l] = rtruncnorm(1, a = 0, mean = qnorm(mu_z[i, l]))
          }
        }
      }
    }
    
    # Z_1 = latent variable for treatment level 1
    # building intermediate values
    pesi_1 = t(omega_1)
    mu_z = cbind((pesi_1[, 1]), (pesi_1[, 2] / (1 - pesi_1[, 1])))
    
    if (L_1 > 3) {
      mu_z = cbind(mu_z, sapply(3:(L_1 - 1), function(l) (pesi_1[, l] / (1 - apply(pesi_1[, 1:(l - 1)], 1, sum)))))
    }
    
    # Ensure that the value of mu_z is within a reasonable range
    mu_z[is.nan(mu_z)] = 1
    mu_z[mu_z > 1] = 1
    mu_z = mu_z - 9.9e-15 * (mu_z > (1 - 1e-16))
    
    # Add smoothing to ensure there is no zero probability
    mu_z = pmax(pmin(mu_z, 1 - 1e-10), 1e-10)
    
    # updating Z_1
    for (i in 1:n1) {
      for (l in 1:(min(xi_1[i], L_1 - 1))) {
        if (l > 1) {
          if (l < xi_1[i]) {
            Z_1[[i]][l] = rtruncnorm(1, b = 0, mean = qnorm(mu_z[i, l]))
          } else {
            Z_1[[i]][l] = rtruncnorm(1, a = 0, mean = qnorm(mu_z[i, l]))
          }
        } else {
          if (l < xi_1[i]) {
            Z_1[[i]][l] = rtruncnorm(1, b = 0, mean = qnorm(mu_z[i, l]))
          } else {
            Z_1[[i]][l] = rtruncnorm(1, a = 0, mean = qnorm(mu_z[i, l]))
          }
        }
      }
    }
    
    
    ########################################################
    # ----------- Confounder-Dependent Weights  ------------
    ########################################################
    
    # beta_0 = beta for treatment level 0
    clusters_temp = which(units_for_clusters_0 != 0)
    if (max(clusters_temp) == L_0) {
      clusters_temp = clusters_temp[-length(clusters_temp)]
    }
    
    for (l in clusters_temp) {
      val = which(xi_0 >= l)
      
      if (length(val) > n_X) {
        z_tilde = unlist(sapply(val, function(i) Z_0[[i]][l]))
        x_tilde = cbind(rep(1, length(val)), matrix(X[T0[val], ], ncol = n_X - 1))
        
        XtX = t(x_tilde) %*% x_tilde
        epsilon <- 1e-6
        A = diag(n_X) / p_beta[2] + XtX + diag(epsilon, n_X)  # 加入扰动项确保数值稳定
        
        if (kappa(A) < 1e10) {  # Determine whether it is stable and reversible
          V = solve(A)
          beta_0[(l * n_X - (n_X - 1)):(l * n_X)] = rmvnorm(
            1,
            V %*% (1 / p_beta[2] * (diag(n_X) %*% rep(p_beta[1], n_X)) + t(x_tilde) %*% z_tilde),
            V
          )[, 1:n_X]
        } else {
         
          beta_0[(l * n_X - (n_X - 1)):(l * n_X)] = rmvnorm(
            1,
            rep(p_beta[1], n_X),
            diag(n_X) / p_beta[2]
          )[, 1:n_X]
        }
      } else {
        beta_0[(l * n_X - (n_X - 1)):(l * n_X)] = rmvnorm(
          1,
          rep(p_beta[1], n_X),
          diag(n_X) / p_beta[2]
        )[, 1:n_X]
      }
    }
    
    # beta_1 = beta for treatment level 1
    clusters_temp = which(units_for_clusters_1 != 0)
    if (max(clusters_temp) == L_1) {
      clusters_temp = clusters_temp[-length(clusters_temp)]
    }
    
    for (l in clusters_temp) {
      val = which(xi_1 >= l)
      
      if (length(val) > n_X) {
        z_tilde = unlist(sapply(val, function(i) Z_1[[i]][l]))
        x_tilde = cbind(rep(1, length(val)), matrix(X[T1[val], ], ncol = n_X - 1))
        
        XtX = t(x_tilde) %*% x_tilde
        epsilon <- 1e-6
        A = diag(n_X) / p_beta[2] + XtX + diag(epsilon, n_X)
        
        if (kappa(A) < 1e10) {
          V = solve(A)
          beta_1[(l * n_X - (n_X - 1)):(l * n_X)] = rmvnorm(
            1,
            V %*% (1 / p_beta[2] * (diag(n_X) %*% rep(p_beta[1], n_X)) + t(x_tilde) %*% z_tilde),
            V
          )[, 1:n_X]
        } else {
          beta_1[(l * n_X - (n_X - 1)):(l * n_X)] = rmvnorm(
            1,
            rep(p_beta[1], n_X),
            diag(n_X) / p_beta[2]
          )[, 1:n_X]
        }
      } else {
        beta_1[(l * n_X - (n_X - 1)):(l * n_X)] = rmvnorm(
          1,
          rep(p_beta[1], n_X),
          diag(n_X) / p_beta[2]
        )[, 1:n_X]
      }
    }
    
    ##############################################################
    # ----------- imputing missing outcomes allocation -----------
    ##############################################################
    
    # -----   Y_MIS: missing outcome: Y(1-t)   -----
    # here we are just imputing the cluster allocation (no the Y values)
    
    # estimate cluster allocation for each iteration
    if(r>R_burnin){
      
      # level t=0 observed --> level T=1 missing
      om_0=sapply(T0, function(i) omega(X=X[i,],beta=beta_1))
      om_0[which(is.nan(om_0))]=0
      clusters_Ymiss_0=sapply(1:n0, function(i) (1:L_1)%*%(rmultinom(1,1,om_0[,i])))
      
      # level t=1 observed --> level T=0 missing
      om_1=sapply(T1, function(i) omega(X=X[i,],beta=beta_0))
      om_1[which(is.nan(om_1))]=0
      clusters_Ymiss_1=sapply(1:n1, function(i) (1:L_0)%*%(rmultinom(1,1,om_1[,i])))
    }
    
    # -----   saving information   -----
    # parameters
    post_eta[,r]=c(eta_0,eta_1)
    post_var[,r]=c(sigma_0,sigma_1)
    post_beta[,r]=c(beta_0,beta_1)
    # cluster allocation for each iteration
    if(r>R_burnin){
      cluster_allocation_0[T0,r-R_burnin]=xi_0
      cluster_allocation_0[T1,r-R_burnin]=clusters_Ymiss_1
      cluster_allocation_1[T1,r-R_burnin]=xi_1
      cluster_allocation_1[T0,r-R_burnin]=clusters_Ymiss_0
      
      Y0_imp[,r-R_burnin]=rnorm(n,eta_0[cluster_allocation_0[,r-R_burnin]],
                                sigma_0[cluster_allocation_0[,r-R_burnin]])
      Y1_imp[,r-R_burnin]=rnorm(n,eta_1[cluster_allocation_1[,r-R_burnin]],
                                sigma_1[cluster_allocation_1[,r-R_burnin]])
    }
    tau_samples[r - R_burnin] <- mean(Y1_imp[,r - R_burnin] - Y0_imp[,r - R_burnin])
  }
  
    #
    if (r%%500==0) {
      print(paste0(r,"/",R," iterations"))
      flush.console()
    }
  }
  
  ############################################
  # -----   point estimation partition   -----
  ############################################
  
  # partition: cluster allocation
  partition=estimation_partition(cluster_allocation_0,cluster_allocation_1)
  # posterior mean of tau=Y(1)-Y(0)
  atoms=posterior_atoms(partition, R, R_burnin)
  Y0_imp_median=apply(Y0_imp,1,median)
  Y1_imp_median=apply(Y1_imp,1,median)
  
  tau=Y1_imp_median-Y0_imp_median
  tau_=rep(NA,n)
  tau_[T0]=Y1_imp_median[T0]-Y_obs[T0]
  tau_[T1]=Y_obs[T1]-Y0_imp_median[T1]
  
  
  print(paste0("sample ",c, " done"))
  
  return(list(# remove the "#" if you want the chains of these parameters
    #post_eta=post_eta,  
    #post_var=post_var,
    #post_beta=post_beta,
    partition=partition,
    #atoms=atoms,
    tau=tau,
    tau_samples = tau_samples,
    Y0_imp=Y0_imp,
    Y1_imp=Y1_imp
    #tau_=tau_
    ))



