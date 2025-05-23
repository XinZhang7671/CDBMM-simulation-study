load("results_5cov_3g.RData")

#libraries
library(RColorBrewer)
library(knitr)
library(dplyr)
library(gridExtra)
library(grid)
library(ggplot2)
library(mcclust)
library(parallel)


scenario_1 <- results_5cov_3g$Scenario_1
scenario_2 <- results_5cov_3g$Scenario_2
scenario_3 <- results_5cov_3g$Scenario_3
scenario_4 <- results_5cov_3g$Scenario_4
scenario_5 <- results_5cov_3g$Scenario_5
scenario_6 <- results_5cov_3g$Scenario_6
scenario_7 <- results_5cov_3g$Scenario_7
scenario_8 <- results_5cov_3g$Scenario_8

# True tau (ITE)
true_tau_1=scenario_1$input_data$Y[2,]-scenario_1$input_data$Y[1,]
true_tau_2=scenario_2$input_data$Y[2,]-scenario_2$input_data$Y[1,]
true_tau_3=scenario_3$input_data$Y[2,]-scenario_3$input_data$Y[1,]
true_tau_4=scenario_4$input_data$Y[2,]-scenario_4$input_data$Y[1,]
true_tau_5=scenario_5$input_data$Y[2,]-scenario_5$input_data$Y[1,]
true_tau_6=scenario_6$input_data$Y[2,]-scenario_6$input_data$Y[1,]
true_tau_7=scenario_7$input_data$Y[2,]-scenario_7$input_data$Y[1,]
true_tau_8=scenario_8$input_data$Y[2,]-scenario_8$input_data$Y[1,]

# True ATE
true_ATE_1=mean(true_tau_1)
true_ATE_2=mean(true_tau_2)
true_ATE_3=mean(true_tau_3)
true_ATE_4=mean(true_tau_4)
true_ATE_5=mean(true_tau_5)
true_ATE_6=mean(true_tau_6)
true_ATE_7=mean(true_tau_7)
true_ATE_8=mean(true_tau_8)

# All samples of CDBMM ITE (500X100)
tau_1_matrix=sapply(scenario_1[["CDBMM_results"]], function(res) res[["tau"]])
tau_2_matrix=sapply(scenario_2[["CDBMM_results"]], function(res) res[["tau"]])
tau_3_matrix=sapply(scenario_3[["CDBMM_results"]], function(res) res[["tau"]])
tau_4_matrix=sapply(scenario_4[["CDBMM_results"]], function(res) res[["tau"]])
tau_5_matrix=sapply(scenario_5[["CDBMM_results"]], function(res) res[["tau"]])
tau_6_matrix=sapply(scenario_6[["CDBMM_results"]], function(res) res[["tau"]])
tau_7_matrix=sapply(scenario_7[["CDBMM_results"]], function(res) res[["tau"]])
tau_8_matrix=sapply(scenario_8[["CDBMM_results"]], function(res) res[["tau"]])

# CDBMM ATE (1:100)
CDBMM_ATE_1 <- colMeans(tau_1_matrix)
CDBMM_ATE_2 <- colMeans(tau_2_matrix)
CDBMM_ATE_3 <- colMeans(tau_3_matrix)
CDBMM_ATE_4 <- colMeans(tau_4_matrix)
CDBMM_ATE_5 <- colMeans(tau_5_matrix)
CDBMM_ATE_6 <- colMeans(tau_6_matrix)
CDBMM_ATE_7 <- colMeans(tau_7_matrix)
CDBMM_ATE_8 <- colMeans(tau_8_matrix)

# bias_ATE 
bias_ATE_1=CDBMM_ATE_1 - true_ATE_1
bias_ATE_2=CDBMM_ATE_2 - true_ATE_2
bias_ATE_3=CDBMM_ATE_3 - true_ATE_3
bias_ATE_4=CDBMM_ATE_4 - true_ATE_4
bias_ATE_5=CDBMM_ATE_5 - true_ATE_5
bias_ATE_6=CDBMM_ATE_6 - true_ATE_6
bias_ATE_7=CDBMM_ATE_7 - true_ATE_7
bias_ATE_8=CDBMM_ATE_8 - true_ATE_8

# mse_ATE 
mse_ATE_1=(CDBMM_ATE_1 - true_ATE_1)^2
mse_ATE_2=(CDBMM_ATE_2 - true_ATE_2)^2
mse_ATE_3=(CDBMM_ATE_3 - true_ATE_3)^2
mse_ATE_4=(CDBMM_ATE_4 - true_ATE_4)^2
mse_ATE_5=(CDBMM_ATE_5 - true_ATE_5)^2
mse_ATE_6=(CDBMM_ATE_6 - true_ATE_6)^2
mse_ATE_7=(CDBMM_ATE_7 - true_ATE_7)^2
mse_ATE_8=(CDBMM_ATE_8 - true_ATE_8)^2




# RAND index for clusters obtained with CDBMM
true_clusters <- list(
  rep(1, length(scenario_1$S_cluster)),         # Scenario 1: homogeneous
  rep(1, length(scenario_2$S_cluster)),         # Scenario 2: homogeneous
  scenario_3$S_cluster,                # start with Scenario 3: heterogeneous
  scenario_4$S_cluster,
  scenario_5$S_cluster,
  scenario_6$S_cluster,
  scenario_7$S_cluster,
  scenario_8$S_cluster
)

samples=100
L_0=12
rand_CDBMM_all <- lapply(1:8, function(i) {
  true_cluster <- true_clusters[[i]]
  CDBMM_result <- results_5cov_3g[[i]]$CDBMM_results
  
  sapply(1:samples, function(s) {
    pred_cluster <- CDBMM_result[[s]]$partition
    mcclust::arandi(true_cluster, pred_cluster[,1] * (pred_cluster[,2] + L_0))
  })
})

for (i in 1:2) {
  all_single_group <- all(sapply(results_5cov_3g[[i]]$CDBMM_results, function(res) {
    length(unique(res$partition[,1] * (res$partition[,2] + L_0))) == 1
  }))
  
  if (all_single_group) {
    rand_CDBMM_all[[i]] <- rep(1, samples)
  }
}

names(rand_CDBMM_all) <- paste0("scenario_", 1:8) 

mean_bias_CDBMM <- sapply(1:8, function(i) mean(get(paste0("bias_ATE_", i))))
mean_bias_CDBMM

mean_mse_CDBMM <- sapply(1:8, function(i) mean(get(paste0("mse_ATE_", i))))
mean_mse_CDBMM


#########################################################################
#        ---  CDBMM ATE boxplot   ----
#########################################################################
bias_data <- data.frame(
  BIAS = c(bias_ATE_1, bias_ATE_2, bias_ATE_3,
          bias_ATE_4, bias_ATE_5, bias_ATE_6,
          bias_ATE_7, bias_ATE_8),
  Scenario = factor(rep(paste0("Scenario ", 1:8), each = length(bias_ATE_1)))
)

b <- ggplot(bias_data, aes(x = Scenario, y = BIAS, fill = Scenario)) +
  geom_boxplot(color = "gray30", outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Scenario), width = 0.2, alpha = 0.4, size = 1) +
  scale_fill_brewer(palette = "Pastel2") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Scenario",
       y = "Bias") +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major.y = element_line(color = "grey90", size = 0.4),  # background horizontal lines
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

ggsave("Bias_Boxplot.pdf", plot = b, width = 8, height = 5.5, dpi = 300)

mse_data <- data.frame(
  MSE = c(mse_ATE_1, mse_ATE_2, mse_ATE_3,
          mse_ATE_4, mse_ATE_5, mse_ATE_6,
          mse_ATE_7, mse_ATE_8),
  Scenario = factor(rep(paste0("Scenario ", 1:8), each = length(mse_ATE_1)))
)

p <- ggplot(mse_data, aes(x = Scenario, y = MSE, fill = Scenario)) +
  geom_boxplot(color = "gray30", outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Scenario), width = 0.2, alpha = 0.4, size = 1) +
  scale_fill_brewer(palette = "Pastel2") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = "Scenario",
       y = "Mean Squared Error (MSE)") +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major.y = element_line(color = "grey90", size = 0.4),  # background horizontal lines
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  )

ggsave("MSE_Boxplot.pdf", plot = p, width = 8, height = 5.5, dpi = 300)



#########################################################################
#        ---  Relative Error Ratio   ----
#########################################################################

relative_error_ratios <- sapply(1:8, function(i) {
  mse_i <- get(paste0("mse_ATE_", i))           # Posterior ATE MSE
  ATE_i <- get(paste0("true_ATE_", i))         # True ITE
  mean(mse_i) / (ATE_i^2 + 1e-8)            # Relative error ratio
})


result_table <- data.frame(
  Scenario = paste0("Scenario ", 1:8),
  Relative_Error_Ratio = round(relative_error_ratios, 4),
  Error_Level = ifelse(relative_error_ratios < 0.05, "Excellent",
                       ifelse(relative_error_ratios < 0.10, "Acceptable", "Needs Improvement"))
)

kable(result_table, caption = "Relative Error Ratios for ATE Estimation")

table_plot <- tableGrob(result_table, rows = NULL)
ggsave("relative_error_ratio.pdf", table_plot, width = 8, height = 4)


#########################################################################
#        ---  CDBMM and CART rand table   ----
#########################################################################

rand_means <- sapply(rand_CDBMM_all, mean)

rand_table <- data.frame(
  Scenario = names(rand_CDBMM_all),
  Rand_Index_Mean = round(rand_means, 4)
)

kable(rand_table, caption = "Adjusted Rand Index (ARI) across Scenarios")
table_plot <- tableGrob(rand_table, rows = NULL)
ggsave("rand_index_CDBMM.pdf", table_plot, width = 8, height = 4)