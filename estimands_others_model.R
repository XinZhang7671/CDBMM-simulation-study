load("BART_results.RData")

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


# All samples of BART ITE (500X100)
BART_tau_1_matrix=sapply(BART_results[["BART_scenario_1"]], function(res) res[["tau"]])
BART_tau_2_matrix=sapply(BART_results[["BART_scenario_2"]], function(res) res[["tau"]])
BART_tau_3_matrix=sapply(BART_results[["BART_scenario_3"]], function(res) res[["tau"]])
BART_tau_4_matrix=sapply(BART_results[["BART_scenario_4"]], function(res) res[["tau"]])
BART_tau_5_matrix=sapply(BART_results[["BART_scenario_5"]], function(res) res[["tau"]])
BART_tau_6_matrix=sapply(BART_results[["BART_scenario_6"]], function(res) res[["tau"]])
BART_tau_7_matrix=sapply(BART_results[["BART_scenario_7"]], function(res) res[["tau"]])
BART_tau_8_matrix=sapply(BART_results[["BART_scenario_8"]], function(res) res[["tau"]])

# BART ATE (1:100)
BART_ATE_1 <- colMeans(BART_tau_1_matrix)
BART_ATE_2 <- colMeans(BART_tau_2_matrix)
BART_ATE_3 <- colMeans(BART_tau_3_matrix)
BART_ATE_4 <- colMeans(BART_tau_4_matrix)
BART_ATE_5 <- colMeans(BART_tau_5_matrix)
BART_ATE_6 <- colMeans(BART_tau_6_matrix)
BART_ATE_7 <- colMeans(BART_tau_7_matrix)
BART_ATE_8 <- colMeans(BART_tau_8_matrix)

# bias_ATE 
BART_bias_ATE_1=BART_ATE_1 - true_ATE_1
BART_bias_ATE_2=BART_ATE_2 - true_ATE_2
BART_bias_ATE_3=BART_ATE_3 - true_ATE_3
BART_bias_ATE_4=BART_ATE_4 - true_ATE_4
BART_bias_ATE_5=BART_ATE_5 - true_ATE_5
BART_bias_ATE_6=BART_ATE_6 - true_ATE_6
BART_bias_ATE_7=BART_ATE_7 - true_ATE_7
BART_bias_ATE_8=BART_ATE_8 - true_ATE_8

# mse_ATE 
BART_mse_ATE_1=(BART_ATE_1 - true_ATE_1)^2
BART_mse_ATE_2=(BART_ATE_2 - true_ATE_2)^2
BART_mse_ATE_3=(BART_ATE_3 - true_ATE_3)^2
BART_mse_ATE_4=(BART_ATE_4 - true_ATE_4)^2
BART_mse_ATE_5=(BART_ATE_5 - true_ATE_5)^2
BART_mse_ATE_6=(BART_ATE_6 - true_ATE_6)^2
BART_mse_ATE_7=(BART_ATE_7 - true_ATE_7)^2
BART_mse_ATE_8=(BART_ATE_8 - true_ATE_8)^2


mean_mse_BART <- sapply(1:8, function(i) mean(get(paste0("BART_mse_ATE_", i))))
mean_mse_BART

#########################################################################
#        ---  comparison of Relative Error Ratio   ----
#########################################################################

BART_relative_error_ratios <- sapply(1:8, function(i) {
  BART_mse_i <- get(paste0("BART_mse_ATE_", i))           # Posterior ATE MSE
  ATE_i <- get(paste0("true_ATE_", i))         # True ITE
  mean(BART_mse_i) / (ATE_i^2 + 1e-8)            # Relative error ratio
})


BART_result_table <- data.frame(
  Scenario = paste0("Scenario ", 1:8),
  Relative_Error_Ratio = round(BART_relative_error_ratios, 4),
  Error_Level = ifelse(BART_relative_error_ratios < 0.05, "Excellent",
                       ifelse(BART_relative_error_ratios < 0.10, "Acceptable", "Needs Improvement"))
)

kable(BART_result_table, caption = "Relative Error Ratios for ATE Estimation")

table_plot <- tableGrob(BART_result_table, rows = NULL)
ggsave("BART_relative_error_ratio.pdf", table_plot, width = 8, height = 4)


# comparison
# 合并为一个对比数据框
compare_table <- data.frame(
  Scenario = paste0("Scenario ", 1:8),
  CDBMM_Relative_Error = round(relative_error_ratios, 4),
  BART_Relative_Error = round(BART_relative_error_ratios, 4)
)

# 添加哪一个更好的判断列
compare_table$Better_Method <- ifelse(
  compare_table$CDBMM_Relative_Error < compare_table$BART_Relative_Error,
  "CDBMM", "BART"
)

kable(compare_table, caption = "Comparison of Relative Error Ratios Between CDBMM and BART")
table_plot <- tableGrob(compare_table, rows = NULL)
ggsave("Relative_Error_Comparison.pdf", table_plot, width = 8, height = 4)


#########################################################################
#        ---  comparison of boxplot  ----
#########################################################################

# 合并 CDBMM 的 mse
mse_CDBMM <- list(mse_ATE_1, mse_ATE_2, mse_ATE_3, mse_ATE_4,
                  mse_ATE_5, mse_ATE_6, mse_ATE_7, mse_ATE_8)
mse_CDBMM_df <- data.frame(
  MSE = unlist(mse_CDBMM),
  Scenario = rep(paste0("Scenario_", 1:8), each = length(mse_ATE_1)),
  Method = "CDBMM"
)

# 合并 BART 的 mse
mse_BART <- list(BART_mse_ATE_1, BART_mse_ATE_2, BART_mse_ATE_3, BART_mse_ATE_4,
                 BART_mse_ATE_5, BART_mse_ATE_6, BART_mse_ATE_7, BART_mse_ATE_8)
mse_BART_df <- data.frame(
  MSE = unlist(mse_BART),
  Scenario = rep(paste0("Scenario_", 1:8), each = length(BART_mse_ATE_1)),
  Method = "BART"
)

# 合并为一个长数据框
mse_data <- rbind(mse_CDBMM_df, mse_BART_df)


ggplot(mse_data, aes(x = Scenario, y = MSE, fill = Method)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, position = position_dodge(0.8)) +
  scale_fill_manual(values = c("CDBMM" = "#1f77b4", "BART" = "#ff7f0e")) +  # 可调整颜色
  theme_minimal(base_size = 14) +
  labs(y = "MSE of ATE",
       x = "Scenario") +
  theme_classic(base_size = 13) +
  theme(
    panel.grid.major.y = element_line(color = "grey90", size = 0.4),  # background horizontal lines
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right"
  )




#########################################################################
#        ---  comparison of ARI  ----
#########################################################################

# RAND index for clusters obtained with CART
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

rand_CART_all <- lapply(1:8, function(i) {
  true_cluster <- true_clusters[[i]]
  CART_result <- CART_results[[i]]  # 注意不加 $CART_scenario
  
  sapply(1:samples, function(s) {
    pred_cluster <- CART_result[[s]]$partition  # 是向量
    mcclust::arandi(true_cluster, pred_cluster)
  })
})

names(rand_CART_all) <- paste0("scenario_", 1:8)

rand_CART_means <- sapply(rand_CART_all, mean)
rand_CDBMM_means <- sapply(rand_CDBMM_all, mean)

# 构建对比表格
rand_compare_table <- data.frame(
  Scenario = paste0("Scenario ", 1:8),
  CDBMM_Rand = round(rand_CDBMM_means, 4),
  CART_Rand = round(rand_CART_means, 4),
  Better_Method = ifelse(rand_CDBMM_means >= rand_CART_means, "CDBMM", "CART")
)

# 以 R 控制台中展示（可选）
kable(rand_compare_table, caption = "Adjusted Rand Index Comparison: CDBMM vs CART")

# 生成 PDF 图像形式的表格
table_plot <- tableGrob(rand_compare_table, rows = NULL)
ggsave("rand_index_comparison.pdf", table_plot, width = 8, height = 4)