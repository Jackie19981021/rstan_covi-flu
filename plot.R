
# 安装并加载必要的R包
library(rstan)
library(rstantools)
library(dplyr)
library(tidyr)
library(ggplot2)
# 读取Stan模型代码
seir_model_code <- readLines("/Users/chenjiaqi/Desktop/rstan/test_new.stan")

# 编译Stan模型
seir_model <- stan_model(model_code = seir_model_code)

# Step 1: Read the CSV file
covid_data <- read.csv("/Users/chenjiaqi/Downloads/RAT2/covid_data.csv")

# Step 2: Extract the specific rows from the 'daily_cases' column
# Note: Adjust the row indices if your CSV file includes a header
dc <- covid_data$daily_cases[1028:1118] 
npi <- covid_data$NPI[1028:1118] 
fc <- covid_data$FLU[1028:1118] 
# 准备数据
data <- data.frame(week = 91,
                   covid_cases = dc,
                   flu_cases = fc)

data_list <- list(T = nrow(data), 
                  covid_cases = data$covid_cases,
                  flu_cases = data$flu_cases,
                  N = 7340000,
                  phi_covid = 1,
                  phi_flu = 1,
                  a=npi    )

# 优化参数的初始值
# opt <- optimizing(
#   seir_model,
#   data = data_list,
#   hessian = TRUE
# )

# 使用优化后的初始值运行MCMC采样
fit <- sampling(
  seir_model,
  data = data_list,
  chains = 4,
  iter = 2000,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  # init = list(
  #   list(
  #     beta1 = opt$par[1], 
  #     beta2 = opt$par[2], 
  #     sigma1 = opt$par[3], 
  #     sigma2 = opt$par[4], 
  #     gamma1 = opt$par[5], 
  #     gamma2 = opt$par[6], 
  #     rho12 = opt$par[7], 
  #     rho21 = opt$par[8], 
  #     S0 = opt$par[9], 
  #     E10 = opt$par[10], 
  #     I10 = opt$par[11], 
  #     R10 = opt$par[12], 
  #     E20 = opt$par[13], 
  #     I20 = opt$par[14], 
  #     R20 = opt$par[15]
  #   ),
  #   list(
  #     beta1 = opt$par[1], 
  #     beta2 = opt$par[2], 
  #     sigma1 = opt$par[3], 
  #     sigma2 = opt$par[4], 
  #     gamma1 = opt$par[5], 
  #     gamma2 = opt$par[6], 
  #     rho12 = opt$par[7], 
  #     rho21 = opt$par[8], 
  #     S0 = opt$par[9], 
  #     E10 = opt$par[10], 
  #     I10 = opt$par[11], 
  #     R10 = opt$par[12], 
  #     E20 = opt$par[13], 
  #     I20 = opt$par[14], 
  #     R20 = opt$par[15]
  #   ),
  #   list(
  #     beta1 = opt$par[1], 
  #     beta2 = opt$par[2], 
  #     sigma1 = opt$par[3], 
  #     sigma2 = opt$par[4], 
  #     gamma1 = opt$par[5], 
  #     gamma2 = opt$par[6], 
  #     rho12 = opt$par[7], 
  #     rho21 = opt$par[8], 
  #     S0 = opt$par[9], 
  #     E10 = opt$par[10], 
  #     I10 = opt$par[11], 
  #     R10 = opt$par[12], 
  #     E20 = opt$par[13], 
  #     I20 = opt$par[14], 
  #     R20 = opt$par[15]
  #   ),
  #   list(
  #     beta1 = opt$par[1], 
  #     beta2 = opt$par[2], 
  #     sigma1 = opt$par[3], 
  #     sigma2 = opt$par[4], 
  #     gamma1 = opt$par[5], 
  #     gamma2 = opt$par[6], 
  #     rho12 = opt$par[7], 
  #     rho21 = opt$par[8], 
  #     S0 = opt$par[9], 
  #     E10 = opt$par[10], 
  #     I10 = opt$par[11], 
  #     R10 = opt$par[12], 
  #     E20 = opt$par[13], 
  #     I20 = opt$par[14], 
  #     R20 = opt$par[15]
   # )
  #)
)


posterior_samples <- rstan::extract(fit)
print(names(posterior_samples))
print(fit, probs = c(0.025, 0.5, 0.975))

# 绘制诊断图
rstan::traceplot(fit)

# # 创建数据框
beta_data <- data.frame(
  Beta1 = posterior_samples$beta1,
  Beta2 = posterior_samples$beta2
)

# 绘制 beta 的分布图
ggplot(beta_data, aes(x = Beta1)) +
  geom_density(fill = "blue", alpha = 0.5) +
  geom_density(aes(x = Beta2), fill = "red", alpha = 0.5) +
  labs(title = "Distribution of Beta1(COVID-19) and Beta2(FLU)",
       x = "Value",
       y = "Density") +
  theme_minimal()




# 计算COVID-19和流感的预测平均值
covid_pred_mean <- apply(posterior_samples$covid_cases_pred, 2, mean)
flu_pred_mean <- apply(posterior_samples$flu_cases_pred, 2, mean)

library(ggplot2)

# 准备数据用于绘图
plot_data <- data.frame(Week = data$week, Actual_Covid = data$covid_cases, Predicted_Covid = covid_pred_mean)


# 计算COVID-19预测值的95%置信区间，忽略NA和NaN值
covid_lower_ci <- apply(posterior_samples$covid_cases_pred, 2, quantile, probs = 0.025,na.rm = TRUE)


covid_upper_ci <- apply(posterior_samples$covid_cases_pred, 2, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

# 计算流感预测值的95%置信区间
flu_lower_ci <- apply(posterior_samples$flu_cases_pred, 2, quantile, probs = 0.025, na.rm = TRUE)
flu_upper_ci <- apply(posterior_samples$flu_cases_pred, 2, quantile, probs = 0.975, na.rm = TRUE)









library(ggplot2)

# 准备绘图数据，合并COVID和流感数据
plot_data <- data.frame(
  Week = rep(data$week, times = 2),
  Actual = c(data$covid_cases, data$flu_cases),
  Predicted = c(covid_pred_mean, flu_pred_mean),
  Lower_CI = c(covid_lower_ci, flu_lower_ci),
  Upper_CI = c(covid_upper_ci, flu_upper_ci),
  Disease = rep(c("COVID-19", "Influenza"), each = nrow(data))
)

# 绘制实际病例、预测病例和置信区间
ggplot(plot_data, aes(x = Week)) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI, fill = Disease), alpha = 0.2) +
  geom_line(aes(y = Actual, colour = "Actual"), size = 1.2) +
  geom_line(aes(y = Predicted, colour = "Predicted"), size = 1.2, linetype = "dashed") +
  scale_fill_manual(values = c("COVID-19" = "orange", "Influenza" = "lightblue")) +
  scale_colour_manual("", values = c("Actual" = "black", "Predicted" = "red")) +
  facet_wrap(~Disease, scales = "free_y") +  # 使用facet_wrap分面显示不同疾病
  labs(title = "Case Comparison: COVID-19 vs. Influenza", y = "Cases", x = "Week") +
  theme_minimal() +
  theme(legend.position = "bottom")

# 假设你的数据有 T 个时间点
# 假设 T 是时间步长，用于生成周数据
# weeks <- 1:T
# 
# # 创建数据框
# seir_data <- data.frame(
#   Week = rep(weeks, times = 8),
#   Value = c(posterior_samples$S_covid, posterior_samples$E_covid, posterior_samples$I_covid, posterior_samples$R_covid, 
#             posterior_samples$S_flu, posterior_samples$E_flu, posterior_samples$I_flu, posterior_samples$R_flu),
#   State = rep(c("Susceptible", "Exposed", "Infected", "Recovered"), each = length(weeks)),
#   Disease = rep(c("COVID-19", "Influenza"), each = 4 * length(weeks))
# )
# 
# 
# library(ggplot2)
# 
# # 绘制SEIR模型图像
# ggplot(seir_data, aes(x = Week, y = Value, color = State, group = State)) +
#   geom_line() +
#   facet_wrap(~Disease, scales = "free_y") +
#   labs(title = "SEIR Model for COVID-19 and Influenza",
#        x = "Week",
#        y = "Number of Individuals") +
#   scale_color_manual(values = c("Susceptible" = "blue", "Exposed" = "orange", "Infected" = "red", "Recovered" = "green")) +
#   theme_minimal() +
#   theme(legend.position = "right")





library(ggplot2)
#weeks <- 1:length(Rt_covid)  # 假设 Rt_covid 和 Rt_flu 长度相同

# # 准备绘图数据
# rt_data <- data.frame(
#   Week = rep(weeks, times = 2),
#   Rt = c(Rt_covid, Rt_flu),
#   Disease = rep(c("COVID-19", "Influenza"), each = length(Rt_covid))
# )
# 
# 
# # 绘制R_t演化曲线
# ggplot(rt_data, aes(x = Week, y = Rt, color = Disease)) +
#   geom_line() +
#   geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # 添加R_t=1的参考线
#   labs(title = "Evolution of Rt for COVID-19 and Influenza",
#        x = "Week",
#        y = "Rt (Effective Reproduction Number)") +
#   scale_color_manual(values = c("COVID-19" = "blue", "Influenza" = "green")) +
#   theme_minimal() +
#   theme(legend.position = "bottom")
# 打印出后验样本中的所有变量名


# 提取各个部分的后验样本
S_covid <- posterior_samples$S_covid
E_covid <- posterior_samples$E_covid
I_covid <- posterior_samples$I_covid
R_covid <- posterior_samples$R_covid

# 时间点
T <- ncol(S_covid)  # 使用列数来表示时间点数
time_points <- 1:T
# 计算均值
S_mean <- colMeans(S_covid)
E_mean <- colMeans(E_covid)
I_mean <- colMeans(I_covid)
R_mean <- colMeans(R_covid)

# 计算95%置信区间
S_ci <- apply(S_covid, 2, quantile, probs = c(0.025, 0.975))
E_ci <- apply(E_covid, 2, quantile, probs = c(0.025, 0.975))
I_ci <- apply(I_covid, 2, quantile, probs = c(0.025, 0.975))
R_ci <- apply(R_covid, 2, quantile, probs = c(0.025, 0.975))
# 创建数据框
covid_data <- data.frame(
  Day = rep(time_points, times = 4),
  Population = c(S_mean, E_mean, I_mean, R_mean),
  Type = rep(c("Susceptible", "Exposed", "Infected", "Recovered"), each = T),
  Lower = c(S_ci[1,], E_ci[1,], I_ci[1,], R_ci[1,]),
  Upper = c(S_ci[2,], E_ci[2,], I_ci[2,], R_ci[2,])
)

# 绘制
ggplot(covid_data, aes(x = Day, y = Population, color = Type)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Type), alpha = 0.2) +
  geom_line() +
  labs(title = "COVID-19 SEIR Model with Confidence Intervals", x = "Day", y = "Population") +
  theme_minimal() +
  scale_color_manual(values = c("blue", "green", "red", "purple")) +
  scale_fill_manual(values = c("blue", "green", "red", "purple")) +
  theme(legend.position = "bottom")


library(ggplot2)
library(dplyr)

# 提取 R_t 数据
rt_covid <- posterior_samples$R_t_covid
rt_flu <- posterior_samples$R_t_flu

# 计算每个时间点的 R_t 平均值
rt_covid_mean <- apply(rt_covid, 2, mean)
rt_flu_mean <- apply(rt_flu, 2, mean)

# 创建一个数据框，包括时间点和 R_t 的平均值
rt_data <- data.frame(
  Week = 1:nrow(data),
  Rt_Covid = rt_covid_mean,
  Rt_Flu = rt_flu_mean
)

# 绘制 R_t
p <- ggplot(rt_data, aes(x = Week)) +
  geom_line(aes(y = Rt_Covid, colour = "COVID-19"), size = 1.2) +
  geom_line(aes(y = Rt_Flu, colour = "Influenza"), size = 1.2) +
  labs(title = "Effective Reproductive Number (R_t) Over Time",
       x = "Week",
       y = "Effective Reproductive Number (R_t)") +
  scale_colour_manual(values = c("COVID-19" = "red", "Influenza" = "blue")) +
  theme_minimal()

print(p)






