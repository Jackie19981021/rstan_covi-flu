# 安装并加载必要的R包
library(rstan)
library(rstantools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# 读取流感数据
seir_model_code <- readLines("/Users/chenjiaqi/Desktop/rstan/week_data.stan")
# 编译Stan模型
seir_model <- stan_model(model_code = seir_model_code)
flu_data <- read_csv("/Users/chenjiaqi/Downloads/RAT2/weekly_covid_cases.csv")

# 读取新冠数据
covid_data <- read_csv("/Users/chenjiaqi/Downloads/RAT2/covid_data.csv")

# 限制新冠数据到112天
dc <- covid_data$daily_cases[755:852] 
# 正确提取流感数据的12个周数据
fc <- flu_data$flu_cases1[131:144]
npi <- covid_data$NPI[755:852] 
data_list <- list(  T = length(dc),  # 新冠数据天数
                    W = length(fc),  # 流感数据周数
                  covid_cases = dc,
                  flu_cases = fc,
                  N = 7340000,
                  #phi_covid = 1,
                  #phi_flu = 1,
                  a=npi    )
# 使用优化找到初始参数估计
# 优化参数的初始值
opt <- optimizing(
   seir_model,
   data = data_list,
   hessian = TRUE
 )


# 准备MCMC采样的初始参数
# initial_values <- function() {
#   list(
#     beta1 = opt$par["beta1"],
#     beta2 = opt$par["beta2"],
#     sigma1 = opt$par["sigma1"],
#     sigma2 = opt$par["sigma2"],
#     gamma1 = opt$par["gamma1"],
#     gamma2 = opt$par["gamma2"],
#     rho1 = opt$par["rho1"],
#     rho2 = opt$par["rho2"],
#     S0 = opt$par["S0"],
#     E10 = opt$par["E10"],
#     I10 = opt$par["I10"],
#     R10 = opt$par["R10"],
#     E20 = opt$par["E20"],
#     I20 = opt$par["I20"],
#     R20 = opt$par["R20"],
#     phi_covid_inv= opt$par["phi_covid_inv"],
#     phi_flu_inv= opt$par["phi_flu_inv"]
#     # 其他需要的参数也要设置
#   )
# }

fit <- sampling(
  seir_model,
  data = data_list,
  chains = 4,
  iter = 2000,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  #init = initial_values
)

posterior_samples <- rstan::extract(fit)
print(names(posterior_samples))
print(fit, probs = c(0.025, 0.5, 0.975))
# # 绘制诊断图
# rstan::traceplot(fit)
rstan::traceplot(fit, pars = c("beta1", "beta2", "sigma1", "sigma2", "gamma1", "gamma2", "rho1", "rho2", "S0", "S120", "S210", "E10", "I10", "R10", "E20", "I20", "R20", "susc1", "susc2",  "adjustment_parameter1", "adjustment_parameter2"))

# 创建数据框
beta_data <- data.frame(
  Value = c(posterior_samples$beta1, posterior_samples$beta2),
  Disease = rep(c("COVID-19", "FLU"), each = length(posterior_samples$beta1))
)

# 绘制 beta 的分布图
ggplot(beta_data, aes(x = Value, fill = Disease)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("COVID-19" = "blue", "FLU" = "red"),
                    name = "Disease Type",
                    labels = c("Beta1 (COVID-19)", "Beta2 (FLU)")) +
  labs(title = "Distribution of Beta1 (COVID-19) and Beta2 (FLU)",
       x = "Value",
       y = "Density") +
  theme_minimal()

posterior_samples <- rstan::extract(fit)
# 计算COVID-19和流感的预测平均值
covid_pred_mean <- apply(posterior_samples$covid_cases_pred, 2, mean)
flu_pred_mean <- apply(posterior_samples$flu_cases_pred_weekly, 2, mean)

# 计算预测的95%置信区间
covid_lower_ci <- apply(posterior_samples$covid_cases_pred, 2, quantile, probs = 0.025)
covid_upper_ci <- apply(posterior_samples$covid_cases_pred, 2, quantile, probs = 0.975)
flu_lower_ci <- apply(posterior_samples$flu_cases_pred_weekly, 2, quantile, probs = 0.025)
flu_upper_ci <- apply(posterior_samples$flu_cases_pred_weekly, 2, quantile, probs = 0.975)

# 准备数据框用于绘图
covid_data_plot <- data.frame(
  Day = 1:length(dc),
  Actual = dc,
  Predicted = covid_pred_mean,
  Lower_CI = covid_lower_ci,
  Upper_CI = covid_upper_ci
)

flu_data_plot <- data.frame(
  Week = 1:length(fc),
  Actual = fc,
  Predicted = flu_pred_mean,
  Lower_CI = flu_lower_ci,
  Upper_CI = flu_upper_ci
)

# 绘制COVID-19和流感的后验拟合图
library(ggplot2)

# COVID-19图
p1 <- ggplot(covid_data_plot, aes(x = Day)) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = Predicted), color = "blue", linetype = "dashed") +
  geom_point(aes(y = Actual), color = "red") +
  labs(title = "Posterior Fit for COVID-19 Cases", x = "Day", y = "Number of Cases")

# 流感图
p2 <- ggplot(flu_data_plot, aes(x = Week)) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), fill = "green", alpha = 0.2) +
  geom_line(aes(y = Predicted), color = "green", linetype = "dashed") +
  geom_point(aes(y = Actual), color = "orange") +
  labs(title = "Posterior Fit for Weekly Flu Cases", x = "Week", y = "Number of Cases")

# 打印图形
print(p1)
print(p2)
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

# 提取 R_t 数据
# 提取 R_t 数据
rt_covid <- posterior_samples$R_t_covid
rt_flu <- posterior_samples$R_t_flu

# 计算每个时间点的 R_t 平均值和标准差
rt_covid_mean <- apply(rt_covid, 2, mean)
rt_flu_mean <- apply(rt_flu, 2, mean)
rt_covid_sd <- apply(rt_covid, 2, sd)
rt_flu_sd <- apply(rt_flu, 2, sd)

# 创建一个数据框，包括时间点和 R_t 的平均值和标准差
# 使用 rt_covid 或 rt_flu 的列数来确定时间点的数量
rt_data <- data.frame(
  Week = 1:ncol(rt_covid),  # 使用 ncol() 确定时间点数量
  Rt_Covid = rt_covid_mean,
  Rt_Flu = rt_flu_mean,
  Rt_Covid_SD = rt_covid_sd,
  Rt_Flu_SD = rt_flu_sd
)


# 绘制 R_t 与置信区间
p <- ggplot(rt_data, aes(x = Week)) +
  geom_line(aes(y = Rt_Covid, colour = "COVID-19"), size = 1.2) +
  geom_ribbon(aes(ymin = Rt_Covid - Rt_Covid_SD, ymax = Rt_Covid + Rt_Covid_SD, fill = "COVID-19"), alpha = 0.2) +
  geom_line(aes(y = Rt_Flu, colour = "Influenza"), size = 1.2) +
  geom_ribbon(aes(ymin = Rt_Flu - Rt_Flu_SD, ymax = Rt_Flu + Rt_Flu_SD, fill = "Influenza"), alpha = 0.2) +
  labs(title = "Effective Reproductive Number (R_t) Over Time",
       x = "Day",
       y = "Effective Reproductive Number (R_t)") +
  scale_colour_manual(values = c("COVID-19" = "red", "Influenza" = "blue")) +
  scale_fill_manual(values = c("COVID-19" = "red", "Influenza" = "blue"), guide = FALSE) +
  theme_minimal()

print(p)



# 提取phi_covid和phi_flu的样本
phi_covid_samples <- posterior_samples$phi_covid
phi_flu_samples <- posterior_samples$phi_flu

# 使用ggplot2绘制后验分布
df <- data.frame(Phi_Covid = phi_covid_samples, Phi_Flu = phi_flu_samples)
p1 <- ggplot(df, aes(x = Phi_Covid)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Posterior Distribution of phi_covid")

p2 <- ggplot(df, aes(x = Phi_Flu)) +
  geom_density(fill = "red", alpha = 0.5) +
  labs(title = "Posterior Distribution of phi_flu")

# 显示图形
print(p1)
print(p2)





