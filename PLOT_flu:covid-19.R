# 安装并加载必要的R包
library(rstan)
library(rstantools)
library(dplyr)
library(tidyr)
library(ggplot2)
#library(gridExtra)
# 读取Stan模型代码
#seir_model_code <- readLines("/Users/chenjiaqi/Desktop/rstan/test_new.stan")
seir_model_code <- readLines("/Users/chenjiaqi/Desktop/rstan/test_new.stan")
# 编译Stan模型
seir_model <- stan_model(model_code = seir_model_code)

# Step 1: Read the CSV file
covid_data <- read.csv("/Users/chenjiaqi/Downloads/RAT2/covid_data.csv")

# Step 2: Extract the specific rows from the 'daily_cases' column
# Note: Adjust the row indices if your CSV file includes a header
dc <- covid_data$daily_cases[755:880] 
npi <- covid_data$NPI[755:880] 
fc <- covid_data$FLU1[755:880] 
#dc <- covid_data$daily_cases[755:880] 
# 准备数据
data <- data.frame(week = 1:126,
                   covid_cases = dc,
                   flu_cases = fc)

data_list <- list(T = nrow(data), 
                  covid_cases = data$covid_cases,
                  flu_cases = data$flu_cases,
                  N = 7340000,
                  #phi_covid = 1,
                  #phi_flu = 1,
                  a=npi    )

# 优化参数的初始值 
 opt <- optimizing(
   seir_model,
   data = data_list,
   hessian = TRUE
 )
 
 # 准备MCMC采样的初始参数
 initial_values <- function() {
   list(
     beta1 = opt$par["beta1"],
     beta2 = opt$par["beta2"],
     #sigma1 = opt$par["sigma1"],
     sigma2 = opt$par["sigma2"],
     #gamma1 = opt$par["gamma1"],
     #gamma2 = opt$par["gamma2"],
     rho1 = opt$par["rho1"],
     rho2 = opt$par["rho2"],
     S0 = opt$par["S0"],
     E10 = opt$par["E10"],
     I10 = opt$par["I10"],
     R10 = opt$par["R10"],
     E20 = opt$par["E20"],
     I20 = opt$par["I20"],
     R20 = opt$par["R20"],
     phi_covid_inv= opt$par["phi_covid_inv"],
     phi_flu_inv= opt$par["phi_flu_inv"]
     # 其他需要的参数也要设置
   )
 }
 

# 使用优化后的初始值运行MCMC采样
fit <- sampling(
  seir_model,
  data = data_list,
  chains = 4,
  iter = 2000,
  control = list(adapt_delta = 0.98, max_treedepth = 17),
  init = initial_values
  # init = list(
  #    list(
  #      beta1 = opt$par[1], 
  #      beta2 = opt$par[2], 
  #      sigma1 = opt$par[3], 
  #      sigma2 = opt$par[4], 
  #      gamma1 = opt$par[5], 
  #      gamma2 = opt$par[6], 
  #      rho12 = opt$par[7], 
  #      rho21 = opt$par[8], 
  #      S0 = opt$par[9], 
  #      E10 = opt$par[10], 
  #      I10 = opt$par[11], 
  #      R10 = opt$par[12], 
  #      E20 = opt$par[13], 
  #      I20 = opt$par[14], 
  #      R20 = opt$par[15]
  #    ),
  #    list(
  #      beta1 = opt$par[1], 
  #      beta2 = opt$par[2], 
  #      sigma1 = opt$par[3], 
  #      sigma2 = opt$par[4], 
  #      gamma1 = opt$par[5], 
  #      gamma2 = opt$par[6], 
  #      rho12 = opt$par[7], 
  #      rho21 = opt$par[8], 
  #      S0 = opt$par[9], 
  #      E10 = opt$par[10], 
  #      I10 = opt$par[11], 
  #      R10 = opt$par[12], 
  #      E20 = opt$par[13], 
  #      I20 = opt$par[14], 
  #      R20 = opt$par[15]
  #    ),
  #    list(
  #      beta1 = opt$par[1], 
  #      beta2 = opt$par[2], 
  #      sigma1 = opt$par[3], 
  #      sigma2 = opt$par[4], 
  #      gamma1 = opt$par[5], 
  #      gamma2 = opt$par[6], 
  #      rho12 = opt$par[7], 
  #      rho21 = opt$par[8], 
  #      S0 = opt$par[9], 
  #      E10 = opt$par[10], 
  #      I10 = opt$par[11], 
  #      R10 = opt$par[12], 
  #      E20 = opt$par[13], 
  #      I20 = opt$par[14], 
  #      R20 = opt$par[15]
  #    ),
  #    list(
  #      beta1 = opt$par[1], 
  #      beta2 = opt$par[2], 
  #      sigma1 = opt$par[3], 
  #      sigma2 = opt$par[4], 
  #      gamma1 = opt$par[5], 
  #      gamma2 = opt$par[6], 
  #      rho12 = opt$par[7], 
  #      rho21 = opt$par[8], 
  #      S0 = opt$par[9], 
  #      E10 = opt$par[10], 
  #      I10 = opt$par[11], 
  #      R10 = opt$par[12], 
  #      E20 = opt$par[13], 
  #      I20 = opt$par[14], 
  #      R20 = opt$par[15]
  #  )
  # )
)


posterior_samples <- rstan::extract(fit)
print(names(posterior_samples))
print(fit, probs = c(0.025, 0.5, 0.975))

# # 绘制诊断图
# rstan::traceplot(fit)
rstan::traceplot(fit, pars = c("beta1", "beta2",  "sigma2", "gamma1", "gamma2", "rho1", "rho2", "eta1", "eta2","S0", "S120", "S210", "E10", "I10", "R10", "E20", "I20", "R20", "susc1", "susc2",  "adjustment_parameter1", "adjustment_parameter2"))
#rstan::traceplot(fit, pars = c("beta1", "beta2",  "sigma2",  "rho1", "rho2", "eta1", "eta2","S0", "S120", "S210", "E10", "I10", "R10", "E20", "I20", "R20", "susc1", "susc2",  "adjustment_parameter1", "adjustment_parameter2"))

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
  labs(title = "Fit Cases: COVID-19 vs. Influenza", y = "Cases", x = "Day") +
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

# 提取 R_t 数据
rt_covid <- posterior_samples$R_t_covid
rt_flu <- posterior_samples$R_t_flu
library(ggplot2)

# 计算每个时间点的 R_t 平均值和标准差
rt_covid_mean <- apply(rt_covid, 2, mean)
rt_flu_mean <- apply(rt_flu, 2, mean)
# 计算每个时间点的 R_t 的 2.5% 和 97.5% 分位数作为置信区间
rt_covid_lower <- apply(rt_covid, 2, quantile, probs = 0.025)
rt_covid_upper <- apply(rt_covid, 2, quantile, probs = 0.975)
rt_flu_lower <- apply(rt_flu, 2, quantile, probs = 0.025)
rt_flu_upper <- apply(rt_flu, 2, quantile, probs = 0.975)

# 创建一个数据框，包括时间点和 R_t 的平均值以及置信区间
rt_data <- data.frame(
  Week = 1:ncol(rt_covid),  # 确保 Week 与数据列数匹配
  Rt_Covid = rt_covid_mean,
  Rt_Flu = rt_flu_mean,
  Rt_Covid_Lower = rt_covid_lower,
  Rt_Covid_Upper = rt_covid_upper,
  Rt_Flu_Lower = rt_flu_lower,
  Rt_Flu_Upper = rt_flu_upper
)
# 绘制 R_t 与置信区间
p <- ggplot(rt_data, aes(x = Week)) +
  geom_ribbon(aes(ymin = Rt_Covid_Lower, ymax = Rt_Covid_Upper, fill = "COVID-19"), alpha = 0.2) +
  geom_line(aes(y = Rt_Covid, colour = "COVID-19"), size = 1.2) +
  geom_ribbon(aes(ymin = Rt_Flu_Lower, ymax = Rt_Flu_Upper, fill = "Influenza"), alpha = 0.2) +
  geom_line(aes(y = Rt_Flu, colour = "Influenza"), size = 1.2) +
  labs(title = "Effective Reproductive Number (R_t) Over Time with Confidence Intervals",
       x = "Week",
       y = "Effective Reproductive Number (R_t)") +
  scale_colour_manual(values = c("COVID-19" = "red", "Influenza" = "blue")) +
  scale_fill_manual(values = c("COVID-19" = "red", "Influenza" = "blue"), guide = FALSE) +
  theme_minimal()

print(p)

# 假设你已经有了effect_data数据框
# 添加周信息


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


#library(rstan)

# 假设你已经有了从 Stan 模型运行后的拟合对象 `fit`
posterior_samples <- rstan::extract(fit)

# 提取每天对 COVID-19 和流感 R_t 的 NPI 影响
npi_effect_covid <- posterior_samples$npi_effect_on_rt_covid
npi_effect_flu <- posterior_samples$npi_effect_on_rt_flu

# 计算每天的平均值和置信区间
npi_effect_covid_mean <- apply(npi_effect_covid, 2, mean)
npi_effect_flu_mean <- apply(npi_effect_flu, 2, mean)
npi_effect_covid_ci <- apply(npi_effect_covid, 2, quantile, probs = c(0.025, 0.975))
npi_effect_flu_ci <- apply(npi_effect_flu, 2, quantile, probs = c(0.025, 0.975))

# 创建数据框架
effect_data <- data.frame(
  Day = 1:ncol(npi_effect_covid),
  Effect_Covid = npi_effect_covid_mean,
  Lower_CI_Covid = npi_effect_covid_ci[1,],
  Upper_CI_Covid = npi_effect_covid_ci[2,],
  Effect_Flu = npi_effect_flu_mean,
  Lower_CI_Flu = npi_effect_flu_ci[1,],
  Upper_CI_Flu = npi_effect_flu_ci[2,]
)


# 绘制 NPI 对 R_t 的影响
p <- ggplot(effect_data, aes(x = Day)) +
  geom_ribbon(aes(ymin = Lower_CI_Covid, ymax = Upper_CI_Covid, fill = "COVID-19"), alpha = 0.3) +
  geom_line(aes(y = Effect_Covid, colour = "COVID-19"), size = 1) +
  geom_ribbon(aes(ymin = Lower_CI_Flu, ymax = Upper_CI_Flu, fill = "Influenza"), alpha = 0.3) +
  geom_line(aes(y = Effect_Flu, colour = "Influenza"), size = 1) +
  labs(title = "Daily NPI Impact on R_t for COVID-19 and Influenza",
       x = "Day",
       y = "NPI Effect on R_t") +
  scale_colour_manual(values = c("COVID-19" = "red", "Influenza" = "blue")) +
  scale_fill_manual(values = c("COVID-19" = "red", "Influenza" = "lightblue")) +
  theme_minimal()

# 显示图形
print(p)

effect_data$Week <- ceiling(effect_data$Day / 7)


library(gridExtra)

# 计算每周的平均值和置信区间
weekly_data <- effect_data %>%
  mutate(Week = ceiling(Day / 7)) %>%
  group_by(Week) %>%
  summarise(
    Effect_Covid_Mean = mean(Effect_Covid),
    Lower_CI_Covid = mean(Lower_CI_Covid),
    Upper_CI_Covid = mean(Upper_CI_Covid),
    Effect_Flu_Mean = mean(Effect_Flu),
    Lower_CI_Flu = mean(Lower_CI_Flu),
    Upper_CI_Flu = mean(Upper_CI_Flu)
  )

# 绘制每周的箱线图
p1 <- ggplot(weekly_data, aes(x = as.factor(Week))) +
  geom_boxplot(aes(ymin = Lower_CI_Covid, lower = Lower_CI_Covid, middle = Effect_Covid_Mean, upper = Upper_CI_Covid, ymax = Upper_CI_Covid, fill = "COVID-19"), stat = "identity", alpha = 0.3) +
  labs(title = "Weekly NPI Impact on Rt for COVID-19",
       x = "Week",
       y = "Reduction in Rt (%)") +
  scale_fill_manual(values = c("COVID-19" = "red")) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent)

p2 <- ggplot(weekly_data, aes(x = as.factor(Week))) +
  geom_boxplot(aes(ymin = Lower_CI_Flu, lower = Lower_CI_Flu, middle = Effect_Flu_Mean, upper = Upper_CI_Flu, ymax = Upper_CI_Flu, fill = "Influenza"), stat = "identity", alpha = 0.3) +
  labs(title = "Weekly NPI Impact on R_t for Influenza",
       x = "Week",
       y = "Reduction in Rt (%)") +
  scale_fill_manual(values = c("Influenza" = "blue")) +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent)

# 创建两张子图
grid.arrange(p1, p2, ncol = 1)


# 假设已经有了后验分布的样本
# 提取 alpha1 和 alpha2 的样本
alpha1_samples <- posterior_samples$alpha1
alpha2_samples <- posterior_samples$alpha2

# 计算 NPI 对 R_t 的影响百分比
alpha1_effect <- 100 * (1 - exp(-alpha1_samples))
alpha2_effect <- 100 * (1 - exp(-alpha2_samples))

# 将数据整合到一个数据框中
effects_df <- data.frame(
  Effect = c(alpha1_effect, alpha2_effect),
  Virus = factor(rep(c("COVID-19", "Influenza"), each = length(alpha1_effect)))
)

# 使用 ggplot2 绘制箱线图
library(ggplot2)
p <- ggplot(effects_df, aes(x = Virus, y = Effect, fill = Virus)) +
  geom_boxplot() +
  labs(title = "NPI Effect on R_t for COVID-19 and Influenza",
       x = "",
       y = "Percent Reduction in R_t") +
  scale_fill_manual(values = c("COVID-19" = "red", "Influenza" = "blue")) +
  theme_minimal()

print(p)



