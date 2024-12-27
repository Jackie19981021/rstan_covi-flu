data {
  int<lower=1> T;       // 时间点数量
  int covid_cases[T];   // COVID-19病例数
  int<lower=1> W;                 // 周数
  int flu_cases[W];               // 每周流感病例数
  real N;               // 总人口数
  //real<lower=0> phi_covid;  // COVID-19负二项分布分散参数
  //real<lower=0> phi_flu;    // 流感负二项分布分散参数
  real a[T];
}

transformed data{
  //real sigma1=0.2;
  //real I10=0.2;
}
parameters {
  real<lower=0,upper=5> beta1;       // COVID-19传播率
  real<lower=0,upper=5> beta2;       // 流感传播率
  real<lower=0,upper=1> sigma1;      // COVID-19潜伏期率
  real<lower=0,upper=1> sigma2;      // 流感潜伏期率
  real<lower=0,upper=1> gamma1;      // COVID-19恢复率
  real<lower=0,upper=1> gamma2;      // 流感恢复率
  real<lower=0,upper=0.1> rho1;       // COVID-19患者对流感的暂时免疫率
  real<lower=0,upper=0.1> rho2;       // 流感患者对COVID-19的暂时免疫率
  real<lower=0,upper=1> eta1;       // COVID-19患者对流感的暂时免疫率
  real<lower=0,upper=1> eta2;
  real<lower=0,upper=7000000> S0;          // 初始易感人群
  real<lower=0,upper=7000000> S120;          // 2菌株易感人群
  real<lower=0,upper=7000000> S210;          // 1菌株易感人群
  real<lower=-100,upper=-5> E10;         // 初始COVID-19潜伏人群
  real<lower=-100,upper=-5> I10;         // 初始COVID-19感染人群
  real<lower=0> R10;         // 初始COVID-19恢复人群
  real<lower=-100,upper=-5> E20;         // 初始流感潜伏人群 
  real<lower=-100,upper=-5> I20;         // 初始流感感染人群
  real<lower=0> R20;         // 初始流感恢复人群
  real<lower=0,upper=1> susc1;
  real<lower=0,upper=1> susc2;
  real<lower=0,upper=1> adjustment_parameter1;
  real<lower=0,upper=1> adjustment_parameter2; // 调整参数
  real<lower=0,upper=1> alpha;//npi参数
  real<lower=0> phi_covid_inv;
  real<lower=0> phi_flu_inv;
  //real beta_base;      // 基本传播率
  //real<lower=0> beta_amp;  // 传播率的季节性幅度
  //real<lower=0> beta_freq; // 传播率的频率
  // 其他参数
  
  
  
  
}



transformed parameters {
  //real beta[T];  // 时间变化的beta

  //for (t in 1:T) {
    //beta[t] = beta_base + beta_amp * sin(beta_freq * t + 0);  // 假设季节性变化是正弦波形
  //}
  real S[T];        // 易感人群
  real S12[T];        // 易感人群
  real S21[T];        // 易感人群
  real E1[T];       // COVID-19潜伏人群
  real I1[T];       // COVID-19感染人群
  real R1[T];       // COVID-19恢复人群
  real E2[T];       // 流感潜伏人群
  real I2[T];       // 流感感染人群 
  real R2[T];       // 流感恢复人群
  //real fitcov_cases[T];
  //real fitflu_cases[T];
  real phi_covid = 1. / phi_covid_inv;
  real phi_flu = 1. / phi_flu_inv;

  // S[1] = S0;
  // E1[1] = E10;
  // I1[1] = I10; 
  // R1[1] = R10;
  // E2[1] = E20;
  // I2[1] = I20;
  // R2[1] = R20;
  
  
  
  E1[1] =exp(E10)*susc1*N;
  I1[1] =exp(I10)*susc1*N;
  S[1] =N-E10-E20-I10-I20-R10-R20-S120-S210;
  S12[1] =susc2*N-I20-E20;
  S21[1] =susc1*N-I10-E10;
  R1[1] =(1-susc1)*N;
  E2[1] =exp(E20)*susc2*N;
  I2[1] =exp(I20)*susc2*N;
  R2[1] =(1-susc2)*N;
  
  
  
  
  for (t in 2:T) {
    // real dS = -(beta1 * I1[t-1] + beta2 * I2[t-1]) * S[t-1] / N + 
    //           rho12 * R1[t-1] * I2[t-1] / N + rho21 * R2[t-1] * I1[t-1] / N;
    // real dE1 = beta1 * I1[t-1] * ((S[t-1] - rho12 * R1[t-1] * I2[t-1] / N) / N);
    // real dI1 = sigma1 * E1[t-1] - gamma1 * I1[t-1];
    // real dR1 = gamma1 * I1[t-1] - rho12 * R1[t-1] * I2[t-1] / N;
    // real dE2 = beta2 * I2[t-1] * ((S[t-1] - rho21 * R2[t-1] * I1[t-1] / N) / N);
    // real dI2 = sigma2 * E2[t-1] - gamma2 * I2[t-1];
    // real dR2 = gamma2 * I2[t-1] - rho21 * R2[t-1] * I1[t-1] / N;
    
    // real dS = -(beta1 * (1-a[t-1]) * I1[t-1] + beta2 * (1-a[t-1]) * I2[t-1]) * S[t-1] / N + 
    //        rho12* (1-a[t-1]) * R1[t-1] * I2[t-1] / N + rho21* (1-a[t-1]) * R2[t-1] * I1[t-1] / N;
    // real dE1 = beta1 * (1-a[t-1]) * I1[t-1] * ((S[t-1] - rho12* (1-a[t-1]) * R1[t-1] * I2[t-1] / N) / N);
    // real dE2 = beta2 * (1-a[t-1]) * I2[t-1] * ((S[t-1] - rho21* (1-a[t-1]) * R2[t-1] * I1[t-1] / N) / N);
    // real dI1 = sigma1 * E1[t-1] - gamma1 * I1[t-1];
    // real dR1 = gamma1 * I1[t-1] - rho12 * (1-a[t-1])* R1[t-1] * I2[t-1] / N;
    // real dI2 = sigma2 * E2[t-1] - gamma2 * I2[t-1];
    // real dR2 = gamma2 * I2[t-1] - rho21 * (1-a[t-1])* R2[t-1] * I1[t-1] / N;
    
    
     real dS = -(beta1* exp(-alpha*a[t-1]) * I1[t-1] + beta2* exp(-alpha*a[t-1]) * I2[t-1]) * S[t-1] / N;
     real dE1 = beta1* exp(-alpha*a[t-1])*S21[t-1]*I1[t-1]/N+beta1* exp(-alpha*a[t-1])*S[t-1]*I1[t-1]/N-sigma1*E1[t-1];
     //real dI1[t-1] * ((S[t-1] - rho12 * R1[t-1] * I2[t-1] / N) / N);
     real dI1 = sigma1 * E1[t-1] - gamma1 * I1[t-1];
     real dR1 = gamma1 * I1[t-1] - rho1 * R1[t-1];
     real dE2 = beta2* exp(-alpha*a[t-1])*S12[t-1]*I2[t-1]/N+beta2* exp(-alpha*a[t-1])*S[t-1]*I2[t-1]/N-sigma2*E2[t-1];
    // real dE2 = beta2 * I2[t-1] * ((S[t-1] - rho21 * R2[t-1] * I1[t-1] / N) / N);
     real dI2 = sigma2 * E2[t-1] - gamma2 * I2[t-1];
     real dR2 = gamma2 * I2[t-1] - rho2 * R2[t-1];
     real dS21 = -beta1* exp(-alpha*a[t-1]) * I1[t-1] * S21[t-1] / N+ rho2 * R2[t-1];
     real dS12 = -beta2* exp(-alpha*a[t-1]) * I2[t-1] * S12[t-1] / N+ rho1 * R1[t-1];
    
    S[t] = S[t-1] + dS;
    S12[t] = S12[t-1] + dS12;
    S21[t] = S21[t-1] + dS21;
    E1[t] = E1[t-1] + dE1;
    I1[t] = I1[t-1] + dI1;
    R1[t] = R1[t-1] + dR1;
    E2[t] = E2[t-1] + dE2; 
    I2[t] = I2[t-1] + dI2;
    R2[t] = R2[t-1] + dR2;

    //fitcov_cases[t]=(I1[t]-I1[t-1])*adjustment_parameter1;
    //fitflu_cases[t]=(I2[t]-I2[t-1])*adjustment_parameter2;
    

  }
}

model {
  // 先验分布
  beta1 ~ normal(1, 0.5);
  beta2 ~ normal(0.4, 1); 
  sigma1~ normal(0.5, 0.5); 
  gamma1~ normal(0.2, 0.5);
  gamma2~ normal(0.2, 0.5);
  rho2 ~ normal(0.01, 0.01);
  rho1 ~ normal(0.05, 0.05);
  phi_covid_inv ~ exponential(5);
  phi_flu_inv ~ exponential(5);
  //I2 ~ normal(5, 5);
  // 似然函数
  for (t in 1:T) {
    covid_cases[t] ~ neg_binomial_2(E1[t]* sigma1*adjustment_parameter1, phi_covid);
    //flu_cases[t] ~ neg_binomial_2(E2[t]*sigma2*adjustment_parameter2, phi_flu);
    //covid_cases[t] ~ neg_binomial_2(fitcov_cases[t], phi_covid);
    //flu_cases[t] ~ neg_binomial_2(fitflu_cases[t], phi_flu);
    }
    
  // 每周流感数据的处理
  int flu_index = 1;  // 周索引
  int flu_day_count = 0;  // 累计到每周的天数
  real weekly_flu_pred = 0;  // 每周预测病例的累计

  for (t in 1:T) {
    weekly_flu_pred += E2[t] * sigma2 * adjustment_parameter2;  // 累加预测病例

    if ((t % 7 == 0) || (t == T)) {  // 每周的最后一天或数据的最后一天
      flu_cases[flu_index] ~ neg_binomial_2(weekly_flu_pred, phi_flu);  // 周数据模型
      flu_index += 1;  // 移动到下一周
      weekly_flu_pred = 0;  // 重置周预测累计器
    }
  }    
    
    
}
generated quantities {
    real S_sum[T];
    real S_covid[T];
    real E_covid[T];
    real I_covid[T];
    real R_covid[T];
    real S_flu[T];
    real E_flu[T];
    real I_flu[T];
    real R_flu[T];
    int covid_cases_pred[T];  // 每日COVID-19病例预测
    int flu_cases_pred_weekly[W];  // 每周流感病例预测
    real weekly_flu_pred = 0;  // 用于累计每周的流感病例预测值
    int flu_index = 1;  // 流感周索引
    real R_t_covid[T];
    real R_t_flu[T];

    for (t in 1:T) {
        S_sum[t] = S[t];
        S_covid[t] = S21[t];
        E_covid[t] = E1[t];
        I_covid[t] = I1[t];
        R_covid[t] = R1[t];

        S_flu[t] = S12[t];
        E_flu[t] = E2[t];
        I_flu[t] = I2[t];
        R_flu[t] = R2[t];

        covid_cases_pred[t] = neg_binomial_2_rng(E1[t]* sigma1*adjustment_parameter1+1e-5, phi_covid);
        weekly_flu_pred += E2[t] * sigma2 * adjustment_parameter2;  // 累计每日流感病例预测

        if ((t % 7 == 0) || (t == T)) {  // 每周的最后一天或数据的最后一天
            flu_cases_pred_weekly[flu_index] = neg_binomial_2_rng(weekly_flu_pred + 1e-5, phi_flu);
            flu_index += 1;
            weekly_flu_pred = 0;  // 重置累计器
        }

        // Calculate effective reproduction number (R_t) for each disease
        //R_t_covid[t] = (beta1 * S_covid[t]) / (gamma1 * N);
        //R_t_flu[t] = (beta2 * S_flu[t]) / (gamma2 * N);
        //R_t_covid[t] = beta1*susc1/gamma1;
        //R_t_flu[t] = beta2*susc2/gamma2;
        R_t_covid[t] = (beta1* exp(-alpha*a[1]) /gamma1) *(S_covid[t] /N);
        R_t_flu[t] = (beta2*exp(-alpha*a[1])/gamma2) *(S_flu[t] /N);
    }
}
