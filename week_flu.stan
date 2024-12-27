data {
  int<lower=1> T;             // 时间点数量
  int covid_cases[T];         // COVID-19病例数
  int flu_cases[T];           // 流感病例数（现在可能不再直接使用）
  real N;                    // 总人口数
  real<lower=0> phi_covid;   // COVID-19负二项分布分散参数
  real<lower=0> phi_flu;     // 流感负二项分布分散参数
  real a[T];
  int<lower=0> flu_week_cases[(T + 6) / 7];  // 每七天累计的流感病例数
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
  //real beta_base;      // 基本传播率
  //real<lower=0> beta_amp;  // 传播率的季节性幅度
  //real<lower=0> beta_freq; // 传播率的频率
  // 其他参数
  
  
  
  
}



transformed parameters {
  real S[T];        // 易感人群
  real S12[T];      // 易感人群
  real S21[T];      // 易感人群
  real E1[T];       // COVID-19潜伏人群
  real I1[T];       // COVID-19感染人群
  real R1[T];       // COVID-19恢复人群
  real E2[T];       // 流感潜伏人群
  real I2[T];       // 流感感染人群 
  real R2[T];       // 流感恢复人群
  real flu_cases_weekly_sum[(T + 6) / 7];
  real week_index;

  // 初始化周累计病例数数组
  for (i in 1:((T + 6) / 7)) {
    flu_cases_weekly_sum[i] = 0;
  }

  // 首天的初始状态
  E1[1] = exp(E10) * susc1 * N;
  I1[1] = exp(I10) * susc1 * N;
  S[1] = N - E10 - E20 - I10 - I20 - R10 - R20 - S120 - S210;
  S12[1] = susc2 * N - I20 - E20;
  S21[1] = susc1 * N - I10 - E10;
  R1[1] = (1 - susc1) * N;
  E2[1] = exp(E20) * susc2 * N;
  I2[1] = exp(I20) * susc2 * N;
  R2[1] = (1 - susc2) * N;

  for (t in 2:T) {
    real dS = -(beta1 * (1 - a[t-1]) * I1[t-1] + beta2 * (1 - a[t-1]) * I2[t-1]) * S[t-1] / N;
    real dE1 = beta1 * (1 - eta1) * (1 - a[t-1]) * S21[t-1] * I1[t-1] / N + beta1 * (1 - a[t-1]) * S[t-1] * I1[t-1] / N - sigma1 * E1[t-1];
    real dI1 = sigma1 * E1[t-1] - gamma1 * I1[t-1];
    real dR1 = gamma1 * I1[t-1] - rho1 * R1[t-1];
    real dE2 = beta2 * (1 - eta2) * (1 - a[t-1]) * S12[t-1] * I2[t-1] / N + beta2 * (1 - a[t-1]) * S[t-1] * I2[t-1] / N - sigma2 * E2[t-1];
    real dI2 = sigma2 * E2[t-1] - gamma2 * I2[t-1];
    real dR2 = gamma2 * I2[t-1] - rho2 * R2[t-1];
    real dS21 = -(1 - eta1) * beta1 * (1 - a[t-1]) * I1[t-1] * S21[t-1] / N + rho2 * R2[t-1];
    real dS12 = -(1 - eta2) * beta2 * (1 - a[t-1]) * I2[t-1] * S12[t-1] / N + rho1 * R1[t-1];

    S[t] = S[t-1] + dS;
    S12[t] = S12[t-1] + dS12;
    S21[t] = S21[t-1] + dS21;
    E1[t] = E1[t-1] + dE1;
    I1[t] = I1[t-1] + dI1;
    R1[t] = R1[t-1] + dR1;
    E2[t] = E2[t-1] + dE2;
    I2[t] = I2[t-1] + dI2;
    R2[t] = R2[t-1] + dR2;

    // 累积每七天的流感病例数
    week_index = (t + 6) / 7;
    flu_cases_weekly_sum[int_step(week_index)] += sigma2*E2[t-1]; // 使用感染人数来累计
  }
}


model {
  // 先验分布
  beta1 ~ normal(2, 1);
  beta2 ~ normal(0.4, 1); 
  sigma1~ normal(0.5, 0.5); 
  gamma1~ normal(0.2, 0.5);
  gamma2~ normal(0.2, 0.5);
  //rho2 ~ normal(0.01, 0.01);
  //rho1 ~ normal(0.05, 0.05);
  //I2 ~ normal(5, 5);
  // 似然函数
  for (i in 1:((T + 6) / 7)) {
    flu_week_cases[i] ~ neg_binomial_2(flu_cases_weekly_sum[i] * adjustment_parameter2, phi_flu);
  }
  for (t in 1:T) {
    covid_cases[t] ~ neg_binomial_2(E1[t]* sigma1*adjustment_parameter1, phi_covid);
    //flu_cases[t] ~ neg_binomial_2(E2[t]*sigma2*adjustment_parameter2, phi_flu);
    //covid_cases[t] ~ neg_binomial_2(fitcov_cases[t], phi_covid);
    //flu_cases[t] ~ neg_binomial_2(fitflu_cases[t], phi_flu);
    
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
    int covid_cases_pred[T];
    int flu_week_cases_pred[(T + 6) / 7]; // Weekly flu case predictions
    real flu_weekly_infections[(T + 6) / 7]; // Weekly sum of flu infections
    real R_t_covid[T];
    real R_t_flu[T];
    real weekly_flu_sum = 0; // Initialize weekly flu sum
    int current_week = (1 + 6) / 7; // Initialize the current week based on the first day

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
        
        covid_cases_pred[t] = neg_binomial_2_rng(E1[t] * sigma1 * adjustment_parameter1 + 1e-5, phi_covid);
        
        weekly_flu_sum += E_flu[t] * sigma2;
        if (t % 7 == 0 || t == T) {
            current_week = (t + 6) / 7;
            flu_week_cases_pred[current_week] = neg_binomial_2_rng(weekly_flu_sum* adjustment_parameter2 + 1e-5, phi_flu);
            flu_weekly_infections[current_week] = weekly_flu_sum;
            weekly_flu_sum = 0;
        }

        R_t_covid[t] = (beta1 * 0.33 / gamma1) * (S_covid[t] / N);
        R_t_flu[t] = (beta2 * 0.33 / gamma2) * (S_flu[t] / N);
    }
}
