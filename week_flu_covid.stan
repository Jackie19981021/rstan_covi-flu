data {
  int<lower=1> T;       // 时间点数量
  //int covid_cases[T];   // COVID-19病例数
  //int flu_cases[T];     // 流感病例数
  int<lower=1> W;       // 周数量
  int covid_cases_week[W]; // 每周COVID-19病例数
  int flu_cases_week[W];   // 每周流感病例数
  real N;               // 总人口数
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
  //real<lower=-100,upper=-5> E10;         // 初始COVID-19潜伏人群
  //real<lower=-100,upper=-5> I10;         // 初始COVID-19感染人群
  real<lower=0,upper=7000000> E10;         // 初始COVID-19潜伏人群
  real<lower=0,upper=7000000> I10; // 初始COVID-19感染人群
  real<lower=0> R10;         // 初始COVID-19恢复人群
  //real<lower=-100,upper=-5> E20;         // 初始流感潜伏人群 
  //real<lower=-100,upper=-5> I20;         // 初始流感感染人群
  real<lower=0,upper=7000000> E20;         // 初始流感潜伏人群 
  real<lower=0,upper=7000000> I20; // 初始流感感染人群
  real<lower=0> R20;         // 初始流感恢复人群
  real<lower=0,upper=1> susc1;
  real<lower=0,upper=1> susc2;
  real<lower=0,upper=1> adjustment_parameter1;
  real<lower=0,upper=1> adjustment_parameter2; // 调整参数
  real<lower=0,upper=1> alpha1;//npi参数_covid
  real<lower=0,upper=1> alpha2;//npi参数
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
  real E1_week[W];    // COVID-19每周潜伏人群
  real I1_week[W];    // COVID-19每周感染人群 
  real E2_week[W];    // 流感每周潜伏人群
  real I2_week[W];    // 流感每周感染人群
  
    // 初始化每周累积量为0
  for (w in 1:W) {
    E1_week[w] = 0;
    I1_week[w] = 0;
    E2_week[w] = 0;
    I2_week[w] = 0;
  }
  // S[1] = S0;
  // E1[1] = E10;
  // I1[1] = I10; 
  // R1[1] = R10;
  // E2[1] = E20;
  // I2[1] = I20;
  // R2[1] = R20;
  
  
  // 
  // E1[1] =exp(E10)*susc1*N;
  // I1[1] =exp(I10)*susc1*N;
  // S[1] =N-E10-E20-I10-I20-R10-R20-S120-S210;
  // S12[1] =susc2*N-I20-E20;
  // S21[1] =susc1*N-I10-E10;
  // R1[1] =(1-susc1)*N;
  // E2[1] =exp(E20)*susc2*N;
  // I2[1] =exp(I20)*susc2*N;
  // R2[1] =(1-susc2)*N;
  
  
    S[1] =N*(1-susc1-susc2)-E10-E20-I10-I20-R10-R20;
    S12[1] =N*susc1;
    S21[1] =N*susc2;
  // // 
    E1[1] =E10;
    I1[1] =I10;
    R1[1] =R10;
  // // 
    E2[1] =E20;
    I2[1] =I20;
    R2[1] =R20;
  
  
  
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
    
    
     real dS = -(beta1* exp(-alpha1*a[t-1]) * I1[t-1] + beta2* exp(-alpha2*a[t-1]) * I2[t-1]) * S[t-1] / N;
     real dE1 = beta1*(1-eta1)* exp(-alpha1*a[t-1])*S21[t-1]*I1[t-1]/N+beta1* exp(-alpha1*a[t-1])*S[t-1]*I1[t-1]/N-sigma1*E1[t-1];
     //real dI1[t-1] * ((S[t-1] - rho12 * R1[t-1] * I2[t-1] / N) / N);
     real dI1 = sigma1 * E1[t-1] - gamma1 * I1[t-1];
     real dR1 = gamma1 * I1[t-1] - rho1 * R1[t-1];
     real dE2 = beta2*(1-eta2)* exp(-alpha2*a[t-1])*S12[t-1]*I2[t-1]/N+beta2* exp(-alpha2*a[t-1])*S[t-1]*I2[t-1]/N-sigma2*E2[t-1];
    // real dE2 = beta2 * I2[t-1] * ((S[t-1] - rho21 * R2[t-1] * I1[t-1] / N) / N);
     real dI2 = sigma2 * E2[t-1] - gamma2 * I2[t-1];
     real dR2 = gamma2 * I2[t-1] - rho2 * R2[t-1];
     real dS21 = -(1-eta1)*beta1* exp(-alpha1*a[t-1]) * I1[t-1] * S21[t-1] / N+ rho2 * R2[t-1];
     real dS12 = -(1-eta2)*beta2* exp(-alpha2*a[t-1]) * I2[t-1] * S12[t-1] / N+ rho1 * R1[t-1];
    
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
// 计算每周的累积量
for (w in 1:W) {
  for (t in ((w-1)*7+1):(w*7)) {
    E1_week[w] += E1[t];
    I1_week[w] += I1[t];
    E2_week[w] += E2[t];
    I2_week[w] += I2[t];
  }
}  
  
  
  
}

model {
  // 先验分布
  // beta1 ~ normal(1, 0.5);
  // beta2 ~ normal(0.4, 1); 
  // sigma1~ normal(0.5, 0.5); 
  // sigma2~ normal(0.2, 0.5);
  // gamma1~ normal(0.2, 0.5);
  // gamma2~ normal(0.2, 0.5);
  // rho2 ~ normal(0.01, 0.01);
  // rho1 ~ normal(0.05, 0.05);
  phi_covid_inv ~ exponential(5);
  phi_flu_inv ~ exponential(5);
  //I2 ~ normal(5, 5);
  // 似然函数
  // 使用周数据拟合模型
  for (w in 1:W) {
    covid_cases_week[w] ~ neg_binomial_2(E1_week[w] * sigma1 * adjustment_parameter1, phi_covid);
    flu_cases_week[w] ~ neg_binomial_2(E2_week[w] * sigma2 * (100000/N)* adjustment_parameter2, phi_flu);
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
    real R_t_covid[T];
    real R_t_flu[T];
    
    real npi_effect_on_rt_covid[T]; // NPI对COVID-19 R_t的每日影响
    real npi_effect_on_rt_flu[T];   // NPI对流感 R_t的每日影响
    
    real E1_week_pred[W]; // COVID-19每周预测潜伏人群
    real I1_week_pred[W]; // COVID-19每周预测感染人群
    real E2_week_pred[W]; // 流感每周预测潜伏人群
    real I2_week_pred[W]; // 流感每周预测感染人群
    
    int covid_cases_week_pred[W]; // COVID-19每周预测病例数
    int flu_cases_week_pred[W];   // 流感每周预测病例数

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

        real R_t_covid_no_npi = (beta1 / gamma1) * (S_covid[t] / N);
        real R_t_flu_no_npi = (beta2 / gamma2) * (S_flu[t] / N);
    
        R_t_covid[t] = (beta1 * exp(-alpha1 * a[t]) / gamma1) * (S_covid[t] / N);
        R_t_flu[t] = (beta2 * exp(-alpha2 * a[t]) / gamma2) * (S_flu[t] / N);
    
        npi_effect_on_rt_covid[t] = R_t_covid_no_npi - R_t_covid[t];
        npi_effect_on_rt_flu[t] = R_t_flu_no_npi - R_t_flu[t];
    }

    for (w in 1:W) {
        E1_week_pred[w] = 0;
        I1_week_pred[w] = 0;
        E2_week_pred[w] = 0;
        I2_week_pred[w] = 0;
        covid_cases_week_pred[w] = 0;
        flu_cases_week_pred[w] = 0;

        for (t in ((w-1)*7+1):(w*7)) {
            E1_week_pred[w] += E1[t];
            I1_week_pred[w] += I1[t];
            E2_week_pred[w] += E2[t];
            I2_week_pred[w] += I2[t];
        }
        covid_cases_week_pred[w] = neg_binomial_2_rng(E1_week_pred[w] * sigma1 * adjustment_parameter1 + 1e-5, phi_covid);
        flu_cases_week_pred[w] = neg_binomial_2_rng(E2_week_pred[w] * sigma2* (100000/N) * adjustment_parameter2 + 1e-5, phi_flu);
    }
}
