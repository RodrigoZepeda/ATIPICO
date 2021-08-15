data {
  int<lower=0> N;
  vector[N] Tobs;
  vector[N] Edad;
  
  //Edades a interpolar en el resultado para reporte
  int<lower=0> M; 
  vector[M] EdadesResultado; 
}

parameters {
  vector[2] gamma;
  vector<lower=0>[2] nu;
  vector<lower=0>[2] beta;
  
  //Probability of outlier
  real<lower=0, upper=1> theta;
  
  //Roundong and outlier error as params (latent variables)
  vector<lower=-3.5, upper=3.5>[N] error1;
  vector<lower=0, upper=1000>[N] error2;
}

transformed parameters {
  vector[N] Treal;
  vector[N] Outlier;
  vector[N] mu_real;
  vector[N] mu_outlier;
  
  //Parámetro de medias de la regresión
  mu_real    = nu[1] + gamma[1]*Edad;
  mu_outlier = nu[2] + gamma[2]*Edad;
  
  //Datos con redondeo en función de los observados
  Treal   = Tobs + error1; 
  Outlier = Tobs + error2; 
}

model {
  nu       ~ cauchy(0, 2.5);
  beta     ~ cauchy(0, 2.5);
  gamma    ~ normal(0, 1000);
  error1   ~ uniform(-3.5, 3.5);
  error2   ~ uniform(0, 1000);
  theta    ~ beta(0.25,1);
  for (n in 1:N)
  target += log_mix(theta,
                    gamma_lpdf(Treal[n] | mu_real[n] / beta[1], beta[1]),
                    gamma_lpdf(Outlier[n] | mu_outlier[n] / beta[2], beta[2]));
}

generated quantities {
  
  vector[M] Tpred_0;
  vector[M] Outlierpred_0;
  
  real mu_real_pred_0;
  real mu_outlier_pred_0;

  
  for (m in 1:M){
    //Sims
    mu_real_pred_0    = nu[1] + gamma[1]*EdadesResultado[m];
    mu_outlier_pred_0 = nu[2] + gamma[2]*EdadesResultado[m];
  
    Tpred_0[m]         = gamma_rng(mu_real_pred_0 ./ beta[1], beta[1]);
    Outlierpred_0[m]   = gamma_rng(mu_outlier_pred_0 ./ beta[2], beta[2]);
  }
}
