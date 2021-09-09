data {
  int<lower=0> N;
  vector[N] Tobs;
  vector[N] Edad;
  
  //Edades a interpolar en el resultado para reporte
  int<lower=0> M; 
  vector[M] EdadesResultado; 
}

parameters {
  vector[1] gamma;
  vector<lower=0>[1] nu;
  vector<lower=0>[1] beta;
  
  //Roundong and outlier error as params (latent variables)
  vector<lower=-3.5, upper=3.5>[N] error1;
}

transformed parameters {
  vector[N] Treal;
  vector[N] Outlier;
  vector[N] mu_real;
  vector[N] mu_outlier;
  
  //Parámetro de medias de la regresión
  mu_real    = nu[1] + gamma[1]*log(Edad);
  
  //Datos con redondeo en función de los observados
  Treal   = Tobs + error1; 
  
}

model {
  nu       ~ cauchy(0, 2.5);
  beta     ~ cauchy(0, 2.5);
  gamma    ~ normal(0, 1000);
  error1   ~ uniform(-3.5, 3.5);
  for (n in 1:N)
  target += gamma_lpdf(Treal[n] | mu_real[n] / beta[1], beta[1]);
}

generated quantities {
  
  vector[M] Tpred_0;
  real mu_real_pred_0;
  
  
  for (m in 1:M){
    //Sims
    mu_real_pred_0    = nu[1] + gamma[1]*log(EdadesResultado[m]);
    Tpred_0[m]         = gamma_rng(mu_real_pred_0 ./ beta[1], beta[1]);
  }
}
