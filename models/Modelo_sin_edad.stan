data {
  int<lower=0> N;
  vector[N] Tobs;
  vector[N] Sexo;

  //Edades a interpolar en el resultado para reporte
  int<lower=0> M; 
  vector[M] EdadesResultado; 
}

parameters {
  vector[2] eta;
  vector<lower=0>[2] nu;
  vector<lower=0>[2] beta;
  
  //Probability of outlier
  real<lower=0, upper=1> theta;
  
  //Roundong and outlier error as params (latent variables)
  vector[N] error1;
  vector[N] error2;
}

transformed parameters {
  vector[N] Treal;
  vector[N] Outlier;
  vector[N] mu_real;
  vector[N] mu_outlier;
  
  //Parámetro de medias de la regresión
  mu_real    = nu[1]  + eta[1]* Sexo;
  mu_outlier = nu[2]  + eta[2]* Sexo;
  
  //Datos con redondeo en función de los observados
  Treal   = Tobs + error1; 
  Outlier = Tobs + error2; 
}

model {
  nu       ~ cauchy(0, 2.5);
  beta     ~ cauchy(0, 2.5);
  eta      ~ normal(0, 1000);
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
  vector[M] Tpred_1;
  
  //Sims
  for (m in 1:M){
    Tpred_0[m]         = gamma_rng( nu[1] ./ beta[1], beta[1]);
    Tpred_1[m]         = gamma_rng((nu[1] + eta[1])./ beta[1], beta[1]);
  }

}
