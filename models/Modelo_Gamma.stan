data {
  int<lower=0> N;
  vector[N] Tobs;
  
  //Edades a interpolar en el resultado para reporte
  int<lower=0> M; 
  vector[M] EdadesResultado; 
}

parameters {
  real<lower=0> alpha[1];
  real<lower=0> beta[1];
  vector<lower=-3.5, upper=3.5>[N] error1;
}

transformed parameters {
  vector[N] Treal;
  Treal   = Tobs + error1; 
}

model {
  alpha    ~ cauchy(0, 2.5);
  beta     ~ cauchy(0, 2.5);
  error1   ~ uniform(-3.5, 3.5);
  for (n in 1:N)
  target += gamma_lpdf(Treal[n] | alpha[1], beta[1]);
}

generated quantities {
  vector[M] Tpred_0;
  vector[M] Tpred_1;
  
  for (m in 1:M){
    Tpred_0[m]     = gamma_rng(alpha[1], beta[1]);
    Tpred_1[m]     = gamma_rng(alpha[1], beta[1]);
  }
}
