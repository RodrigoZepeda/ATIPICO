//Modelo de redondeo semanal
// ------------------------------------------------------------------------
// Modelo para ajuste de tiempos de recuperación. 
// Treal = Tiempo de recuperación real (días)
// Tobs  = Tiempo de recuperación observado con fuerte énfasis en redondeos
// de caracter semanal (7, 14, 28, etc). 
// theta = Probabilidad de redondeo (no todos los datos son redondeados)
// El modelo está dado por:
//        Treal ~ Gamma(alpha, beta)
//        error ~ Uniforme(-0.5, 0.5)
// Con probabilidad theta:
//    Tobs = redondeo_7(Treal) = 7[Treal/7 + error]
// con probabilidad 1 - theta:
//    Tobs = Treal
// Lo cual implica que:
//    Treal = Tobs - 7*error

data {
  int<lower=0> N;
  vector[N] Tobs;
}

parameters {
  real<lower=0> alpha[2];
  real<lower=0> beta[2];
  real<lower=0, upper=1> theta;
  vector<lower=-3.5, upper=3.5>[N] error1;
  vector<lower=0, upper=1000>[N] error2;
}

transformed parameters {
  vector[N] Treal;
  vector[N] Outlier;
  Treal   = Tobs + error1; 
  Outlier = Tobs + error2; 
}

model {
  alpha    ~ cauchy(0, 2.5);
  beta     ~ cauchy(0, 2.5);
  error1   ~ uniform(-3.5, 3.5);
  error2   ~ uniform(0, 1000);
  theta   ~ beta(0.25,1);
  for (n in 1:N)
  target += log_mix(theta,
                    gamma_lpdf(Treal[n] | alpha[1], beta[1]),
                    gamma_lpdf(Outlier[n] | alpha[2], beta[2]));
}

generated quantities {
  real Tpred       = gamma_rng(alpha[1], beta[1]);
  real Outlierpred = gamma_rng(alpha[2], beta[2]);
}
