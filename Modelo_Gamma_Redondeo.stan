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
  real<lower=0, upper=1> lambda;
  vector<lower=-3.5, upper=3.5>[N] error;
}

transformed parameters {
  vector[N] Treal;
  Treal = Tobs + error; //The Jacobian es 1 there is no need
}

model {
  alpha[1] ~ cauchy(0, 2.5);
  beta[1]  ~ cauchy(0, 2.5);
  alpha[2] ~ cauchy(100, 2.5);
  beta[2]  ~ cauchy(100, 2.5);
  error    ~ uniform(-3.5, 3.5);
  lambda   ~ beta(0.25,1);
  for (n in 1:N)
  target += log_mix(lambda,
                    gamma_lpdf(Treal[n] | alpha, beta),
                    gamma_lpdf(Tobs[n] | alpha[, beta[2]));
}

generated quantities {
  real Tpred = gamma_rng(alpha[1], beta[1]);
  real Tout  = gamma_rng(alpha[2], beta[2]);
}
