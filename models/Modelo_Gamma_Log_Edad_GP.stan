data {
  int<lower=0> Ninput;
  vector[Ninput] Tobs;
  real Edad[Ninput];
  
  //Edades a interpolar en el resultado para reporte
  int<lower=0> Npredict; 
  real EdadesPredict[Npredict]; 
}

transformed data {
  //To remain positive definite
  real delta = 1e-9;
  
  //Vector of predictors see
  //https://mc-stan.org/docs/2_27/functions-reference/covariance.html
  int<lower = 0> N = Npredict + Ninput;
  real EdadModel[N];
  EdadModel[1:Ninput]   = Edad[1:Ninput];
  EdadModel[1:Npredict] = EdadesPredict[1:Npredict];
}

parameters {
  
  //Probability of outlier
  real<lower=0, upper=1> theta;
  
  //Intercept and variance
  vector<lower=0>[2] nu;
  vector<lower=0>[2] beta;

  //Parameters for GP
  real<lower=0> rho;
  real<lower=0> alpha;
  real<lower=0> sigma;
  vector[N] eta;

  //Roundong and outlier error as params (latent variables)
  vector<lower=-3.5, upper=3.5>[N] error1;
  vector<lower=0, upper=1000>[N]   error2;

}

transformed parameters {
  vector[N] Treal;
  vector[N] Outlier;
  
  //Gaussian process
  vector[N] f;
  {
    matrix[N, N] L_K;
    matrix[N, N] K = gp_exp_quad_cov(EdadModel, alpha, rho);

    L_K = cholesky_decompose(add_diag(K, delta));
    f = L_K * eta;
  }
  
  //Error for outliers
  Treal   = Tobs + error1; 
  Outlier = Tobs + error2; 
  
}

model {
  //Priors to be defined later
  nu       ~ cauchy(0, 2.5);
  beta     ~ cauchy(0, 2.5);
  rho      ~ inv_gamma(5, 5);
  alpha    ~ std_normal();
  sigma    ~ std_normal();
  eta      ~ std_normal();
  error1   ~ uniform(-3.5, 3.5);
  error2   ~ uniform(0, 1000);
  theta    ~ beta(0.25,1);
  for (n in 1:Ninput)
  target += log_mix(theta,
                    gamma_lpdf(Treal[n] | (nu[1] + f[n]) / beta[1], beta[1]),
                    gamma_lpdf(Outlier[n] | nu[2] / beta[2], beta[2]));
}

generated quantities {
  
  vector[Npredict] Tpred;
  for (n2 in 1:Npredict)
    Tpred[n2] = gamma_rng( (nu[1] + f[Ninput + n2])/beta[1], beta[1]);
  
}
