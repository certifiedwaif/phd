data {
  int<lower=1> N; // Number of observations
  int<lower=1> P; // Number of fixed effects covariates
  int<lower=1> B; // Block size
  int<lower=1> M; // Number of subjects
  int<lower=0> y[N]; // Estimated treatment effects
  vector[P] X[N]; // Fixed effects covariate matrix
  vector[M*B] Z[N]; // Random effects covariate matrix
}

transformed data {
  vector[P] zeros_beta;
  vector[B] zeros_u;

  for (i in 1:P)
    zeros_beta[i] <- 0.0;

  for (i in 1:B)
    zeros_u[i] <- 0.0;
}

parameters {
  vector[P] vbeta; 
  vector[B] u[M];
  cov_matrix[B] S;
  cov_matrix[P] BetaPriorInv;
  real<lower=0, upper=1> rho;
  cov_matrix[B] sigma_u;
}

model {
  real eta;

  rho ~ beta(.01, .01);
  sigma_u ~ inv_wishart(B+1, S);

  vbeta ~ multi_normal_prec(zeros_beta, BetaPriorInv);
  
  for (m in 1:M)
    u[m] ~ multi_normal(zeros_u, sigma_u);
  
  for (n in 1:N) {    
    eta <- dot_product(X[n], vbeta);
    for (m in 1:M)
      for (b in 1:B)
        eta <- eta + Z[n,(m-1)*B+b]*u[m][b];

    if (y[n] == 0)
      increment_log_prob(log_sum_exp(bernoulli_log(0, rho),
        bernoulli_log(1, rho) + poisson_log(y[n], exp(eta))));
    else
      increment_log_prob(bernoulli_log(1, rho) + poisson_log(y[n], exp(eta)));
  };
}