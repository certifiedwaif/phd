require(Matrix)

trace = function(X)
{
  return(sum(diag(X)))
}

eye = function(p)
{
  return(diag(rep(1.0, p)))
}

fit_linear_mixed_model = function(Y, X, Z, A_epsilon, B_epsilon, A_u, B_u, K, sigma2_beta, p, n)
{
  calcLogLik = function()
  {
    
    logLik = 0
    logLik = logLik + 0.5*(p + sum(K))- 0.5*n*log(2*pi) - 0.5*p * log(sigma2_beta)
    logLik = logLik + 0.5 * log(det(Sigma_q_beta_u)) - (sum(mu_q_beta_u^2) + trace(Sigma_q_beta_u)) / (2*sigma2_beta)
    logLik = logLik + A_epsilon * log(B_epsilon) - (A_epsilon + 0.5 * n) * log(B_q_sigma2_epsilon)
    logLik = logLik + lgamma(A_epsilon + 0.5*n) - lgamma(A_epsilon)
    logLik = logLik + sum(A_u * log(B_u) - (A_u + 0.5*K) * log(B_q_sigma2_u) + lgamma(A_u + 0.5 * K) - lgamma(A_u))
    
    return(logLik)
  }
  
  C = cbind(X, Z)
  # Initialise parameters ----
  B_q_sigma2_epsilon = 1
  B_q_sigma2_u = 1
  
  # Cycle ----
  lastLogLik = -1e6
  logLik = -1e6 + 1
  while (logLik > lastLogLik)
  {
    lastLogLik = logLik
    Sigma_q_beta_u_inv = (A_epsilon + 0.5*n)/(B_q_sigma2_epsilon) * t(C) %*% C
    #browser()
    Sigma_q_beta_u_inv = Sigma_q_beta_u_inv + bdiag(eye(p) / sigma2_beta, (A_u + 0.5 * K)/(B_q_sigma2_u) * eye(K))
    Sigma_q_beta_u = solve(Sigma_q_beta_u_inv)
    mu_q_beta_u = (A_epsilon + 0.5 * n)/(B_q_sigma2_epsilon) * Sigma_q_beta_u %*% t(C) %*% Y
    B_q_sigma2_epsilon = B_epsilon + 0.5*(sum((Y - C %*% mu_q_beta_u)^2) + trace(t(C) %*% C %*% Sigma_q_beta_u))
    # Number of columns in the Z matrix
    p_x = dim(X)[2]
    p_z = dim(Z)[2]
    #print(p)
    random_int_idx = (p_x+1):(p_x+p_z)
    #print(random_int_idx)
    B_q_sigma2_u = B_u + 0.5 * (sum(mu_q_beta_u[random_int_idx]^2) + trace(Sigma_q_beta_u[random_int_idx, random_int_idx]))
    logLik = calcLogLik()
    
    params = list(Sigma_q_beta_u, mu_q_beta_u, B_q_sigma2_epsilon, B_q_sigma2_u, logLik)
  }
  return(params)
}
