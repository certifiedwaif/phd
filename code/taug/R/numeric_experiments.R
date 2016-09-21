# numeric_experiments.R
#' @import purrr

#' @export
compare_means <- function(n, p, R2, beta_hat_LS)
{
  # Compare E(g/(1+g) | y) beta_hat_LS with (1 + tau_g)^{-1} beta_hat_LS
  # beta_hat_LS appears in both expressions, so we can just compare
  # E(g/(1+g) | y) with (1 + tau_g)^{-1}
  G_1 <- E_g_one_plus_g(n, p, R2)
  exact_mean <- G_1 * beta_hat_LS
  tau_g <- tau_g(n, p, R2)
  approx_mean <- (1 + tau_g)^(-1) * beta_hat_LS
  diff_mean <- abs(exact_mean - approx_mean)

  return(diff_mean)
}

#' @export
compare_covs <- function(n, p, R2, beta_hat_LS, mXTX_inv)
{
  # Compare (G_2 - G_1^2) beta_hat_LS beta_hat_LS^T + (n / (n - 2)) (G_1 - G_2 R^2) (X^T X)^{-1}
  # with tau_sigma^2 (1 + tau_g)^{-1} (X^T X)^{-1}
  G_1 <- E_g_one_plus_g(n, p, R2)
  G_2 <- E_g_one_plus_g_squared(n, p, R2)
  exact_cov <- (G_2 - G_1^2) * tcrossprod(beta_hat_LS) + (n / (n - 2)) * (G_1 - G_2 * R2) * mXTX_inv
  tau_g <- tau_g(n, p, R2)
  tau_sigma2 <- (1 - (1 + tau_g)^(-1) * R2)^(-1)
  approx_cov <- tau_sigma2^(-1) * (1 + tau_g)^(-1) * mXTX_inv
  diff_cov <- abs(exact_cov - approx_cov)

  return(diff_cov)
}

#' @export
compare_approx_exact <- function()
{
  vp <- c(10, 20, 50, 100, 500, 1000)
  vn <- c(1.1, 1.25, 1.5, 2, 5, 10)
  vR2 <- seq(from=0.00, to=1.00, by = 1e-2)
  df <- map(vp, function(p) {
          map(vn, function(n) {
            map(vR2, function(R2) {c(p*n, p, R2)})
          })
        }) %>% as.data.frame %>% t
  colnames(df) <- c("vn", "vp", "vR2")
  df <- as.data.frame(df)
  exact_ints <- 1:nrow(df) %>% map(function(x) {
    cat(x, "\n")
    with(df[x,], taug::g_ints(vn, vp, vR2))
  })
}
