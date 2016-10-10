# numeric_experiments.R
#' @import purrr
#' @import dplyr
#' @import ggplot2

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
create_param_df <- function()
{
  #vp <- c(10, 20, 50, 100)
  vp <- c(10, 20, 50, 100, 500, 1000)
  vn <- c(1.1, 1.25, 1.5, 2, 5, 10)
  vR2 <- seq(from=0.00, to=1.00, by = 1e-2)

  param_df <- map(vp, function(p) {
    map(vn, function(n) {
      map(vR2, function(R2) {c(p*n, p, R2)})
    })
  }) %>% as.data.frame %>% t
  rownames(param_df) <- NULL
  colnames(param_df) <- c("vn", "vp", "vR2")
  return(as.data.frame(param_df))
}

#' @export
create_approx_exact_df <- function()
{
  param_df <- create_param_df()
  vp <- param_df$vp
  vn <- param_df$vn
  vR2 <- param_df$vR2

  exact_ints_df <- param_df %>% with(taug::g_ints(vn, vp, vR2)) %>% as.data.frame

  vtau_g <- param_df %>% with(pmap_dbl(list(vn, vp, vR2), tau_g))
  vlog_p <- param_df %>% with(pmap_dbl(list(vn, vp, vR2), log_p))
  vtau_sigma2 <- (1 - (1 + vtau_g)^(-1) * vR2)^(-1)
  vs <- (0.5 * (vn + vp) * (1 - (1 + vtau_g)^(-1)*vR2))
  vc <- 0.5 * ((1 + (1 + vtau_g)^(-1))^(-1) * (1 + vtau_g)^(-2) * vn * vR2 + (1 + vtau_g)^(-1) * vp)
  det_mXTX <- rep(1, length(vtau_g))
  det_mSigma <- vtau_sigma2 * (1 + vtau_g)^(-1) / det_mXTX
  velbo <- param_df %>% with(pmap_dbl(list(vn, vp, vc, vs, vtau_g, det_mXTX, det_mSigma), elbo))
  combined_df <- cbind(param_df, exact_ints_df) %>%
                  mutate(vtau_g = vtau_g, vtau_sigma2 = vtau_sigma2) %>%
                  mutate(approx_mean = (1 + vtau_g)^(-1), exact_mean = G_1) %>%
                  mutate(approx_var = vtau_sigma2^(-1) * (1 + vtau_g)^(-1),
                          exact_var = (vn / (vn - 2)) * (G_1 - G_2 * vR2)) %>%
                  mutate(vlog_p = vlog_p, velbo = velbo)

  return(combined_df)
}

#' @export
plot_graphs <- function() {
  # Split by n, p and plot
  combined_df <- create_approx_exact_df()
  label_fn <- function(df) {
    map(df, ~ str_c("p = ", .))
  }
  approx_means_df <- combined_df %>%
                      select(one_of(c("vp", "vn", "vR2", "approx_mean"))) %>%
                      rename(mean = approx_mean) %>%
                      mutate(mean_type = "Approximate")
  exact_means_df <- combined_df %>%
                    select(one_of(c("vp", "vn", "vR2", "exact_mean"))) %>%
                    rename(mean = exact_mean) %>%
                    mutate(mean_type = "Exact")
  means_df <- bind_rows(approx_means_df, exact_means_df) %>%
              mutate(mean_type = factor(mean_type))
  means_df %>%
    ggplot(aes(x=vR2, y=mean)) +
      geom_line(aes(group=interaction(vn, mean_type), color=vn, linetype=mean_type)) +
      facet_wrap(~vp, labeller=label_fn) +
      xlab("Correlation co-efficient") + ylab("Shrinkage")
}
