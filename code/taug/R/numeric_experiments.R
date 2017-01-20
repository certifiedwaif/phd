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

#' The exact variance of g/(1 + g)
#'
#' @param n The number of observations
#' @param p The number of covariates
#' @param R2 The correlation co-efficient squared
#' @return The exact variance of g/(1 + g)
#' @export
exact_var_g <- function(n, p, R2)
{
  G_1 <- E_g_one_plus_g(n, p, R2)
  G_2 <- E_g_one_plus_g_squared(n, p, R2)
  exact_var_g <- G_2 - G_1^2
  return(exact_var_g)
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
      map(vR2, function(R2) {c(n, p*n, p, R2)})
    })
  }) %>% as.data.frame %>% t
  rownames(param_df) <- NULL
  colnames(param_df) <- c("vkappa", "vn", "vp", "vR2")
  return(as.data.frame(param_df))
}

E_q_g <- function(n, p, R2, c)
{
  a <- -3. / 4.
  b <- (n - p) / 2. - 2. - a
  result <- exp(-c) * c^(n/2. - p - a - 1) / 2. * gamma((n - p) / 2. - n / 2. + p + a) * Re(whittakerW(c, -3. * n / 4. - a / 2. + 1. / 2., -n / 4. + p / 2. + a / 2.))
  cat("result ", result, "\n")
  return(result)
}

#' @export
create_approx_exact_df <- function()
{
  param_df <- create_param_df() # %>% filter(vp == 10. && vkappa == 1.1)
  vp <- param_df$vp
  vn <- param_df$vn
  vkappa <- param_df$vkappa
  vR2 <- param_df$vR2

  exact_ints_df <- g_ints(vn, vp, vR2) %>% as.data.frame

  vtau_g <- pmap_dbl(list(vn, vp, vR2), tau_g)
  vlog_p <- pmap_dbl(list(vn, vp, vR2), log_p)
  vtau_sigma2 <- (1 - (1 + vtau_g)^(-1) * vR2)^(-1)
  vexact_precision <- pmap_dbl(list(vn, vp, vR2), exact_precision)
  vexact_precision <- 1 / vexact_precision
  vr <- (vn + vp) / 2
  vs <- (0.5 * (vn + vp) * (1 - (1 + vtau_g)^(-1)*vR2))
  E_q_sigma2_inv <- vs / (vr - 1)
  #vc <- 0.5 * ((1 + (1 + vtau_g)^(-1))^(-1) * (1 + vtau_g)^(-2) * vn * vR2 + (1 + vtau_g)^(-1) * vp)
  vc <- 0.5 * vn * vR2 / ((1 + vtau_g) * (1 - vR2 + vtau_g)) + 0.5 * vp / (1 + vtau_g)
  det_mXTX <- rep(1, length(vtau_g))
  log_det_mXTX <- log(det_mXTX)
  det_mSigma <- log(vtau_sigma2)^vp * (1 + vtau_g)^(-vp) / det_mXTX
  log_det_mSigma <- vp * log(vtau_sigma2) - vp * log(1 + vtau_g) - log_det_mXTX
  velbo <- pmap_dbl(list(vn, vp, vc, vs, vtau_g, log_det_mXTX, log_det_mSigma), elbo)

  exact_var_g <- pmap_dbl(list(vn, vp, vR2), exact_var_g)
  approx_var_g <- pmap_dbl(list(vn, vp, vR2), var_g_over_one_plus_g)

  vE_g_y <- pmap_dbl(list(vn, vp, vR2), E_g_y)
  vE_q_g <- pmap_dbl(list(vn, vp, vR2, vc), E_q_g_c)

  vsigma2 <- seq(from=.01, to=10, length.out=length(vn))
  # vp_sigma2_y <- pmap_dbl(list(vn, vp, vR2, vsigma2), p_sigma2_y)
  # vq_sigma2 <- pmap_dbl(list(vr, vs, vsigma2), q_sigma2)
  vaccuracy_sigma2 <- pmap_dbl(list(vn, vp, vR2, vr, vs), accuracy_sigma2)
  vaccuracy_g <- pmap_dbl(list(vn, vp, vR2, vc), accuracy_g)

  combined_df <- cbind(param_df, exact_ints_df) %>%
                  mutate(vtau_g = vtau_g, vtau_sigma2 = vtau_sigma2) %>%
                  mutate(approx_mean = (1 + vtau_g)^(-1), exact_mean = G_1) %>%
                  mutate(approx_var = vtau_sigma2^(-1) * (1 + vtau_g)^(-1),
                          exact_var = (vn / (vn - 2)) * (G_1 - G_2 * vR2)) %>%
                  mutate(vlog_p = vlog_p, velbo = velbo) %>%
                  mutate(exact_var_g = exact_var_g, approx_var_g = approx_var_g) %>%
                  mutate(E_q_sigma2_inv = E_q_sigma2_inv) %>%
                  mutate(vexact_precision = vexact_precision) %>%
                  mutate(vE_g_y = vE_g_y) %>%
                  mutate(vE_q_g = vE_q_g) %>%
                  mutate(vsigma2 = vsigma2) %>%
                  # mutate(vp_sigma2_y = vp_sigma2_y) %>%
                  # mutate(vq_sigma2 = vq_sigma2) %>%
                  mutate(vaccuracy_sigma2 = vaccuracy_sigma2) %>%
                  mutate(vaccuracy_g = vaccuracy_g)

  return(combined_df)
}

#' @export
plot_graphs <- function() {
  # Split by n, p and plot
  combined_df <- create_approx_exact_df()
  label_fn <- function(df) {
    map(df, ~ str_c("p = ", .))
  }

  pdf("Shrinkage.pdf")
  approx_means_df <- combined_df %>%
                      select(one_of(c("vp", "vkappa", "vR2", "approx_mean"))) %>%
                      rename(mean = approx_mean, Kappa = vkappa) %>%
                      mutate(Type = "Approximate")
  exact_means_df <- combined_df %>%
                    select(one_of(c("vp", "vkappa", "vR2", "exact_mean"))) %>%
                    rename(mean = exact_mean, Kappa = vkappa) %>%
                    mutate(Type = "Exact")
  means_df <- bind_rows(approx_means_df, exact_means_df) %>%
              mutate(Type = factor(Type))
  means_df %>%
    ggplot(aes(x=vR2, y=mean)) +
      geom_line(aes(group=interaction(Kappa, Type), color=Kappa, linetype=Type)) +
      facet_wrap(~vp, labeller=label_fn) +
      xlab(TeX("$R^2$")) + ylab("Shrinkage") + theme_bw() + scale_colour_gradient(low="#777777", high="#ffffff")
  dev.off()

  # approx_log_p_df <- combined_df %>%
  #                     select(one_of(c("vp", "vkappa", "vR2", "velbo"))) %>%
  #                     rename(log_p = velbo) %>%
  #                     mutate(log_p_type = "Approximate")
  # exact_log_p_df <- combined_df %>%
  #                     select(one_of(c("vp", "vkappa", "vR2", "vlog_p"))) %>%
  #                     rename(log_p = vlog_p) %>%
  #                     mutate(log_p_type = "Exact")
  # log_p_df <- bind_rows(approx_log_p_df, exact_log_p_df) %>%
  #             mutate(log_p_type = factor(log_p_type))
  pdf("Relative_error_log_p.pdf")
  combined_df %>% filter(vR2 <= 0.99) %>%
    mutate(relative_error=(vlog_p - velbo)/vlog_p) %>%
    rename(Kappa = vkappa) %>%
    ggplot(aes(x=vR2, y=relative_error)) +
      geom_line(aes(group=Kappa, color=Kappa)) +
      facet_wrap(~vp, labeller=label_fn) +
      xlab(TeX("$R^2$")) + ylab("Relative error") + ylim(-.15, 2) + theme_bw() + scale_colour_gradient(low="#777777", high="#ffffff")
  dev.off()

  pdf("Log_of_Relative_error_of_Variance_of_g.pdf")
  combined_df %>% filter(vR2 <= 0.99) %>%
    mutate(relative_error=log10((approx_var_g - exact_var_g)/exact_var_g)) %>%
    rename(Kappa = vkappa) %>%
    ggplot(aes(x=vR2, y=relative_error)) +
      geom_line(aes(group=Kappa, color=Kappa)) +
      facet_wrap(~vp, labeller=label_fn) +
      xlab(TeX("$R^2$")) + ylab(TeX("$\\log_{10}$ of relative error")) + theme_bw() + scale_colour_gradient(low="#777777", high="#ffffff")
  dev.off()

  pdf("Variance.pdf")
  approx_var_df <- combined_df %>%
                      select(one_of(c("vp", "vkappa", "vR2", "approx_var"))) %>%
                      rename(var = approx_var, Kappa = vkappa) %>%
                      mutate(Type = "Approximate")
  exact_var_df <- combined_df %>%
                    select(one_of(c("vp", "vkappa", "vR2", "exact_var"))) %>%
                    rename(var = exact_var, Kappa = vkappa) %>%
                    mutate(Type = "Exact")
  var_df <- bind_rows(approx_var_df, exact_var_df) %>%
              mutate(Type = factor(Type))
  var_df %>%
    ggplot(aes(x=vR2, y=var)) +
      geom_line(aes(group=interaction(Kappa, Type), color=Kappa, linetype=Type)) +
      facet_wrap(~vp, labeller=label_fn) +
      xlab(TeX("$R^2$")) + ylab("Variance") + theme_bw() + scale_colour_gradient(low="#777777", high="#ffffff")
  dev.off()

  pdf("Precision.pdf")
  approx_precision_df <- combined_df %>%
                      select(one_of(c("vp", "vkappa", "vR2", "E_q_sigma2_inv"))) %>%
                      rename(precision = E_q_sigma2_inv, Kappa = vkappa) %>%
                      mutate(Type = "Approximate")
  exact_precision_df <- combined_df %>%
                    select(one_of(c("vp", "vkappa", "vR2", "vexact_precision"))) %>%
                    rename(precision = vexact_precision, Kappa = vkappa) %>%
                    mutate(Type = "Exact")
  precision_df <- bind_rows(approx_precision_df, exact_precision_df) %>%
              mutate(Type = factor(Type))
  precision_df %>%
    ggplot(aes(x=vR2, y=precision)) +
      geom_line(aes(group=interaction(Kappa, Type), color=Kappa, linetype=Type)) +
      facet_wrap(~vp, labeller=label_fn) +
      xlab(TeX("$R^2$")) + ylab("Precision") + theme_bw() + scale_colour_gradient(low="#777777", high="#ffffff")
  dev.off()

  pdf("Accuracy_sigma2.pdf")
  combined_df %>% dplyr::filter(vaccuracy_sigma2 >= 0.0 & vaccuracy_sigma2 <= 1.0) %>%
    rename(Kappa = vkappa) %>%
    ggplot(aes(x=vR2, y=vaccuracy_sigma2)) +
    geom_line(aes(group=Kappa, color=Kappa)) +
    facet_wrap(~vp, labeller=label_fn) +
    xlab(TeX("$R^2$")) + ylab(TeX("Accuracy $\\sigma^2$")) + theme_bw() + scale_colour_gradient(low="#777777", high="#ffffff")
  dev.off()

  pdf("Accuracy_g.pdf")
  combined_df %>%
    rename(Kappa = vkappa) %>%
    ggplot(aes(x=vR2, y=vaccuracy_g)) +
    geom_line(aes(group=Kappa, color=Kappa)) +
    facet_wrap(~vp, labeller=label_fn) +
    xlab(TeX("$R^2$")) + ylab("Accuracy g") + theme_bw() + scale_colour_gradient(low="#777777", high="#ffffff")
  dev.off()
}

plot_vw1_vw2 <- function() {
  bodyfat_vw1 <- read.csv(file = "bodyfat_vw1.csv", header=FALSE)
  GradRate_vw1 <- read.csv(file = "GradRate_vw1.csv", header=FALSE)
  Hitters_vw1 <- read.csv(file = "Hitters_vw1.csv", header=FALSE)
  USCrime_vw1 <- read.csv(file = "USCrime_vw1.csv", header=FALSE)
  Wage_vw1 <- read.csv(file = "Wage_vw1.csv", header=FALSE)
  bodyfat_vw2 <- read.csv(file = "bodyfat_vw2.csv", header=FALSE)
  GradRate_vw2 <- read.csv(file = "GradRate_vw2.csv", header=FALSE)
  Hitters_vw2 <- read.csv(file = "Hitters_vw2.csv", header=FALSE)
  USCrime_vw2 <- read.csv(file = "USCrime_vw2.csv", header=FALSE)
  Wage_vw2 <- read.csv(file = "Wage_vw2.csv", header=FALSE)
}
