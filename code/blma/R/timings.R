#' Time various model selection methods, and produce a table of the results
#'
#' @param y The vector y of the response variable
#' @param X The matrix X of covariates
#'
#' @return A tibble of timing results with the columns: package, prior and time
#'
#' @export
timings <- function(y, X)
{
  # Package, prior
  packages <- c("BLMA",
                "BLMA",
                "BAS",
                "BAS",
                "BVS",
                "BMS",
                "BLMA",
                "BLMA",
                "BAS",
                "BLMA",
                "BLMA",
                "BLMA",
                "BVS",
                "BLMA",
                "BLMA")

  priors <- c("BIC",
              "ZE",
              "hyper-g",
              "hyper-g-laplace",
              "Liangetal",
              "g",
              "liang_g1",
              "liang_g2",
              "hyper-g-n",
              "liang_g_n_appell",
              "liang_g_n_quad",
              "liang_g_n_approx",
              "Robust",
              "robust_bayarri1",
              "robust_bayarri2")

  tbl <- tibble(package=packages, prior=priors)
  result <- map2(tbl$package, tbl$prior, function(package, prior.val) {
    cat("package", package, "prior.val", prior.val, "\n")
    start_time <- proc.time()[3]
    if (package == "BAS") {
      library(BAS)
      bas_result <- bas.lm(y~X, prior=prior.val, model=uniform())
      vinclusion_prob <- bas_result$probne0
    }
    if (package == "BVS") {
      library(BayesVarSel)
      bvs_result <- Bvs(formula="y~.", data=data.frame(y=y, X=X), prior.betas=prior.val, prior.models="Constant",
                        time.test=FALSE, n.keep=50000)
      vinclusion_prob <- bvs_result$inclprob
    }
    if (package == "BMS") {
      library(BMS)
      bms_result <- bms(cbind(y, X), nmodel=50000, mcmc="enumerate", g="hyper=3", mprior="uniform")
      vinclusion_prob <- coef(bms_result, order.by.pip=FALSE)[,1]
    }
    if (package == "BLMA") {
      library(correlation)
      library(appell)
      blma_result <- blma(y, X, prior.val)
      vinclusion_prob <- blma_result$vinclusion_prob
    }
    end_time <- proc.time()[3]
    time <- end_time - start_time
    return(list(time=time, vinclusion_prob=vinclusion_prob))
  })
  tbl$time <- map_dbl(result, ~.x$time)
  tbl$vinclusion_prob <- map(result, ~.x$vinclusion_prob)
  tbl
}
