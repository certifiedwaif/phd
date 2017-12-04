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
  tbl$time <- map2_dbl(tbl$package, tbl$prior, function(package, prior.val) {
    cat("package", package, "prior.val", prior.val, "\n")
    start_time <- proc.time()[3]
    if (package == "BAS") {
      library(BAS)
      bas.lm(y~X, prior=prior.val, model=uniform())
    }
    if (package == "BVS") {
      library(BayesVarSel)
      Bvs(formula="y~.", data=data.frame(y=y, X=X), prior.betas=prior.val, prior.models="Constant",
          time.test=FALSE, n.keep=50000)
    }
    if (package == "BMS") {
      library(BMS)
      bms(cbind(y, X), nmodel=50000, mcmc="enumerate", g="hyper=3", mprior="uniform")
    }
    if (package == "BLMA") {
      library(correlation)
      all_correlations_mX(vy, mX, prior.val)
    }
    end_time <- proc.time()[3]
    end_time - start_time
  })
  tbl
}
