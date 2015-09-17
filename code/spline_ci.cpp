/*

  spline_ci.cpp

*/

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "eigenmvn.h"

using std::vector;
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;

// [[Rcpp::export]]
Rcpp::List spline_ci(int num_reps, MapMatd mC_x, MapVecd vmu_vb, MapMatd mLambda_vb, MapMatd vbeta, MapMatd vu)
{
  VectorXd f_hat_vb(num_reps);
  VectorXd f_hat_mcmc(num_reps);
  
  int sample_idx;
  Eigen::EigenMultivariateNormal<double> vb_nu(mC_x * vmu_vb, mC_x * mLambda_vb * mC_x.transpose());

  // for (j in 1:NUM_REPS) {
  for (int j = 0; j < num_reps; j++) {
    //   # Generate samples - simulate vbeta, vu
    //   f_hat_vb[j] <- t(rmvnorm(1, mC_x %*% vmu_vb, mC_x %*% mLambda_vb %*% t(mC_x)))
    f_hat_vb[j] = vb_nu.samples(1)(0, 0);
    //   sample_idx <- sample(dim(fit$vbeta)[1], 1)
    sample_idx = rand() % vbeta.rows();
    //   vbeta_mcmc <- fit$vbeta[sample_idx, ]
    VectorXd vbeta_mcmc = vbeta.row(sample_idx);
    //   vu_mcmc <- fit$vu[sample_idx, ]
    VectorXd vu_mcmc = vu.row(sample_idx);
    //   vnu_mcmc <- c(vbeta_mcmc, vu_mcmc)
    VectorXd vnu_mcmc(vbeta_mcmc.rows() + vu_mcmc.rows());
    vnu_mcmc << vbeta_mcmc, vu_mcmc;
    //   f_hat_mcmc[j] <- t(mC_x %*% vnu_mcmc)
    f_hat_mcmc[j] = (mC_x * vnu_mcmc)(0, 0);
  // }
  }

  return Rcpp::List::create(Rcpp::Named("f_hat_vb") = f_hat_vb,
                            Rcpp::Named("f_hat_mcmc") = f_hat_mcmc);
}
                                                                                                                                                                                                  