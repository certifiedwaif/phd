#include "taug.h"

#include <Rcpp.h>
#include <RcppEigen.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]

using namespace Rcpp;
using namespace Eigen;

//' @importFrom Rcpp evalCpp
//' @useDynLib taug

//' Calculate tau_g
//'
//' @param n The number of samples
//' @param mGraycode The graycode matrix
//' @param vR2 The vector of correlations
//' @param vlog_ZE The vector of log(ZE)
//' @export
// [[Rcpp::export]]
List rcpp_taug(unsigned int n, NumericMatrix mGraycode, NumericVector vR2,
                        NumericVector vlog_ZE) {
  Map<MatrixXd> mGraycode_m = as< Map<MatrixXd> >(mGraycode);
  Map<VectorXd> vR2_m = as< Map<VectorXd> >(vR2);
  Map<VectorXd> vlog_ZE_m = as< Map<VectorXd> >(vlog_ZE);
  VectorXd vp_m;
  VectorXd vq_m;
  tau_g(n, mGraycode_m, vR2_m, vlog_ZE_m, vp_m, vq_m);
  NumericVector vp;
  NumericVector vq;
  List result = List::create(Named("vp", vp),
                              Named("vq",  vq));
  return result;
}

//' Calculate the first integral for the posterior variance of beta_hat
//'
//' @param n The number of samples
//' @param p The number of covariates
//' @param R2 The correlation
//' @param upper_limit The upper limit of integration
//' @return The value of the integral
//' @export
// [[Rcpp::export]]
double var_int1(int n, double p, double R2, double upper_limit)
{
  auto a = -3./4.;
  auto b = (n - p) / 2. - 2. - a;
  return trapint([=](int i)->double {
    return static_cast<double>(upper_limit * (i + 1)) / static_cast<double>(GRID_POINTS);
  },
  [=](double g)->double {
    #ifndef _OPENMP
      Rcout << "Called with " << g << std::endl;
    #endif
    return exp((b + 1.) * log(g) - log(1.+g) - n/2. * log(1. + (1.-R2)*g));
  });
}


//' Calculate the second integral for the posterior variance of beta_hat
//'
//' @param n The number of samples
//' @param p The number of covariates
//' @param R2 The correlation
//' @param upper_limit The upper limit of integration
//' @return The value of the integral
//' @export
// [[Rcpp::export]]
double var_int2(int n, int p, double R2, double upper_limit)
{
  auto a = -3./4.;
  auto b = (n - p) / 2. - 2. - a;
  return trapint([=](int i)->double {
    return static_cast<double>(upper_limit * (i + 1)) / static_cast<double>(GRID_POINTS);
  },
  [=](double g)->double {
    auto result = (b + 2.) * log(g) - 2. * log(1. + g) - n / 2. * log(1. + (1.-R2) * g);
    Rcout << g << " " << result << std::endl;
    return result;
  });
}
