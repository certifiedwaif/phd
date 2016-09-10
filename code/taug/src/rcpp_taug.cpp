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
List taug(unsigned int n, NumericMatrix mGraycode, NumericVector vR2,
                        NumericVector vlog_ZE) {
  Map<MatrixXd> mGraycode_m = as< Map<MatrixXd> >(mGraycode);
  Map<VectorXd> vR2_m = as< Map<VectorXd> >(vR2);
  Map<VectorXd> vlog_ZE_m = as< Map<VectorXd> >(vlog_ZE);
  VectorXd vp_m;
  VectorXd vq_m;
  tau_g(n, mGraycode_m, vR2_m, vlog_ZE_m, vp_m, vq_m);
  NumericVector vp(wrap(vp_m));
  NumericVector vq(wrap(vq_m));
  List result = List::create(Named("vp", vp),
                              Named("vq",  vq));
  return result;
}


//' Calculate the first integral for the posterior variance of beta_hat
//'
//' @param n The number of samples
//' @param p The number of covariates
//' @param R2 The correlation
//' @return The value of the integral
//' @export
// [[Rcpp::export]]
double var_int1(int n, double p, double R2)
{
  auto a = -3./4.;
  auto b = (n - p) / 2. - 2. - a;

  return (gsl_sf_beta(p / 2. + a + 1., b + 2.) / gsl_sf_beta(p / 2. + a + 1., b + 1.)) * gsl_sf_hyperg_2F1(p / 2. + a + 1., 1., n / 2. + 1., R2);
}


//' Calculate the second integral for the posterior variance of beta_hat
//'
//' @param n The number of samples
//' @param p The number of covariates
//' @param R2 The correlation
//' @return The value of the integral
//' @export
// [[Rcpp::export]]
double var_int2(int n, int p, double R2)
{
  auto a = -3./4.;
  auto b = (n - p) / 2. - 2. - a;

  return (gsl_sf_beta(p / 2. + a + 1., b + 3.) / gsl_sf_beta(p / 2. + a + 1., b + 1.)) * gsl_sf_hyperg_2F1(p / 2. + a + 1., 2, n / 2. + 2., R2);
}


//' Calculate posterior variance integrals based on vectors of n, p and R2 values
//'
//' @param vn A numeric vector of n values
//' @param vp A numeric vector of p values
//' @param vR2 A numeric value of 
//' @return List of the values of the posterior variance integrals. var_int1 contains the result of the first
//' integral, while var_int2 contains the second.
//' @export
// [[Rcpp::export]]
List var_ints(NumericVector vn, NumericVector vp, NumericVector vR2)
{
  if (vn.size() != vp.size()) throw std::range_error("vn and vp are different sizes");
  if (vn.size() != vR2.size()) throw std::range_error("vn and vR2 are different sizes");

  VectorXd var_ints_1(vn.size());
  VectorXd var_ints_2(vn.size());

  #pragma omp parallel for
  for (auto i = 0; i < vn.size(); i++) {
    var_ints_1(i) = var_int1(vn[i], vp[i], vR2[i]);
    var_ints_2(i) = var_int2(vn[i], vp[i], vR2[i]);
  }

  NumericVector result1(Rcpp::wrap(var_ints_1));
  NumericVector result2(Rcpp::wrap(var_ints_2));

  return List::create(Named("var_int1") = result1,
                      Named("var_int2") =  result2);
}
