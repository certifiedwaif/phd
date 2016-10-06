#include "taug.h"

#include <Rcpp.h>
#include <RcppEigen.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>
#include <sstream>

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


void handler (const char * reason,
              const char * file,
              int line,
              int gsl_errno)
{
  std::stringstream ss;
  ss << reason << " in file " << file << " on " << line << ", errorno " << gsl_errno << std::endl;
  throw std::range_error(ss.str());
}


//' Calculate the first integral for the posterior variance of beta_hat
//'
//' @param n The number of samples
//' @param p The number of covariates
//' @param R2 The correlation
//' @return The value of the integral
//' @export
// [[Rcpp::export]]
double E_g_one_plus_g(int n, double p, double R2)
{
  gsl_set_error_handler_off();
  auto a = -3./4.;
  auto b = (n - p) / 2. - 2. - a;
  gsl_sf_result result;
  auto alpha_1 = p / 2. + a + 1.;
  auto beta_1 = b + 2.;
  gsl_sf_beta_e(alpha_1, beta_1, &result);
  auto num = result.val;
  auto alpha_2 = p / 2. + a + 1.;
  auto beta_2 = b + 1.;
  gsl_sf_beta_e(p / 2. + a + 1., beta_2, &result);
  auto den = result.val;
  gsl_sf_hyperg_2F1_e(p / 2. + a + 1., 1., n / 2. + 1., R2, &result);
  auto conf = result.val;

  if (den == 0.0)
    return NA_REAL;
  else
    return (num / den) * conf;
}


//' Calculate the second integral for the posterior variance of beta_hat
//'
//' @param n The number of samples
//' @param p The number of covariates
//' @param R2 The correlation
//' @return The value of the integral
//' @export
// [[Rcpp::export]]
double E_g_one_plus_g_squared(int n, int p, double R2)
{
  gsl_set_error_handler_off();

  auto a = -3./4.;
  auto b = (n - p) / 2. - 2. - a;
  gsl_sf_result result;
  gsl_sf_beta_e(p / 2. + a + 1., b + 3., &result);
  auto num = result.val;
  gsl_sf_beta_e(p / 2. + a + 1., b + 1., &result);
  auto den = result.val;
  gsl_sf_hyperg_2F1_e(p / 2. + a + 1., 2, n / 2. + 2., R2, &result);
  auto conf = result.val;

  if (den == 0.0)
    return NA_REAL;
  else
    return (num / den) * conf;
}


//' Calculate posterior variance integrals based on vectors of n, p and R2 values
//'
//' @param vn A numeric vector of n values
//' @param vp A numeric vector of p values
//' @param vR2 A numeric value of
//' @return List of the values of the posterior variance integrals. E_g_one_plus_g contains the result of the first
//' integral, while E_g_one_plus_g_squared contains the second.
//' @export
// [[Rcpp::export]]
List g_ints(NumericVector vn, NumericVector vp, NumericVector vR2)
{
  if (vn.size() != vp.size()) throw std::range_error("vn and vp are different sizes");
  if (vn.size() != vR2.size()) throw std::range_error("vn and vR2 are different sizes");

  VectorXd var_ints_1(vn.size());
  VectorXd var_ints_2(vn.size());

  // #pragma omp parallel for
  for (auto i = 0; i < vn.size(); i++) {
    var_ints_1(i) = E_g_one_plus_g(vn[i], vp[i], vR2[i]);
    var_ints_2(i) = E_g_one_plus_g_squared(vn[i], vp[i], vR2[i]);
  }

  NumericVector result1(Rcpp::wrap(var_ints_1));
  NumericVector result2(Rcpp::wrap(var_ints_2));

  return List::create(Named("G_1") = result1,
                      Named("G_2") =  result2);
}


//' Calculate marginal log-likelihood log p(y)
//'
//' @param n The number of observations
//' @param p The number of covariates
//' @param R2 The correlation co-efficient, squared
//' @return The value of the marginal log-likelihood log p(y)
//' @export
// [[Rcpp::export]]
double log_p(int n, int p, double R2)
{
  auto a = -3./4.;
  auto b = (n - p) / 2. - 2. - a;

  auto result = gsl_sf_lngamma(p/2. + a + 1.) - n/2. * log(n * PI) - gsl_sf_lngamma((n - p) / 2.);
  result = result - gsl_sf_lngamma(a + 1.) - ((n - p) / 2. - a - 1.) * log(1 - R2);

  return result;
}

//' Calculate the variational lower bound
//'
//' @param n The number of observations
//' @param p The number of covariates
//' @param c The variational parameter c
//' @param s The variational parameter s
//' @param tau_g The variational parameter tau_g
//' @param det_XTX The determinant of X^T X
//' @param det_mSigma The determinant of mSigma
//' @return The variational lower bound
//' @export
// [[Rcpp::export]]
double elbo(int n, int p, double c, double s, double tau_g, double det_XTX, double det_mSigma)
{
  auto a = -3./4.;
  auto b = (n - p) / 2. - 2. - a;
  auto beta = c;
  auto nu = n / 2. - p - a - 1.;
  auto mu = -(n - p/2.) + 1.;
  auto Z = pow(beta, (nu - 1.) / 2.) * gsl_sf_gamma(1. - mu - nu) * exp(- beta/(2.*nu)) * whittakerW(beta, (nu - 1.)/2 + mu, -nu / 2.);

  return p / 2. - n / 2. * log(2 * PI) + 0.5 * log(det_XTX) - gsl_sf_lnbeta(a + 1., b + 1.) \
          - (n + p) / 2. * log(s) + gsl_sf_lngamma((n + p) / 2.) + c * tau_g + log(Z) + 0.5 * log(det_mSigma);
}