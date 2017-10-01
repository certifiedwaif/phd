#include "taug.h"

#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppEigen)]]

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
