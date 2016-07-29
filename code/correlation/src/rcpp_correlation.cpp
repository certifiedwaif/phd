#include <Rcpp.h>
#include "correlation.hpp"

using namespace Eigen;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rcpp_all_correlations_mX(NumericVector vy, NumericMatrix mX, int intercept_col,
										bool bIntercept, bool bCentre) {
	return all_correlations_mX_cpp(vy, mX, intercept_col, bIntercept, bCentre);
}

// [[Rcpp::export]]
NumericVector rcpp_all_correlations_mX_mZ(NumericVector vy, NumericMatrix mX, NumericMatrix mZ,
					  int intercept_col, bool bIntercept, bool bCentre) {
	return all_correlations_mX_mZ_cpp(vy, mX, mZ, intercept_col, bIntercept, bCentre);
}
