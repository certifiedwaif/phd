#include <Rcpp.h>
#include "correlation.hpp"

using namespace Eigen;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rcpp_all_correlations_mX(NumericVector vy, NumericMatrix mX, int intercept_col,
										bool bIntercept, bool bCentre) {
	Map<VectorXd> vy_m = as< Map<VectorXd> >(vy);
	Map<MatrixXd> mX_m = as< Map<MatrixXd> >(mX);
	VectorXd result = all_correlations_mX_cpp(vy_m, mX_m, intercept_col, bIntercept, bCentre);
	NumericVector result(wrap(result));
	return result;
}

// [[Rcpp::export]]
NumericVector rcpp_all_correlations_mX_mZ(NumericVector vy, NumericMatrix mX, NumericMatrix mZ,
					  int intercept_col, bool bIntercept, bool bCentre) {
	Map<VectorXd> vy_m = as< Map<VectorXd> >(vy);
	Map<MatrixXd> mX_m = as< Map<MatrixXd> >(mX);
	Map<MatrixXd> mZ_m = as< Map<MatrixXd> >(mZ);
	VectorXd result = all_correlations_mX_mZ_cpp(vy_m, mX_m, mZ_m, intercept_col, bIntercept, bCentre);
	NumericVector result(wrap(result));
	return result;
}
