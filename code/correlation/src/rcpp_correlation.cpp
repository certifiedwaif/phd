#include <Rcpp.h>
#if defined(_OPENMP)
  #include <omp.h>
#endif;

#include "graycode.h"
#include "correlation.h"

using namespace Eigen;
using namespace Rcpp;

//' @importFrom Rcpp evalCpp
//' @useDynLib correlation

//' Calculate the correlation of all the sub-models of mX and vy
//'
//' @param vy Vector of responses
//' @param mX Covariate matrix
//' @param intercept_col The index of the column in mX containing the intercept, if any
//' @param bNatural_Order Whether to return the results in natural order or graycode order. Defaults to graycode order.
//' @param bIntercept Logical value indicating whether there is an intercept column or not
//' @param bCentre Logical value indicating whether to centre the response vector and covariance matrix or not
//' @param cores The number of cores to use
//' @return A list containing vR2, the vector of correlations for each model and vp_gamma, the vector of
//' number of covariates for each model
//' @export
// [[Rcpp::export]]
List all_correlations_mX(NumericVector vy, NumericMatrix mX, int intercept_col = 1,
													bool bNatural_Order = false, bool bIntercept = false, bool bCentre = false,
													int cores = 1) {
	Map<VectorXd> vy_m = as< Map<VectorXd> >(vy);
	Map<MatrixXd> mX_m = as< Map<MatrixXd> >(mX);
	#if defined(_OPENMP)
		omp_set_num_threads(cores);
	#endif;
	List result = all_correlations_mX_cpp(vy_m, mX_m, intercept_col - 1, bNatural_Order, bIntercept, bCentre);
	return result;
}

//' Calculate the correlation of all the sub-models mX/mZ and vy, where mX is fixed in every model and the sub-models of mZ are included
//'
//' @param vy Vector of responses
//' @param mX Fixed covariate matrix
//' @param mZ Varying covariate matrix
//' @param intercept_col The index of the column in mX containing the intercept, if any
//' @param bNatural_Order Whether to return the results in natural order or graycode order. Defaults to graycode order.
//' @param bIntercept Logical value indicating whether there is an intercept column or not
//' @param bCentre Logical value indicating whether to centre the response vector and covariance matrix or not
//' @param cores The number of cores to use
//' @return A list containing vR2, the vector of correlations and vp_gamma, the vector of number of covariates
//' for each model
//' @export
// [[Rcpp::export]]
List all_correlations_mX_mZ(NumericVector vy, NumericMatrix mX, NumericMatrix mZ,
                            int intercept_col = 1,
                            bool bNatural_Order = false, bool bIntercept = false, bool bCentre = false,
                            int cores = 1) {
	Map<VectorXd> vy_m = as< Map<VectorXd> >(vy);
	Map<MatrixXd> mX_m = as< Map<MatrixXd> >(mX);
	Map<MatrixXd> mZ_m = as< Map<MatrixXd> >(mZ);
	#if defined(_OPENMP)
		omp_set_num_threads(cores);
	#endif;
	List result = all_correlations_mX_mZ_cpp(vy_m, mX_m, mZ_m, intercept_col - 1, bNatural_Order, bIntercept, bCentre);
	return result;
}

//' Return the graycode matrix
//'
//' @param varying The number of covariates varying in the graycode matrix
//' @param fixed The number of fixed covariates in the graycode matrix. These covariates will always be included
//' @return The graycode matrix. The number of fixed columns will be included in the lower indexed columns
//' as 1s, while the higher indexed columns will varying depending on whether each covariate in the varying
//' set of covariates is included or not.
//' @export
// [[Rcpp::export]]
IntegerMatrix graycode(unsigned int varying, unsigned int fixed = 0) {
	Graycode gray(fixed, varying);
	MatrixXi result = gray.to_MatrixXi();
	return wrap(result);
}
