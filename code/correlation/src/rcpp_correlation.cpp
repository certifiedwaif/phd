#include <Rcpp.h>
#if defined(_OPENMP)
  #include <omp.h>
#endif;

#include "graycode.h"
#include "correlation.h"
#include <string>

using namespace Eigen;
using namespace Rcpp;
using namespace std;

//' @importFrom Rcpp evalCpp
//' @useDynLib correlation

//' Calculate the correlation of all the sub-models of mX and vy
//'
//' @param vy Vector of responses
//' @param mX Covariate matrix
//' @param prior The g-prior to use. The choices of g-prior available are "maruyama", "BIC", "ZE",
//' "liang_g1", "liang_g2", "liang_g_n_appell", "liang_g_approx", "liang_g_n_quad",
//' "robust_bayarri1" and "robust_bayarri2"
//' @param intercept_col The index of the column in mX containing the intercept, if any
//' @param bNatural_Order Whether to return the results in natural order or graycode order. Defaults to graycode order.
//' @param bIntercept Logical value indicating whether there is an intercept column or not
//' @param bCentre Logical value indicating whether to centre the response vector and covariance matrix or not
//' @param cores The number of cores to use
//' @return A list containing 
//' vR2, the vector of correlations for each model
//' vp_gamma, the vector of number of covariates for each model
//' vlogp, the vector of logs of the likelihoods of each model
//' vinclusion_prob, the vector of inclusion probabilities for each of the covariates
//' @examples
//' library(MASS)
//'
//' mD <- UScrime
//' notlog <- c(2,ncol(UScrime))
//' mD[,-notlog] <- log(mD[,-notlog])
//'
//' for (j in 1:ncol(mD)) {
//'   mD[,j] <- (mD[,j] - mean(mD[,j]))/sd(mD[,j])
//' }
//'
//' varnames <- c(
//'   "log(AGE)",
//'   "S",
//'   "log(ED)",
//'   "log(Ex0)",
//'   "log(Ex1)",
//'   "log(LF)",
//'   "log(M)",
//'   "log(N)",
//'   "log(NW)",
//'   "log(U1)",
//'   "log(U2)",
//'   "log(W)",
//'   "log(X)",
//'   "log(prison)",
//'   "log(time)")
//'
//' y.t <- mD$y
//' X.f <- data.matrix(cbind(mD[1:15]))
//' colnames(X.f) <- varnames 
//' corr_result <- all_correlations_mX(y.t, X.f, "maruyama")
//' > str(corr_result)
//' List of 4
//'  $ vR2            : num [1:32768] 0 0.00759 0.01 0.00822 0.13921 ...
//'  $ vp_gamma       : int [1:32768] 0 1 2 1 2 3 2 1 2 3 ...
//'  $ vlogp          : num [1:32768] 6.92e-310 -8.51 -1.30e+01 -8.50 -9.74 ...
//'  $ vinclusion_prob: num [1:15] 0.284 0.054 0.525 0.679 0.344 ...
//' @export
// [[Rcpp::export]]
List all_correlations_mX(NumericVector vy, NumericMatrix mX, std::string prior, int intercept_col = 1,
													bool bNatural_Order = false, bool bIntercept = false, bool bCentre = false,
													int cores = 1) {
	Map<VectorXd> vy_m = as< Map<VectorXd> >(vy);
	Map<MatrixXd> mX_m = as< Map<MatrixXd> >(mX);
	#if defined(_OPENMP)
		omp_set_num_threads(cores);
	#endif;
	List result = all_correlations_mX_cpp(vy_m, mX_m, prior, intercept_col - 1, bNatural_Order, bIntercept,
								bCentre);
	return result;
}

//' Calculate the correlation of all the sub-models mX/mZ and vy, where mX is fixed in every model and the sub-models of mZ are included
//'
//' @param vy Vector of responses
//' @param mX Fixed covariate matrix
//' @param mZ Varying covariate matrix
//' @param prior The g-prior to use. The choices of g-prior available are "maruyama", "BIC", "ZE",
//' "liang_g1", "liang_g2", "liang_g_n_appell", "liang_g_approx", "liang_g_n_quad",
//' "robust_bayarri1" and "robust_bayarri2"
//' @param intercept_col The index of the column in mX containing the intercept, if any
//' @param bNatural_Order Whether to return the results in natural order or graycode order. Defaults to graycode order.
//' @param bIntercept Logical value indicating whether there is an intercept column or not
//' @param bCentre Logical value indicating whether to centre the response vector and covariance matrix or not
//' @param cores The number of cores to use
//' @return A list containing 
//' vR2, the vector of correlations for each model
//' vp_gamma, the vector of number of covariates for each model
//' vlogp, the vector of logs of the likelihoods of each model
//' vinclusion_prob, the vector of inclusion probabilities for each of the covariates
//' @examples
//' library(MASS)
//'
//' mD <- UScrime
//' notlog <- c(2,ncol(UScrime))
//' mD[,-notlog] <- log(mD[,-notlog])
//'
//' for (j in 1:ncol(mD)) {
//'   mD[,j] <- (mD[,j] - mean(mD[,j]))/sd(mD[,j])
//' }
//'
//' varnames <- c(
//'   "log(AGE)",
//'   "S",
//'   "log(ED)",
//'   "log(Ex0)",
//'   "log(Ex1)",
//'   "log(LF)",
//'   "log(M)",
//'   "log(N)",
//'   "log(NW)",
//'   "log(U1)",
//'   "log(U2)",
//'   "log(W)",
//'   "log(X)",
//'   "log(prison)",
//'   "log(time)")
//'
//' y.t <- mD$y
//' X.f <- data.matrix(cbind(mD[, 1:10]))
//' colnames(X.f) <- varnames 
//' Z.f <- data.matrix(cbind(mD[, 11:15]))
//' corr_result <- all_correlations_mX_mZ(y.t, X.f, Z.f, "maruyama")
//' > str(corr_result)
//' List of 4
//'  $ vR2            : num [1:32] 0 0.719 0.724 0.688 0.771 ...
//'  $ vp_gamma       : int [1:32] 0 11 12 11 12 13 12 11 12 13 ...
//'  $ vlogp          : num [1:32] 9.56e-316 -1.21e+01 -1.40e+01 -1.46e+01 -9.60 ...
//'  $ vinclusion_prob: num [1:15] 1 1 1 1 1 1 1 1 1 1 ...
//' 
//' @export
// [[Rcpp::export]]
List all_correlations_mX_mZ(NumericVector vy, NumericMatrix mX, NumericMatrix mZ, std::string prior,
                            int intercept_col = 1,
                            bool bNatural_Order = false, bool bIntercept = false, bool bCentre = false,
                            int cores = 1) {
	Map<VectorXd> vy_m = as< Map<VectorXd> >(vy);
	Map<MatrixXd> mX_m = as< Map<MatrixXd> >(mX);
	Map<MatrixXd> mZ_m = as< Map<MatrixXd> >(mZ);
	#if defined(_OPENMP)
		omp_set_num_threads(cores);
	#endif;
	List result = all_correlations_mX_mZ_cpp(vy_m, mX_m, mZ_m, prior, intercept_col - 1, bNatural_Order, bIntercept, bCentre);
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
