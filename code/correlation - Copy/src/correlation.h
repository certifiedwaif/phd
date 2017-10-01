// correlation.hpp

#include <string>
#include <Rcpp.h>
#include <RcppEigen.h>
typedef unsigned int uint;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

// "If you want work with matrices, you should use C++ with Eigen or Armadillo. It's pretty fast." - Hadley Wickham,
// completely unprompted.

#pragma once

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::string;

MatrixXd parseCSVfile_double(string infilename);
VectorXd all_correlations_mX_mZ_cpp(VectorXd vy, MatrixXd mX, MatrixXd mZ, const uint intercept_col,
																		const bool bIntercept = false, const bool bCentre = true,
																		int cores = 1);
VectorXd all_correlations_mX_cpp(VectorXd vy, MatrixXd mX, const uint intercept_col,
																 const bool bIntercept = false, const bool bCentre = true,
																 int cores = 1);
