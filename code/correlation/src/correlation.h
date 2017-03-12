// correlation.hpp

#include <string>
#include <Rcpp.h>
#include <RcppEigen.h>
#include "graycode.h"
typedef unsigned int uint;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]

// "If you want work with matrices, you should use C++ with Eigen or Armadillo. It's pretty fast." - Hadley Wickham,
// completely unprompted.

#pragma once

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixBase;
using std::string;

MatrixXd parseCSVfile_double(string infilename);
VectorXd all_correlations_mX_mZ_cpp(VectorXd vy, MatrixXd mX, MatrixXd mZ, const uint intercept_col,
																		const bool bIntercept = false, const bool bCentre = true);
VectorXd all_correlations_mX_cpp(VectorXd vy, MatrixXd mX, const uint intercept_col,
																 const bool bIntercept = false, const bool bCentre = true);
template <typename Derived1, typename Derived2>
MatrixBase<Derived2>& get_cols(const MatrixBase<Derived1>& m1, const dbitset& gamma, MatrixBase<Derived2>& m2);
template <typename Derived1, typename Derived2>
MatrixBase<Derived2>& get_rows(const MatrixBase<Derived1>& m1, const dbitset& gamma, MatrixBase<Derived2>& m2);
