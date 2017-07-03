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

using namespace Rcpp;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixBase;
using std::string;

MatrixXd parseCSVfile_double(string infilename);
List all_correlations_mX_mZ_cpp(VectorXd vy, MatrixXd mX, MatrixXd mZ, const uint intercept_col,
																		const bool bNatural_Order = false, const bool bIntercept = false, const bool bCentre = true);
List all_correlations_mX_cpp(VectorXd vy, MatrixXd mX, const uint intercept_col,
																		const bool bNatural_Order = false, const bool bIntercept = false, const bool bCentre = true);
template <typename Derived1, typename Derived2>
Eigen::MatrixBase<Derived2>& get_cols(const Eigen::MatrixBase<Derived1>& m1, const dbitset& gamma, Eigen::MatrixBase<Derived2>& m2);
template <typename Derived1, typename Derived2>
Eigen::MatrixBase<Derived2>& get_rows(const Eigen::MatrixBase<Derived1>& m1, const dbitset& gamma, Eigen::MatrixBase<Derived2>& m2);
template <typename Derived1, typename Derived2>
Eigen::MatrixBase<Derived2>& rank_one_update(const dbitset& gamma, const uint col_abs, const uint min_idx,
	const uint fixed,
const Eigen::MatrixBase<Derived1>& mXTX, const Eigen::MatrixBase<Derived1>& mA, Eigen::MatrixBase<Derived2>& mA_prime, bool& bLow);
template <typename Derived1, typename Derived2>
Eigen::MatrixBase<Derived2>& rank_one_downdate(const uint col_abs, const uint min_idx, const uint fixed,
const Eigen::MatrixBase<Derived1>& mA, Eigen::MatrixBase<Derived2>& mA_prime);
void centre(VectorXd& v);
