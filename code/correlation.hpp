// correlation.hpp

#pragma once

// "If you want work with matrices, you should use C++ with Eigen or Armadillo. It's pretty fast." - Hadley Wickham,
// completely unprompted.

#include <string>
#include <Eigen/Dense>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::string;

MatrixXd parseCSVfile_double(string infilename);
VectorXd all_correlations_mX_mZ_cpp(VectorXd vy, MatrixXd mX, MatrixXd mZ, const uint intercept_col,
																		const bool bIntercept = false, const bool bCentre = true);
VectorXd all_correlations_mX_cpp(VectorXd vy, MatrixXd mX, const uint intercept_col,
																 const bool bIntercept = false, const bool bCentre = true);