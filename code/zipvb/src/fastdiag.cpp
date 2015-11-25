/*

  fastdiag.cpp

  Provide a fast version of the calculation mC mLambda mC^T in C++ using Eigen

*/

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include "zipvb_types.hpp"

// [[Rcpp::export]]
VectorXd fastdiag(MapMatd lambda, MapMatd C)
{

  if (C.cols() != lambda.cols()) {
    stop("mC and mLambda do not have the same numbers of columns");
  }

  VectorXd result(C.rows());

  for (int i = 0; i < C.rows(); i++) {
    result[i] = C.row(i) * lambda * C.row(i).transpose();
  }

  return result;
}

// [[Rcpp::export]]
VectorXd fastdiag2(MapMatd R, MapMatd C)
{

  if (C.cols() != R.cols()) {
    stop("mR and mR do not have the same numbers of columns");
  }

  // for (int i = 0; i < C.rows(); i++) {
  //   result[i] = (C.row(i) * R.triangularView<Eigen::Lower>()).squaredNorm();
  // }

  VectorXd result = (C * R.triangularView<Eigen::Lower>()).rowwise().squaredNorm();

  return result;
}

// [[Rcpp::export]]
VectorXd fastsolve(MapMatd R, MapMatd C)
{

  if (C.cols() != R.cols()) {
    stop("mC and mR do not have the same numbers of columns");
  }

  // for (int i = 0; i < C.rows(); i++) {
  //   VectorXd a = R.triangularView<Eigen::Lower>().solve(C.row(i).transpose());
  //   result[i] = a.dot(a);
  // }

  MatrixXd sol = R.triangularView<Eigen::Lower>().solve(C.transpose());
  VectorXd result = sol.colwise().squaredNorm();

  return result;
}
