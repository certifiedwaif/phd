/*

  fastdiag.cpp

  Provide a fast version of the calculation mC mLambda mC^T in C++ using Eigen

*/

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using std::vector;
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef Eigen::Map<Eigen::MatrixXd> MapMatd;

// [[Rcpp::export]]
VectorXd fastdiag(MapMatd C, MapMatd lambda)
{
  VectorXd result(C.rows());

  if (C.cols() != lambda.cols()) {
    stop("mC and mLambda do not have the same numbers of columns");
  }  

  for (int i = 0; i < C.rows(); i++) {
    result[i] = C.row(i) * lambda * C.row(i).transpose();
  }

  return result;
}

// [[Rcpp::export]]
VectorXd fastsolve(MapMatd C, MapMatd R)
{
  VectorXd result(C.rows());

  if (C.cols() != R.cols()) {
    stop("mC and mR do not have the same numbers of columns");
  }  

  // for (int i = 0; i < C.rows(); i++) {
  //   VectorXd a = R.triangularView<Eigen::Lower>().solve(C.row(i).transpose());
  //   result[i] = a.dot(a);
  // }

  MatrixXd sol = R.triangularView<Eigen::Lower>().solve(C.transpose());
  result = sol.colwise().squaredNorm();

  return result;
}
