/*

  fastcorrelation.cpp

*/

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;

// [[Rcpp::export]]
VectorXd correlation(MapVecd vy, MapMatd mX, MapMatd mZ)
{
  MatrixXd mC(mX.rows(), mX.cols() + mZ.cols());
  mC << mX, mZ;
  
  MatrixXd m1(mX.cols() + mZ.cols(), mX.cols() + mZ.cols()); // Symmetric
  MatrixXd m2(mX.cols() + mZ.cols(), 1);
  m1 << mX.transpose() * mX, mX.transpose() * mZ,
        mZ.transpose() * mX, mZ.transpose() * mZ;
  m2 << mX.transpose() * vy,
        mZ.transpose() * vy;
  VectorXd R2 = (vy.transpose() * mC * (m1).inverse() * m2) / vy.squaredNorm();

  return R2;
}
