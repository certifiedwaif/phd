/*

  fastcorrelation.cpp

*/

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Rcpp;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::MatrixXd;

typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;

// From the Wikipedia page on Gray code
/*
        The purpose of this function is to convert an unsigned
        binary number to reflected binary Gray code.

        The operator >> is shift right. The operator ^ is exclusive or.
*/
unsigned int binaryToGrey(unsigned int num)
{
  return (num >> 1) ^ num;
}

/*
        The purpose of this function is to convert a reflected binary
        Gray code number to a binary number.
*/
unsigned int greyToBinary(unsigned int num)
{
  unsigned int mask;
  for (mask = num >> 1; mask != 0; mask = mask >> 1)
  {
    num = num ^ mask;
  }
  return num;
}

VectorXd binaryToVec(unsigned int num, unsigned int p)
{
  VectorXd result(p);
  for (unsigned int i = 0; i < p; i++) {
    result[(p - 1) - i] = num & 1;
    num >>= 1;
  }
  return(result);
}

MatrixXd greycode(unsigned int p)
{
  unsigned int rows = 1 << p;
  MatrixXd result(rows, p);
  for (unsigned int i = 0; i < rows; i++) {
    result.row(i) = binaryToVec(binaryToGrey(i), p).transpose();
  }
  return(result);
}

// [[Rcpp::export]]
VectorXd all_correlations(const MapVecd vy, const MapMatd mX, const MapMatd mZ)
{
  const unsigned int n = mX.rows();
  const unsigned int p = mX.cols();
  const unsigned int m = mZ.cols();

  // Generate greycode matrix
  MatrixXd mGrey = greycode(m);
  VectorXd vR2_all(mGrey.rows());
  
  // Loop through models, updating and downdating m1_inverse as necessary
  for (unsigned int row = 1; row < mGrey.rows(); row++) {
    // Construct mZ_gamma
    RowVectorXd vGreycodeRow = mGrey.row(row);
    // cout << vGreycodeRow << endl;
    unsigned int one_count = vGreycodeRow.sum();

    MatrixXd mZ_gamma(n, one_count);
    unsigned int mZ_gamma_col = 0;
    for (unsigned int mZ_col = 0; mZ_col < m; mZ_col++) {
      if (vGreycodeRow(mZ_col) == 1) {
        mZ_gamma.col(mZ_gamma_col) = mZ.col(mZ_col);
        mZ_gamma_col++;
      }
    }
    
    const unsigned int m_gamma = mZ_gamma.cols();
    MatrixXd m1(p + m_gamma, p + m_gamma); // Symmetric
    VectorXd m2(p + m_gamma);
    m1 << mX.transpose() * mX, mX.transpose() * mZ_gamma,
          mZ_gamma.transpose() * mX, mZ_gamma.transpose() * mZ_gamma;
    MatrixXd m1_inv = m1.inverse();
  
    m2 << mX.transpose() * vy,
          mZ_gamma.transpose() * vy;
  
    MatrixXd mC(n, p + m_gamma);
    mC << mX, mZ_gamma;
    VectorXd vR2 = (vy.transpose() * mC * m1_inv * m2) / vy.squaredNorm();
    vR2_all(row) = vR2(0);
  }

  return vR2_all;
}

// [[Rcpp::export]]
VectorXd one_correlation(const MapVecd vy, const MapMatd mX, const MapMatd mZ)
{
  const unsigned int n = mX.rows();
  const unsigned int p = mX.cols();
  const unsigned int m = mZ.cols();

  MatrixXd mC(n, p + m);
  mC << mX, mZ;
  
  MatrixXd m1(p + m, p + m); // Symmetric
  VectorXd m2(p + m);
  m1 << mX.transpose() * mX, mX.transpose() * mZ,
        mZ.transpose() * mX, mZ.transpose() * mZ;
  m2 << mX.transpose() * vy,
        mZ.transpose() * vy;
  VectorXd R2 = (vy.transpose() * mC * m1.inverse() * m2) / vy.squaredNorm();

  return R2;
}
