// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Eigen/IterativeLinearSolvers>

using std::vector;
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::ConjugateGradient;

typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::SparseMatrix<double> SpMat;

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
SpMat fastinv(const MapMatd Rd, const int p, const int m, const int blocksize, const int spline_dim)
{
  // Construct sparse version of R
  
  typedef Eigen::Triplet<double> T;
  std:vector<T> triplets;
  triplets.reserve(m*blocksize/2 + p*p/2 + p*m*blocksize);
  // Insert entries for random effects
  for (int m_idx = 0; m_idx < m; m_idx++) {
    // 1 for the first row
    // ...
    // b for the bth row
    for (int b_idx = 0; b_idx < blocksize; b_idx++) {
      const int offset = m_idx*blocksize;
      for (int row_idx = offset; row_idx < offset+b_idx; row_idx++) {
        int i = offset+b_idx;
        int j = offset+row_idx;
        triplets.push_back(T(i, j, Rd(i, j)));
      }
    }
  }
  
  // Insert entries for fixed effects
  for (int p_idx = m*blocksize; p_idx < m*blocksize+p; p_idx++) {
    for (int row_idx = 0; row_idx < p_idx; row_idx++) {
      int i = p_idx;
      int j = row_idx;
      triplets.push_back(T(i, j, Rd(i, j)));
    }
  }

  SpMat Rsp;
  Rsp.setFromTriplets(triplets.begin(), triplets.end());
  
  // Solve for RHS = I
  ConjugateGradient<SparseMatrix<double> > solver;
  solver.compute(Rsp);
  SparseMatrix<double> I(Rd.rows(),Rd.cols());
  I.setIdentity();
  
  return solver.solve(I);
}