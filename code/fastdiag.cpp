// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>

using std::vector;
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::SimplicialCholesky;
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
  
  //std::cout << "Built triplets list" << std::endl;
  triplets.reserve(m*blocksize/2 + p*p/2 + p*m*blocksize);
  // Insert entries for random effects
  for (int m_idx = 0; m_idx < m; m_idx++) {
    // 1 for the first row
    // ...
    // b for the bth row
    for (int b_idx = 0; b_idx < blocksize; b_idx++) {
      const int offset = m_idx*blocksize;
      //std::cout << "offset " << offset << std::endl;
      for (int row_idx = offset; row_idx < offset+b_idx+1; row_idx++) {
        //std::cout << "row_idx" << row_idx << std::endl;
        int i = offset+b_idx;
        int j = offset+row_idx;
        //std::cout << "loop: i " << i << " j " << j << " Rd(i, j) " << Rd(i, j) << std::endl;

        triplets.push_back(T(i, j, Rd(i, j)));
      }
    }
  }
  
  //std::cout << "loop2" << std::endl;
  // Insert entries for fixed effects
  for (int p_idx = m*blocksize; p_idx < m*blocksize+p; p_idx++) {
    //std::cout << "p_idx " << p_idx << std::endl;
    for (int row_idx = 0; row_idx < p + m*blocksize; row_idx++) {
      int i = p_idx;
      int j = row_idx;
      
      //std::cout << "row_idx " << row_idx << std::endl;
      
      //std::cout << "loop2: i " << i << " j " << j << " Rd(i, j) " << Rd(i, j) << std::endl;

    triplets.push_back(T(i, j, Rd(i, j)));
    }
  }
  
  //std::cout << "Constructing sparse matrix" << std::endl;
  SpMat Rsp(Rd.rows(), Rd.cols());
  Rsp.setFromTriplets(triplets.begin(), triplets.end());
  
  //std::cout << Rsp << std::endl;
  
  // Solve for RHS = I
  //std::cout << "Solving" << std::endl;
  Eigen::SimplicialCholesky<SpMat, Eigen::Lower, Eigen::ColMajor> solver;
  //std::cout << "maxIterations" << solver.maxIterations() << std::endl;
  //solver.setTolerance(1e-99);
  solver.compute(Rsp);
  SparseMatrix<double> I(Rd.rows(),Rd.cols());
  I.setIdentity();
  //std::cout << "I " << I << std::endl;
  //std::cout << "Returning" << std::endl;
  return solver.solve(I);
}

// [[Rcpp::export]]
VectorXd fastsolve(SpMat R, SpMat C)
{
  Eigen::SimplicialCholesky<SpMat, Eigen::Lower, Eigen::ColMajor> solver;
  solver.compute(R);
  return solver.solve(C);
}