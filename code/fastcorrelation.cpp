/*

  fastcorrelation.cpp

*/

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
#include <omp.h>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <boost/tokenizer.hpp>
#include <boost/dynamic_bitset.hpp>
#include "correlation.hpp"
  
using namespace Rcpp;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::RowVectorXi;
using Eigen::MatrixXi;
using namespace std;

typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;


// [[Rcpp::export]]
VectorXd all_correlations_mX(MapVecd vy, MapMatd mX, unsigned int intercept_col,
                              bool bIntercept, bool bCentre, unsigned int num_threads=1)
{
  omp_set_num_threads(num_threads);
  return all_correlations_mX_cpp(vy, mX, intercept_col, bIntercept, bCentre);
}

// [[Rcpp::export]]
VectorXd all_correlations_mX_mZ(MapVecd vy, MapMatd mX, MapMatd mZ, unsigned int intercept_col,
                                bool bIntercept, bool bCentre, unsigned int num_threads=1)
{
  omp_set_num_threads(num_threads);
  return all_correlations_mX_mZ_cpp(vy, mX, mZ, intercept_col, bIntercept, bCentre);
}

// [[Rcpp::export]]
MatrixXi graycode(unsigned int varying, unsigned int fixed, unsigned int num_threads=1)
{
  omp_set_num_threads(num_threads);
  Graycode graycode(fixed, varying);
  return graycode.to_MatrixXi();
}
