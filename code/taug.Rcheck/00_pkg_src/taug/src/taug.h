#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppGSL)]]
using Eigen::DenseBase;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixXi;


void tau_g(int n, const MatrixXd& mGraycode, const VectorXd& vR2, const VectorXd& vlog_ZE,
           VectorXd& vp, VectorXd& vq);
