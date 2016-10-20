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

const int GRID_POINTS = 10000;

void tau_g(int n, const MatrixXd& mGraycode, const VectorXd& vR2, const VectorXd& vlog_ZE,
           VectorXd& vp, VectorXd& vq);
double trapint(std::function<double(int)> x, std::function<double(double)> f);
double whittakerW(double z, double kappa, double mu);
double log_p(int n, int p, double R2);
double elbo(int n, int p, double c, double s, double tau_g, double log_det_XTX, double log_det_mSigma);
