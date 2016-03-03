/*

  fastcorrelation.cpp

*/

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(BH)]]
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


typedef boost::dynamic_bitset<> dbitset;
const bool NUMERIC_FIX = false;

// From the Wikipedia page on Gray code
/*
  The purpose of this function is to convert an unsigned
  binary number to reflected binary Gray code.

  The operator >> is shift right. The operator ^ is exclusive or.
*/
unsigned int binary_to_grey(unsigned int num)
{
  return (num >> 1) ^ num;
}


/*
  The purpose of this function is to convert a reflected binary
  Gray code number to a binary number.
*/
unsigned int grey_to_binary(unsigned int num)
{
  unsigned int mask;
  for (mask = num >> 1; mask != 0; mask = mask >> 1) {
    num = num ^ mask;
  }
  return num;
}


VectorXd binary_to_vec(unsigned int num, unsigned int p)
{
  VectorXd result(p);
  for (unsigned int i = 0; i < p; i++) {
    result[(p - 1) - i] = num & 1;
    num >>= 1;
  }
  return(result);
}


VectorXd grey_vec(unsigned int i, unsigned int p)
{
  return binary_to_vec(binary_to_grey(i), p).transpose();
}


MatrixXd greycode(unsigned int p)
{
  unsigned int rows = 1 << p;
  MatrixXd result(rows, p);
  for (unsigned int i = 0; i < rows; i++) {
    result.row(i) = grey_vec(i, p);
  }
  return(result);
}


MatrixXd& sherman_morrison(MatrixXd& mA_inv, const VectorXd vu, const VectorXd vv)
{
  // MatrixXd mA_inv_tilde(mA_inv.rows(), mA_inv.cols());
  const double lambda = (vv.transpose() * mA_inv * vu).value();

  mA_inv += -((mA_inv * vu) * (vv.transpose() * mA_inv)) / (1 + lambda);

  return mA_inv;
}


dbitset& greycode(const unsigned int idx, const unsigned int p, dbitset& bs_ret)
{
  dbitset bs(p, idx);
  bs =  bs ^ (bs >> 1);
  bs_ret = bs;

  return bs_ret;
}


dbitset& greycode_change(const unsigned int idx, const unsigned int p, dbitset& gamma, bool& update,
unsigned int& min_idx, unsigned int& diff_idx,
unsigned int& bits_set)
{
  dbitset bs_curr(p);
  dbitset bs_prev(p);

  bs_curr = greycode(idx, p, bs_curr);
  bs_prev = greycode(idx - 1, p, bs_prev);
  #ifdef DEBUG
  cout << "Previous gamma: " << bs_prev << endl;
  cout << "Current gamma:  " << bs_curr << endl;
  #endif

  // Find the LSB.
  min_idx = min(bs_prev.find_first(), bs_curr.find_first());

  // Find bit that has changed.
  diff_idx = (bs_curr ^ bs_prev).find_first();

  // Has it been set, or unset?
  update = bs_curr[diff_idx];

  bits_set = bs_curr.count();

  gamma = bs_curr;
  return gamma;
}


MatrixXd& get_mX_gamma(const MatrixXd& mX, const dbitset& gamma, MatrixXd& mX_gamma)
{
  vector<unsigned int> gamma_columns;
  unsigned int p_gamma = 0;

  // Find which columns are set in gamma, and how many
  for (unsigned int i = 0; i < gamma.size(); i++) {
    if (gamma[i]) {
      gamma_columns.push_back(i);
      p_gamma++;
    }
  }

  // Extract columns from mX, put into mX_gamma
  for (unsigned int i = 0; i < p_gamma; i++) {
    mX_gamma.col(i) = mX.col(gamma_columns[i]);
  }

  return mX_gamma;
}


double sd(const VectorXd& x)
{
  const unsigned int n = x.size();
  return (x.array() - x.mean()).square().sum() / (n - 1);
}


void centre(VectorXd& v)
{
  const unsigned int n = v.size();
  VectorXd v_bar(n);
  v_bar.fill(v.mean());
  double sd_v = sd(v);

  if (sd_v > 0.0) {
    VectorXd centred_v = v - v_bar;
    v = centred_v / sd_v;
  }
  else {
    // If sd_v is zero, the data is constant. So leave it alone.
  }
}


void show_matrix_difference(ostream& os, const MatrixXd& m1, const MatrixXd& m2, const double epsilon = 1e-8)
{
  // Check that m1 and m2 have the same dimensions.
  if (!(m1.rows() == m2.rows() && m1.cols() == m2.cols())) {
    os << "Dimensions of m1 and m2 do not match, cannot display difference." << endl;
    return;
  }

  // Iterate through the elements of m1 and m2, looking for and reporting differences
  for (unsigned int i = 0; i < m1.rows(); i++) {
    for (unsigned int j = 0; j < m1.cols(); j++) {
      if (abs(m1(i, j) - m2(i, j)) > epsilon) {
        os << "Row " << i << ", column " << j << " m1 " << m1(i, j) << " m2 " << m2(i, j);
        os <<  " difference " << m1(i, j) - m2(i, j);
        os << " relative difference " << (m1(i, j) - m2(i, j)) / m1(i, j) << endl;
      }
    }
  }
}


MatrixXd& reorder_last_row_column_to_ith(MatrixXd& m, const unsigned int i, const unsigned int p_gamma)
{
  // Construct a permutation matrix mPerm which interchanges the p_gamma-th and i-th rows/columns,
  // because at the moment, the p_gamma-th row/column of mA_prime contains the inverse row/column for
  // the i-th row/column of mX_gamma_prime.transpose() * mX_gamma_prime.
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> mPerm(p_gamma);
  // mPerm.setIdentity();
  // Get p_gamma index
  int *ind = mPerm.indices().data();
  for (unsigned int idx = 0; idx < p_gamma; idx++) {
    if (idx < i) {
      ind[idx] = idx;
    }
    else if (idx == i) {
      ind[idx] = p_gamma - 1;
    }
    else {
      ind[idx] = idx - 1;
    }
  }

  #ifdef DEBUG
  cout << endl << "mPerm " << mPerm.toDenseMatrix() << endl;
  cout << "m before permutation " << endl << m << endl;
  #endif

  // Permute the rows and columns of m_prime.
  m = m * mPerm;
  m = mPerm.transpose() * m;

  #ifdef DEBUG
  cout << "m after permutation " << endl << m << endl;
  #endif

  return m;
}


MatrixXd& reorder_ith_row_column_to_last(MatrixXd& m, const unsigned int i, const unsigned int p_gamma)
{
  // Construct a permutation matrix mPerm which interchanges the p_gamma-th and i-th rows/columns,
  // because at the moment, the p_gamma-th row/column of mA_prime contains the inverse row/column for
  // the i-th row/column of mX_gamma_prime.transpose() * mX_gamma_prime.
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> mPerm(p_gamma);
  // mPerm.setIdentity();

  // Get p_gamma index
  int *ind = mPerm.indices().data();
  for (unsigned int idx = 0; idx < (p_gamma - 1); idx++) {
    if (idx < i) {
      ind[idx] = idx;
    }
    else {
      ind[idx] = idx + 1;
    }
  }
  ind[p_gamma - 1] = i;

  #ifdef DEBUG
  cout << "mPerm" << endl << mPerm.toDenseMatrix() << endl;
  cout << "m before permutation " << endl << m << endl;
  #endif

  // Permute the rows and columns of m_prime.
  m = m * mPerm;
  m = mPerm.transpose() * m;

  #ifdef DEBUG
  cout << "m after permutation " << endl << m << endl;
  #endif

  return m;
}


//' Perform the rank one update on (X_gamma^T X_gamma)^{-1}
//'
//' @param[in]  mX_gamma       The previous iteration's matrix of covariates
//' @param[in]  min_idx        The minimum column number, relative to mX
//' @param[in]  i              The column number of the column we're adding to mX_gamma_prime
//' @param[in]  vx             The new column vector to be added to mX_gamma
//' @param[in]  mX_gamma_prime The current iteration's matrix of covariates
//' @param[in]  mA             The previous iteration's inverse i.e. (X_gamma^T X_gamma)^{-1}
//' @param[out] mA_prime       The current iteration's inverse i.e. (X_gamma_prime^T X_gamma_prime)^{-1}
//'                            which we are calculating using this function.
//' @return                    The new inverse (mX_gamma_prime^T mX_gamma_prime)^{-1}
MatrixXd& rank_one_update(const MatrixXd& mX_gamma, unsigned int min_idx, unsigned int i, const VectorXd& vx,
const MatrixXd& mX_gamma_prime, const MatrixXd& mA, MatrixXd& mA_prime)
{
  const unsigned int p_gamma_prime = mX_gamma_prime.cols();

  #ifdef DEBUG
  const unsigned int n = mX_gamma_prime.rows();
  cout << "mX_gamma " << mX_gamma.topRows(min(n, 6u)) << endl;
  cout << "p_gamma_prime " << p_gamma_prime << endl;
  cout << "i " << i << endl;
  cout << "vx " << vx.head(min(n, 6u)) << endl;
  cout << "mA " << mA << endl;
  #endif

  // Construct mA_prime
  VectorXd v1 = mX_gamma.transpose() * vx;
  const double b = 1 / (vx.dot(vx) - v1.dot(mA * v1));
  // b is supposed to be positive definite.
  #ifdef DEBUG
  cout << "b " << b << endl;
  #endif
  const double epsilon = 1e-4;
  if (b > epsilon) {
    // Do rank one update
    MatrixXd m1 = vx.transpose() * mX_gamma * mA;
    VectorXd v2 = -mA * mX_gamma.transpose() * vx * b;
    mA_prime << mA + b * m1.transpose() * m1, v2,
      v2.transpose(), b;

    if ((i - min_idx) < p_gamma_prime) {
      // Construct a permutation matrix mPerm which interchanges the p_gamma_prime-th and i-th rows/columns,
      // because at the moment, the p_gamma_prime-th row/column of mA_prime contains the inverse row/column for
      // the i-th row/column of mX_gamma_prime.transpose() * mX_gamma_prime.
      mA_prime = reorder_last_row_column_to_ith(mA_prime, (i - min_idx), p_gamma_prime);
    }
    #ifdef DEBUG
    // Check that mA_prime is really an inverse for mX_gamma_prime.transpose() * mX_gamma_prime
    MatrixXd identity_prime = (mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime;
    if (!identity_prime.isApprox(MatrixXd::Identity(p_gamma_prime, p_gamma_prime))) {
      // This inverse is nonsense. Do a full inversion.
      MatrixXd mA_prime_full = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();

      show_matrix_difference(cout, mA_prime, mA_prime_full);
      mA_prime = mA_prime_full;
      // cout << "(mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime " << identity_prime << endl;
    }
    #endif
  }
  else {
    // Perform full inverse
    mA_prime = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
  }

  #ifdef DEBUG
  // Check that mA_prime is really an inverse for mX_gamma_prime.transpose() * mX_gamma_prime
  MatrixXd identity_prime = (mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime;
  cout << "(mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime" << endl;
  cout << identity_prime << endl;
  #endif

  return mA_prime;
}


//' Perform the rank one downdate on (X_gamma^T X_gamma)^{-1}
//'
//' @param[in]  mX_gamma_prime The current iteration's matrix of covariates
//' @param[in]  min_idx        The minimum column number, relative to mX
//' @param[in]  i              The index of the column we are downdating.
//' @param[in/out]  mA         The previous iteration's inverse i.e. (X_gamma^T X_gamma)^{-1}
//' @param[out] mA_prime       The current iteration's inverse i.e. (X_gamma_prime^T X_gamma_prime)^{-1}
//' @return                    The new inverse (mX_gamma_prime^T mX_gamma_prime)^{-1}
MatrixXd& rank_one_downdate(const MatrixXd& mX_gamma, const MatrixXd& mX_gamma_prime,
const unsigned int min_idx, const unsigned int i,
MatrixXd& mA, MatrixXd& mA_prime)
{
  const unsigned int p_gamma_prime = mX_gamma_prime.cols();
  const unsigned int p_gamma = mX_gamma.cols();

  // Move i-th row/column to the end, and the (i+1)-th to p_gamma_prime-th rows/columns up/left by one.
  if ((i - min_idx) < p_gamma_prime) {
    // Construct a permutation matrix mPerm which interchanges the p_gamma_prime-th and i-th rows/columns,
    // because at the moment, the p_gamma_prime-th row/column of mA_prime contains the inverse row/column for
    // the i-th row/column of mX_gamma_prime.transpose() * mX_gamma_prime.
    mA = reorder_ith_row_column_to_last(mA, (i - min_idx), p_gamma);
    #ifdef DEBUG
    // Check that this is really an inverse for a re-ordered X_gamma^T X_gamma
    MatrixXd mXTX = mX_gamma.transpose() * mX_gamma;
    mXTX = reorder_ith_row_column_to_last(mXTX, (i - min_idx), p_gamma);
    MatrixXd identity_mXTX = mXTX * mA;
    // Assert identity_mXTX.isApprox(MatrixXd::Identity(p_gamma, p_gamma));
    if (!identity_mXTX.isApprox(MatrixXd::Identity(p_gamma, p_gamma)) && NUMERIC_FIX) {
      cout << "mXTX * mA" << endl;
      cout << identity_mXTX << endl;
    }
    #endif
  }

  MatrixXd mA_11 = mA.topLeftCorner(p_gamma_prime, p_gamma_prime);
  VectorXd va_12;
  RowVectorXd va_21;
  const double a_22 = mA(p_gamma - 1, p_gamma - 1);

  va_12 = mA.col(p_gamma - 1).head(p_gamma - 1);
  va_21 = va_12.transpose();
  mA_prime = mA_11 - (va_12 * va_21) / a_22;

  #ifdef DEBUG
  // Check that mA_prime is really an inverse for mX_gamma_prime.transpose() * mX_gamma_prime
  MatrixXd identity_prime = (mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime;
  if (!identity_prime.isApprox(MatrixXd::Identity(p_gamma_prime, p_gamma_prime)) && NUMERIC_FIX) {
    cout << "(mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime" << endl;
    cout << identity_prime << endl;
    cout << "Iterative calculation of inverse is wrong, recalculating ..." << endl;
    // This inverse is nonsense. Do a full inversion.
    MatrixXd mA_prime_full = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
    show_matrix_difference(cout, mA_prime, mA_prime_full);
    // TODO: Find the differences between mA_prime_full and mA_prime
    identity_prime = (mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime_full;
    mA_prime = mA_prime_full;
  }

  // Check that mA_prime is really an inverse for mX_gamma_prime.transpose() * mX_gamma_prime
  cout << "(mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime" << endl;
  cout << identity_prime << endl;
  #endif

  return mA_prime;
}


//' Calculate the correlations for every subset of the covariates in mX
//'
//' @param[in] vy            A vector of responses of length n
//' @param[in] mX            A matrix of covariates of size n x p
//' @param[in] intercept_col The column index of the intercept column, if there is one.
//' @param[in] bIntercept    Boolean whether there's an intercept or not
//' @param[in] bCentre       Boolean whether to centre vy and mX or not
//' @return                  A vector of length 2^p of the correlations
//' @export
VectorXd all_correlations_cpp(VectorXd vy, MatrixXd mX, const unsigned int intercept_col,
const unsigned int max_iterations, const bool bIntercept = false, const bool bCentre = true)
{
  // unsigned int idx;       // The main loop index
  // The number of observations
  const unsigned int n = mX.rows();
  // The number of covariates
  const unsigned int p = mX.cols();
  // Vector of correlations for all models
  VectorXd vR2_all(max_iterations);
  // MatrixXd mA;        // The inverse of (X^T X) for the previous iteration
  bool bmA_set = false;    // Whether mA has been set yet
  bool bUpdate;        // True for an update, false for a downdate
  unsigned int min_idx;    // The minimum index of covariate that's included in the model
  unsigned int diff_idx;     // The covariate which is changing
  double numerator;      // The numerator in the correlation calculation
  double R2;           // Correlation
  dbitset gamma(p);      // The model gamma
  unsigned int p_gamma_prime;  // The number of columns in the matrix mX_gamma_prime
  unsigned int p_gamma;    // The number of columns in the matrix mX_gamma
  VectorXd vx;         // The column vector for the current covariate
  vector< MatrixXd > vec_mA(p);
  vector< MatrixXd > vec_mX_gamma(p);

  // Pre-allocated memory
  for (unsigned int i = 0; i < p; i++) {
    vec_mA[i].resize(i + 1, i + 1);
    vec_mX_gamma[i].resize(n, i + 1);
  }

  if (bCentre) {
    // Centre vy
    // centre(vy);

    // Centre non-intercept columns of mX
    for (unsigned int i = 0; i < mX.cols(); i++) {
      // Skip intercept column if there is one.
      if (bIntercept && i == intercept_col)
        continue;

      vx = mX.col(i);
      centre(vx);
      mX.col(i) = vx;
    }
  }

  // Loop through models, updating and downdating mA as necessary
  #pragma omp parallel for\
    firstprivate(bmA_set, gamma, vx, vec_mX_gamma, vec_mA)\
    private(numerator, R2, min_idx, diff_idx, p_gamma_prime, p_gamma, bUpdate)\
      shared(cout, vy, mX, vR2_all)\
      default(none)
  for (unsigned int idx = 1; idx < max_iterations; idx++) {
    #ifdef DEBUG
    cout << endl << "Iteration " << idx << endl;
    #endif
    // By properties of Greycode, only one element can be different. And it's either one higher or
    // one lower.
    // Check if update or downdate, and for which variable
    // gamma = greycode(idx, p, gamma);
    gamma = greycode_change(idx, p, gamma, bUpdate, min_idx, diff_idx, p_gamma_prime);

    MatrixXd& mA = vec_mA[p_gamma - 1];
    MatrixXd& mX_gamma = vec_mX_gamma[p_gamma - 1];
    // Get mX matrix for gamma
    // MatrixXd mX_gamma_prime(n, p_gamma_prime);
    MatrixXd& mX_gamma_prime = vec_mX_gamma[p_gamma_prime - 1];
    mX_gamma_prime = get_mX_gamma(mX, gamma, mX_gamma_prime);
    // MatrixXd mA_prime(p_gamma_prime, p_gamma_prime);
    MatrixXd& mA_prime = vec_mA[p_gamma_prime - 1];
    vx = mX.col(diff_idx);

    // If we haven't previously calculated this inverse, calculate it the first time.
    if (!bmA_set) {
      // Calculate full inverse mA, O(p^3)
      mA_prime = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
                 // [y^T X, y^T x]^T
      VectorXd v1(p_gamma_prime);
      v1 = vy.transpose() * mX_gamma_prime;
      numerator = (v1.transpose() * mA_prime * v1).value();
      R2 = numerator / vy.squaredNorm();
      #ifdef DEBUG
      cout << "Numerator " << numerator << " denominator " << vy.squaredNorm();
      cout << " R2 " << R2 << endl;
      #endif
      vR2_all(idx) = R2;

      bmA_set = true;
    }
    else {
      if (bUpdate) {
        // Rank one update of mA_prime from mA
        #ifdef DEBUG
        cout << "Updating " << diff_idx << endl;
        cout << "p_gamma_prime " << p_gamma_prime << endl;
        #endif

        // Calculate [y^T X, y^T x]^T
        VectorXd v1(p_gamma_prime);
        // VectorXd v2, v3;
        // v2 = vy.transpose() * vx;
        // v3 = vy.transpose() * mX_gamma;
        // v1 << v2, v3;
        v1 = vy.transpose() * mX_gamma_prime;

        mA_prime = rank_one_update(mX_gamma, min_idx, diff_idx, vx, mX_gamma_prime, mA, mA_prime);

        numerator = (v1.transpose() * mA_prime * v1).value();
        R2 = numerator / vy.squaredNorm();
        vR2_all(idx) = R2;
        #ifdef DEBUG
        cout << "v1 " << v1 << endl;
        cout << "mA_prime " << mA_prime << endl;
        cout << "Numerator " << numerator << " denominator " << vy.squaredNorm();
        cout << " R2 " << R2 << endl;
        cout << "vR2_all(" << idx << ") = " << R2 << endl;
        #endif
      }
      else {
        // Rank one downdate
        #ifdef DEBUG
        cout << "Downdating " << diff_idx << endl;
        #endif
        MatrixXd& mA_prime = vec_mA[p_gamma_prime - 1];
        mA_prime = rank_one_downdate(mX_gamma, mX_gamma_prime, min_idx, diff_idx, mA, mA_prime);

        // Calculate correlation
        numerator = (vy.transpose() * mX_gamma_prime * mA_prime * mX_gamma_prime.transpose() * vy).value();
        R2 = numerator / vy.squaredNorm();
        vR2_all(idx) = R2;
        #ifdef DEBUG
        cout << "Numerator " << numerator << " denominator " << vy.squaredNorm();
        cout << " R2 " << R2 << endl;
        cout << "vR2_all(" << idx << ") = " << R2 << endl;
        #endif
      }
    }

    p_gamma = p_gamma_prime;
  }
  return vR2_all;
}

// [[Rcpp::export]]
VectorXd all_correlations(MapVecd vy, MapMatd mX, unsigned int intercept_col,
                          bool bIntercept, bool bCentre)
{
  VectorXd vR2_all(max_iterations);
  const unsigned int p = mX.cols();
  const unsigned int max_iterations = 1 << p;

  vR2_all = all_correlations_cpp(vy, mX, intercept_col, max_iterations, bIntercept, bCentre);

  return vR2_all;
}
