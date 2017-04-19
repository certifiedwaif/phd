// correlation.cpp
#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

#include "correlation.h"
#include "graycode.h"

#include <sys/time.h>
#include <unistd.h>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/tokenizer.hpp>

using namespace boost;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::MatrixXd;
using Eigen::MatrixBase;
using Eigen::VectorXi;
using Eigen::RowVectorXi;
using Eigen::MatrixXi;
using namespace std;
using namespace Rcpp;

const bool NUMERIC_FIX = false;
// #define DEBUG

// Code copied from here: https://gist.github.com/stephenjbarr/2266900
MatrixXd parseCSVfile_double(string infilename)
{
	ifstream in(infilename.c_str());
	if (!in.is_open()) {
		throw runtime_error("Could not open " + infilename);
	}

	typedef tokenizer< escaped_list_separator<char> > Tokenizer;

	vector< string > vec;
	string line;
	vector< vector< string > > matrows;

	while (getline(in, line)) {
		Tokenizer tok(line);
		vec.assign(tok.begin(),tok.end());

		// // Print each row
		// copy(vec.begin(), vec.end(),
		//      ostream_iterator<string>(cout, "|"));
		// cout << "\n----------------------" << endl;

		matrows.push_back(vec);
	}
	in.close();

	// FIGURE OUT HOW MANY OF THE ROWS HAVE THE RIGHT NUMBER
	// OF COLUMNS
	uint Nrows = matrows.size();
	uint Ncols = matrows[0].size();
	uint Ngoodrows = 0;
	for(uint i = 0; i < Nrows; i++) {
		if(matrows[i].size() == Ncols) {
			Ngoodrows++;
		}
	}

	// TRANSFORM THE VECTOR OF ROWS INTO AN EIGEN INTEGER MATRIX
	MatrixXd xmat = MatrixXd(Ngoodrows, Ncols);
	// cout << "INPUT MATRIX: " << Nrows << "x" << Ncols << endl;

	int rc = 0;

	for(uint i = 0; i < Nrows; i++) {
		uint rowsize = matrows[i].size();

		if(rowsize != Ncols) {
			// cout << "Row " << i << " has bad column count" << endl;
			continue;
		}

		for(uint j = 0; j < Ncols; j++) {
			xmat(rc,j) = strtod(matrows[i][j].c_str(), NULL);
		}
		rc++;
	}

	return(xmat);
}


vector<uint>& get_indices_from_dbitset(const dbitset& gamma, vector<uint>& v)
{
	for (auto i = 0; i < gamma.size(); i++) {
		if (gamma[i]) {
			v.push_back(i);
		}
	}

	return v;
}


// Get columns
// Need a corresponding get rows
// And get rows and columns
template <typename Derived1, typename Derived2>
Eigen::MatrixBase<Derived2>& get_cols(const Eigen::MatrixBase<Derived1>& m1, const dbitset& gamma, Eigen::MatrixBase<Derived2>& m2)
{
	// Special case of get_rows_and_cols
	vector<uint> columns;
	columns = get_indices_from_dbitset(gamma, columns);

	for (auto i = 0; i < columns.size(); i++) {
		m2.col(i) = m1.col(columns[i]);
	}

	return m2;
}


template <typename Derived1, typename Derived2>
Eigen::MatrixBase<Derived2>& get_rows(const Eigen::MatrixBase<Derived1>& m1, const dbitset& gamma, Eigen::MatrixBase<Derived2>& m2)
{
	// Special case of get_rows_and_cols
	vector<uint> rows;
	rows = get_indices_from_dbitset(gamma, rows);
	// MatrixXd m2(rows.size(), m1.cols());

	for (auto i = 0; i < rows.size(); i++) {
		m2.row(i) = m1.row(rows[i]);
	}

	return m2;
}


template <typename Derived1, typename Derived2>
Eigen::MatrixBase<Derived2>& get_rows_and_cols(const Eigen::MatrixBase<Derived1>& m1, const dbitset& rows_bs,
const dbitset& cols_bs, Eigen::MatrixBase<Derived2>& m2)
{
	vector<uint> row_indices;
	row_indices = get_indices_from_dbitset(rows_bs, row_indices);
	vector<uint> col_indices;
	col_indices = get_indices_from_dbitset(cols_bs, col_indices);

	// Matrices are stored in column major order, so the intent is to access contiguous memory in
	// sequence by looping down each row in the inner loop.
	for (auto j = 0; j < col_indices.size(); j++) {
		for (auto i = 0; i < row_indices.size(); i++) {
			m2(i, j) = m1(row_indices[i], col_indices[j]);
		}
	}

	return m2;
}


double sd(const VectorXd& x)
{
	const auto n = x.size();
	return (x.array() - x.mean()).square().sum() / (n - 1);
}


void centre(VectorXd& v)
{
	const auto n = v.size();
	VectorXd v_bar(n);
	v_bar.fill(v.mean());
	auto sd_v = sd(v);

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
	for (auto i = 0; i < m1.rows(); i++) {
		for (auto j = 0; j < m1.cols(); j++) {
			if (abs(m1(i, j) - m2(i, j)) > epsilon) {
				os << "Row " << i << ", column " << j << " m1 " << m1(i, j) << " m2 " << m2(i, j);
				os <<  " difference " << m1(i, j) - m2(i, j);
				os << " relative difference " << (m1(i, j) - m2(i, j)) / m1(i, j) << endl;
			}
		}
	}
}


// Perform the rank one update on (X_gamma^T X_gamma)^{-1}
template <typename Derived1, typename Derived2>
Eigen::MatrixBase<Derived2>& rank_one_update(const dbitset& gamma, const uint col_abs, const uint min_idx,
	const uint fixed,
const Eigen::MatrixBase<Derived1>& mXTX, const Eigen::MatrixBase<Derived1>& mA, Eigen::MatrixBase<Derived2>& mA_prime, bool& bLow)
{
	auto p = mXTX.cols();
	auto p_gamma_prime = mA_prime.cols();
	auto p_gamma = mA.cols();

	// Construct mA_prime
	// b = 1 / (x^T x - x^T X_gamma A X_gamma^T x)
	auto xTx = mXTX(col_abs, col_abs);
	dbitset col_bs(p);
	col_bs[col_abs] = true;
	MatrixXd X_gamma_T_x(p_gamma, 1);
	#ifdef DEBUG
	Rcpp::Rcout << "gamma " << gamma << endl;
	Rcpp::Rcout << "col_bs " << col_bs << endl;
	#endif
	X_gamma_T_x = get_rows_and_cols(mXTX, gamma, col_bs, X_gamma_T_x);

	// const double b = 1 / (x^T x - x^T X_gamma A X_gamma^T x).value();
	const auto b = 1 / (xTx - (X_gamma_T_x.transpose() * mA * X_gamma_T_x).value());
	// b is supposed to be positive definite.
	#ifdef DEBUG
	Rcpp::Rcout << "b " << b << endl;
	#endif
	const auto epsilon = 1e-12;
	if (b > epsilon) {
		// Do rank one update
		// Matrix m1 = A + b A X_gamma^T x x^T X_gamma A
		// The relative column index
		auto col = col_abs - min_idx + fixed;
		MatrixXd X_gamma_x(p_gamma, p_gamma);
		MatrixXd A_X_gamma_T_x = mA * X_gamma_T_x;
		// Re-arrange.
		MatrixXd b_X_gamma_T_x_A = b * A_X_gamma_T_x;
		if (col == 0) {
			mA_prime << b, -b * A_X_gamma_T_x.transpose(),
				-b * A_X_gamma_T_x, mA + b * A_X_gamma_T_x * A_X_gamma_T_x.transpose();
		}
		else if (0 < col && col < p_gamma_prime - 1) {
			MatrixXd m1(p_gamma, p_gamma);
			m1 = mA + b * A_X_gamma_T_x * A_X_gamma_T_x.transpose();
			// mA_prime << 1, 2, 3,
			// 						4, 5, 6,
			// 						7, 8, 9;
			mA_prime.topLeftCorner(col, col) = m1.topLeftCorner(col, col);
			mA_prime.row(col).leftCols(col) = -b_X_gamma_T_x_A.topRows(col).transpose();
			mA_prime.col(col).bottomRows(p_gamma_prime - (col + 1)) = -b_X_gamma_T_x_A.bottomRows(p_gamma_prime - (col + 1));
			mA_prime(col, col) = b;
			mA_prime.bottomLeftCorner(p_gamma_prime - (col + 1), col) = m1.bottomLeftCorner(p_gamma_prime - (col + 1), col);
			mA_prime.bottomRightCorner(p_gamma_prime - (col + 1), p_gamma_prime - (col + 1)) = m1.bottomRightCorner(p_gamma_prime - (col + 1), p_gamma_prime - (col + 1));

			// Should take advantage of the symmetry of mA_prime. For now, just fill in upper triangular entries.
			#ifdef DEBUG
			Rcpp::Rcout << "mA_prime " << endl << mA_prime << endl;
			#endif
			for (auto j = 0; j < p_gamma_prime; j++) {
				for (auto i = 0; i < j; i++) {
					mA_prime(i, j) = mA_prime(j, i);
				}
			}
			#ifdef DEBUG
			Rcpp::Rcout << "mA_prime " << mA_prime << endl;
			#endif
		} else																	 // col == p_gamma_prime
		{
			mA_prime << mA + b * A_X_gamma_T_x * A_X_gamma_T_x.transpose(), -b_X_gamma_T_x_A,
				-b_X_gamma_T_x_A.transpose(), b;
		}
		bLow = false;
	}
	else {
		// Signal that a rank one update was impossible so that the calling code can perform a full inversion.
		bLow = true;
		// bLow = false;
	}

	return mA_prime;
}


// Perform the rank one downdate on (X_gamma^T X_gamma)^{-1}
template <typename Derived1, typename Derived2>
Eigen::MatrixBase<Derived2>& rank_one_downdate(const uint col_abs, const uint min_idx, const uint fixed,
const Eigen::MatrixBase<Derived1>& mA, Eigen::MatrixBase<Derived2>& mA_prime)
{
	auto p_gamma_prime = mA_prime.cols();
	auto p_gamma = mA.cols();
	// The relative column index
	auto col = col_abs - min_idx + fixed;

	// Need to deal with three cases
	if (col == 0) {
		const MatrixXd mA_11 = mA.bottomRightCorner(p_gamma_prime, p_gamma_prime);
		const VectorXd va_12 = mA.col(0).tail(p_gamma_prime);
		// const MatrixXd va_12 = mA.block(0, 0, p_gamma_prime, 1);
		auto a_22 = mA(0, 0);
		mA_prime = mA_11 - (va_12 * va_12.transpose()) / a_22;
	}
	else if (1 <= col && col <= p_gamma - 1) {
		// 1 2 3
		// 4 5 6
		// 7 8 9
		MatrixXd mA_11(p_gamma_prime, p_gamma_prime);
		VectorXd va_12(p_gamma_prime);
		mA_11.topLeftCorner(col, col) = mA.topLeftCorner(col, col);
		mA_11.bottomRows(p_gamma_prime - col).leftCols(col) = mA.bottomRows(p_gamma_prime - col).leftCols(col);
		mA_11.bottomRightCorner(p_gamma_prime - col, p_gamma_prime - col) = mA.bottomRightCorner(p_gamma_prime - col, p_gamma_prime - col);
		va_12.head(col) = mA.col(col).head(col);
		va_12.tail(p_gamma_prime - col) = mA.col(col).tail(p_gamma_prime - col);
		auto a_22 = mA(col, col);
		mA_prime = mA_11 - (va_12 * va_12.transpose()) / a_22;
	} else																		 // col == p_gamma
	{
		const MatrixXd mA_11 = mA.topLeftCorner(p_gamma_prime, p_gamma_prime);
		const VectorXd va_12 = mA.col(p_gamma - 1).head(p_gamma - 1);
		auto a_22 = mA(p_gamma - 1, p_gamma - 1);
		mA_prime = mA_11 - (va_12 * va_12.transpose()) / a_22;
	}

	// Should take advantage of the symmetry of mA_prime. For now, just fill in upper triangular entries.
	for (auto j = 0; j < p_gamma_prime; j++) {
		for (auto i = 0; i < j; i++) {
			mA_prime(i, j) = mA_prime(j, i);
		}
	}

	return mA_prime;
}


template <typename Derived1, typename Derived2>
void update_mA_prime(bool bUpdate, const dbitset& gamma,
const uint col, const uint min_idx, const uint fixed,
const Eigen::MatrixBase<Derived1>& mXTX, const Eigen::MatrixBase<Derived1>& mA, Eigen::MatrixBase<Derived2>& mA_prime,
bool& bLow)
{
	if (bUpdate) {
		// Rank one update of mA_prime from mA
		#ifdef DEBUG
		Rcpp::Rcout << "Updating " << col << endl;
		#endif
		mA_prime = rank_one_update(gamma, col, min_idx, fixed, mXTX, mA, mA_prime, bLow);
	}
	else {
		// Rank one downdate
		#ifdef DEBUG
		Rcpp::Rcout << "Downdating " << col << endl;
		#endif
		mA_prime = rank_one_downdate(col, min_idx, fixed, mA, mA_prime);
	}
}


// Calculate the correlations for every subset of the covariates in mX
List all_correlations_main(const Graycode& graycode, VectorXd vy, MatrixXd mX, const uint fixed,
	const uint intercept_col, const uint max_iterations, const bool bIntercept = false,
	const bool bCentre = true)
{
	const uint n = mX.rows();									 // The number of observations
	const uint p = mX.cols();									 // The number of covariates
	VectorXd vR2_all(max_iterations);					 // Vector of correlations for all models
	VectorXi vpgamma_all(max_iterations);			 // Vector of number of covariates included in each model
	bool bmA_set = false;											 // Whether mA has been set yet
	bool bUpdate;															 // True for an update, false for a downdate
	uint diff_idx;														 // The covariate which is changing
	uint min_idx;															 // The minimum bit which is set in gamma_prime
	dbitset gamma(p);													 // The model gamma
	dbitset gamma_prime(p);										 // The model gamma
	uint p_gamma_prime;												 // The number of columns in the matrix mX_gamma_prime
	uint p_gamma;															 // The number of columns in the matrix mX_gamma
	vector<MatrixXd> vec_mA(p);
	vector<MatrixXd> vec_mX_gamma(p);
	vector<MatrixXd> vec_m1(p);
	const MatrixXd mXTX = mX.transpose() * mX;
	const MatrixXd mXTy = mX.transpose() * vy;
	const double yTy = vy.squaredNorm();

	// Pre-allocate memory
	for (uint i = 0; i < p; i++) {
		vec_mA[i].resize(i + 1, i + 1);
		vec_mX_gamma[i].resize(n, i + 1);
		vec_m1[i].resize(i + 1, 1);
	}

	if (bCentre) {
		// Centre vy
		// centre(vy);

		// Centre non-intercept columns of mX
		for (uint i = 0; i < mX.cols(); i++) {
			// Skip intercept column if there is one.
			if (bIntercept && i == intercept_col)
				continue;

			VectorXd vcol = mX.col(i);
			centre(vcol);
			mX.col(i) = vcol;
		}
	}

	// Loop through models, updating and downdating mA as necessary
	#pragma omp parallel for\
		firstprivate(gamma, gamma_prime, bmA_set, vec_mX_gamma, vec_mA, vec_m1)\
		private(diff_idx, min_idx, p_gamma_prime, p_gamma, bUpdate)\
			shared(mX, vR2_all, graycode)\
			default(none)
	for (uint idx = 1; idx < max_iterations; idx++) {
		#ifdef DEBUG
		Rcpp::Rcout << endl << "Iteration " << idx << endl;
		#endif
		// By properties of Greycode, only one element can be different. And it's either one higher or
		// one lower.
		// Check if update or downdate, and for which variable
		gamma = gamma_prime;
		gamma_prime = graycode[idx];

		#ifdef DEBUG
		Rcpp::Rcout << "Previous gamma: " << gamma << endl;
		Rcpp::Rcout << "Current gamma:  " << gamma_prime << endl;
		#endif

		graycode.change(gamma_prime, gamma, bUpdate, diff_idx, min_idx, p_gamma_prime);

		#ifdef DEBUG
		Rcpp::Rcout << "min_idx " << min_idx << " diff_idx " << diff_idx;
		Rcpp::Rcout << " bUpdate " << (bUpdate ? "true" : "false") << endl;
		#endif

		// Get mX matrix for gamma
		MatrixXd& mA = vec_mA[p_gamma - 1];
		MatrixXd& mA_prime = vec_mA[p_gamma_prime - 1];
		MatrixXd& mX_gamma_prime = vec_mX_gamma[p_gamma_prime - 1];

		// Only needed when bmA_set is false.
		// If we haven't previously calculated this inverse, calculate it the first time.
		if (!bmA_set) {
			// Calculate full inverse mA, O(p^3)
			mX_gamma_prime = get_cols(mX, gamma_prime, mX_gamma_prime);
			mA_prime = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
			bmA_set = true;
		}
		else {
			bool bLow;
			#ifdef DEBUG
			Rcpp::Rcout << "mA_prime before update " << mA_prime << endl;
			#endif
			update_mA_prime(bUpdate, gamma,
				diff_idx, min_idx, fixed,
				mXTX,   mA, mA_prime,
				bLow);
			if (bLow) {
				mX_gamma_prime = get_cols(mX, gamma_prime, mX_gamma_prime);
				mA_prime = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
			}
			#ifdef DEBUG
			Rcpp::Rcout << "mA_prime after update " << mA_prime << endl;
			#endif

			#ifdef DEBUG
			// Check that mA_prime is really an inverse for mX_gamma_prime.transpose() * mX_gamma_prime
			mX_gamma_prime = get_cols(mX, gamma_prime, mX_gamma_prime);
			MatrixXd identity_prime = (mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime;
			if (!identity_prime.isApprox(MatrixXd::Identity(p_gamma_prime, p_gamma_prime)) && NUMERIC_FIX) {
				Rcpp::Rcout << "(mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime" << endl;
				Rcpp::Rcout << identity_prime << endl;
				Rcpp::Rcout << "Iterative calculation of inverse is wrong, recalculating ..." << endl;
				// This inverse is nonsense. Do a full inversion.
				MatrixXd mA_prime_full = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
				show_matrix_difference(Rcpp::Rcout, mA_prime, mA_prime_full);
				// TODO: Find the differences between mA_prime_full and mA_prime
				identity_prime = (mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime_full;
				mA_prime = mA_prime_full;
			}

			// Check that mA_prime is really an inverse for mX_gamma_prime.transpose() * mX_gamma_prime
			Rcpp::Rcout << "(mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime" << endl;
			Rcpp::Rcout << identity_prime << endl;
			#endif
		}

		double R2;
		double numerator;
		// Can pre-compute, using vXTy
		// VectorXd v1(p_gamma_prime);
		// v1 = vy.transpose() * mX_gamma_prime;
		MatrixXd& m1 = vec_m1[p_gamma_prime - 1];
		m1 = get_rows(mXTy, gamma_prime, m1);
		numerator = (m1.transpose() * mA_prime * m1).value();
		R2 = numerator / yTy;
		#ifdef DEBUG
		Rcpp::Rcout << "m1 " << m1 << endl;
		Rcpp::Rcout << "mA_prime " << mA_prime << endl;
		Rcpp::Rcout << "Numerator " << numerator << " denominator " << yTy;
		Rcpp::Rcout << " R2 " << R2 << endl;
		#endif
		vR2_all(idx) = R2;											 // Calculate correlation
		vpgamma_all(idx) = p_gamma_prime;

		p_gamma = p_gamma_prime;
		// FIXME: How do you do this in a thread-safe way?
		// Rcpp::checkUserInterrupt();
	}
	VectorXd vR2(max_iterations);
	VectorXi vp_gamma(max_iterations);
	for (auto i = 1; i < max_iterations; i++) {
		// if (i <= 10)
		// 	Rcpp::Rcout << i << " to " << graycode.gray_to_binary(i) << std::endl;
		vR2(i) = vR2_all(graycode.gray_to_binary(i));
		vp_gamma(i) = vpgamma_all(graycode.gray_to_binary(i));
	}
	return List::create(Named("vR2") = vR2,
											Named("vp_gamma") = vp_gamma);
}

// [[Rcpp:export]]
List all_correlations_mX_cpp(VectorXd vy, MatrixXd mX, const uint intercept_col,
const bool bIntercept, const bool bCentre)
{
	const uint p = mX.cols();
	const uint fixed = 0;
	const uint max_iterations = 1 << p;

	Graycode graycode(p);
	return all_correlations_main(graycode, vy, mX, fixed, intercept_col, max_iterations, bIntercept, bCentre);
}


// Calculate the correlations for every subset of the covariates in mX
// [[Rcpp:export]]
List all_correlations_mX_mZ_cpp(VectorXd vy, MatrixXd mX, MatrixXd mZ, const uint intercept_col,
																		const bool bIntercept, const bool bCentre)
{
	const uint n = mX.rows();
	const uint p1 = mX.cols();
	const uint p2 = mZ.cols();
	MatrixXd mC(n, p1 + p2);
	const uint max_iterations = 1 << p2;

	// mC << mX, mZ;
	mC.leftCols(p1) = mX;
	mC.rightCols(p2) = mZ;
	Graycode graycode(p1, p2);
	return all_correlations_main(graycode, vy, mC, p1, intercept_col, max_iterations, bIntercept, bCentre);
}


VectorXd one_correlation(VectorXd vy, MatrixXd mX, MatrixXd mZ)
{
	const uint n = mX.rows();
	const uint p = mX.cols();
	const uint m = mZ.cols();

	MatrixXd mC(n, p + m);
	mC << mX, mZ;

	MatrixXd m1(p + m, p + m);								 // Symmetric
	VectorXd m2(p + m);
	m1 << mX.transpose() * mX, mX.transpose() * mZ,
		mZ.transpose() * mX, mZ.transpose() * mZ;
	m2 << mX.transpose() * vy,
		mZ.transpose() * vy;
	VectorXd R2 = (vy.transpose() * mC * m1.inverse() * m2) / vy.squaredNorm();

	return R2;
}
