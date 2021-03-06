// correlation.cpp

// "If you want work with matrices, you should use C++ with Eigen or Armadillo. It's pretty fast." - Hadley Wickham,
// completely unprompted.

#include <sys/time.h>
#include <unistd.h>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <boost/tokenizer.hpp>
#include <boost/dynamic_bitset.hpp>

using namespace boost;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::RowVectorXi;
using Eigen::MatrixXi;
using namespace std;

typedef dynamic_bitset<> dbitset;
const bool NUMERIC_FIX = false;

// Code copied from here: https://gist.github.com/stephenjbarr/2266900
MatrixXd parseCSVfile_double(string infilename)
{
	ifstream in(infilename.c_str());
	if (!in.is_open()) return MatrixXd(1,1);

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
	unsigned int Nrows = matrows.size();
	unsigned int Ncols = matrows[0].size();
	unsigned int Ngoodrows = 0;
	for(unsigned int i = 0; i < Nrows; i++) {
		if(matrows[i].size() == Ncols) {
			Ngoodrows++;
		}
	}

	// TRANSFORM THE VECTOR OF ROWS INTO AN EIGEN INTEGER MATRIX
	MatrixXd xmat = MatrixXd(Ngoodrows, Ncols);
	// cout << "INPUT MATRIX: " << Nrows << "x" << Ncols << endl;

	int rc = 0;

	for(unsigned int i = 0; i < Nrows; i++) {
		unsigned int rowsize = matrows[i].size();

		if(rowsize != Ncols) {
			// cout << "Row " << i << " has bad column count" << endl;
			continue;
		}

		for(unsigned int j = 0; j < Ncols; j++) {
			xmat(rc,j) = strtod(matrows[i][j].c_str(), NULL);
		}
		rc++;
	}

	return(xmat);
}


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


dbitset& greycode(const unsigned int idx, const unsigned int p, dbitset& bs_ret)
{
	dbitset bs(p, idx);
	bs =  bs ^ (bs >> 1);
	bs_ret = bs;

	return bs_ret;
}


dbitset& greycode_change(const unsigned int idx, const unsigned int p, dbitset& gamma_prime, const dbitset& gamma, bool& update,
unsigned int& diff_idx, unsigned int& min_idx,
unsigned int& bits_set)
{
	gamma_prime = greycode(idx, p, gamma_prime);
	#ifdef DEBUG
	cout << "Previous gamma: " << gamma << endl;
	cout << "Current gamma:  " << gamma_prime << endl;
	#endif

	// Find the LSB.
	min_idx = min(gamma.find_first(), gamma_prime.find_first());

	// Find bit that has changed.
	diff_idx = (gamma_prime ^ gamma).find_first();

	// Has it been set, or unset?
	update = gamma_prime[diff_idx];

	bits_set = gamma_prime.count();

	return gamma_prime;
}


vector<unsigned int>& get_indices_from_dbitset(const dbitset& gamma, vector<unsigned int>& v)
{
	for (unsigned int i = 0; i < gamma.size(); i++) {
		if (gamma[i]) {
			v.push_back(i);
		}
	}

	return v;
}


// Get columns
// Need a corresponding get rows
// And get rows and columns
MatrixXd& get_cols(const MatrixXd& m1, const dbitset& gamma, MatrixXd& m2)
{
	// Special case of get_rows_and_cols
	vector<unsigned int> columns;
	columns = get_indices_from_dbitset(gamma, columns);

	for (unsigned int i = 0; i < columns.size(); i++) {
		m2.col(i) = m1.col(columns[i]);
	}

	return m2;
}


MatrixXd& get_rows(const MatrixXd& m1, const dbitset& gamma, MatrixXd& m2)
{
	// Special case of get_rows_and_cols
	vector<unsigned int> rows;
	rows = get_indices_from_dbitset(gamma, rows);
	// MatrixXd m2(rows.size(), m1.cols());

	for (unsigned int i = 0; i < rows.size(); i++) {
		m2.row(i) = m1.row(rows[i]);
	}

	return m2;
}


MatrixXd& get_rows_and_cols(const MatrixXd& m1, const dbitset& rows_bs, const dbitset& cols_bs, MatrixXd& m2)
{
	vector<unsigned int> row_indices;
	row_indices = get_indices_from_dbitset(rows_bs, row_indices);
	vector<unsigned int> col_indices;
	col_indices = get_indices_from_dbitset(cols_bs, col_indices);

	// Matrices are stored in column major order, so the intent is to access contiguous memory in
	// sequence by looping down each row in the inner loop.
	for (unsigned int j = 0; j < col_indices.size(); j++) {
		for (unsigned int i = 0; i < row_indices.size(); i++) {
			m2(i, j) = m1(row_indices[i], col_indices[j]);
		}
	}

	return m2;
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


//' Perform the rank one update on (X_gamma^T X_gamma)^{-1}
//'
//' @param[in]  gamma    The previous iteration's matrix of covariates
//' @param[in]  i              The column number of the column we're adding to mX_gamma_prime
//' @param[in]  min_idx             The new column vector to be added to mX_gamma
//' @param[in]  mXTX 		 The current iteration's matrix of covariates
//' @param[in]  mA             The previous iteration's inverse i.e. (X_gamma^T X_gamma)^{-1}
//' @param[out] mA_prime       The current iteration's inverse i.e. (X_gamma_prime^T X_gamma_prime)^{-1}
//'                            which we are calculating using this function.
//' @param[out] bLow       The current iteration's inverse i.e. (X_gamma_prime^T X_gamma_prime)^{-1}
//' @return                    The new inverse (mX_gamma_prime^T mX_gamma_prime)^{-1}
MatrixXd& rank_one_update(const dbitset& gamma, const unsigned int col_abs, const unsigned int min_idx,
const MatrixXd& mXTX, MatrixXd& mA, MatrixXd& mA_prime, bool& bLow)
{
	const unsigned int p = mXTX.cols();
	const unsigned int p_gamma_prime = mA_prime.cols();
	const unsigned int p_gamma = mA.cols();

	// Construct mA_prime
	// b = 1 / (x^T x - x^T X_gamma A X_gamma^T x)
	double xTx = mXTX(col_abs, col_abs);
	dbitset col_bs(p);
	col_bs[col_abs] = true;
	MatrixXd X_gamma_T_x(p_gamma, 1);
	#ifdef DEBUG
	cout << "gamma " << gamma << endl;
	cout << "col_bs " << col_bs << endl;
	#endif
	X_gamma_T_x = get_rows_and_cols(mXTX, gamma, col_bs, X_gamma_T_x);

	// const double b = 1 / (x^T x - x^T X_gamma A X_gamma^T x).value();
	const double b = 1 / (xTx - (X_gamma_T_x.transpose() * mA * X_gamma_T_x).value());
	// b is supposed to be positive definite.
	#ifdef DEBUG
	cout << "b " << b << endl;
	#endif
	const double epsilon = 1e-12;
	if (b > epsilon) {
		// Do rank one update
		// Matrix m1 = A + b A X_gamma^T x x^T X_gamma A
		// The relative column index
		const unsigned int col = col_abs - min_idx;
		MatrixXd X_gamma_x(p_gamma, p_gamma);
		MatrixXd A_X_gamma_T_x = mA * X_gamma_T_x;
		// Re-arrange.
		MatrixXd b_X_gamma_T_x_A = mA * X_gamma_T_x * b;
		if (col == 0) {
			mA_prime << b, -b_X_gamma_T_x_A.transpose(),
				-b_X_gamma_T_x_A, mA + b * A_X_gamma_T_x * A_X_gamma_T_x.transpose();
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
			cout << "mA_prime " << endl << mA_prime << endl;
			#endif
			for (unsigned int j = 0; j < p_gamma_prime; j++) {
				for (unsigned int i = 0; i < j; i++) {
					mA_prime(i, j) = mA_prime(j, i);
				}
			}
			#ifdef DEBUG
			cout << "mA_prime " << mA_prime << endl;
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


//' Perform the rank one downdate on (X_gamma^T X_gamma)^{-1}
//'
//' @param[in]  mX_gamma_prime The current iteration's matrix of covariates
//' @param[in]  min_idx        The minimum column number, relative to mX
//' @param[in]  i              The index of the column we are downdating.
//' @param[in/out]  mA         The previous iteration's inverse i.e. (X_gamma^T X_gamma)^{-1}
//' @param[out] mA_prime       The current iteration's inverse i.e. (X_gamma_prime^T X_gamma_prime)^{-1}
//' @return                    The new inverse (mX_gamma_prime^T mX_gamma_prime)^{-1}
MatrixXd& rank_one_downdate(const unsigned int col_abs, const unsigned int min_idx, MatrixXd& mA, MatrixXd& mA_prime)
{
	const unsigned int p_gamma_prime = mA_prime.cols();
	const unsigned int p_gamma = mA.cols();
	// The relative column index
	const unsigned int col = col_abs - min_idx;

	MatrixXd mA_11(p_gamma_prime, p_gamma_prime);
	VectorXd va_12(p_gamma_prime);
	double a_22;
	// Need to deal with three cases
	if (col == 0) {
		mA_11 = mA.bottomRightCorner(p_gamma_prime, p_gamma_prime);
		va_12 = mA.col(0).tail(p_gamma_prime);
		a_22 = mA(0, 0);
	}
	else if (1 <= col && col <= p_gamma - 1) {
		// 1 2 3
		// 4 5 6
		// 7 8 9
		mA_11.topLeftCorner(col, col) = mA.topLeftCorner(col, col);
		mA_11.bottomRows(p_gamma_prime - col).leftCols(col) = mA.bottomRows(p_gamma_prime - col).leftCols(col);
		mA_11.bottomRightCorner(p_gamma_prime - col, p_gamma_prime - col) = mA.bottomRightCorner(p_gamma_prime - col, p_gamma_prime - col);
		va_12.head(col) = mA.col(col).head(col);
		va_12.tail(p_gamma_prime - col) = mA.col(col).tail(p_gamma_prime - col);
		a_22 = mA(col, col);

	} else																		 // col == p_gamma
	{
		mA_11 = mA.topLeftCorner(p_gamma_prime, p_gamma_prime);
		va_12 = mA.col(p_gamma - 1).head(p_gamma - 1);
		a_22 = mA(p_gamma - 1, p_gamma - 1);
	}
	mA_prime = mA_11 - (va_12 * va_12.transpose()) / a_22;

	// Should take advantage of the symmetry of mA_prime. For now, just fill in upper triangular entries.
	for (unsigned int j = 0; j < p_gamma_prime; j++) {
		for (unsigned int i = 0; i < j; i++) {
			mA_prime(i, j) = mA_prime(j, i);
		}
	}

	return mA_prime;
}


void update_mA_prime(bool bUpdate, const dbitset& gamma,
const unsigned int col, const unsigned int min_idx,
const MatrixXd& mXTX, MatrixXd& mA, MatrixXd& mA_prime,
bool& bLow)
{
	if (bUpdate) {
		// Rank one update of mA_prime from mA
		#ifdef DEBUG
		cout << "Updating " << col << endl;
		#endif
		mA_prime = rank_one_update(gamma, col, min_idx, mXTX, mA, mA_prime, bLow);
	}
	else {
		// Rank one downdate
		#ifdef DEBUG
		cout << "Downdating " << col << endl;
		#endif
		mA_prime = rank_one_downdate(col, min_idx, mA, mA_prime);
	}
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
// [[Rcpp::export]]
VectorXd all_correlations_mX_mZ(VectorXd vy, MatrixXd mX, MatrixXd mZ, const unsigned int intercept_col,
const unsigned int max_iterations, const bool bIntercept = false, const bool bCentre = true)
{
	VectorXd vR2_all(max_iterations);
	return vR2_all;
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
// [[Rcpp::export]]
VectorXd all_correlations_mX(VectorXd vy, MatrixXd mX, const unsigned int intercept_col,
const unsigned int max_iterations, const bool bIntercept = false, const bool bCentre = true)
{
	const unsigned int n = mX.rows();					 // The number of observations
	const unsigned int p = mX.cols();					 // The number of covariates
	VectorXd vR2_all(max_iterations);					 // Vector of correlations for all models
	bool bmA_set = false;											 // Whether mA has been set yet
	bool bUpdate;															 // True for an update, false for a downdate
	unsigned int diff_idx;										 // The covariate which is changing
	unsigned int min_idx;											 // The minimum bit which is set in gamma_prime
	dbitset gamma(p);													 // The model gamma
	dbitset gamma_prime(p);										 // The model gamma
	unsigned int p_gamma_prime;								 // The number of columns in the matrix mX_gamma_prime
	unsigned int p_gamma;											 // The number of columns in the matrix mX_gamma
	vector< MatrixXd > vec_mA(p);
	vector< MatrixXd > vec_mX_gamma(p);
	vector< MatrixXd > vec_m1(p);
	const MatrixXd mXTX = mX.transpose() * mX;
	const MatrixXd mXTy = mX.transpose() * vy;
	const double yTy = vy.squaredNorm();

	// Pre-allocate memory
	for (unsigned int i = 0; i < p; i++) {
		vec_mA[i].resize(i + 1, i + 1);
		vec_mX_gamma[i].resize(n, i + 1);
		vec_m1[i].resize(i + 1, 1);
	}

	if (bCentre) {
		// Centre vy
		// centre(vy);

		// Centre non-intercept columns of mX
		for (unsigned int i = 0; i < mX.cols(); i++) {
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
		firstprivate(bmA_set, gamma, gamma_prime, vec_mX_gamma, vec_mA, vec_m1)\
		private(diff_idx, min_idx, p_gamma_prime, p_gamma, bUpdate)\
			shared(cout, mX, vR2_all)\
			default(none)
	for (unsigned int idx = 1; idx < max_iterations; idx++) {
		#ifdef DEBUG
		cout << endl << "Iteration " << idx << endl;
		#endif
		// By properties of Greycode, only one element can be different. And it's either one higher or
		// one lower.
		// Check if update or downdate, and for which variable
		gamma = gamma_prime;
		gamma_prime = greycode_change(idx, p, gamma_prime, gamma, bUpdate, diff_idx, min_idx, p_gamma_prime);

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
			update_mA_prime(bUpdate, gamma,
				diff_idx, min_idx,
				mXTX,   mA, mA_prime,
				bLow);
			if (bLow) {
				mX_gamma_prime = get_cols(mX, gamma_prime, mX_gamma_prime);
				mA_prime = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
			}

			#ifdef DEBUG
			// Check that mA_prime is really an inverse for mX_gamma_prime.transpose() * mX_gamma_prime
			mX_gamma_prime = get_cols(mX, gamma_prime, mX_gamma_prime);
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
		cout << "m1 " << m1 << endl;
		cout << "mA_prime " << mA_prime << endl;
		cout << "Numerator " << numerator << " denominator " << yTy;
		cout << " R2 " << R2 << endl;
		#endif
		vR2_all(idx) = R2;											 // Calculate correlation

		p_gamma = p_gamma_prime;
	}
	return vR2_all;
}


VectorXd one_correlation(VectorXd vy, MatrixXd mX, MatrixXd mZ)
{
	const unsigned int n = mX.rows();
	const unsigned int p = mX.cols();
	const unsigned int m = mZ.cols();

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


MatrixXd anscombe()
{
	Eigen::Matrix<double, 11, 8, Eigen::RowMajor> mAnscombe;

	mAnscombe <<
		10.0,  8.04, 10.0, 9.14, 10.0, 7.46, 8.0,  6.58,
		8.0,  6.95, 8.0,  8.14, 8.0,  6.77, 8.0,  5.76,
		13.0, 7.58, 13.0, 8.74, 13.0, 12.74,  8.0,  7.71,
		9.0,  8.81, 9.0,  8.77, 9.0,  7.11, 8.0,  8.84,
		11.0, 8.33, 11.0, 9.26, 11.0, 7.81, 8.0,  8.47,
		14.0, 9.96, 14.0, 8.10, 14.0, 8.84, 8.0,  7.04,
		6.0,  7.24, 6.0,  6.13, 6.0,  6.08, 8.0,  5.25,
		4.0,  4.26, 4.0,  3.10, 4.0,  5.39, 19.0, 12.50,
		12.0, 10.84,  12.0, 9.13, 12.0, 8.15, 8.0,  5.56,
		7.0,  4.82, 7.0,  7.26, 7.0,  6.42, 8.0,  7.91,
		5.0,  5.68, 5.0,  4.74, 5.0,  5.73, 8.0,  6.89;

	return mAnscombe;
}


void check_greycode()
{
	cout << "Unflipped" << endl;
	MatrixXd mGreycode_R = parseCSVfile_double("greycode.csv");
	unsigned int n = mGreycode_R.rows();
	unsigned int p = mGreycode_R.cols();
	MatrixXd mGreycode_Cpp(n, p);
	for (unsigned int i = 0; i < mGreycode_Cpp.rows(); i++) {
		mGreycode_Cpp.row(i) = grey_vec(i, p);
	}

	for (unsigned int i = 0; i < 10; i++) {
		cout << "R   " << i << ": " << mGreycode_R.row(i) << endl;
		cout << "C++ " << i << ": " << mGreycode_Cpp.row(i) << endl;
		cout << endl;
	}
	show_matrix_difference(cout, mGreycode_R, mGreycode_Cpp);
}


void check_anscombe()
{
	// Simpler test case - Anscombe's quartet
	const bool intercept = false, centre = true;
	MatrixXd mAnscombe = anscombe();

	VectorXd vy = mAnscombe.col(0);
	MatrixXd mX = mAnscombe.middleCols(1, 3);
	#ifdef DEBUG
	cout << mAnscombe << endl;
	cout << vy << endl;
	cout << mX << endl;
	#endif
	VectorXd vR2_all = all_correlations_mX(vy, mX, 0, intercept, centre);

	// Test case
	VectorXd expected_correlations(8);
	expected_correlations << 0, 0.7615888, 0.83919, 0.9218939, 0.9075042, 0.666324;
}


void check_update()
{
	// First column
	// Last column
	// In the middle
}


void check_downdate()
{
	// First column
	// Last column
	// In the middle
}


int main_test()
{
	check_downdate();
	check_update();

	return 0;
}


int main()
{
	//VectorXd R2_one = one_correlation(vy, mX, mZ);
	// cout << R2_one << endl;

	// cout << greycode(19) << endl;

	// return 0;

	const bool intercept = false, centre = true;

	VectorXd vy = parseCSVfile_double("vy.csv");
	MatrixXd mX = parseCSVfile_double("mX.csv");
	const unsigned int p = mX.cols();
	// The number of greycode combinations, 2^p
	const unsigned int greycode_rows = (1 << p);
	const unsigned int max_iterations = greycode_rows;
	// const unsigned int max_iterations = 100000;
	// unsigned int p = mX.cols();

	struct timeval start, end;
	long mtime, seconds, useconds;
	gettimeofday(&start, NULL);
	VectorXd vR2_all = all_correlations_mX(vy, mX, 0, max_iterations, intercept, centre);
	gettimeofday(&end, NULL);

	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;

	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
	cout << "Elapsed time in milliseconds: " << mtime << endl;

	VectorXd vExpected_correlations = parseCSVfile_double("Hitters_exact2.csv");
	cout << "i,R2" << endl;
	for (int i = 1; i < vR2_all.size(); i++) {
		double diff = vR2_all(i) - vExpected_correlations(i);
		// const double epsilon = 1e-8;
		// if (abs(diff) > epsilon) {
		// cout << grey_vec(i - 1, p) << " to " << grey_vec(i, p) << endl;
		cout << i << ", C++ R2 " << vR2_all(i) << " R R2 " << vExpected_correlations(i);
		cout << " difference " << diff << endl;
		// }
	}

	return 0;
}
