// correlation.cpp

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


void greycode_change(const unsigned int idx, const unsigned int p, bool& update, 
											unsigned int& min_idx, unsigned int& diff_idx)
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
}


MatrixXd& get_mX_gamma(MatrixXd mX, dbitset gamma, MatrixXd& mX_gamma)
{
	vector<unsigned int> gamma_columns;
	unsigned int p_gamma = 0;
	const unsigned int n = mX.rows();

	// Find which columns are set in gamma, and how many
	for (unsigned int i = 0; i < gamma.size(); i++) {
		if (gamma[i]) {
			gamma_columns.push_back(i);
			p_gamma++;
		}
	}

	// Construct a matrix of size n by p_gamma, where p_gamma is the number of columns set in gamma
	MatrixXd mX_gamma_prime(n, p_gamma);
	mX_gamma = mX_gamma_prime;

	// Extract columns from mX, put into mX_gamma
	for (unsigned int i = 0; i < p_gamma; i++) {
		mX_gamma.col(i) = mX.col(gamma_columns[i]);
	}

	return mX_gamma;
}


double sd(VectorXd x)
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


void show_matrix_difference(ostream& os, const MatrixXd m1, const MatrixXd m2, const double epsilon = 1e-8)
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
	mPerm.setIdentity();
	// Get p_gamma index
	int *ind = mPerm.indices().data();
	int tmp = ind[p_gamma - 1];
	// Move every index from idx-1 to idx
	for (int idx = p_gamma - 1; idx > ((int) i) ; idx--) {
		ind[idx] = ind[idx - 1];
	}
	// Set index i to previous p_gamma index
	ind[i] = tmp;

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
	mPerm.setIdentity();
	// Get p_gamma index
	int *ind = mPerm.indices().data();
	int tmp = ind[i];
	// Move every index from idx-1 to idx
	for (int idx = i; idx < ((int)p_gamma - 1); idx++) {
		ind[idx] = ind[idx + 1];
	}
	// Set index i to previous p_gamma index
	ind[p_gamma - 1] = tmp;

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
//' @return                    The new inverse (mX_gamma_prime^T mX_gamma_prime)^{-1}
MatrixXd& rank_one_update(MatrixXd mX_gamma, unsigned int min_idx, unsigned int i, VectorXd vx,
													MatrixXd mX_gamma_prime, MatrixXd mA, MatrixXd& mA_prime)
{
	const unsigned int p_gamma = mX_gamma.cols();
	const unsigned int p_gamma_prime = mX_gamma_prime.cols();
	const unsigned int n = mX_gamma_prime.rows();
	#ifdef DEBUG
	cout << "mX_gamma " << mX_gamma.topRows(min(n, 6u)) << endl;
	cout << "p_gamma_prime " << p_gamma_prime << endl;
	cout << "i " << i << endl;
	cout << "vx " << vx.head(min(n, 6u)) << endl;
	cout << "mA " << mA << endl;
	#endif

	// Construct mA_prime
	const double b = 1 / (vx.transpose() * vx - vx.transpose() * mX_gamma * mA * mX_gamma.transpose() * vx).value();
	// b is supposed to be positive definite.
	#ifdef DEBUG
	cout << "b " << b << endl;
	#endif
	const double epsilon = 1e-4;
	if (b > epsilon) {
		// Do rank one update
		mA_prime << mA + b * mA * mX_gamma.transpose() * vx * vx.transpose() * mX_gamma * mA, -mA * mX_gamma.transpose() * vx * b,
			-b * vx.transpose() * mX_gamma * mA, b;

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
//' @param[in]  mA             The previous iteration's inverse i.e. (X_gamma^T X_gamma)^{-1}
//' @param[out] mA_prime       The current iteration's inverse i.e. (X_gamma_prime^T X_gamma_prime)^{-1}
//' @return                    The new inverse (mX_gamma_prime^T mX_gamma_prime)^{-1}
MatrixXd& rank_one_downdate(MatrixXd mX_gamma, MatrixXd mX_gamma_prime, unsigned int min_idx, unsigned int i,
														MatrixXd mA, MatrixXd& mA_prime)
{
	const unsigned int p_gamma_prime = mX_gamma_prime.cols();
	const unsigned int p_gamma = p_gamma_prime + 1;

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
// [[Rcpp::export]]
VectorXd all_correlations(VectorXd vy, MatrixXd mX, const unsigned int intercept_col,
const bool bIntercept = false, const bool bCentre = true)
{
	unsigned int idx;			 // The main loop index
	// The number of covariates
	const unsigned int p = mX.cols();
	// The number of greycode combinations, 2^p
	const unsigned int greycode_rows = (1 << p);
	// Vector of correlations for all models
	VectorXd vR2_all(greycode_rows);
	MatrixXd mA;				 // The inverse of (X^T X) for the previous iteration
	bool bmA_set = false;		 // Whether mA has been set yet
	bool bUpdate;				 // True for an update, false for a downdate
	unsigned int min_idx;		   // The minimum index of covariate that's included in the model
	unsigned int diff_idx;		 // The covariate which is changing
	double numerator;			 // The numerator in the correlation calculation
	double R2;					 // Correlation
	dbitset gamma;				 // The model gamma
	MatrixXd mX_gamma;			 // The matrix of covariates for the previous gamma
	MatrixXd mX_gamma_prime;	 // The matrix of covariates for the current gamma
	unsigned int p_gamma_prime;	 // The number of columns in the matrix mX_gamma
	VectorXd vx;				 // The column vector for the current covariate

	if (bCentre) {
		// Centre vy
		centre(vy);

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

	// Private: mA, bMa_set, bUpdate, diff_idx, numerator, R2, gamma, mX_gamma, mX_gamma_prime, p_gamma_prime, vx
	// Shared: p, greycode_rows, vR2_all, vy, mX, intercept_col, bIntercept, bCentre
	// Loop through models, updating and downdating mA as necessary
	// #pragma omp parallel for firstprivate(bmA_set, mA, gamma, mX_gamma, mX_gamma_prime, vx)\ 
	// 	private(numerator, R2, diff_idx, p_gamma_prime, bUpdate)\ 
	// 		shared(cout, vy, mX, vR2_all)\ 
	// 		default(none)
	for (idx = 1; idx < greycode_rows; idx++) {
		#ifdef DEBUG
		cout << endl << "Iteration " << idx << endl;
		#endif
		// By properties of Greycode, only one element can be different. And it's either one higher or
		// one lower.
		// Check if update or downdate, and for which variable
		gamma = greycode(idx, p, gamma);
		greycode_change(idx, p, bUpdate, min_idx, diff_idx);

		// Get mX matrix for gamma
		mX_gamma_prime = get_mX_gamma(mX, gamma, mX_gamma_prime);
		p_gamma_prime = mX_gamma_prime.cols();
		MatrixXd mA_prime(p_gamma_prime, p_gamma_prime);
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
			mA = mA_prime;

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

				// Save mA
				mA = mA_prime;
			}
			else {
				// Rank one downdate
				#ifdef DEBUG
				cout << "Downdating " << diff_idx << endl;
				#endif
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

				// Save mA
				mA = mA_prime;
			}
		}

		// Save mX_gamma
		mX_gamma = mX_gamma_prime;
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

	MatrixXd m1(p + m, p + m);	 // Symmetric
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

	mAnscombe << 10.0,  8.04, 10.0, 9.14, 10.0, 7.46, 8.0,  6.58,
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
	VectorXd vR2_all = all_correlations(vy, mX, 0, intercept, centre);

	// Test case
	VectorXd expected_correlations(8);
	expected_correlations << 0, 0.7615888, 0.83919, 0.9218939, 0.9075042, 0.666324;
}


void check_update()
{
	MatrixXd mX(5, 3);
	mX << 7, 2, 3,
		4, 5, 6,
		7, 12, 9,
		8, 11, 12,
		17, 14, 15;

	VectorXd vx(5);
	vx << 12, -7, 6, 5, -9;

	MatrixXd mA = (mX.transpose() * mX).inverse();
	MatrixXd expected_mA_prime, actual_mA_prime(4, 4);

	// Check add last column
	MatrixXd mX_last_col(5, 4);
	mX_last_col << mX, vx;
	actual_mA_prime = rank_one_update(mX, 0, 3, vx, mX_last_col, mA, actual_mA_prime);
	expected_mA_prime = (mX_last_col.transpose() * mX_last_col).inverse();
	assert(expected_mA_prime.isApprox(actual_mA_prime));

	// Check add first column
	MatrixXd mX_first_col(5, 4);
	mX_first_col << vx, mX;
	actual_mA_prime = rank_one_update(mX, 0, 0, vx, mX_first_col, mA, actual_mA_prime);
	expected_mA_prime = (mX_first_col.transpose() * mX_first_col).inverse();
	assert(expected_mA_prime.isApprox(actual_mA_prime));

	// Check add the middle column
	MatrixXd mX_middle_col(5, 4);
	mX_middle_col << mX.leftCols(1), vx, mX.rightCols(2);
	actual_mA_prime = rank_one_update(mX, 0, 1, vx, mX_middle_col, mA, actual_mA_prime);
	expected_mA_prime = (mX_middle_col.transpose() * mX_middle_col).inverse();
	assert(expected_mA_prime.isApprox(actual_mA_prime));

	// TODO: Add test cases for when min != 0
}


void check_downdate()
{
	MatrixXd mX(5, 3);
	mX << 7, 2, 3,
		4, 5, 6,
		7, 12, 9,
		8, 11, 12,
		17, 14, 15;
	MatrixXd mA = (mX.transpose() * mX).inverse();
	MatrixXd expected_mA_prime, actual_mA_prime;

	// Check removing last column
	MatrixXd mX_no_last_col = mX.leftCols(2);
	actual_mA_prime = rank_one_downdate(mX, mX_no_last_col, 0, 2, mA, actual_mA_prime);
	expected_mA_prime = (mX_no_last_col.transpose() * mX_no_last_col).inverse();
	assert(expected_mA_prime.isApprox(actual_mA_prime));

	// Check removing first column
	MatrixXd mX_no_first_col = mX.rightCols(2);
	actual_mA_prime = rank_one_downdate(mX, mX_no_first_col, 0, 0, mA, actual_mA_prime);
	expected_mA_prime = (mX_no_first_col.transpose() * mX_no_first_col).inverse();
	assert(expected_mA_prime.isApprox(actual_mA_prime));

	// Check removing the middle column
	MatrixXd mX_no_middle_col(5, 2);
	mX_no_middle_col << mX.col(0), mX.col(2);
	actual_mA_prime = rank_one_downdate(mX, mX_no_middle_col, 0, 1, mA, actual_mA_prime);
	expected_mA_prime = (mX_no_middle_col.transpose() * mX_no_middle_col).inverse();
	assert(expected_mA_prime.isApprox(actual_mA_prime));

	// TODO: Add test cases for when min != 0
}


int main_test()
{
	check_downdate();
	check_update();

	return 0;
}


int main()
{
	const bool intercept = false, centre = true;
	//VectorXd R2_one = one_correlation(vy, mX, mZ);
	// cout << R2_one << endl;

	VectorXd vy = parseCSVfile_double("vy.csv");
	MatrixXd mX = parseCSVfile_double("mX.csv");
	// unsigned int p = mX.cols();
	VectorXd vR2_all = all_correlations(vy, mX, 0, intercept, centre);

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
