// correlation.cpp

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
			xmat(rc,j) = int(round(strtod(matrows[i][j].c_str(), NULL)));
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
	for (mask = num >> 1; mask != 0; mask = mask >> 1)
	{
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

MatrixXd greycode(unsigned int p)
{
	unsigned int rows = 1 << p;
	MatrixXd result(rows, p);
	for (unsigned int i = 0; i < rows; i++) {
		result.row(i) = binary_to_vec(binary_to_grey(i), p).transpose();
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
	bs = (bs >> 1) ^ bs;
	bs_ret = bs;

	return bs_ret;
}

void greycode_change(const unsigned int idx, const unsigned int p, bool& update, unsigned int& diff_idx,
										 bool bDebug)
{
	dbitset bs_curr(p);
	dbitset bs_prev(p);

	bs_curr = greycode(idx, p, bs_curr);
	bs_prev = greycode(idx - 1, p, bs_prev);
	if (bDebug) {
		cout << "Current gamma:  " << bs_curr << endl;
		cout << "Previous gamma: " << bs_prev << endl;
	}

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
	} else {
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
				os << "Row " << i << ", column " << j << " " << m1(i, j) << " " << m2(i, j) << " ";
				os << m1(i, j) - m2(i, j) << " " << (m1(i, j) - m2(i, j)) / m1(i, j) << endl;
			}
		}
	}
}

//' Calculate the correlations for every subset of the covariates in mX
//'
//' @param vy            A vector of responses of length n
//' @param mX            A matrix of covariates of size n x p
//' @param intercept_col The column index of the intercept column, if there is one.
//' @param bIntercept    Boolean whether there's an intercept or not
//' @param bCentre       Boolean whether to centre vy and mX or not
//' @param bDebug        Boolean whether to issue debugging prints or not
//' @return              A vector of length 2^p of the correlations
//' @export
// [[Rcpp::export]]
VectorXd all_correlations(VectorXd vy, MatrixXd mX, unsigned int intercept_col, bool bIntercept = false, 
													bool bCentre = true, bool bDebug = false)
{
	const unsigned int p = mX.cols();            // The number of covariates
	const unsigned int greycode_rows = (1 << p); // The number of greycode combinations, 2^p
	VectorXd vR2_all(greycode_rows);             // Vector of correlations for all models
	MatrixXd mA;                                 // The inverse of (X^T X) for the previous iteration
	bool bmA_set = false;                        // Whether mA has been set yet
	bool bUpdate;																 // True for an update, false for a downdate
	unsigned int diff_idx;                       // The covariate which is changing
	double numerator;														 // The numerator in the correlation calculation
	double R2;                                	 // Correlation
	dbitset gamma;											 				 // The model gamma
	MatrixXd mX_gamma;													 // The matrix of covariates for the previous gamma
	MatrixXd mX_gamma_prime;									 	 // The matrix of covariates for the current gamma
	unsigned int p_gamma;												 // The number of columns in the matrix mX_gamma
	VectorXd vx;                                 // The column vector for the current covariate
	
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

	// Loop through models, updating and downdating mA as necessary
	for (unsigned int idx = 1; idx < greycode_rows; idx++) {
		if (bDebug) cout << "Iteration " << idx << endl;
		// By properties of Greycode, only one element can be different. And it's either one higher or
		// one lower.
		// Check if update or downdate, and for which variable
		gamma = greycode(idx, p, gamma);
		greycode_change(idx, p, bUpdate, diff_idx, bDebug);

		// Get mX matrix for gamma
		mX_gamma_prime = get_mX_gamma(mX, gamma, mX_gamma_prime);
		p_gamma = mX_gamma_prime.cols();
		vx = mX.col(diff_idx);
		
		// If we haven't previously calculated this inverse, calculate it the first time.
		if (!bmA_set) {
			// Calculate full inverse mA, O(p^3)
			MatrixXd mA_prime(p_gamma, p_gamma);
			mA_prime = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
			VectorXd v1(p_gamma); // [y^T X, y^T x]^T
			v1 << vy.transpose() * mX_gamma_prime;
			if (bDebug) {
				cout << v1.size() << endl;
				cout << mA_prime.cols() << endl;
			}
			numerator = (v1 * mA_prime * v1.transpose()).value();
			R2 = numerator / vy.squaredNorm();
			vR2_all(idx) = R2;
			mA = mA_prime;

			bmA_set = true;
		} else {
			if (bUpdate) {
				// Rank one update
				// Construct mA
				if (bDebug) {
					cout << "Updating " << diff_idx << endl;
	
					cout << "vy.size() " << vy.size() << endl;
					cout << "vx.size() " << vx.size() << endl;
					cout << "mX_gamma.cols() " << mX_gamma.cols() << endl;
					cout << "mX_gamma.rows() " << mX_gamma.rows() << endl;
					cout << "mX_gamma_prime.cols() " << mX_gamma_prime.cols() << endl;
					cout << "mX_gamma_prime.rows() " << mX_gamma_prime.rows() << endl;
					cout << "p_gamma " << p_gamma << endl;
					cout << "mA.cols() " << mA.cols() << endl;
					cout << "mA.rows() " << mA.rows() << endl;
				}
				VectorXd v1(p_gamma); // [y^T X, y^T x]^T
				VectorXd v2, v3;
				v2 = vy.transpose() * mX_gamma;
				v3 = vy.transpose() * vx;
				if (bDebug) {
					cout << "v1.size() " << v1.size() << endl;
					cout << "v2.size() " << v2.size() << endl;
					cout << "v3.size() " << v3.size() << endl;
				}
				v1 << v2, v3;
				MatrixXd mA_prime(p_gamma, p_gamma);
				const double b = 1 / (vx.transpose() * vx - vx.transpose() * mX_gamma * mA * mX_gamma.transpose() * vx).value();
				// b is supposed to be positive definite.
				if (bDebug) {
					cout << "b " << b << endl;
				}
				const double epsilon = 1e-4;
				if (b > epsilon) {
					// Do rank one update
					mA_prime << mA + b * mA * mX_gamma.transpose() * vx * vx.transpose() * mX_gamma * mA, -mA * mX_gamma.transpose() * vx * b,
											-b * vx.transpose() * mX_gamma * mA, b;
					// Check that mA_prime is really an inverse for mX_gamma_prime.transpose() * mX_gamma_prime
					MatrixXd identity_prime = (mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime;
					if (!identity_prime.isApprox(MatrixXd::Identity(p_gamma, p_gamma))) {
						// This inverse is nonsense. Do a full inversion.
						MatrixXd mA_prime_full = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
						// TODO: Find the differences between mA_prime_full and mA_prime
						if (bDebug) show_matrix_difference(cout, mA_prime, mA_prime_full);
						mA_prime = mA_prime_full;					
					}
					if (bDebug) {
						cout << "mA_prime.cols() " << mA_prime.cols() << endl;
						// cout << "(mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime " << identity_prime << endl;
					}
				} else {
					// Perform full inverse
					mA_prime = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
				}
				if (bDebug) {
					// Check that mA_prime is really an inverse for mX_gamma_prime.transpose() * mX_gamma_prime
					MatrixXd identity_prime = (mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime;
					cout << "(mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime" << endl;
					cout << identity_prime << endl;
				}

				numerator = (v1.transpose() * mA_prime * v1).value();
				R2 = numerator / vy.squaredNorm();
				if (bDebug) {
					cout << "Numerator " << numerator << " denominator " << vy.squaredNorm();
					cout << " R2 " << R2 << endl;
				}
				vR2_all(idx) = R2;

				// Save mA
				mA = mA_prime;
			} else {
				// Rank one downdate
				if (bDebug) cout << "Downdating " << diff_idx << endl;

				// Calculate correlation
				MatrixXd mA_11 = mA.topLeftCorner(p_gamma, p_gamma);
				VectorXd va_12;
				RowVectorXd va_21;
				const double a_22 = mA(p_gamma, p_gamma);
				
				// Remember that Eigen's indexing is zero-based i.e. from 0 to n - 1, so mA.col(p_gamma) is actually
				// accessing the p_gamma + 1 th column.
				va_12 = mA.col(p_gamma).head(p_gamma);
				va_21 = va_12.transpose();
				if (bDebug) {
					cout << "mA_11.cols() " << mA_11.cols() << endl;
					cout << "va_12.size() " << va_12.size() << endl;
					cout << "va_12.size() " << va_12.size() << endl;
				}
				MatrixXd mA_prime(p_gamma, p_gamma);
				mA_prime = mA_11 - (va_12 * va_21) / a_22;

				// Check that mA_prime is really an inverse for mX_gamma_prime.transpose() * mX_gamma_prime
				MatrixXd identity_prime = (mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime;
				if (!identity_prime.isApprox(MatrixXd::Identity(p_gamma, p_gamma))) {
					// This inverse is nonsense. Do a full inversion.
					MatrixXd mA_prime_full = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
					if (bDebug) show_matrix_difference(cout, mA_prime, mA_prime_full);
					// TODO: Find the differences between mA_prime_full and mA_prime
					mA_prime = mA_prime_full;					
				}

				if (bDebug) {
					cout << "mA_prime.cols() " << mA_prime.cols() << endl;
					// Check that mA_prime is really an inverse for mX_gamma_prime.transpose() * mX_gamma_prime
					MatrixXd identity_prime = (mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime;
					cout << "(mX_gamma_prime.transpose() * mX_gamma_prime) * mA_prime" << endl;
					cout << identity_prime << endl;
				}

				VectorXd numerator;
				numerator = vy.transpose() * mX_gamma_prime * mA_prime * mX_gamma_prime.transpose() * vy;;
				R2 = (numerator / vy.squaredNorm()).value();
				if (bDebug) {
					cout << "Numerator " << numerator << " denominator " << vy.squaredNorm();
					cout << " R2 " << R2 << endl;
				}
				vR2_all(idx) = R2;

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
	
	MatrixXd m1(p + m, p + m); // Symmetric
	VectorXd m2(p + m);
	m1 << mX.transpose() * mX, mX.transpose() * mZ,
				mZ.transpose() * mX, mZ.transpose() * mZ;
	m2 << mX.transpose() * vy,
				mZ.transpose() * vy;
	VectorXd R2 = (vy.transpose() * mC * m1.inverse() * m2) / vy.squaredNorm();

	return R2;
}

int main()
{
	VectorXd vy = parseCSVfile_double("vy.csv");
	MatrixXd mX = MatrixXd::Ones(263, 1);
	MatrixXd mZ = parseCSVfile_double("mX.csv");
	MatrixXd mC(263, mZ.cols() + 1);
	mC << mX, mZ;

	Eigen::initParallel();
	Eigen::setNbThreads(1);

	//VectorXd R2_one = one_correlation(vy, mX, mZ);
	// cout << R2_one << endl;

	const bool intercept = true, centre = true, debug = true;
	VectorXd vR2_all = all_correlations(vy, mC, 0, intercept, centre, debug);

	cout << "i,R2" << endl;
	for (int i = 0; i < vR2_all.size(); i++) {
		if (vR2_all(i) > 1.0 || vR2_all(i) < 0.0) {
			cout << i << ", " << vR2_all(i) << endl;
		}
	}
	
	return 0;
}
