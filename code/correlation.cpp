// correlation.cpp

#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <boost/tokenizer.hpp>

using namespace boost;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::RowVectorXi;
using Eigen::MatrixXi;
using namespace std;

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
  int Nrows = matrows.size();
  int Ncols = matrows[0].size();
  int Ngoodrows = 0;
  for(int i = 0; i < Nrows; i++) {
    if(matrows[i].size() == Ncols) {
			Ngoodrows++;
    }
  }

  // TRANSFORM THE VECTOR OF ROWS INTO AN EIGEN INTEGER MATRIX
  MatrixXd xmat = MatrixXd(Ngoodrows, Ncols);
  cout << "INPUT MATRIX: " << Nrows << "x" << Ncols << endl;

  int rc = 0;

  for(int i = 0; i < Nrows; i++) {
    int rowsize = matrows[i].size();

    if(rowsize != Ncols) {
			cout << "Row " << i << " has bad column count" << endl;
			continue;
    } 

    for(int j = 0; j < Ncols; j++) {
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

VectorXd make_ve(const unsigned int size, const unsigned int p)
{
	VectorXd ve(size);
	ve.setZero(size);
	ve(p - 1) = 1.; // ve_{p}
	
	return ve;
}

VectorXd all_correlations_sherman_morrison(VectorXd vy, MatrixXd mX, MatrixXd mZ)
{
	const unsigned int n = mX.rows();
	const unsigned int p = mX.cols();
	const unsigned int m = mZ.cols();
	const unsigned int greycode_rows = (1 << m);

	// Generate greycode matrix
	// MatrixXd mGrey = greycode(m);
	VectorXd vR2_all(greycode_rows);
	MatrixXd m1_inv_last;
	MatrixXd mC_last;
	
	bool bm1_inv_last_set = false;
	// Loop through models, updating and downdating m1_inverse as necessary
	for (unsigned int row = 1; row < greycode_rows; row++) {
		// Construct mZ_gamma
		RowVectorXd v_greycode_row = binary_to_vec(binary_to_grey(row), m);
		unsigned int one_count = v_greycode_row.sum();
		MatrixXd mZ_gamma(n, one_count);
    // cout << v_greycode_row << endl;
		
		// MatrixXd m1_inv = m1.inverse();
		// TODO: Rank one updates and downdates
		// By properties of greycode, only one element can be different. And it's either one higher or
		// one lower.
		// Check if update or downdate, and for which variable
		RowVectorXd v_diff = binary_to_vec(binary_to_grey(row), m) - binary_to_vec(binary_to_grey(row - 1), m);
		unsigned int diff_idx;
		for (diff_idx = 0; v_diff(diff_idx) == 0; diff_idx++);
		VectorXd vz = mZ.col(diff_idx);
		unsigned int mZ_gamma_col = 0;
		for (unsigned int mZ_col_idx = 0; mZ_col_idx < m; mZ_col_idx++) {
			if (v_greycode_row(mZ_col_idx) == 1.) {
				mZ_gamma.col(mZ_gamma_col) = mZ.col(mZ_col_idx);
				mZ_gamma_col++;
			}
		}

    const unsigned int m_gamma = mZ_gamma.cols();
    MatrixXd m1(p + m_gamma, p + m_gamma); // Symmetric
    VectorXd m2(p + m_gamma);
		m1 << mX.transpose() * mX, mX.transpose() * mZ_gamma,
					mZ_gamma.transpose() * mX, mZ_gamma.transpose() * mZ_gamma;

		m2 << mX.transpose() * vy,
					mZ_gamma.transpose() * vy;
	
    MatrixXd mC(n, p + m_gamma);
    mC << mX, mZ_gamma;
    VectorXd vR2;
		// If we haven't previously calculated this inverse, calculate it the first time.
		if (!bm1_inv_last_set) {
			m1_inv_last = m1.inverse();
			vR2 = (vy.transpose() * mC * m1_inv_last * m2) / vy.squaredNorm();
			vR2_all(row) = vR2.value();
			bm1_inv_last_set = true;
		} else {
			if (v_diff(diff_idx) == 1.) {
				// Update
				// Construct m1_inv
				MatrixXd m1_inv(m1_inv_last.rows() + 1, m1_inv_last.cols() + 1);
				cout << "Updating " << diff_idx << endl;
				cout << "m1_inv_last.cols() " << m1_inv_last.cols() << endl;
				cout << "mC_last.cols() " << mC_last.cols() << endl;
				cout << "mC.cols() " << mC.cols() << endl;

				// const double b = 1 / (vx.transpose() * vx - vx.transpose() * mX * m1_inv_last * mX.transpose() * vx)(0);
				// m1_inv << m1_inv_last + b * m1_inv_last * mX * vx * vx.transpose() * mX * m1_inv_last, -m1_inv_last * mX.transpose() * vx * b,
				// 					-b * vx.transpose() * mX * m1_inv_last, b;
				m1_inv << m1_inv_last, VectorXd::Zero(m1_inv_last.rows()),
									VectorXd::Zero(m1_inv_last.rows()).transpose(), 1.;

				// Figure out update outer products
				VectorXd ve = make_ve(m1_inv.rows(), m1_inv.rows());

				VectorXd vv(m1_inv.rows()); // Form the vector [C^T z, 0]^T
				vv << mC_last.transpose() * vz,
							0.;

				RowVectorXd vv2(m1_inv.rows()); // Form the vector [C^T z, z^T z]
				vv2 << (mC_last.transpose() * vz).transpose(), vz.squaredNorm();
				
				// Idea: Make m1_inv symmetric. Then just do one Sherman-Morrison update.
				// Sherman-Morrison
				m1_inv = sherman_morrison(m1_inv, vv, ve);
				m1_inv = sherman_morrison(m1_inv, ve, vv2);
				
				// Calculate correlation
				vR2 = (vy.transpose() * mC * m1_inv * m2) / vy.squaredNorm();
				vR2_all(row) = vR2.value();

				// Save m1_inv
				m1_inv_last = m1_inv;
			} else {
				// Downdate
				cout << "Downdating " << diff_idx << endl;
				// const double c = 1 / vx.squaredNorm();
				// VectorXd vb = -mX.transpose() * vx; // FIXME: I don't think this expression is right
				// m1_inv = m1_inv_last - vb * vb.transpose() / c;
				// Construct m1_inv
				MatrixXd m1_inv = m1_inv_last;
				
				// Figure out downdate outer products
				VectorXd ve = make_ve(m1_inv.rows(), m1_inv.rows());
				
				VectorXd vv(m1_inv.rows()); // Form the vector [C^T z, 0]^T
				vv << mC_last.transpose() * vz,
							0.;

				RowVectorXd vv2(m1_inv.rows()); // Form the vector [C^T z, z^T z]
				vv2 << (mC_last.transpose() * vz).transpose(), vz.squaredNorm();

				// Sherman-Morrison
				m1_inv = sherman_morrison(m1_inv, -vv, ve);
				m1_inv = sherman_morrison(m1_inv, ve, -vv2);

				// Take the upper left block of the resulting m1_inv to be our new inverse.
				m1_inv = m1_inv.topLeftCorner(m1_inv.rows() - 1, m1_inv.cols() - 1);

				// Calculate correlation
				vR2 = (vy.transpose() * mC * m1_inv * m2) / vy.squaredNorm();
				vR2_all(row) = vR2.value();

				// Save m1_inv
				m1_inv_last = m1_inv;
			}
		}
		mC_last = mC;
	}

	return vR2_all;
}

VectorXd all_correlations_block_inverse(VectorXd vy, MatrixXd mX, MatrixXd mZ)
{
	const unsigned int n = mX.rows();
	const unsigned int p = mX.cols();
	const unsigned int m = mZ.cols();
	const unsigned int greycode_rows = (1 << m);
	VectorXd vR2_all(greycode_rows);
	bool bmA_set = false;
	
	// Loop through models, updating and downdating m1_inverse as necessary
	for (unsigned int row = 1; row < greycode_rows; row++) {
		// By properties of Greycode, only one element can be different. And it's either one higher or
		// one lower.
		// Check if update or downdate, and for which variable
		VectorXd v_diff(p + m);
    VectorXd vR2;
    unsigned int diff_idx; // The covariate which is changing
    VectorXd vx = mX.col(diff_idx);

		// If we haven't previously calculated this inverse, calculate it the first time.
		if (!bmA_set) {
			// Calculate full inverse
			bmA_set = true;
		} else {
			if (v_diff(diff_idx) == 1.) {
				// Rank one update
				// Construct mA
				cout << "Updating " << diff_idx << endl;

				const double b = 1 / (vx.transpose() * vx - vx.transpose() * mX * m1_inv_last * mX.transpose() * vx).value();
				vR2 = (vy.transpose() * mC * m1_inv * m2) / vy.squaredNorm();
				vR2_all(row) = vR2.value();

				// Save mA
			} else {
				// Rank one downdate
				cout << "Downdating " << diff_idx << endl;
				// Construct m1_inv
				MatrixXd m1_inv = m1_inv_last;
				

				// Calculate correlation
				vR2 = ;
				vR2_all(row) = vR2.value();

				// Save mA
				
			}
		}
		mC_last = mC;
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

int main(int argc, char **argv)
{
	// Eigen::initParallel();
	// Eigen::setNbThreads(4);

	VectorXd vy = parseCSVfile_double("vy.csv");
	MatrixXd mX = MatrixXd::Ones(263, 1);
	MatrixXd mZ = parseCSVfile_double("mX.csv");

	VectorXd R2_one = one_correlation(vy, mX, mZ);
	cout << R2_one << endl;

	VectorXd R2_all = all_correlations(vy, mX, mZ);
	cout << R2_all.head(10) << endl;
	
	return 0;
}