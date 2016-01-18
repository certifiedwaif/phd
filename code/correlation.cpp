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
unsigned int binaryToGrey(unsigned int num)
{
  return (num >> 1) ^ num;
}

/*
        The purpose of this function is to convert a reflected binary
        Gray code number to a binary number.
*/
unsigned int greyToBinary(unsigned int num)
{
  unsigned int mask;
  for (mask = num >> 1; mask != 0; mask = mask >> 1)
  {
    num = num ^ mask;
  }
  return num;
}

VectorXd binaryToVec(unsigned int num, unsigned int p)
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
    result.row(i) = binaryToVec(binaryToGrey(i), p).transpose();
  }
  return(result);
}

VectorXd all_correlations(VectorXd vy, MatrixXd mX, MatrixXd mZ)
{
	const unsigned int n = mX.rows();
	const unsigned int p = mX.cols();
  const unsigned int m = mZ.cols();

	// Generate greycode matrix
	MatrixXd mGrey = greycode(m);
	VectorXd vR2_all(mGrey.rows());
	MatrixXd mA;
	// Loop through models, updating and downdating m1_inverse as necessary
	for (unsigned int row = 1; row < mGrey.rows(); row++) {
		// Construct mZ_gamma
		RowVectorXd vGreycodeRow = mGrey.row(row);
    // cout << vGreycodeRow << endl;
		unsigned int one_count = vGreycodeRow.sum();

		MatrixXd mZ_gamma(n, one_count);
		unsigned int mZ_gamma_col = 0;
		for (unsigned int mZ_col = 0; mZ_col < m; mZ_col++) {
			if (vGreycodeRow(mZ_col) == 1) {
				mZ_gamma.col(mZ_gamma_col) = mZ.col(mZ_col);
				mZ_gamma_col++;
			}
		}
		
    const unsigned int m_gamma = mZ_gamma.cols();
    MatrixXd m1(p + m_gamma, p + m_gamma); // Symmetric
    VectorXd m2(p + m_gamma);
		m1 << mX.transpose() * mX, mX.transpose() * mZ_gamma,
					mZ_gamma.transpose() * mX, mZ_gamma.transpose() * mZ_gamma;
		
		// Idea: Allocate matrices and so on outside of the loop. Then resize them inside.
		// MatrixXd m1_inv = m1.inverse();
		// TODO: Rank one updates and downdates
		// By properties of greycode, only one element can be different. And it's either one higher or
		// one lower.
		// Check if update or downdate, and for which variable
		RowVectorXd vDiff = mGrey.row(row) - mGrey.row(row - 1);
		unsigned int diff_idx;
		for (diff_idx = 0; vDiff(diff_idx) == 0; diff_idx++);
		MatrixXd m1_inv(p + m_gamma, p + m_gamma);
		VectorXd vx = mZ.col(diff_idx);
		mA = (mX.transpose() * mX).inverse(); // Re-use from last time
    
		if (vDiff(diff_idx) == 1) {
			// Update
			const double b = 1 / (vx.transpose() * vx - vx.transpose() * mX * mA * mX.transpose() * vx)(0);
			m1_inv << mA + b * mA * mX * vx * vx.transpose() * mX * mA, -mA * mX.transpose() * vx * b,
								-b * vx.transpose() * mX * mA, b;
		} else {
			// Downdate
			const double c = 1 / vx.squaredNorm();
			VectorXd vb = -mX.transpose() * vx; // FIXME: I don't think this expression is right
			m1_inv = mA - vb * vb.transpose() / c;
		}
	
		m2 << mX.transpose() * vy,
					mZ_gamma.transpose() * vy;
	
    MatrixXd mC(n, p + m_gamma);
    mC << mX, mZ_gamma;
		VectorXd vR2 = (vy.transpose() * mC * m1_inv * m2) / vy.squaredNorm();
		vR2_all(row) = vR2(0);
		mA = m1_inv;
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
	VectorXd vy = parseCSVfile_double("vy.csv");
	MatrixXd mX = MatrixXd::Ones(263, 1);
  MatrixXd mZ = parseCSVfile_double("mX.csv");

	VectorXd R2_one = one_correlation(vy, mX, mZ);
	cout << R2_one << endl;

	VectorXd R2_all = all_correlations(vy, mX, mZ);
	cout << R2_all.head(10) << endl;
	
	return 0;
}