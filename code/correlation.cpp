// correlation.cpp

#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <boost/tokenizer.hpp>

using namespace boost;
using Eigen::VectorXd;
using Eigen::MatrixXd;
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

int main(int argc, char **argv)
{
	VectorXd vy = parseCSVfile_double("vy.csv");
	MatrixXd mX = parseCSVfile_double("mX.csv");
	MatrixXd mZ = MatrixXd::Ones(263, 1); // Sparse

	MatrixXd mC(mX.rows(), mX.cols() + mZ.cols());
	mC << mX, mZ;
	
	MatrixXd m1(mX.cols() + mZ.cols(), mX.cols() + mZ.cols()); // Symmetric
	MatrixXd m2(mX.cols() + mZ.cols(), 1);
	m1 << mX.transpose() * mX, mX.transpose() * mZ,
				mZ.transpose() * mX, mZ.transpose() * mZ;
	m2 << mX.transpose() * vy,
				mZ.transpose() * vy;
	VectorXd R2 = (vy.transpose() * mC * (m1).inverse() * m2) / vy.squaredNorm();

	cout << R2 << endl;
	return 0;
}