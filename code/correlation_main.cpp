// correlation_main.cpp

#include <sys/time.h>
#include <unistd.h>
#include <iostream>
#include <Eigen/Dense>
#include "correlation.hpp"

using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::MatrixXd;
using namespace std;

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


// void check_graycode()
// {
// 	cout << "Unflipped" << endl;
// 	MatrixXd mGreycode_R = parseCSVfile_double("graycode.csv");
// 	unsigned int n = mGreycode_R.rows();
// 	unsigned int p = mGreycode_R.cols();
// 	MatrixXd mGreycode_Cpp(n, p);
// 	for (unsigned int i = 0; i < mGreycode_Cpp.rows(); i++) {
// 		mGreycode_Cpp.row(i) = gray_vec(i, p);
// 	}

// 	for (unsigned int i = 0; i < 10; i++) {
// 		cout << "R   " << i << ": " << mGreycode_R.row(i) << endl;
// 		cout << "C++ " << i << ": " << mGreycode_Cpp.row(i) << endl;
// 		cout << endl;
// 	}
// 	show_matrix_difference(cout, mGreycode_R, mGreycode_Cpp);
// }


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
	VectorXd vR2_all = all_correlations_mX_cpp(vy, mX, 0, intercept, centre);

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
	const bool intercept = false, centre = true;
	VectorXd vy = parseCSVfile_double("vy.csv");
	MatrixXd mC = parseCSVfile_double("mX.csv");
	
	MatrixXd mX = mC.leftCols(10);
	MatrixXd mZ = mC.rightCols(9);
	VectorXd vR2_all_mX_mZ = all_correlations_mX_mZ_cpp(vy, mX, mZ, 0, intercept, centre);
	cout << "i,R2" << endl;
	for (uint i = 1; i < vR2_all_mX_mZ.size(); i++) {
		// const double epsilon = 1e-8;
		// if (abs(diff) > epsilon) {
		// cout << grey_vec(i - 1, p) << " to " << grey_vec(i, p) << endl;
		cout << i << ", C++ R2 " << vR2_all_mX_mZ(i) << endl;
		// }
	}

	struct timeval start, end;
	long mtime, seconds, useconds;
	gettimeofday(&start, NULL);
	VectorXd vR2_all_mX = all_correlations_mX_cpp(vy, mC, 0, intercept, centre);
	gettimeofday(&end, NULL);

	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec;

	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;
	cout << "Elapsed time in milliseconds: " << mtime << endl;

	VectorXd vExpected_correlations = parseCSVfile_double("Hitters_exact2.csv");
	cout << "i,R2" << endl;
	for (uint i = 1; i < vR2_all_mX.size(); i++) {
		double diff = vR2_all_mX(i) - vExpected_correlations(i);
		// const double epsilon = 1e-8;
		// if (abs(diff) > epsilon) {
		// cout << grey_vec(i - 1, p) << " to " << grey_vec(i, p) << endl;
		cout << i << ", C++ R2 " << vR2_all_mX(i) << " R R2 " << vExpected_correlations(i);
		cout << " difference " << diff << endl;
		// }
	}

	return 0;
}
