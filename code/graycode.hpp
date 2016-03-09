// graycode.h

#ifndef GRAYCODE_HPP
#define GRAYCODE_HPP

#include <Eigen/Dense>

using namespace boost;
using Eigen::VectorXd;
using Eigen::RowVectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXi;
using Eigen::RowVectorXi;
using Eigen::MatrixXi;
using namespace std;

typedef dynamic_bitset<> dbitset;

struct Graycode {
	Graycode(unsigned int _p);
	Graycode(unsigned int _fixed, unsigned int _varying);
	unsigned int fixed;
	unsigned int varying;
	unsigned int size;

	unsigned int binary_to_gray(const unsigned int num);
	unsigned int gray_to_binary(const unsigned int num);
	VectorXd binary_to_vec(const unsigned int num);
	VectorXd gray_vec(const unsigned int i);
	// operator MatrixXd();
	dbitset operator[](const unsigned int idx);
	void change(const dbitset& gamma_prime, const dbitset& gamma,
						 	bool& update, unsigned int& diff_idx, unsigned int& min_idx,
							unsigned int& bits_set);
};


Graycode::Graycode(unsigned int _p)
{ 
	fixed = 0; 
	varying = _p; 
	size = fixed + varying;
}

Graycode::Graycode(unsigned int _fixed, unsigned int _varying)
{
	fixed = _fixed; 
	varying = _varying; 
	size = fixed + varying;
}


unsigned int Graycode::binary_to_gray(unsigned int num)
{
	return (num >> 1) ^ num;
}


/*
	The purpose of this function is to convert a reflected binary
	Gray code number to a binary number.
*/
unsigned int Graycode::gray_to_binary(unsigned int num)
{
	unsigned int mask;
	for (mask = num >> 1; mask != 0; mask = mask >> 1) {
		num = num ^ mask;
	}
	return num;
}


VectorXd Graycode::binary_to_vec(unsigned int num)
{
	VectorXd result(size);
	for (unsigned int i = 0; i < size; i++) {
		result[(size - 1) - i] = num & 1;
		num >>= 1;
	}
	return(result);
}


VectorXd Graycode::gray_vec(unsigned int i)
{
	return binary_to_vec(binary_to_gray(i)).transpose();
}


// Graycode::MatrixXd operator()
// {
// 	unsigned int rows = 1 << size;
// 	MatrixXd result(rows, size);
// 	for (unsigned int i = 0; i < rows; i++) {
// 		result.row(i) = gray_vec(i);
// 	}
// 	return(result);
// }


dbitset Graycode::operator[](const unsigned int idx)
{
	dbitset bs(size, idx);
	bs =  bs ^ (bs >> 1);

	return bs;
}


void Graycode::change(const dbitset& gamma_prime, const dbitset& gamma,
													bool& update, unsigned int& diff_idx, unsigned int& min_idx, unsigned int& bits_set)
{
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
}

#endif
