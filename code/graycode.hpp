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
	const unsigned int fixed;
	const unsigned int varying;
	const unsigned int size;

	unsigned int binary_to_gray(const unsigned int num);
	unsigned int gray_to_binary(const unsigned int num);
	VectorXd binary_to_vec(const unsigned int num);
	VectorXd gray_vec(const unsigned int i);
	MatrixXd to_MatrixXd();
	dbitset operator[](const unsigned int idx) const;
	void change(const dbitset& gamma_prime, const dbitset& gamma,
						 	bool& update, unsigned int& diff_idx, unsigned int& min_idx,
							unsigned int& bits_set) const;
};


Graycode::Graycode(unsigned int _p) : fixed(0), varying(_p), size(fixed + varying)
{ 
}

Graycode::Graycode(unsigned int _fixed, unsigned int _varying) : fixed(_fixed), varying(_varying), size(fixed + varying)
{
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


MatrixXd Graycode::to_MatrixXd()
{
	unsigned int rows = 1 << varying;
	MatrixXd result(rows, size);
	#pragma omp parallel for
	for (unsigned int i = 0; i < rows; i++) {
		dbitset bs = (*this)[i];
		for (unsigned int j = 0; j < size; j++) {
			result(i, j) = bs[j] ? 1 : 0;
		}
	}
	return(result);
}


dbitset Graycode::operator[](const unsigned int idx) const
{
	dbitset bs_varying(varying, idx);
	bs_varying =  bs_varying ^ (bs_varying >> 1);
	if (fixed != 0) {
		dbitset bs(size);

		for (unsigned int i = 0; i < fixed; i++) {
			bs[i] = true;
		}

		for (unsigned int i = 0; i < varying; i++) {
			bs[i + fixed] = bs_varying[i];
		}

		return bs;
	} else {
		return bs_varying;
	}
}


void Graycode::change(const dbitset& gamma_prime, const dbitset& gamma,
													bool& update, unsigned int& diff_idx, unsigned int& min_idx, unsigned int& bits_set) const
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
