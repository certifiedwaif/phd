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
typedef unsigned int uint;

struct Graycode {
	Graycode(uint _p);
	Graycode(uint _fixed, uint _varying);
	const uint fixed;
	const uint varying;
	const uint size;

	uint binary_to_gray(const uint num);
	uint gray_to_binary(const uint num);
	VectorXd binary_to_vec(const uint num);
	VectorXd gray_vec(const uint i);
	MatrixXi to_MatrixXi();
	dbitset operator[](const uint idx) const;
	void change(const dbitset& gamma_prime, const dbitset& gamma,
						 	bool& update, uint& diff_idx, uint& min_idx,
							uint& bits_set) const;
};


Graycode::Graycode(uint _p) : fixed(0), varying(_p), size(fixed + varying)
{ 
}

Graycode::Graycode(uint _fixed, uint _varying) : fixed(_fixed), varying(_varying), size(fixed + varying)
{
}


uint Graycode::binary_to_gray(uint num)
{
	return (num >> 1) ^ num;
}


/*
	The purpose of this function is to convert a reflected binary
	Gray code number to a binary number.
*/
uint Graycode::gray_to_binary(uint num)
{
	uint mask;
	for (mask = num >> 1; mask != 0; mask = mask >> 1) {
		num = num ^ mask;
	}
	return num;
}


VectorXd Graycode::binary_to_vec(uint num)
{
	VectorXd result(size);
	for (uint i = 0; i < size; i++) {
		result[(size - 1) - i] = num & 1;
		num >>= 1;
	}
	return(result);
}


VectorXd Graycode::gray_vec(uint i)
{
	return binary_to_vec(binary_to_gray(i)).transpose();
}


MatrixXi Graycode::to_MatrixXi()
{
	uint rows = 1 << varying;
	MatrixXi result(rows, size);
	#pragma omp parallel for
	for (uint i = 0; i < rows; i++) {
		dbitset bs = (*this)[i];
		for (uint j = 0; j < size; j++) {
			result(i, j) = bs[j] ? 1 : 0;
		}
	}
	return(result);
}


dbitset Graycode::operator[](const uint idx) const
{
	dbitset bs_varying(varying, idx);
	bs_varying =  bs_varying ^ (bs_varying >> 1);
	if (fixed != 0) {
		dbitset bs(size);

		for (uint i = 0; i < fixed; i++) {
			bs[i] = true;
		}

		for (uint i = 0; i < varying; i++) {
			bs[i + fixed] = bs_varying[i];
		}

		return bs;
	} else {
		return bs_varying;
	}
}


void Graycode::change(const dbitset& gamma_prime, const dbitset& gamma,
													bool& update, uint& diff_idx, uint& min_idx, uint& bits_set) const
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
