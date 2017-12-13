// graycode.h

// [[Rcpp::depends(BH)]]

#pragma once

#include <Rcpp.h>
#include <boost/dynamic_bitset.hpp>
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

	uint binary_to_gray(const uint num) const;
	uint gray_to_binary(const uint num) const;
	VectorXd binary_to_vec(const uint num);
	VectorXd gray_vec(const uint i);
	MatrixXi to_MatrixXi() const;
	dbitset operator[](const uint idx) const;
	void change(const dbitset& gamma_prime, const dbitset& gamma,
						 	bool& update, uint& diff_idx, uint& min_idx,
							uint& bits_set) const;
};
