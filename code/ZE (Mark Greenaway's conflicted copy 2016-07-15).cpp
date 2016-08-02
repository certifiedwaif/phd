/*
 * main.cpp
 *
 *  Created on: 20/07/2015
 *      Author: markg
 */

// [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
#include <iostream>
#include <list>
#include <utility>
#include <armadillo>

// using namespace Rcpp;
using namespace std;
using namespace arma;

typedef uint uint;

// From the Wikipedia page on Gray code
/*
        The purpose of this function is to convert an unsigned
        binary number to reflected binary Gray code.

        The operator >> is shift right. The operator ^ is exclusive or.
*/
uint binary_to_gray(uint num)
{
  return (num >> 1) ^ num;
}

/*
        The purpose of this function is to convert a reflected binary
        Gray code number to a binary number.
*/
uint gray_to_binary(uint num)
{
  uint mask;
  for (mask = num >> 1; mask != 0; mask = mask >> 1)
  {
    num = num ^ mask;
  }
  return num;
}

vec binary_to_vec(uint num, uint p)
{
  vec result(p);
  for (uint i = 0; i < p; i++) {
    result[(p - 1) - i] = num & 1;
    num >>= 1;
  }
  return(result);
}

mat greycode(int p)
{
  uint rows = 1 << p;
  mat result(rows, p);
  for (uint i = 0; i < rows; i++) {
    result.row(i) = binary_to_vec(binary_to_gray(i), p).t();
  }
  return(result);
}

mat diff(mat x)
{
  mat d(x.n_rows - 1, x.n_cols);
  // For each column, calculate the difference between the current row and the previous row
  for (uint i = 0; i < x.n_rows - 1; i++) {
    for (uint j = 0; j < x.n_cols; j++) {
      d(i, j) = x(i + 1, j) - x(i, j);
    }
  }
  return(d);
}

uvec get_rows(vector<pair<int, int> > pairs)
{
  uvec row_idx(pairs.size());
  for (uint i = 0; i < pairs.size(); i++) {
    row_idx[i] = pairs[i].first;
  }
  row_idx = unique(row_idx);
  // cout << "get_rows" << row_idx;
  return(row_idx);
}

uvec get_cols(vector<pair<int, int> > pairs)
{
  uvec col_idx(pairs.size());
  for (uint i = 0; i < pairs.size(); i++) {
    col_idx[i] = pairs[i].second;
  }
  col_idx = unique(col_idx);
  // cout << "get_cols" << col_idx;
  return(col_idx);
}

vec all_but_k_vec(vec v, uint k)
{
	vec result(v.n_rows - 1);
	uint i_idx = 0;
	for (uint i = 0; i < v.n_rows; i++) {
		if (i != k) {
			result[i_idx] = v[i];
			i_idx++;
		}
	}
	return result;
}

mat all_but_k_mat(mat m, uint k)
// Omit the kth row and column from m
{
	mat result(m.n_rows - 1, m.n_cols - 1);
	uint i_idx = 0, j_idx = 0;
	for (uint i = 0; i < m.n_rows; i++) {
		j_idx = 0;
		for (uint j = 0; j < m.n_cols; j++) {
			if (j != k && i != k) {
				result(i_idx, j_idx) = m(i, j);
				j_idx++;
			}
		}
		if (i != k) {
			i_idx++;
		}
	}
	return result;
}

vec ZE_exact_fast_cpp(vec vy, mat mX, bool LARGEP)
{
  int n = vy.n_elem;
  int p = mX.n_cols;
  int q = 0;
  mat mA = greycode(p);
  mat mD = diff(mA);
  vec vonep(p);
  vonep.ones();
  vector<bool> vs(mD.n_rows);
  vec result = mD * vonep;
  for (uint i = 0; i < mD.n_rows; i++) {
    vs[i] = result[i] == 1 ? true : false;
  }
  vec vq = mA * vonep;
  vec one_to_p(p);
  for (uint i = 0; i < p; i++) {
    one_to_p[i] = i + 1;
  }
  vec vw = abs(mD * one_to_p);

  vector<mat> lmZ(p);
  vector<vector<pair<uint, uint> > > linds11(p);
  vector<vector<pair<uint, uint> > > linds12(p);
  vector<vector<pair<uint, uint> > > linds21(p);
  vector<vector<pair<uint, uint> > > linds22(p);
  for (uint i=0; i < p; i++) {
    lmZ[i] = zeros<mat>(i + 1, i + 1);
    for (uint inds = 0; inds < i; inds++) {
      linds11[i].push_back(make_pair(inds, inds));
      linds12[i].push_back(make_pair(inds, i));
      linds21[i].push_back(make_pair(i, inds));
      linds22[i].push_back(make_pair(i, i));
    }
  }

  mat XTX = mX.t() * mX;
  vec XTy = mX.t() * vy;
  vec yTy = vy.t() * vy;
  vec vR2(mA.n_rows);
  vR2.zeros();
  vec Zb;
  uvec inds;
  uvec indsNew;
  uint k;
  for (uint j = 1; j < mA.n_rows; j++) {
    uvec newVar(1);
    newVar[0] = vw[j-1] - 1;
    if (vs[j-1]) {
      // Adding variable
      indsNew = join_vert(inds, newVar);
    } else {
      // Removing variable
      uvec kvec = find(inds == (vw[j-1] - 1));
	  k = kvec[0];
      indsNew = inds.rows(find(inds != newVar[0]));
    }
    vec b = XTy.elem(indsNew);
    q = vq[j] - 1;
    if (q == 0) {
      lmZ[0] = 1/XTX.submat(indsNew, indsNew);
      Zb = lmZ[0] * b;
    } else {
      if (vs[j-1]) {
        uvec newVar2(1);
        newVar2[0] = vw[j-1] - 1;
        vec v = XTX.submat(inds, newVar2);
        vec Zv = lmZ[q - 1] * v;
        vec d(Zv.n_rows);
        vec val = 1/(XTX.submat(newVar2, newVar2) - sum(v % Zv));
        d.fill(val[0]);
        lmZ[q].submat(get_rows(linds11[q]), get_cols(linds11[q])) = lmZ[q - 1] + d % Zv*Zv.t();
        vec val2 = -d % Zv;
        uvec rows = get_rows(linds12[q]);
        uvec cols = get_cols(linds12[q]);
        mat mat2(rows.n_rows, cols.n_rows);
        for (uint i2 = 0; i2 < cols.n_rows; i2++) {
        	mat2.col(i2) = val2;
        }

        lmZ[q].submat(rows, cols) = mat2;
        // lmZ[q][linds21[q]] <- -d*Zv;
        lmZ[q].submat(get_rows(linds21[q]), get_cols(linds21[q])) = mat2.t();
        uvec rows2 = get_rows(linds22[q]);
        uvec cols2 = get_cols(linds22[q]);
        mat d2(rows2.n_rows, cols2.n_cols);
        d2.fill(val[0]);
        lmZ[q].submat(rows2, cols2) = d2;
      } else {
        vec Z12;
        Z12 = all_but_k_vec(lmZ[q + 1].col(k), k);
        mat ZO = Z12 * Z12.t();
        mat m = all_but_k_mat(lmZ[q + 1], k);
        lmZ[q] = m - ZO / lmZ[q + 1](k,k);
      }
      Zb = lmZ[q] * b;
    }
    vR2[j] = sum(b % Zb);
    inds = indsNew;
  }
  vR2 = vR2 / yTy(0, 0);
  return(vR2);
}

// [[Rcpp::export]]
// List ZE_exact_fast(NumericVector vy_R, NumericMatrix mX_R)
// {
// 	vec vy = Rcpp::as<vec>(vy_R);
// 	mat mX = Rcpp::as<mat>(mX_R);
// 	vec vR2 = ZE_exact_fast_cpp(vy, mX, false);
// 	return Rcpp::List::create(Rcpp::Named("vR2")=vR2);
// }

int main(int argc, char **argv)
{
  vec vy;
  mat mX;
  vy.load("vy.csv", csv_ascii);
  // cout << "vy" << vy.submat(0, 0, 9, 0);
  mX.load("mX.csv", csv_ascii);
  // cout << "mX" << mX.submat(0, 0, 9, 18);
  //vy << 10 << 20 << endr;
  //mX << 10 << 20 << endr
  //   << 30 << 40 << endr;
  // Load vy
  // Load mX
  vec result = ZE_exact_fast_cpp(vy, mX, false);

  vec expected;
  expected.load("Hitters_exact2.csv", csv_ascii);
  cout << "C++ result, expected result, difference" << endl;
  for (uint i = 0; i < 10; i++) {
    cout << result(i) << " " << expected(i) << " " << result(i) - expected(i) << endl;
  }

  return 0;
}



