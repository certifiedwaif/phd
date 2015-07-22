/*
 * main.cpp
 *
 *  Created on: 20/07/2015
 *      Author: markg
 */

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <list>
#include <utility>

using namespace Rcpp;
using namespace std;
using namespace arma;

// From the Wikipedia page on Gray code
/*
        The purpose of this function is to convert an unsigned
        binary number to reflected binary Gray code.

        The operator >> is shift right. The operator ^ is exclusive or.
*/
unsigned int binaryToGray(unsigned int num)
{
  return (num >> 1) ^ num;
}

/*
        The purpose of this function is to convert a reflected binary
        Gray code number to a binary number.
*/
unsigned int grayToBinary(unsigned int num)
{
  unsigned int mask;
  for (mask = num >> 1; mask != 0; mask = mask >> 1)
  {
    num = num ^ mask;
  }
  return num;
}

vec binaryToVec(unsigned int num, unsigned int p)
{
  vec result(p);
  for (unsigned int i = 0; i < p; i++) {
    result[(p - 1) - i] = num & 1;
    num >>= 1;
  }
  return(result);
}

mat greycode(int p)
{
  unsigned int rows = 1 << p;
  mat result(rows, p);
  for (unsigned int i = 0; i < rows; i++) {
    result.row(i) = binaryToVec(binaryToGray(i), p).t();
  }
  return(result);
}

mat diff(mat x)
{
  mat d(x.n_rows - 1, x.n_cols);
  // For each column, calculate the difference between the current row and the previous row
  for (unsigned int i = 0; i < x.n_rows - 1; i++) {
    for (unsigned int j = 0; j < x.n_cols; j++) {
      d(i, j) = x(i + 1, j) - x(i, j);
    }
  }
  return(d);
}

uvec get_rows(vector<pair<int, int> > pairs)
{
  uvec row_idx(pairs.size());
  for (unsigned int i = 0; i < pairs.size(); i++) {
    row_idx[i] = pairs[i].first;
  }
  row_idx = unique(row_idx);
  // cout << "get_rows" << row_idx;
  return(row_idx);
}

uvec get_cols(vector<pair<int, int> > pairs)
{
  uvec col_idx(pairs.size());
  for (unsigned int i = 0; i < pairs.size(); i++) {
    col_idx[i] = pairs[i].second;
  }
  col_idx = unique(col_idx);
  // cout << "get_cols" << col_idx;
  return(col_idx);
}

vec all_but_k_vec(vec v, unsigned int k)
{
	vec result(v.n_rows - 1);
	unsigned int i_idx = 0;
	for (unsigned int i = 0; i < v.n_rows; i++) {
		if (i != k) {
			result[i_idx] = v[i];
			i_idx++;
		}
	}
	return result;
}

mat all_but_k_mat(mat m, unsigned int k)
// Omit the kth row and column from m
{
	mat result(m.n_rows - 1, m.n_cols - 1);
	unsigned int i_idx = 0, j_idx = 0;
	for (unsigned int i = 0; i < m.n_rows; i++) {
		j_idx = 0;
		for (unsigned int j = 0; j < m.n_cols; j++) {
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
  // n <- length(vy)
  int n = vy.n_elem;
  // p <- ncol(mX)
  int p = mX.n_cols;
  // q <- 0
  int q = 0;
  // mA <- greycode(p)
  mat mA = greycode(p);
  // mD <- diff(mA)
  mat mD = diff(mA);
  // vonep <- matrix(1,p,1)
  vec vonep(p);
  vonep.ones();
  // vs <- (mD%*%vonep)==1
  vector<bool> vs(mD.n_rows);
  vec result = mD * vonep;
  for (unsigned int i = 0; i < mD.n_rows; i++) {
    vs[i] = result[i] == 1 ? true : false;
  }
  // vq <- mA%*%vonep
  vec vq = mA * vonep;
  // cout << "vq " << vq.submat(0, 0, 9, 0);
  // vw <- abs(mD%*%matrix(1:p,p,1))
  vec one_to_p(p);
  for (int i = 0; i < p; i++) {
    one_to_p[i] = i + 1;
  }
  vec vw = abs(mD * one_to_p);
  // cout << "vw " << vw.submat(0, 0, 9, 0);
  // res.con <- ZE.constants(n,p,LARGEP)

  // lmZ <- list()
  std::vector<mat> lmZ(p);
  // linds11 <- list()
  // linds12 <- list()
  // linds21 <- list()
  // linds22 <- list()
  std::vector<std::vector<pair<int, int> > > linds11(p);
  std::vector<std::vector<pair<int, int> > > linds12(p);
  std::vector<std::vector<pair<int, int> > > linds21(p);
  std::vector<std::vector<pair<int, int> > > linds22(p);
  // for (i in 1:p) {
  for (int i=0; i < p; i++) {
  //   lmZ[[i]] <- matrix(0,i,i)
    lmZ[i] = zeros<mat>(i + 1, i + 1);
  //   inds <- 1:(i-1)
  //   mZ11 <- matrix(FALSE,i,i); mZ11[inds,inds] <- TRUE
  //   mZ12 <- matrix(FALSE,i,i); mZ12[inds,i]    <- TRUE
  //   mZ21 <- matrix(FALSE,i,i); mZ21[i,inds]    <- TRUE
  //   mZ22 <- matrix(FALSE,i,i); mZ22[i,i]       <- TRUE
  //   linds11[[i]] <- which(mZ11)
  //   linds12[[i]] <- which(mZ12)
  //   linds21[[i]] <- which(mZ21)
  //   linds22[[i]] <- which(mZ22)
    for (int inds = 0; inds < i; inds++) {
      linds11[i].push_back(make_pair(inds, inds));
      linds12[i].push_back(make_pair(inds, i));
      linds21[i].push_back(make_pair(i, inds));
      linds22[i].push_back(make_pair(i, i));
    }
  }

  // XTX <- t(mX)%*%mX
  mat XTX = mX.t() * mX;
  // XTy <- t(mX)%*%vy
  vec XTy = mX.t() * vy;
  // yTy <- t(vy)%*%vy
  vec yTy = vy.t() * vy;
  // vR2  <- rep(0,nrow(mA))
  vec vR2(mA.n_rows);
  vR2.zeros();
  vec Zb;
  // inds <- c()
  uvec inds;
  uvec indsNew;
  unsigned int k;
  // for (j in 2:nrow(mA))
  for (unsigned int j = 1; j < mA.n_rows; j++) {
    uvec newVar(1);
    newVar[0] = vw[j-1] - 1;
    if (vs[j-1]) {
      // Adding variable
      // indsNew <- c(inds,vw[j-1])
      indsNew = join_vert(inds, newVar);
      // cout << "Add variable" << newVar;
      // cout << "indsNew" << indsNew;
    } else {
      // Removing variable
        // cout << "Remove variable" << newVar;
  //     k <- which(inds==vw[j-1])
      uvec kvec = find(inds == (vw[j-1] - 1));
      // cout << "vw[j-1]" << (vw[j-1] - 1) << endl;
      // cout << "inds " << inds;
      // cout << "kvec " << kvec;
	  k = kvec[0];
      // indsNew <- inds[-k]
      indsNew = inds.rows(find(inds != newVar[0]));
      // cout << "indsNew" << indsNew;
    }
    // b <- XTy[indsNew]
    vec b = XTy.elem(indsNew);
    // cout << "XTy" << XTy;
    // cout << "b" << b;
    // q <- vq[j]
    q = vq[j] - 1;
    // cout << "q" << q << endl;
    // if (q==1) {
    if (q == 0) {
      // lmZ[[1]] <- 1/XTX[indsNew,indsNew]
      lmZ[0] = 1/XTX.submat(indsNew, indsNew);
      // Zb <- lmZ[[1]]*b
      Zb = lmZ[0] * b;
    } else {
      if (vs[j-1]) {
        // v <- XTX[inds,vw[j-1]]
        uvec newVar2(1);
        newVar2[0] = vw[j-1] - 1;
        // cout << "j " << j << " ";
        // cout << "rows" << inds;
        // cout << "cols" << newVar2;
        vec v = XTX.submat(inds, newVar2);
        // Zv <- lmZ[[q-1]]%*%v
        vec Zv = lmZ[q - 1] * v;
        // d <- 1/(XTX[vw[j-1],vw[j-1]] - sum(v*Zv))
        vec d(Zv.n_rows);
        vec val = 1/(XTX.submat(newVar2, newVar2) - sum(v % Zv));
        d.fill(val[0]);
        // cout << "v" << v;
        // cout << "Zv" << Zv;
        // cout << "d" << d;
        // lmZ[[q]][linds11[[q]]] <- lmZ[[q-1]] + d*Zv%*%t(Zv)
        lmZ[q].submat(get_rows(linds11[q]), get_cols(linds11[q])) = lmZ[q - 1] + d % Zv*Zv.t();
        // lmZ[q][linds12[q]] <- -d*Zv;
        vec val2 = -d % Zv;
        uvec rows = get_rows(linds12[q]);
        uvec cols = get_cols(linds12[q]);
        mat mat2(rows.n_rows, cols.n_rows);
        for (unsigned int i2 = 0; i2 < cols.n_rows; i2++) {
        	mat2.col(i2) = val2;
        }

        lmZ[q].submat(rows, cols) = mat2;
        // lmZ[q][linds21[q]] <- -d*Zv;
        lmZ[q].submat(get_rows(linds21[q]), get_cols(linds21[q])) = mat2.t();
        // lmZ[q][linds22[q]] <- d;
        uvec rows2 = get_rows(linds22[q]);
        uvec cols2 = get_cols(linds22[q]);
        mat d2(rows2.n_rows, cols2.n_cols);
        d2.fill(val[0]);
        lmZ[q].submat(rows2, cols2) = d2;
      } else {
        // Z12 <- lmZ[[q+1]][-k,k]

        vec Z12;
        Z12 = all_but_k_vec(lmZ[q + 1].col(k), k);
//        // cout << "Z12" << Z12;
        // ZO <- Z12%*%t(Z12)
        mat ZO = Z12 * Z12.t();
        // lmZ[[q]] <- lmZ[[q+1]][-k,-k] -  ZO/lmZ[[q+1]][k,k]
        mat m = all_but_k_mat(lmZ[q + 1], k);
        lmZ[q] = m - ZO / lmZ[q + 1](k,k);
      }
      // Zb <- lmZ[[q]]%*%b
      Zb = lmZ[q] * b;
    }
    // vR2[j]  <- sum(b*Zb)
    vR2[j] = sum(b % Zb);
    // cout << "vR2[j]" << vR2[j] << endl;
    // inds <- indsNew
    inds = indsNew;
  }
  // vR2 <- vR2/yTy
  vR2 = vR2/yTy(0, 0);
  // cout << "First ten elements of vR2" << vR2.submat(0, 0, 9, 0);
  // Do this in the calling code, where it's easier.
  // vlog.ZE  <-  -res.con$vcon[vq+1]*log(1 - vR2) + res.con$vpen[vq+1]
  // vlog.ZE  = -res_con$vcon[vq+1]*log(1 - vR2) + res_con$vpen[vq+1];
  // return(list(mA=mA,vlog.ZE=vlog.ZE))
  return(vR2);
}

// [[Rcpp::export]]
List ZE_exact_fast(NumericVector vy_R, NumericMatrix mX_R)
{
	vec vy = Rcpp::as<vec>(vy_R);
	mat mX = Rcpp::as<mat>(mX_R);
	vec vR2 = ZE_exact_fast_cpp(vy, mX, false);
	return Rcpp::List::create(Rcpp::Named("vR2")=vR2);
}

// int main(int argc, char **argv)
// {
//   vec vy;
//   mat mX;
//   vy.load("vy.csv", csv_ascii);
//   cout << "vy" << vy.submat(0, 0, 9, 0);
//   mX.load("mX.csv", csv_ascii);
//   cout << "mX" << mX.submat(0, 0, 9, 18);
//   //vy << 10 << 20 << endr;
//   //mX << 10 << 20 << endr
//   //   << 30 << 40 << endr;
//   // Load vy
//   // Load mX
//   vec result = ZE_exact_fast(vy, mX, 19);

//   return 0;
// }



