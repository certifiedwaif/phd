/*
 * main.cpp
 *
 *  Created on: 20/07/2015
 *      Author: markg
 */

#include <iostream>
#include <list>
#include <utility>
#include <armadillo>

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

vec binaryToVec(unsigned int num, int p)
{
  vec result(p);
  for (int i = 0; i < p; i++) {
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
  return(row_idx);
}

uvec get_cols(vector<pair<int, int> > pairs)
{
  uvec col_idx(pairs.size());
  for (unsigned int i = 0; i < pairs.size(); i++) {
    col_idx[i] = pairs[i].second;
  }
  return(col_idx);
}

void ZE_exact_fast(vec vy, mat mX, int LARGEP)
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
  // vw <- abs(mD%*%matrix(1:p,p,1))
  vec one_to_p(p);
  for (int i = 0; i < p; i++) {
    one_to_p[i] = i;
  }
  vec vw = abs(mD * one_to_p);
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
    lmZ[i] = zeros<mat>(i, i);
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
  for (unsigned int j = 1; j < mA.n_rows; j++)
  {
    uvec newVar(1);
    newVar[0] = vw[j-1];
    if (vs[j-1]) {
      // Adding variable
      // indsNew <- c(inds,vw[j-1])
      indsNew = join_vert(inds, newVar);
    } else {
      // Removing varable
  //     k <- which(inds==vw[j-1])
      uvec kvec = find(inds == vw[j-1]);
      k = kvec[1];
      // indsNew <- inds[-k]
      indsNew = find(inds != vw[j-1]);
    }
    // b <- XTy[indsNew]
    vec b = XTy.elem(indsNew);
    // q <- vq[j]
    q = vq[j];
    // if (q==1) {
    if (q == 1) {
      // lmZ[[1]] <- 1/XTX[indsNew,indsNew]
      lmZ[0] = 1/XTX.submat(indsNew, indsNew);
      // Zb <- lmZ[[1]]*b
      Zb = lmZ[0] * b;
    } else {
      if (vs[j-1]) {
        // v <- XTX[inds,vw[j-1]]
        vec newVar2(inds.n_rows);
        newVar2.fill(vw[j-1]);
        vec v = XTX.submat(inds, newVar);
        // Zv <- lmZ[[q-1]]%*%v
        vec Zv = lmZ[q - 2] * v;
        // d <- 1/(XTX[vw[j-1],vw[j-1]] - sum(v*Zv))
        mat d = 1/(XTX.submat(newVar, newVar)) - sum(v*Zv);
        // lmZ[[q]][linds11[[q]]] <- lmZ[[q-1]] + d*Zv%*%t(Zv)
        lmZ[q - 1].submat(get_rows(linds11[q - 1]), get_cols(linds11[q - 1])) = lmZ[q-2] + d*Zv*Zv.t();
        // lmZ[q][linds12[q]] <- -d*Zv;
        lmZ[q - 1].submat(get_rows(linds12[q - 1]), get_cols(linds12[q - 1])) = -d*Zv;
        // lmZ[q][linds21[q]] <- -d*Zv;
        lmZ[q - 1].submat(get_rows(linds21[q - 1]), get_cols(linds21[q - 1])) = -d*Zv;
        // lmZ[q][linds22[q]] <- d;
        lmZ[q - 1].submat(get_rows(linds22[q - 1]), get_cols(linds22[q - 1])) = d;
      } else {
        // Z12 <- lmZ[[q+1]][-k,k]
        mat Z12 = join_horiz(lmZ[q].rows(1, k - 1).col(k), lmZ[q].rows(k + 1, lmZ[q].n_rows).col(k));
        // ZO <- Z12%*%t(Z12)
        mat ZO = Z12 * Z12.t();
        // lmZ[[q]] <- lmZ[[q+1]][-k,-k] -  ZO/lmZ[[q+1]][k,k]
        mat m = join_vert(join_horiz(lmZ[q].submat(0, 0, k-1, k-1), lmZ[q].submat(0, k+1, k-1, k-1)),
                          join_horiz(lmZ[q].submat(k+1, 0, lmZ[q].n_rows, k-1), lmZ[q].submat(k+1, k+1, lmZ[q].n_rows, lmZ[q].n_cols)));
        lmZ[q - 1] = m - ZO/lmZ[q](k,k);
      }
      // Zb <- lmZ[[q]]%*%b
      Zb = lmZ[q - 1] * b;
    }
    // vR2[j]  <- sum(b*Zb)
    vR2[j] = sum(b * Zb);
    // inds <- indsNew
    inds = indsNew;
  }
  // vR2 <- vR2/yTy
  vR2 = vR2/yTy;
  // vlog.ZE  <-  -res.con$vcon[vq+1]*log(1 - vR2) + res.con$vpen[vq+1]
  // vlog.ZE  = -res_con$vcon[vq+1]*log(1 - vR2) + res_con$vpen[vq+1];
  // return(list(mA=mA,vlog.ZE=vlog.ZE))
}

int main(int argc, char **argv)
{
  vec vy;
  mat mX;
  vy.load("vy.csv", csv_ascii);
  mX.load("mX.csv", csv_ascii);
  //vy << 10 << 20 << endr;
  //mX << 10 << 20 << endr
  //   << 30 << 40 << endr;
  // Load vy
  // Load mX
  ZE_exact_fast(vy, mX, 19);

  return 0;
}



