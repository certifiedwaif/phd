#include <iostream>
#include <armadillo>
#include <list>
#include <utility>

using namespace std;
using namespace arma;

mat greycode(int p)
{
	mat result(p, p);
	throw "greycode not implemented yet\n";
	return(result);
}

mat diff(mat x)
{
	mat result(x.n_rows, x.n_cols);
	throw "diff not implemented yet\n";
	return(result);
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
	for (int i = 0; i < mD.n_rows; i++) {
		vs[i] = result[i, 1] == 1 ? true : false;
	}
	// vq <- mA%*%vonep
	vec vq = mA * vonep;
	// vw <- abs(mD%*%matrix(1:p,p,1))
	vec one_to_p;
	for (int i = 0; i < p; i++) {
		one_to_p[i] = i;
	}
	vec vw = abs(mD * one_to_p);
	// res.con <- ZE.constants(n,p,LARGEP) 
	
	// lmZ <- list()
	std::vector<mat> lmZ;
	// linds11 <- list()
	// linds12 <- list()
	// linds21 <- list()
	// linds22 <- list()
	std::vector<std::vector<pair<int, int> > > linds11;
	std::vector<std::vector<pair<int, int> > > linds12;
	std::vector<std::vector<pair<int, int> > > linds21;
	std::vector<std::vector<pair<int, int> > > linds22;
	// for (i in 1:p) {
	for (int i=0; i < p; i++) {
	// 	lmZ[[i]] <- matrix(0,i,i)
		lmZ[i].zeros();
	// 	inds <- 1:(i-1)
	// 	mZ11 <- matrix(FALSE,i,i); mZ11[inds,inds] <- TRUE
	// 	mZ12 <- matrix(FALSE,i,i); mZ12[inds,i]    <- TRUE
	// 	mZ21 <- matrix(FALSE,i,i); mZ21[i,inds]    <- TRUE
	// 	mZ22 <- matrix(FALSE,i,i); mZ22[i,i]       <- TRUE
	// 	linds11[[i]] <- which(mZ11)
	// 	linds12[[i]] <- which(mZ12)
	// 	linds21[[i]] <- which(mZ21)
	// 	linds22[[i]] <- which(mZ22)
		for (int inds = 0; inds < i - 2; inds++) {
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
	for (int j = 1; j < mA.n_rows; j++)
	{
		uvec newVar(1);
		newVar[1] = vw[j-1];
		if (vs[j-1]) {
			// Adding variable
			// indsNew <- c(inds,vw[j-1])
			indsNew = join_horiz(inds, newVar);
		} else {
			// Removing varable
	// 		k <- which(inds==vw[j-1])
			uvec kvec = find(inds == vw[j-1]);
			k = kvec[1];
	// 		indsNew <- inds[-k]
			indsNew = find(inds != vw[j-1]);
		}
	// 	b <- XTy[indsNew]
		vec b = XTy.elem(indsNew);
	// 	q <- vq[j]
		q = vq[j];
	// 	if (q==1) {
		if (q == 1) {
	// 		lmZ[[1]] <- 1/XTX[indsNew,indsNew] 
			lmZ[1] = 1/XTX.submat(indsNew, indsNew);
	// 		Zb <- lmZ[[1]]*b
			vec Zb = lmZ[1] * b;
		} else {
			if (vs[j-1]) {
	// 			v <- XTX[inds,vw[j-1]]
				vec newVar2(inds.n_rows);
				newVar2.fill(vw[j-1]);
				vec v = XTX.submat(inds, newVar);
	// 			Zv <- lmZ[[q-1]]%*%v
				vec Zv = lmZ[q-1] * v;
	// 			d <- 1/(XTX[vw[j-1],vw[j-1]] - sum(v*Zv))
				mat d = 1/(XTX.submat(newVar, newVar)) - sum(v*Zv);
	// 			lmZ[[q]][linds11[[q]]] <- lmZ[[q-1]] + d*Zv%*%t(Zv)
	// 			lmZ[[q]][linds12[[q]]] <- -d*Zv
	// 			lmZ[[q]][linds21[[q]]] <- -d*Zv
	// 			lmZ[[q]][linds22[[q]]] <- d
				uvec row_idx11(linds11.size());
				uvec col_idx11(linds11.size());
				for (int i = 0; i < row_idx11.n_rows; i++) {
					row_idx11[i] = linds11[q][i].first;
					col_idx11[i] = linds11[q][i].second;
				}
				lmZ[q].submat(row_idx11, col_idx11) = lmZ[q-1] + d*Zv*Zv.t();
				// lmZ[q][linds12[q]] <- -d*Zv;
				uvec row_idx12(linds12.size());
				uvec col_idx12(linds12.size());
				for (int i = 0; i < row_idx12.n_rows; i++) {
					row_idx12[i] = linds12[q][i].first;
					col_idx12[i] = linds12[q][i].second;
				}
				lmZ[q].submat(row_idx12, col_idx12) = -d*Zv;
				// lmZ[q][linds21[q]] <- -d*Zv;
				uvec row_idx21(linds21.size());
				uvec col_idx21(linds21.size());
				for (int i = 0; i < row_idx21.n_rows; i++) {
					row_idx21[i] = linds21[q][i].first;
					col_idx21[i] = linds21[q][i].second;
				}
				lmZ[q].submat(row_idx21, col_idx21) = -d*Zv;
				// lmZ[q][linds22[q]] <- d;
				uvec row_idx22(linds22.size());
				uvec col_idx22(linds22.size());
				for (int i = 0; i < row_idx22.n_rows; i++) {
					row_idx22[i] = linds22[q][i].first;
					col_idx22[i] = linds22[q][i].second;
				}
				lmZ[q].submat(row_idx22, col_idx22) = d;
			} else {
	// 			Z12 <- lmZ[[q+1]][-k,k]
				mat Z12 = join_horiz(lmZ[q + 1].rows(1, k - 1).col(k), lmZ[q + 1].rows(k + 1, lmZ[q + 1].n_rows).col(k));
	// 			ZO <- Z12%*%t(Z12)
				mat ZO = Z12 * Z12.t();
	// 			lmZ[[q]] <- lmZ[[q+1]][-k,-k] -  ZO/lmZ[[q+1]][k,k]
				mat m = join_vert(join_horiz(lmZ[q+1].submat(1, 1, k-1, k-1), lmZ[q+1].submat(1, k+1, 1, k-1)),
													join_horiz(lmZ[q+1].submat(k+1, 1, lmZ[q + 1].n_rows, k-1), lmZ[q+1].submat(k+1, k+1, lmZ[q + 1].n_rows, lmZ[q + 1].n_cols)));
				lmZ[q] = m - ZO/lmZ[q+1][k,k];
			}
	// 		Zb <- lmZ[[q]]%*%b
			Zb = lmZ[q] * b;
		} 
	// 	vR2[j]  <- sum(b*Zb) 
		vR2[j] = sum(b * Zb);
	// 	inds <- indsNew
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
	return 0;
}