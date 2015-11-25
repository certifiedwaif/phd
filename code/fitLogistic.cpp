// fitLogistic.R

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::plugins(openmp)]]
#include <omp.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

mat ALL_bohning_cpp(vec vy, mat mX, mat models, double tau, double MAXABSERR);

// [[Rcpp::export]]
List ALL_bohning(NumericVector vy_R, NumericMatrix mX_R, NumericMatrix models_R, double tau, double MAXABSERR)
{
	vec vy = Rcpp::as<vec>(vy_R);
	mat mX = Rcpp::as<mat>(mX_R);
	mat models = Rcpp::as<mat>(models_R);

	mat mBeta = ALL_bohning_cpp(vy, mX, models, tau, MAXABSERR);
	return List::create(Rcpp::Named("mBeta")=mBeta);
}

// ALL.bohning <- function(vy,mX,models,tau=1.0E-5,MAXABSERR=1.0E-3)    
// {
mat ALL_bohning_cpp(const vec vy, const mat mX, const mat models, const double tau, const double MAXABSERR)
{
	// Rcout << "Entering ALL_bohning_cpp" << endl;
// 	MAXITER <- 1000
	const int MAXITER = 1000;

// 	n <- length(vy)
	const int n = vy.n_rows;
// 	p <- ncol(mX)
	const int p = mX.n_cols;
	
// 	mBeta   <- matrix(0,nrow(models),p)
	// Rcout << "Constructing mBeta" << endl;
	mat mBeta(models.n_rows, p);
	
// 	for (j in 1:nrow(models)) 
// 	{	 
	// #pragma omp parallel for shared(mBeta) private(vgamma, inds1, one_to_n, mX1, vmu, mI, mSigma_inv, mA, vmu_old, vxi, vb, vc, err) default(none) schedule(auto)
	// omp_set_num_threads(4);
	// #pragma omp parallel for
	for (int j = 0; j < models.n_rows; j++) {
// 		vgamma <- c(1,models[j,])
		// Rcout << "Constructing vgamma" << endl;
		vec vgamma(models.n_cols + 1);
		vgamma(0) = 1;
		vgamma.rows(1, models.n_cols) = models.row(j).t();
// 		inds1 <- which(vgamma==1)
		// Rcout << "Finding one indices" << endl;
		uvec inds1 = find(vgamma == 1);
// 		q <- length(inds1)
		const int q = inds1.n_rows;
// 		mX1 <- matrix(mX[,inds1],n,q)
		uvec one_to_n = linspace<uvec>(0, mX.n_rows - 1, mX.n_rows);
		// Rcout << "Constructing mX1" << endl;
		// Rcout << "one_to_n" << one_to_n.tail_rows(10) << endl;
		// Rcout << "inds1" << inds1 << endl;
		// Rcout << "mX.n_rows" << mX.n_rows << endl;
		// Rcout << "mX.n_cols" << mX.n_cols << endl;
		mat mX1 = mX.submat(one_to_n, inds1);

		// Rcout << "Constructing vmu" << endl;
// 		vmu <- matrix(0,q,1)
		vec vmu(q, fill::zeros);
// 		mI <- diag(1,q)
		mat mI(q, q, fill::eye);
// 		mSigma.inv <- 0.25*t(mX1)%*%mX1 + tau*mI
		// Rcout << "Constructing mSigma_inv" << endl;
		mat mSigma_inv = 0.25 * mX1.t() * mX1 + tau * mI;
// 		mA <- solve(mSigma.inv,t(mX1),tol=1.0E-99)
		// Rcout << "Constructing mA" << endl;
		mat mA = solve(mSigma_inv, mX1.t());
// 		for (ITER in 1:MAXITER) 
// 		{        
		for (int i = 0; i < MAXITER; i++) {
// 			vmu.old <- vmu
			// Rcout << "Constructing vmu_old" << endl;
			vec vmu_old = vmu;
// 			vxi <- mX1%*%vmu   
			// Rcout << "Constructing vxi" << endl;
			vec vxi = mX1 * vmu;
// 			vb <- 1/(1+exp(-vxi))
			// Rcout << "Constructing vb" << endl;
			vec vb = 1 / (1 + exp(-vxi));
// 			vc <- vy + 0.25*vxi - vb
			// Rcout << "Constructing vc" << endl;
			vec vc = vy + 0.25 * vxi - vb;
// 			vmu <- mA%*%vc
			// Rcout << "Constructing vmu" << endl;
			vmu = mA * vc;
// 			err <- max(abs(vmu.old - vmu))
			// Rcout << "Checking error" << endl;
			double err = max(abs(vmu_old - vmu));
			// Rcout << "Err " << err << endl;
// 			#cat(ITER,err,"\n")
			if (err < MAXABSERR) {

// 			if (err < MAXABSERR) {
// 				break;
// 			} 
// 		}
				// Rcout << "Leaving loop" << endl;
				break;
			}
// 		mBeta[j,inds1] <- vmu
// 	}
		}
		uvec j_vec(1);
		j_vec.fill(j);
		// Rcout << "Filling mBeta" << endl;
		mBeta.submat(j_vec, inds1) = vmu.t();
	}
// 	return(list(mBeta=mBeta))
// }

	return mBeta;
}