

#include <math.h>

/* Global Variables */
#define TINY   1.0E-30 

///////////////////////////////////////////////////////////////////////////////

double gTilde(double z) 
{
    double p;
    double q;
    double zp1;
    double val;
    if (z<0.000001) {
        p = 420.0 + z*(540.0 + z*(180.0 + 12.0*z));
        q = 420.0 + z*(750.0 + z*(380.0 + 47.0*z));
        val = p/q;
    } else {
        zp1 = z + 1.0;
        val = z/(zp1*log(zp1));
    }
    return(val);
}

///////////////////////////////////////////////////////////////////////////////

double hTilde(double z) 
{
    double p;
    double q;
    double zp1;
    if (z<0.000001) {
        p = 18060.0 + z*(14760.0 + z*(2100.0 - 48.0*z));
        q = 18060.0 + z*(41850.0 + z*(30260.0 + 6517.0*z));
        return( p/q );
    } else {
        zp1 = z + 1.0;
        return( z/(zp1*zp1*log(zp1)) );
    }
}

///////////////////////////////////////////////////////////////////////////////

double log_b0_safe(double x) 
{
    if (x > 500.0) {
        return( log(x) );
    }
    if (x>= -35.0) {
        return( log(log(1.0 + exp(x))) );
    }
    if (x>= -500.0) {
        return( x + log1p(-exp(x)) );
    } else {
        return( x );
    }
}

///////////////////////////////////////////////////////////////////////////////

double b0_safe(double x) 
{
    if (x > 500) {
        return( x );
    } else {
        return( log(1.0 + exp(x)) );
    }
}

///////////////////////////////////////////////////////////////////////////////

double b1_safe(double x) 
{
    return( 1.0/(1.0 + exp(-x)) );
}

///////////////////////////////////////////////////////////////////////////////

double b2_safe(double x) 
{
    double b1val = b1_safe(x);
    return( b1val*(1.0 - b1val) );
}

///////////////////////////////////////////////////////////////////////////////

double dsnorm(double x, double mu, double sigma, int ISLOG) 
{
    double z;
    double dmu;
    dmu = x - mu;
    z = dmu/sigma;
    double minusHalfLog2pi = -0.9189385332046727;
    if(ISLOG) {
        return( minusHalfLog2pi - log(sigma) - 0.5*z*z  );
    } else {
        return( exp(minusHalfLog2pi - log(sigma) - 0.5*z*z) ); 
    }
}


///////////////////////////////////////////////////////////////////////////////

void StartingValue(double mu, double sigma, double* x, double* f) 
{
    int i;
    double sigma2 = sigma*sigma;
    double sqrtDisc = sqrt(0.25*mu*mu + sigma2);
    double maxVal = -1.0;
    int maxInd = 0;
        
    double vx[2];
    double vf[2];
    vx[0] = 0.0;
    vx[1] = 0.5*mu + sqrtDisc;
    
    for (i=0;i<2;i++) {
        vf[i] = b0_safe(vx[i])*dsnorm(vx[i],mu,sigma,0);
        if (vf[i]>maxVal) {
            maxVal = vf[i];
            maxInd = i;
        }
    }
    (*x) = vx[maxInd];
    (*f) = maxVal;
}

///////////////////////////////////////////////////////////////////////////////

void ModeAndEffSupp(double mu, double sigma, double *x, double *f, double *L, double *R) 
{
    int i;
    double z; 
    double xVal = (*x);
    double sigma2 = sigma*sigma;
    double sigma2inv = 1.0/sigma2;
    double gTil;
    double hTil;
    double fVal;  
    double fL; 
    double fR;           
    double g;
    double gL;
    double gR;        
    double h;
    double dx;
    double C = -12.0;
    double Lval;
    double Rval;
    double dmu;
    for (i=0;i<100;i++) {
        z = exp(xVal);
        gTil = gTilde(z);
        hTil = hTilde(z);  
        dmu = xVal - mu;
        g = gTil - sigma2inv*dmu;
        h = hTil - gTil*gTil - sigma2inv;
        dx = g/h;
        xVal = xVal - dx;
        if (fabs(dx)<1.0E-2) {
            break;
        } 
    }
    fVal = b0_safe(xVal)*dsnorm(xVal,mu,sigma,0);
    
    dx = 1.0/sqrt(-h);
    
    Lval = xVal - 6.0*dx;
    for (i=0;i<100;i++) {
        fL = log_b0_safe(Lval) + dsnorm(Lval,mu,sigma,1);
        z = exp(Lval);
        gTil = gTilde(z);      
        gL = gTil - (Lval - mu)/sigma2;
        Lval = Lval - (fL - C)/gL;
        if (fabs(fL-C)<1.0E-2) {
            break;
        }
    }    
    
    Rval = xVal + 6.0*dx;
    for (i=0;i<100;i++) {
        fR = log_b0_safe(Rval) + dsnorm(Rval,mu,sigma,1);
        z = exp(Rval);
        gTil = gTilde(z);      
        gR = gTil - (Rval - mu)/sigma2;
        Rval = Rval - (fR - C)/gR;
        if (fabs(fR-C)<1.0E-2) {
            break;
        }
    }       
    
    (*x) = xVal;
    (*f) = fVal;
    (*L) = Lval;
    (*R) = Rval;    
}

///////////////////////////////////////////////////////////////////////////////

double B0_trapint(double mu, double sigma, int N, double *x, double *f, double *L, double *R) 
{
    double Lval = (*L);
    double Rval = (*R);
    int i;
    double dx = (Rval - Lval)/((double) N);
    double val = 0.0;
    double fL;
    double fR;
    fL = b0_safe(Lval)*dsnorm(Lval,mu,sigma,0);
    fR = b0_safe(Rval)*dsnorm(Rval,mu,sigma,0);
    val = 0.5*(fL+fR);
    for (i=0;i<(N-1);i++) {
        Lval = Lval + dx;
        val = val + b0_safe(Lval)*dsnorm(Lval,mu,sigma,0);
    }
    val = val*dx;
    return(val);
}

///////////////////////////////////////////////////////////////////////////////

double B0(double mu, double sigma, int N) 
{    
    double x = 0.0;
    double f = 0.0;
    double L = 0.0;
    double R = 0.0;
    double val = 0.0;
    if (mu<0) {
        return( mu + B0( -mu, sigma, N) ); 
    } else {
        StartingValue(mu, sigma, &x, &f);
        ModeAndEffSupp(mu, sigma, &x, &f, &L, &R); 
        val = B0_trapint(mu, sigma, N, &x, &f, &L, &R);
        return( val );
    }
} 

///////////////////////////////////////////////////////////////////////////////

void R_B0(double *mu, double *sigma, int *N, double *val) 
{
    (*val) = B0( *mu, *sigma, (*N));
} 

///////////////////////////////////////////////////////////////////////////////

void R_vB0(double *mu, double *sigma, int *N, double *val,  int *n) 
{
    int i;
    for (i=0;i<(*n);i++) {
        val[i] = B0( mu[i], sigma[i], (*N));
    }
} 

///////////////////////////////////////////////////////////////////////////////

double B1_trapint(double mu, double sigma, int N, double *x, double *f, double *L, double *R) 
{
    double Lval = (*L);
    double Rval = (*R);
    int i;
    double dx = (Rval - Lval)/((double) N);
    double val = 0.0;
    double fL;
    double fR;
    fL = b1_safe(Lval)*dsnorm(Lval,mu,sigma,0);
    fR = b1_safe(Rval)*dsnorm(Rval,mu,sigma,0);
    val = 0.5*(fL+fR);
    for (i=0;i<(N-1);i++) {
        Lval = Lval + dx;
        val = val + b1_safe(Lval)*dsnorm(Lval,mu,sigma,0);
    }
    val = val*dx;
    return(val);
}

///////////////////////////////////////////////////////////////////////////////

double B1(double mu, double sigma, int N) 
{    
    double x = 0.0;
    double f = 0.0;
    double L = 0.0;
    double R = 0.0;
    double val = 0.0;
    if (mu<0) {
        return( 1.0 - B1( -mu, sigma, N) ); 
    } else {
        StartingValue(mu, sigma, &x, &f);
        ModeAndEffSupp(mu, sigma, &x, &f, &L, &R); 
        val = B1_trapint(mu, sigma, N, &x, &f, &L, &R);
        return( val );
    }
} 

///////////////////////////////////////////////////////////////////////////////

void R_B1(double *mu, double *sigma, int *N, double *val) 
{
    (*val) = B1( *mu, *sigma, (*N));
} 

///////////////////////////////////////////////////////////////////////////////

void R_vB1(double *mu, double *sigma, int *N, double *val,  int *n) 
{
    int i;
    for (i=0;i<(*n);i++) {
        val[i] = B1( mu[i], sigma[i], (*N));
    }
} 

///////////////////////////////////////////////////////////////////////////////

double B2_trapint(double mu, double sigma, int N, double *x, double *f, double *L, double *R) 
{
    double Lval = (*L);
    double Rval = (*R);
    int i;
    double dx = (Rval - Lval)/((double) N);
    double val = 0.0;
    double fL;
    double fR;
    fL = b2_safe(Lval)*dsnorm(Lval,mu,sigma,0);
    fR = b2_safe(Rval)*dsnorm(Rval,mu,sigma,0);
    val = 0.5*(fL+fR);
    for (i=0;i<(N-1);i++) {
        Lval = Lval + dx;
        val = val + b2_safe(Lval)*dsnorm(Lval,mu,sigma,0);
    }
    val = val*dx;
    return(val);
}

///////////////////////////////////////////////////////////////////////////////

double B2(double mu, double sigma, int N) 
{    
    double x = 0.0;
    double f = 0.0;
    double L = 0.0;
    double R = 0.0;
    double val = 0.0;
    if (mu<0) {
        return( B2( -mu, sigma, N) ); 
    } else {
        StartingValue(mu, sigma, &x, &f);
        ModeAndEffSupp(mu, sigma, &x, &f, &L, &R); 
        val = B2_trapint(mu, sigma, N, &x, &f, &L, &R);
        return( val );
    }
} 

///////////////////////////////////////////////////////////////////////////////

void R_B2(double *mu, double *sigma, int *N, double *val) 
{
    (*val) = B2( *mu, *sigma, (*N));
} 

///////////////////////////////////////////////////////////////////////////////

void R_vB2(double *mu, double *sigma, int *N, double *val,  int *n) 
{
    int i;
    for (i=0;i<(*n);i++) {
        val[i] = B2( mu[i], sigma[i], (*N));
    }
} 

///////////////////////////////////////////////////////////////////////////////

void B12_trapint(double mu, double sigma, int N, double *x, double *f, double *L, double *R, double* val1, double* val2) 
{
    double Lval = (*L);
    double Rval = (*R);
    int i;
    double dx = (Rval - Lval)/((double) N);
    double qL  = dsnorm(Lval,mu,sigma,0);
    double fL1 = b1_safe(Lval)*qL;
    double fL2 = b2_safe(Lval)*qL;
    double qR  = dsnorm(Rval,mu,sigma,0);
    double fR1 = b1_safe(Rval)*qR;
    double fR2 = b2_safe(Rval)*qR;
    (*val1) = 0.5*(fL1+fR1);
    (*val2) = 0.5*(fL2+fR2);
    for (i=0;i<(N-1);i++) {
        Lval = Lval + dx;
				qL  = dsnorm(Lval,mu,sigma,0);
        (*val1) = (*val1) + b1_safe(Lval)*qL;
        (*val2) = (*val2) + b2_safe(Lval)*qL;        
    }
    (*val1) = (*val1)*dx;
    (*val2) = (*val2)*dx;
}

///////////////////////////////////////////////////////////////////////////////

void B12(double mu, double sigma, int N, double *val1, double* val2) 
{    
    double x = 0.0;
    double f = 0.0;
    double L = 0.0;
    double R = 0.0;
    if (mu<0) {
        B12( -mu, sigma, N, val1, val2);
        (*val1) = 1.0 - (*val1);
    } else {
        StartingValue(mu, sigma, &x, &f);
        ModeAndEffSupp(mu, sigma, &x, &f, &L, &R); 
        B12_trapint(mu, sigma, N, &x, &f, &L, &R, val1, val2);    
    }
} 

///////////////////////////////////////////////////////////////////////////////

void R_B12(double *mu, double *sigma, double *d, int *N, double *val1, double* val2) 
{
    B12(*mu, *sigma, (*N), val1, val2);
} 

///////////////////////////////////////////////////////////////////////////////

void R_vB12(double *mu, double *sigma, int *N, double *val1, double* val2,  int *n) 
{
    int i;
    double v1;
    double v2;
    for (i=0;i<(*n);i++) {
        B12( mu[i], sigma[i], (*N), &v1, &v2);
        val1[i] = v1;  
        val2[i] = v2;     
    }
} 

///////////////////////////////////////////////////////////////////////////////

void B012_trapint(double mu, double sigma, int N, double *x, double *f, double *L, double *R, double* val0, double* val1, double* val2) 
{
    double Lval = (*L);
    double Rval = (*R);
    int i;
    double dx = (Rval - Lval)/((double) N);
    double qL  = dsnorm(Lval,mu,sigma,0);
    double fL0 = b0_safe(Lval)*qL;
    double fL1 = b1_safe(Lval)*qL;
    double fL2 = b2_safe(Lval)*qL;
    double qR  = dsnorm(Rval,mu,sigma,0);
    double fR0 = b0_safe(Rval)*qR;
    double fR1 = b1_safe(Rval)*qR;
    double fR2 = b2_safe(Rval)*qR;
    (*val0) = 0.5*(fL0+fR0);
    (*val1) = 0.5*(fL1+fR1);
    (*val2) = 0.5*(fL2+fR2);    
    for (i=0;i<(N-1);i++) {
        Lval = Lval + dx;
				qL  = dsnorm(Lval,mu,sigma,0);
        (*val0) = (*val0) + b0_safe(Lval)*qL;
        (*val1) = (*val1) + b1_safe(Lval)*qL;        
        (*val2) = (*val2) + b2_safe(Lval)*qL;    
    }
    (*val0) = (*val0)*dx;    
    (*val1) = (*val1)*dx;
    (*val2) = (*val2)*dx;
}

///////////////////////////////////////////////////////////////////////////////

void B012(double mu, double sigma, int N, double *val0, double* val1, double* val2) 
{    
    double x = 0.0;
    double f = 0.0;
    double L = 0.0;
    double R = 0.0;
    if (mu<0) {
        B012_trapint( -mu, sigma, N, &x, &f, &L, &R, val0, val1, val2);       
        (*val0) = mu + (*val0);
        (*val1) = 1.0 - (*val1);
    } else {
        StartingValue(mu, sigma, &x, &f);
        ModeAndEffSupp(mu, sigma, &x, &f, &L, &R); 
        B012_trapint(mu, sigma, N, &x, &f, &L, &R, val0, val1, val2);    
    }
} 

///////////////////////////////////////////////////////////////////////////////

void R_B012(double *mu, double *sigma, int *N, double *val0, double* val1, double* val2) 
{
    B012(*mu, *sigma, (*N), val0, val1, val2);
} 

///////////////////////////////////////////////////////////////////////////////

void R_vB012(double *mu, double *sigma, int *N, double *val0, double* val1, double* val2,  int *n) 
{
    int i;
    double v0;
    double v1;
    double v2;    
    for (i=0;i<(*n);i++) {
        B012(*mu, *sigma, (*N), &v0, &v1, &v2);
        val0[i] = v0;  
        val1[i] = v1;    
        val2[i] = v2;          
    }
} 

