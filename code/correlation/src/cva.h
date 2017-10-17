// cva.h

#pragma once

double maruyama(const int n, const int p, const double R2, int p_gamma);
double BIC(const int n, const int p, double R2, int vp_gamma);
double ZE(const int n, const int p, double R2, int p_gamma);
double log_hyperg_2F1(double b, double c, double x);
double log_hyperg_2F1_naive(double b, double c, double x);
double liang_g1(const int n, const int p, double R2, int p_gamma);
double liang_g2(const int n, const int p, double R2, int p_gamma);
double liang_g3(const int n, const int p, double R2, int p_gamma);
double robust_bayarri1(const int n, const int p, double R2, int p_gamma);
double robust_bayarri2(const int n, const int p, double R2, int p_gamma);
