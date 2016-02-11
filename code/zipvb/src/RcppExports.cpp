// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "zipvb_types.hpp"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// fastdiag
VectorXd fastdiag(const MapMatd lambda, const MapMatd C);
RcppExport SEXP zipvb_fastdiag(SEXP lambdaSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const MapMatd >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const MapMatd >::type C(CSEXP);
    __result = Rcpp::wrap(fastdiag(lambda, C));
    return __result;
END_RCPP
}
// fastdiag2
VectorXd fastdiag2(const MapMatd R, const MapMatd C);
RcppExport SEXP zipvb_fastdiag2(SEXP RSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const MapMatd >::type R(RSEXP);
    Rcpp::traits::input_parameter< const MapMatd >::type C(CSEXP);
    __result = Rcpp::wrap(fastdiag2(R, C));
    return __result;
END_RCPP
}
// fastsolve
VectorXd fastsolve(const MapMatd R, const MapMatd C);
RcppExport SEXP zipvb_fastsolve(SEXP RSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const MapMatd >::type R(RSEXP);
    Rcpp::traits::input_parameter< const MapMatd >::type C(CSEXP);
    __result = Rcpp::wrap(fastsolve(R, C));
    return __result;
END_RCPP
}