# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @importFrom Rcpp evalCpp
#' @useDynLib correlation
NULL

#' Calculate the correlation of all the sub-models of mX and vy
#'
#' @param vy Vector of responses
#' @param mX Covariate matrix
#' @param intercept_col The index of the column in mX containing the intercept, if any
#' @param bIntercept Logical value indicating whether there is an intercept column or not
#' @param bCentre Logical value indicating whether to centre the response vector and covariance matrix or not
#' @param cores The number of cores to use
#' @return The vector of correlations
#' @export
all_correlations_mX <- function(vy, mX, intercept_col = 1L, bIntercept = FALSE, bCentre = FALSE, cores = 1L) {
    .Call('correlation_all_correlations_mX', PACKAGE = 'correlation', vy, mX, intercept_col, bIntercept, bCentre, cores)
}

#' Calculate the correlation of all the sub-models mX/mZ and vy, where mX is fixed in every model and the sub-models of mZ are included
#'
#' @param vy Vector of responses
#' @param mX Fixed covariate matrix
#' @param mZ Varying covariate matrix
#' @param intercept_col The index of the column in mX containing the intercept, if any
#' @param bIntercept Logical value indicating whether there is an intercept column or not
#' @param bCentre Logical value indicating whether to centre the response vector and covariance matrix or not
#' @param cores The number of cores to use
#' @return The vector of correlations
#' @export
all_correlations_mX_mZ <- function(vy, mX, mZ, intercept_col = 1L, bIntercept = FALSE, bCentre = FALSE, cores = 1L) {
    .Call('correlation_all_correlations_mX_mZ', PACKAGE = 'correlation', vy, mX, mZ, intercept_col, bIntercept, bCentre, cores)
}

#' Return the graycode matrix
#'
#' @param varying The number of covariates varying in the graycode matrix
#' @param fixed The number of fixed covariates in the graycode matrix. These covariates will always be included
#' @return The graycode matrix. The number of fixed columns will be included in the lower indexed columns
#' as 1s, while the higher indexed columns will varying depending on whether each covariate in the varying
#' set of covariates is included or not.
#' @export
graycode <- function(varying, fixed = 0L) {
    .Call('correlation_graycode', PACKAGE = 'correlation', varying, fixed)
}

