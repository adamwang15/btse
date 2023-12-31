# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Bayesian linear model with natural conjugate prior
#'
#' A|Sigma ~ MN(A_l, V_l, Sigma), Sigma ~ IW(S_l, v_l)
#'
NULL

#' Bayesian linear model with independent prior
#'
#' vec_A ~ N(vec_A_l, V_l), Sigma ~ IW(S_l, v_l)
#' y|vec_A,Sigma ~ N(Z*vec_A, I_T \otimes Sigma)
#'
NULL

blm_conjugate_cpp <- function(Y, X, S, prior) {
    .Call(`_btse_blm_conjugate_cpp`, Y, X, S, prior)
}

blm_independent_cpp <- function(Y, X, S, burn, thin, prior) {
    .Call(`_btse_blm_independent_cpp`, Y, X, S, burn, thin, prior)
}

#' Bayesian vector autoregression
#'
NULL

.bvar_cpp <- function(Y, k, model, S, burn, thin, prior) {
    .Call(`_btse_bvar_cpp`, Y, k, model, S, burn, thin, prior)
}

#' Short-run restriction identification
#'
#' Choleski decomposition with ordering
#'
NULL

#' Long-run restriction identification
#'
#' Transform VAR(k) to VMA(infty) then compute Choleski decomposition
#'
NULL

#' Sign restriction identification
#'
#' Draw orthonormal matrices Q until the sign restrictions are satisfied
#'
NULL

identify_shortrun <- function(posterior) {
    .Call(`_btse_identify_shortrun`, posterior)
}

identify_longrun <- function(posterior) {
    .Call(`_btse_identify_longrun`, posterior)
}

identify_sign <- function(posterior, sign) {
    .Call(`_btse_identify_sign`, posterior, sign)
}

#' Estimate impulse response function
#'
NULL

.irf_cpp <- function(posterior, periods, structural = TRUE) {
    .Call(`_btse_irf_cpp`, posterior, periods, structural)
}

