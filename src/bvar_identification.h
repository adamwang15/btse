#ifndef _BVAR_IDENTIFICATION_H_
#define _BVAR_IDENTIFICATION_H_

#include <RcppArmadillo.h>

Rcpp::List identify_shortrun(Rcpp::List posterior);

Rcpp::List identify_longrun(Rcpp::List posterior);

Rcpp::List identify_sign(Rcpp::List posterior, const arma::mat& sign);

#endif  // _BVAR_IDENTIFICATION_H_
