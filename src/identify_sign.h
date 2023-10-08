#ifndef _IDENTIFY_SIGN_H_
#define _IDENTIFY_SIGN_H_

#include <RcppArmadillo.h>

Rcpp::List identify_sign_cpp(Rcpp::List posterior, const arma::mat& sign);

#endif  // _IDENTIFY_SIGN_H_
