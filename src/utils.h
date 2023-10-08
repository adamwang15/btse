#ifndef _UTILS_H_
#define _UTILS_H_

#include <RcppArmadillo.h>

arma::mat matnrnd_cpp(const arma::mat& M, const arma::mat& U, const arma::mat& V);

arma::mat companion_cpp(const arma::mat& B);

arma::mat qr_sign_cpp(const arma::mat& A);

#endif  // _UTILS_H_
