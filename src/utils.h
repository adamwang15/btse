#ifndef _UTILS_H_
#define _UTILS_H_

#include <RcppArmadillo.h>

arma::mat matnrnd_cpp(const arma::mat& M, const arma::mat& U, const arma::mat& V);

arma::mat mvnrnd_inverse_cpp(const arma::mat& mu, const arma::mat& inv_Sigma);

arma::mat wishrnd_inverse_cpp(const arma::mat& inv_S, const int& v);

arma::mat companion_cpp(const arma::mat& B);

arma::mat qr_sign_cpp(const arma::mat& A);

#endif  // _UTILS_H_
