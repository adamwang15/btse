#ifndef _BVAR_H_
#define _BVAR_H_

#include <RcppArmadillo.h>

Rcpp::List bvar_cpp(const arma::mat& Y,
                    const int& k,
                    const int& S,
                    const Rcpp::List& prior);

#endif  // _BVAR_H_
