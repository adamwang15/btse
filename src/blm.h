#ifndef _BLM_H_
#define _BLM_H_

#include <RcppArmadillo.h>

Rcpp::List blm_conjugate_cpp(const arma::mat& Y,
               const arma::mat& X,
               const int& S,
               const Rcpp::List& prior);

Rcpp::List blm_independent_cpp(const arma::mat& Y,
                   const arma::mat& X,
                   const int& S,
                   const int& burn,
                   const int& thin,
                   const Rcpp::List& prior);

#endif  // _BLM_H_
