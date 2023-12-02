#ifndef _KALMAN_FILTER_H_
#define _KALMAN_FILTER_H_

#include <RcppArmadillo.h>

Rcpp::List kalman_filter(const arma::mat& Y,
                         const arma::cube& Z,
                         const arma::cube& T,
                         const arma::cube& SS,
                         const arma::cube& RR,
                         const arma::mat& a_0,
                         const arma::mat& P_0,
                         const bool& smooth = false);

#endif  // _KALMAN_FILTER_H_
