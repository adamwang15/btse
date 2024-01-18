#ifndef _BVAR_ANALYSIS_H_
#define _BVAR_ANALYSIS_H_

#include <RcppArmadillo.h>

arma::cube irf_cpp(const arma::mat& A, const arma::mat& B, const int& periods);

arma::mat historical_decomposition_cpp(arma::cube irf,
                                       arma::mat E,
                                       int var_i,
                                       int start,
                                       int end);

#endif  // _BVAR_ANALYSIS_H_
