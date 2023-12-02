#ifndef _BVAR_ANALYSIS_H_
#define _BVAR_ANALYSIS_H_

#include <RcppArmadillo.h>

Rcpp::List estimate_irf_cpp(Rcpp::List posterior, const int& periods);

#endif  // _BVAR_ANALYSIS_H_
