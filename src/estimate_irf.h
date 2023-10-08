#ifndef _ESTIMATE_IRF_H_
#define _ESTIMATE_IRF_H_

#include <RcppArmadillo.h>

Rcpp::List estimate_irf_cpp(Rcpp::List posterior, const int& periods);

#endif  // _ESTIMATE_IRF_H_
