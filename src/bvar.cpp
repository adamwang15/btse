#include <RcppArmadillo.h>

#include "blm.h"

using namespace arma;

//' Bayesian vector autoregression
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export(.bvar_cpp)]]
Rcpp::List bvar_cpp(const arma::mat& Y,
                    const int& p,
                    const int& S,
                    const Rcpp::List& prior) {
  // p = number of lags
  // construct X from Y
  int T = Y.n_rows;
  int m = Y.n_cols;
  mat X = zeros(T - p, m * p);

  for(int i = 0; i < p; i++) {
    mat Y_lag = Y.rows(i + 1, T - p + i);
    for(int j = 0; j < m; j++) {
      X.col(i * m + j) = Y_lag.col(j);
    }
  }

  return blm_cpp(Y.rows(0, T - p - 1), X, S, prior);
}
