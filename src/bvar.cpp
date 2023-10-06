#include <RcppArmadillo.h>

#include "blm.h"

using namespace arma;

//' Bayesian vector autoregression
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
Rcpp::List bvar_cpp(const arma::mat Y,
                    const int k,
                    const int S,
                    const Rcpp::List prior) {
  // k = number of lags
  // construct X from Y
  int T = Y.n_rows;
  int M = Y.n_cols;
  mat X = zeros(T - k, k * M);

  for(int i = 0; i < k; i++) {
    mat Y_lag = Y.rows(i + 1, T - k + i);
    for(int j = 0; j < M; j++) {
      X.col(i * M + j) = Y_lag.col(j);
    }
  }

  return blm_cpp(Y.rows(0, T - k - 1), X, S, prior);
}
