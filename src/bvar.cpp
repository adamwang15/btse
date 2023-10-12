#include <RcppArmadillo.h>

#include "blm.h"

using namespace arma;

//' Bayesian vector autoregression
//'
//' @name bvar_cpp
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export(.bvar_cpp)]]
Rcpp::List bvar_cpp(const arma::mat& Y,
                    const int& p,
                    const int& S,
                    const Rcpp::List& prior) {
  // p = number of lags
  // construct X from Y
  int n = Y.n_rows;
  int m = Y.n_cols;
  int mp = m * p;
  mat X = zeros(n - p, m * p);

  for(int i = 0; i < p; i++) {
    mat Y_lag = Y.rows(i + 1, n - p + i);
    for(int j = 0; j < m; j++) {
      X.col(i * m + j) = Y_lag.col(j);
    }
  }
  X = join_horiz(X, ones(n - p, 1));

  Rcpp::List posterior = blm_cpp(Y.rows(0, n - p - 1), X, S, prior);
  // cube A = posterior["A"];
  // cube A_noconstant = zeros(mp, m, S);
  // for(int s = 0; s < S; s++) {
  //   A_noconstant.slice(s) = A.slice(s).rows(1, mp);
  // }
  // posterior["A"] = A_noconstant;
  return posterior;
}
