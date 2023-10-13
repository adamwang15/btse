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
                    const int& k,
                    const int& S,
                    const Rcpp::List& prior) {
  // k = number of lags
  // construct X from Y
  int T = Y.n_rows;
  int N = Y.n_cols;
  int Nk = N * k;
  mat X = zeros(T - k, N * k);

  for(int i = 0; i < k; i++) {
    mat Y_lag = Y.rows(i + 1, T - k + i);
    for(int j = 0; j < N; j++) {
      X.col(i * N + j) = Y_lag.col(j);
    }
  }
  X = join_horiz(ones(T - k, 1), X);

  return blm_cpp(Y.rows(0, T - k - 1), X, S, prior);

  Rcpp::List posterior = blm_cpp(Y.rows(0, T - k - 1), X, S, prior);
  cube A = posterior["A"];
  cube A_noconstant = zeros(Nk, N, S);
  for(int s = 0; s < S; s++) {
    A_noconstant.slice(s) = A.slice(s).rows(1, Nk);
  }
  posterior["A"] = A_noconstant;
  return posterior;
}
