#include <RcppArmadillo.h>

#include "blm.h"

using namespace arma;

// Bayesian vector autoregression
//
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export(.bvar_cpp)]]
Rcpp::List bvar_cpp(const arma::mat& Y,
                    const int& k,
                    const std::string& model,
                    const int& S,
                    const int& burn,
                    const int& thin,
                    const Rcpp::List& prior) {
  // k = number of lags
  // construct X from Y
  int T = Y.n_rows;
  int N = Y.n_cols;
  mat X = zeros(T - k, N * k);

  for(int i = 0; i < k; i++) {
    mat Y_lag = Y.rows(k - i - 1, T - i - 2);
    for(int j = 0; j < N; j++) {
      X.col(i * N + j) = Y_lag.col(j);
    }
  }
  X = join_horiz(ones(T - k, 1), X);

  mat Y_head = Y.rows(k, T - 1);
  Rcpp::List posterior;
  if(model == "conjugate") {
    posterior = blm_conjugate_cpp(Y_head, X, S, prior);
  } else {
    posterior = blm_independent_cpp(Y_head, X, S, burn, thin, prior);
  }

  Rcpp::List data;
  data["Y"] = Y_head;
  data["X"] = X;

  posterior["data"] = data;
  posterior["B"] = zeros(1, 1, 1);
  return posterior;
}
