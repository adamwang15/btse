#include <RcppArmadillo.h>

#include "utils.h"

using namespace arma;

//' Estimate impulse response function
//'
//' @name estimate_irf_cpp
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export(.estimate_irf_cpp)]]
Rcpp::List estimate_irf_cpp(Rcpp::List posterior,
                            const int& periods,
                            const bool& structural = true) {
  cube A = posterior["A"];
  cube B = posterior["B"];
  int N = A.n_cols;
  int K = A.n_rows;
  int S = A.n_slices;

  cube irf = cube(N, N * periods, S);
  for(int s = 0; s < A.n_slices; s++) {
    // companion form for MA(infty) representation
    mat A_companion = companion_cpp(A.slice(s).rows(1, K - 1));

    for(int t = 0; t < periods; t++) {
      // slicing starts with 0, both first and last indices included
      // IRF = A_companion^j_[1:N,1:N] * B
      mat A_t = powmat(A_companion, t);
      A_t = A_t.submat(0, 0, N - 1, N - 1);
      mat AB;
      if(structural) {
        AB = A_t * B.slice(s);
      } else {
        AB = A_t;
      }
      irf.slice(s).cols(t * N, (t + 1) * N - 1) = AB;
    }
  }

  posterior["irf"] = irf;
  return posterior;
}
