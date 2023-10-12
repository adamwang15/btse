#include <RcppArmadillo.h>

#include "utils.h"

using namespace arma;

//' Estimate impulse response function
//'
//' @name estimate_irf_cpp
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export(.estimate_irf_cpp)]]
Rcpp::List estimate_irf_cpp(Rcpp::List posterior, const int& periods) {
  cube A = posterior["A"];
  cube B = posterior["B"];
  int m = A.n_cols;
  int S = A.n_slices;

  cube irf = cube(m, m * periods, S);
  for(int s = 0; s < A.n_slices; s++) {
    // companion form for MA(infty) representation
    mat A_companion = companion_cpp(A.slice(s));

    for(int t = 0; t < periods; t++) {
      // slicing starts with 0, both first and last indices included
      // IRF = A_companion^j_[1:m,1:m] * P * Q
      mat A_t = powmat(A_companion, t);
      irf.slice(s).cols(t * m, (t + 1) * m - 1) =
          A_t.submat(0, 0, m - 1, m - 1) * B.slice(s);
    }
  }

  posterior["IRF"] = irf;
  return posterior;
}
