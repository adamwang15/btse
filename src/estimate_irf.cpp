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
  cube PQ = posterior["PQ"];
  int m = A.n_cols;
  int S = A.n_slices;

  cube irf = cube(m, m * periods, S);
  for(int s = 0; s < A.n_slices; s++) {
    // companion form for MA(infty) representation
    mat A_companion = companion_cpp(A.slice(s));
    mat A_products = eye(size(A_companion));

    for(int t = 0; t < periods; t++) {
      // slicing starts with 0, both first and last indices included
      // IRF = A_companion^j_[1:m,1:m] * C * Q
      irf.slice(s).cols(t * m, (t + 1) * m - 1) =
          A_products.submat(0, 0, m - 1, m - 1) * PQ.slice(s);
      A_products = A_products * A_companion;
    }
  }

  posterior["IRF"] = irf;
  return posterior;
}
