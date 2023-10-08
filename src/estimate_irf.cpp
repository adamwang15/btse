#include <RcppArmadillo.h>

#include "utils.h"

using namespace arma;

//' Estimate impulse response function
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export(.estimate_irf_cpp)]]
Rcpp::List estimate_irf_cpp(Rcpp::List posterior, const int& periods) {
  cube B = posterior["B"];
  cube PQ = posterior["PQ"];
  int m = B.n_rows;
  int S = B.n_slices;

  cube irf = cube(m, m * periods, S);
  for(int s = 0; s < B.n_slices; s++) {
    // companion form for MA(infty) representation
    mat A = companion_cpp(B.slice(s));
    mat A_products = eye(size(A));

    for(int t = 0; t < periods; t++) {
      // slicing starts with 0, both first and last indices included
      // IRF = A^j_[1:m,1:m] * C * Q
      irf.slice(s).cols(t * m, (t + 1) * m - 1) =
          A_products.submat(0, 0, m - 1, m - 1) * PQ.slice(s);
      // irf.slice(s).cols(t * m, (t + 1) * m - 1) =
      //     A_products.submat(0, 0, m - 1, m - 1) * eye(m,m);
      A_products = A_products * A;
    }
  }

  posterior["IRF"] = irf;
  return posterior;
}
