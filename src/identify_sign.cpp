#include <RcppArmadillo.h>

#include "utils.h"

using namespace arma;

//' Sign restriction identification
//'
//' Draw orthonormal matrices Q until the sign restrictions are satisfied
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export(.identify_sign_cpp)]]
Rcpp::List identify_sign_cpp(Rcpp::List posterior, const arma::mat& sign) {
  cube Sigma = posterior["Sigma"];
  cube PQ = cube(size(Sigma), fill::zeros);
  int m = Sigma.n_rows;
  int target = accu(abs(sign));

  for(int s = 0; s < Sigma.n_slices; s++) {
    mat P = chol(Sigma.slice(s)).t();  // Choleski decomposition of Sigma

    // draw Q till satisfied
    mat PQ_draw = zeros(m, m);
    while(accu(((PQ_draw % sign) > 0)) != target) {
      mat Q = qr_sign_cpp(mat(m, m, fill::randn));
      PQ_draw = P * Q;
    }
    PQ.slice(s) = PQ_draw;
  }

  posterior["PQ"] = PQ;
  return posterior;
}
