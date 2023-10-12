#include <RcppArmadillo.h>

#include "utils.h"

using namespace arma;

//' Sign restriction identification
//'
//' Draw orthonormal matrices Q until the sign restrictions are satisfied
//'
//' @name identify_sign_cpp
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export(.identify_sign_cpp)]]
Rcpp::List identify_sign_cpp(Rcpp::List posterior, const arma::mat& sign) {
  cube Sigma = posterior["Sigma"];
  cube PQ = cube(size(Sigma), fill::zeros);
  // cube Q = cube(size(Sigma), fill::zeros);
  int m = Sigma.n_rows;
  int q = sign.n_rows;
  int target = accu(abs(sign));

  for(int s = 0; s < Sigma.n_slices; s++) {
    mat P = chol(Sigma.slice(s)).t();  // Choleski decomposition of Sigma

    // draw Q till satisfied
    mat Q_draw;
    mat PQ_draw = zeros(q, q);
    while(accu(((PQ_draw % sign) > 0)) < target) {
      Q_draw = qr_sign_cpp(mat(q, q, fill::randn));
      PQ_draw = P.submat(0, 0, q - 1, q - 1) * Q_draw;
    }
    mat Q = eye(m, m);
    Q.submat(0, 0, q - 1, q - 1) = Q_draw;
    PQ_draw = P * Q;

    PQ.slice(s) = PQ_draw;
    // Q.slice(s) = Q_draw;
  }

  // posterior["Q"] = Q;
  posterior["B"] = PQ;
  return posterior;
}
