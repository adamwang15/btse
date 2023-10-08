#include <RcppArmadillo.h>

using namespace arma;

//' Short-run restriction identification
//'
//' Choleski decomposition with ordering
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export(.identify_shortrun_cpp)]]
Rcpp::List identify_shortrun_cpp(Rcpp::List posterior) {
  cube Sigma = posterior["Sigma"];
  cube P = cube(size(Sigma), fill::zeros);
  for(int s = 0; s < Sigma.n_slices; s++) {
    P.slice(s) = chol(Sigma.slice(s)).t(); // Choleski decomposition of Sigma
  }
  posterior["PQ"] = P;  // Q=I
  return posterior;
}
