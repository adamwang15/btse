#include <RcppArmadillo.h>

using namespace arma;

//' Long-run restriction identification
//'
//' Transform VAR(p) to VMA(infty) then compute Choleski decomposition
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export(.identify_longrun_cpp)]]
Rcpp::List identify_longrun_cpp(Rcpp::List posterior) {
  cube B = posterior["B"];
  cube Sigma = posterior["Sigma"];
  cube P = cube(size(Sigma), fill::zeros);

  int m = B.n_rows;
  int p = B.n_cols / m;
  mat I_m = eye(m, m);
  mat I_mp = repmat(I_m, p, 1); // = [I I ...]'

  mat A = mat(m, m, fill::zeros);
  mat D = mat(m, m, fill::zeros);
  mat inv_D = mat(m, m, fill::zeros);

  for(int s = 0; s < Sigma.n_slices; s++) {
    // D = long run MA(infty) coefficients
    // y^LR = (I-B1-B2-...)^-1*u_t = D*u_t
    inv_D = I_m - B.slice(s) * I_mp;
    D = inv(inv_D);

    // D*u_t = A*e_t => u_t = D^-1*A*e_t
    A = chol(D * Sigma.slice(s) * D.t()).t();
    P.slice(s) = inv_D * A;
  }
  posterior["PQ"] = P;  // Q=I
  return posterior;
}
