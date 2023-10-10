#include <RcppArmadillo.h>

using namespace arma;

//' Long-run restriction identification
//'
//' Transform VAR(p) to VMA(infty) then compute Choleski decomposition
//'
//' @name identify_long_run_cpp
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export(.identify_longrun_cpp)]]
Rcpp::List identify_longrun_cpp(Rcpp::List posterior) {
  cube A = posterior["A"];
  cube Sigma = posterior["Sigma"];
  cube P = cube(size(Sigma), fill::zeros);

  int m = A.n_cols;
  int p = A.n_rows / m;
  mat I_m = eye(m, m);
  mat I_mp = repmat(I_m, 1, p);  // = [I I ...]

  mat C = mat(m, m, fill::zeros);
  mat D = mat(m, m, fill::zeros);
  mat inv_D = mat(m, m, fill::zeros);

  for(int s = 0; s < Sigma.n_slices; s++) {
    // D = long run MA(infty) coefficients
    // y^LR = (I-A1-A2-...)^-1*u_t = D^-1*u_t = C*e_t
    D = I_m - I_mp * A.slice(s);
    inv_D = inv(D);

    // C = chol(D^-1*Sigma*D^-1')
    // D^-1*u_t = C*e_t => u_t = D*C*e_t
    C = chol(inv_D * Sigma.slice(s) * inv_D.t()).t();
    P.slice(s) = D * C;
  }
  posterior["PQ"] = P;  // Q=I
  return posterior;
}
