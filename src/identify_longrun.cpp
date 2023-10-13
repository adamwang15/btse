#include <RcppArmadillo.h>

using namespace arma;

//' Long-run restriction identification
//'
//' Transform VAR(k) to VMA(infty) then compute Choleski decomposition
//'
//' @name identify_long_run_cpp
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export(.identify_longrun_cpp)]]
Rcpp::List identify_longrun_cpp(Rcpp::List posterior) {
  cube A = posterior["A"];
  cube Sigma = posterior["Sigma"];
  cube P = cube(size(Sigma), fill::zeros);

  int N = A.n_cols;
  int K = A.n_rows;
  int k = (K - 1) / N;
  mat I_N = eye(N, N);
  mat I_Nk = repmat(I_N, 1, k);  // = [I I ...]

  mat C = mat(N, N, fill::zeros);
  mat D = mat(N, N, fill::zeros);
  mat inv_D = mat(N, N, fill::zeros);

  for(int s = 0; s < Sigma.n_slices; s++) {
    // D = long run MA(infty) coefficients
    // y^LR = (I-A1-A2-...)^-1*u_t = D^-1*u_t = C*e_t
    D = I_N - I_Nk * A.slice(s).rows(1, K - 1);
    inv_D = inv(D);

    // C = chol(D^-1*Sigma*D^-1')
    // D^-1*u_t = C*e_t => u_t = D*C*e_t i.e. P=D*C
    C = chol(inv_D * Sigma.slice(s) * inv_D.t()).t();
    P.slice(s) = D * C;
  }
  posterior["B"] = P;  // Q=I
  return posterior;
}
