#include <RcppArmadillo.h>

#include "utils.h"

using namespace arma;

//' Short-run restriction identification
//'
//' Choleski decomposition with ordering
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
Rcpp::List identify_shortrun(Rcpp::List posterior) {
  cube Sigma = posterior["Sigma"];
  cube P = cube(size(Sigma), fill::zeros);
  for(int s = 0; s < Sigma.n_slices; s++) {
    P.slice(s) = chol(Sigma.slice(s)).t();  // Choleski decomposition of Sigma
  }
  posterior["B"] = P;  // Q=I
  return posterior;
}

//' Long-run restriction identification
//'
//' Transform VAR(k) to VMA(infty) then compute Choleski decomposition
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
Rcpp::List identify_longrun(Rcpp::List posterior) {
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

//' Sign restriction identification
//'
//' Draw orthonormal matrices Q until the sign restrictions are satisfied
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
Rcpp::List identify_sign(Rcpp::List posterior, const arma::mat& sign) {
  cube Sigma = posterior["Sigma"];
  cube PQ = cube(size(Sigma), fill::zeros);

  int N = Sigma.n_rows;
  // only restricts first q variables, the rest of B are keep as identity matrix
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
    mat Q = eye(N, N);
    Q.submat(0, 0, q - 1, q - 1) = Q_draw;
    PQ_draw = P * Q;

    PQ.slice(s) = PQ_draw;
  }

  posterior["B"] = PQ;
  return posterior;
}
