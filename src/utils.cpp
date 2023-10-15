#include <RcppArmadillo.h>

using namespace arma;

//' Random draw from matrix normal MN(M, U, V)
//'
//' @param M mat n by m, mean matrix
//' @param U mat n by n, covariance matrix of each column
//' @param V mat m by m, covariance matrix of each row
//'
// [[Rcpp:interface(cpp)]]
arma::mat matnrnd_cpp(const arma::mat& M,
                      const arma::mat& U,
                      const arma::mat& V) {
  // vec(A*B*C) = kron(A,C')*vec(B)
  // a = chol(A) => a'*a = A, i.e. a' is lower triangular

  mat X = mat(size(M), fill::randn);
  return M + chol(U).t() * X * chol(V);
}

//' Random draw from multivariate normal N(mu, Sigma) given inverse Sigma
//'
//' @param mu mean vector
//' @param inv_Sigma precision matrix, inverse of Sigma
//'
// [[Rcpp:interface(cpp)]]
arma::mat mvnrnd_inverse_cpp(const arma::mat& mu, const arma::mat& inv_Sigma) {
  mat q = chol(inv_Sigma);
  mat z = mat(size(mu), fill::randn);
  return mu + inv(q) * z;
}

//' Random draw from Wishart distribution W(S, v) given inverse S
//'
//' @param inv_S inverse of scale matrix S
//' @param v degrees of freedom
//'
// [[Rcpp:interface(cpp)]]
arma::mat wishrnd_inverse_cpp(const arma::mat& inv_S, const int& v) {
  mat A = mat(size(inv_S), fill::randn);
  mat L = chol(inv_S).t();
  mat B = L * A;
  return B * B.t() / v;
}

//' Produce companion form
//'
//' @param A vertically stacked autoregressive coefficients
//'
// [[Rcpp:interface(cpp)]]
arma::mat companion_cpp(const arma::mat& A) {
  int mp = A.n_rows;
  int m = A.n_cols;

  return join_vert(A.t(), eye(mp - m, mp));
}

//' QR decomposition, where the diagonal elements of R are positive
//'
// [[Rcpp:interface(cpp)]]
arma::mat qr_sign_cpp(const arma::mat& A) {
  arma::mat Q, R;
  arma::qr_econ(Q, R, A);

  // Check and modify the diagonal elements of R
  for(arma::uword i = 0; i < R.n_cols; ++i) {
    if(R(i, i) < 0) {
      R.col(i) *= -1;  // Change sign of the column
      Q.row(i) *= -1;  // Change sign of the corresponding row in Q
    }
  }

  return Q;
}
