#include <RcppArmadillo.h>

using namespace arma;

// Random draw from matrix normal MN(M, U, V)
//
// @param M mat n by m, mean matrix
// @param U mat n by n, covariance matrix of each column
// @param V mat m by m, covariance matrix of each row
//
// [[Rcpp:interface(cpp)]]
arma::mat matnrnd_cpp(const arma::mat& M,
                      const arma::mat& U,
                      const arma::mat& V) {
  // vec(A*B*C) = kron(A,C')*vec(B)
  // a = chol(A) => a'*a = A, i.e. a' is lower triangular

  mat X = mat(size(M), fill::randn);
  return M + chol(U).t() * X * chol(V);
}

// Random draw from multivariate normal N(mu, Sigma) given inverse Sigma
//
// @param mu mean vector
// @param inv_Sigma precision matrix, inverse of Sigma
//
// [[Rcpp:interface(cpp)]]
arma::mat mvnrnd_inverse_cpp(const arma::mat& mu, const arma::mat& inv_Sigma) {
  mat q = chol(inv_Sigma);
  mat z = mat(size(mu), fill::randn);
  return mu + inv(q) * z;
}

double log_det_cpp(const arma::mat& L) {
  return 2 * sum(log(diagvec(L)));
}

arma::mat inv_chol_cpp(const arma::mat& L) {
  return solve(trimatu(L.t()), solve(trimatl(L), eye(size(L))));
}

double log_mvnpdf_cpp(const arma::mat& x,
                      const arma::mat& mu,
                      const arma::mat& inv_Sigma,
                      const arma::mat& L) {
  double log_det = log_det_cpp(L);
  double log_const = -0.5 * (x.n_rows * log(2 * datum::pi) + log_det);
  mat z = x - mu;
  return log_const - 0.5 * as_scalar(z.t() * inv_Sigma * z);
}

// QR decomposition, where the diagonal elements of R are positive
//
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
