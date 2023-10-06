#include <RcppArmadillo.h>

using namespace arma;

//' Random draw from matrix normal
//'
//' Make a draw from MN(M, U, V)
//'
//' @param M mat n by m, mean matrix
//' @param U mat n by n, covariance matrix of each column
//' @param V mat m by m, covariance matrix of each row
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
arma::mat matnrnd_cpp(const arma::mat& M, const arma::mat& U, const arma::mat& V) {
  // vec(A*B*C) = kron(A,C')*vec(B)
  // a = chol(A) => a'*a = A

  mat X = mat(size(M), fill::randn);
  return M + chol(U).t() * X * chol(V);
}
