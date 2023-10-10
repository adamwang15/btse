#include <RcppArmadillo.h>

#include "utils.h"

using namespace arma;

//' Bayesian linear model
//'
//' Bayesian linear regression with normal-inverse-Wishart conjugate prior
//' prior: A ~ MN(A_0, V_0, Sigma), Sigma ~ IW(Lambda_0, nu_0)
//'
//' @name blm_cpp
//' @param Y dependent variables
//' @param X independent variables
//' @param S sample size
//' @param prior list of priors
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
Rcpp::List blm_cpp(const arma::mat& Y,
                   const arma::mat& X,
                   const int& S,
                   const Rcpp::List& prior) {
  int n = Y.n_rows;
  int m = Y.n_cols;
  int mp = X.n_cols;

  // priors
  mat A_0 = prior["A_0"];
  mat V_0 = prior["V_0"];
  mat Lambda_0 = prior["Lambda_0"];
  int nu_0 = prior["nu_0"];

  // analytic solutions
  mat inv_V_0 = inv_sympd(V_0);
  mat inv_V_n = inv_V_0 + X.t() * X;
  mat V_n = inv_sympd(inv_V_n);
  mat B_n = V_n * (inv_V_0 * A_0 + X.t() * Y);

  // marginal posterior of Sigma
  mat Lambda_n =
      Lambda_0 + Y.t() * Y + A_0.t() * inv_V_0 * A_0 - B_n.t() * inv_V_n * B_n;
  Lambda_n = eye(m, m);
  int nu_n = nu_0 + n;

  // draws
  cube A_draws = zeros(mp, m, S);
  cube Sigma_draws = zeros(m, m, S);
  mat A = zeros(mp, m);
  mat Sigma = zeros(m, m);
  for(int s = 0; s < S; s++) {
    Sigma = iwishrnd(Lambda_n, nu_n);
    A = matnrnd_cpp(B_n, V_n, Sigma);  // conditional posterior of A
    Sigma_draws.slice(s) = Sigma;
    A_draws.slice(s) = A;
  }

  Rcpp::List draws;
  draws["A"] = A_draws;
  draws["Sigma"] = Sigma_draws;
  return draws;
}
