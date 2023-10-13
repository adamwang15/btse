#include <RcppArmadillo.h>

#include "utils.h"

using namespace arma;

//' Bayesian linear model
//'
//' Bayesian linear regression with normal-inverse-Wishart conjugate prior
//' prior: A ~ MN(A_l, V_l, Sigma), Sigma ~ IW(S_l, v_l)
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
  int T = Y.n_rows;
  int N = Y.n_cols;
  int K = X.n_cols;

  // priors
  mat A_l = prior["A_l"];
  mat V_l = prior["V_l"];
  mat S_l = prior["S_l"];
  int v_l = prior["v_l"];

  // analytic solutions
  mat inv_V_l = inv_sympd(V_l);
  mat inv_V_u = inv_V_l + X.t() * X;
  mat V_u = inv_sympd(inv_V_u);
  mat A_u = V_u * (inv_V_l * A_l + X.t() * Y);

  // marginal posterior of Sigma
  mat S_u = S_l + Y.t() * Y + A_l.t() * inv_V_l * A_l - A_u.t() * inv_V_u * A_u;
  S_u = eye(N, N);
  int v_u = v_l + T;

  // draws
  cube A_draws = zeros(K, N, S);
  cube Sigma_draws = zeros(N, N, S);
  mat A = zeros(K, N);
  mat Sigma = zeros(N, N);
  for(int s = 0; s < S; s++) {
    Sigma = iwishrnd(S_u, v_u);
    A = matnrnd_cpp(A_u, V_u, Sigma);  // conditional posterior of A
    Sigma_draws.slice(s) = Sigma;
    A_draws.slice(s) = A;
  }

  Rcpp::List draws;
  draws["A"] = A_draws;
  draws["Sigma"] = Sigma_draws;
  return draws;
}
