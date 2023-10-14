#include <RcppArmadillo.h>

#include "utils.h"

using namespace arma;

//' Bayesian linear model with natural conjugate prior
//'
//' A|Sigma ~ MN(A_l, V_l, Sigma), Sigma ~ IW(S_l, v_l)
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
Rcpp::List blm_conjugate_cpp(const arma::mat& Y,
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
  S_u = symmatu(S_u);
  int v_u = v_l + T;

  // draws
  cube A = zeros(K, N, S);
  cube Sigma = zeros(N, N, S);
  mat A_draw = zeros(K, N);
  mat Sigma_draw = zeros(N, N);

  int progress = S / 10;
  for(int s = 0; s < S; s++) {
    Sigma_draw = iwishrnd(S_u, v_u);
    A_draw = matnrnd_cpp(A_u, V_u, Sigma_draw);  // conditional posterior of A
    Sigma.slice(s) = Sigma_draw;
    A.slice(s) = A_draw;

    if(s % progress == 0) {
      std::cout << "MCMC draws: [" << s + 1 << " / " << S << "]" << std::endl;
    }
  }

  Rcpp::List draws;
  draws["A"] = A;
  draws["Sigma"] = Sigma;
  return draws;
}

//' Bayesian linear model with independent prior
//'
//' vec_A ~ N(vec_A_l, V_l), Sigma ~ IW(S_l, v_l)
//' y|vec_A,Sigma ~ N(Z*vec_A, I_T \otimes Sigma)
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
Rcpp::List blm_independent_cpp(const arma::mat& Y,
                               const arma::mat& X,
                               const int& S,
                               const int& burn,
                               const int& thin,
                               const Rcpp::List& prior) {
  int T = Y.n_rows;
  int N = Y.n_cols;
  int K = X.n_cols;

  // convert representation to y = Z*vec(A) + eps
  mat y = vectorise(Y.t());
  mat Z = zeros(T * N, K * N);
  for(int t = 0; t < T; t++) {
    Z.rows(t * N, (t + 1) * N - 1) = kron(eye(N, N), X.row(t));
  }

  // priors
  mat A_l = prior["A_l"];
  mat V_l = prior["V_l"];
  mat S_l = prior["S_l"];
  int v_l = prior["v_l"];
  mat vec_A_l = vectorise(A_l);
  mat inv_V_l = inv_sympd(V_l);

  // analytic solutions
  mat IVVA = inv_V_l * vec_A_l;
  int v_u = v_l + T;
  mat I_T = eye(T, T);

  // MCMC draws
  int S_thin = (S - burn) / thin;
  int s_thin = 0;
  cube A = zeros(K, N, S_thin);
  cube Sigma = zeros(N, N, S_thin);
  mat vec_A_draw = solve(Z.t() * Z, Z.t() * y);

  // declarations
  mat S_u, y_t, Z_t, u_t, Sigma_draw, inv_Sigma_draw, inv_V_u, L, vec_A_u;
  sp_mat ZKIS, KIS;

  for(int s = 0; s < S; s++) {
    // conditional posterior of Sigma given vec(A)
    S_u = S_l;
    for(int t = 0; t < T; t++) {
      y_t = y.rows(t * N, (t + 1) * N - 1);
      Z_t = Z.rows(t * N, (t + 1) * N - 1);
      u_t = y_t - Z_t * vec_A_draw;
      S_u += u_t * u_t.t();
    }

    Sigma_draw = iwishrnd(S_u, v_u);
    inv_Sigma_draw = inv_sympd(Sigma_draw);

    // conditional posterior of vec(A) given Sigma
    KIS = kron(I_T, inv_Sigma_draw);
    ZKIS = Z.t() * KIS;
    inv_V_u = inv_V_l + ZKIS * Z;
    L = chol(inv_V_u).t();
    vec_A_u = solve(L.t(), solve(L, IVVA + ZKIS * y));
    vec_A_draw = mvnrnd_inverse_cpp(vec_A_u, inv_V_u);

    if(s % thin == 0) {
      std::cout << "MCMC draws: [" << s + 1 << " / " << S << "]" << std::endl;

      if(s >= burn) {
        A.slice(s_thin) = reshape(vec_A_draw, K, N);
        Sigma.slice(s_thin) = Sigma_draw;
        s_thin++;
      }
    }
  }

  Rcpp::List posterior;
  posterior["A"] = A;
  posterior["Sigma"] = Sigma;
  return posterior;
}
