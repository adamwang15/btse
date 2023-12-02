#include <RcppArmadillo.h>

#include "utils.h"

using namespace arma;

//' Kalman filter for linear Gaussian state space model
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
Rcpp::List kalman_filter(const arma::mat& Y,
                         const arma::cube& Z,
                         const arma::cube& T,
                         const arma::cube& SHS,
                         const arma::cube& RQR,
                         const arma::mat& a_0,
                         const arma::mat& P_0,
                         const bool& smooth = false) {
  int n = Y.n_rows;
  int d = a_0.n_rows;
  mat a_tt1 = a_0;
  mat P_tt1 = P_0;

  mat T_t = T.slice(0);
  mat Z_t = Z.slice(0);
  mat SHS_t = SHS.slice(0);
  mat RQR_t = RQR.slice(0);

  mat aa_tt = zeros(d, n);
  mat aa_tt1 = zeros(d, n);
  cube PP_tt = zeros(d, d, n);
  cube PP_tt1 = zeros(d, d, n);

  double log_likelihood = 0;
  mat ZP, y_tt1, F_t, y_t, L, inv_F_t, ZPF, a_tt, P_tt;
  for(int t = 0; t < n; t++) {
    if(T.n_slices > 1) {
      T_t = T.slice(t);
    }
    if(Z.n_slices > 1) {
      Z_t = Z.slice(t);
    }
    if(SHS.n_slices > 1) {
      SHS_t = SHS.slice(t);
    }
    if(RQR.n_slices > 1) {
      RQR_t = RQR.slice(t);
    }

    // step 2: filter
    ZP = Z_t * P_tt1;
    y_tt1 = Z_t * a_tt1;
    F_t = ZP * Z_t.t() + SHS_t;

    // step 3: update
    y_t = Y.row(t).t();
    L = chol(F_t).t();
    inv_F_t = inv_chol_cpp(L);
    ZPF = ZP.t() * inv_F_t;
    a_tt = a_tt1 + ZPF * (y_t - y_tt1);
    P_tt = P_tt1 - ZPF * ZP;

    // store results
    log_likelihood += log_mvnpdf_cpp(y_t, y_tt1, inv_F_t, L);
    aa_tt.col(t) = a_tt;
    aa_tt1.col(t) = a_tt1;
    PP_tt.slice(t) = P_tt;
    PP_tt1.slice(t) = P_tt1;

    // step 1: predict (for next loop)
    a_tt1 = T_t * a_tt;
    P_tt1 = T_t * P_tt * T_t.t() + RQR_t;
  }

  Rcpp::List result;
  result["log_likelihood"] = log_likelihood;

  if(!smooth) {
    return result;
  }

  // Kalman smoother
  mat aa_tT = zeros(d, n);
  aa_tT.col(n - 1) = aa_tt.col(n - 1);
  T_t = T.slice(0);
  for(int t = n - 2; t >= 0; t--) {
    if(T.n_slices > 1) {
      T_t = T.slice(t + 1);
    }
    aa_tT.col(t) =
        aa_tt.col(t) +
        PP_tt.slice(t) * T_t.t() *
            solve(PP_tt1.slice(t + 1), aa_tT.col(t + 1) - aa_tt1.col(t + 1));
  }
  result["a_tT"] = aa_tT.t();

  return result;
}
