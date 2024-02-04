#include <RcppArmadillo.h>

#include "utils.h"

using namespace arma;

// Estimate one particular draw of impulse response function
//
// [[Rcpp:interface(cpp)]]
arma::cube irf_cpp(const arma::mat& A, const arma::mat& B, const int& periods) {
  int N = A.n_cols;
  int K = A.n_rows;
  int k = (K - 1) / N;

  cube irf = zeros(N, N, periods);
  irf.slice(0) = B;

  for(int t = 2; t <= periods; t++) {
    for(int j = 1; j <= std::min(k, t - 1); j++) {
      mat A_j = A.rows(1 + (j - 1) * N, j * N).t();
      irf.slice(t - 1) += A_j * irf.slice(t - j - 1);
    }
  }

  return irf;
}

// Estimate impulse response function
//
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
Rcpp::List irf_cpp(Rcpp::List posterior, const int& periods) {
  cube A = posterior["A"];
  cube B = posterior["B"];

  int N = A.n_cols;
  int S = A.n_slices;

  std::cout << "#############################" << std::endl;
  std::cout << "## Computing IRF for draws ##" << std::endl;
  int progress = floor(S / 10);

  arma::field<arma::cube> irf(S);
  for(int s = 0; s < S; s++) {
    irf(s) = irf_cpp(A.slice(s), B.slice(s), periods);

    if(progress == 0 || (s + 1) % progress == 0) {
      std::cout << "## Progress: [" << s + 1 << "/" << S << "]\r";
    }
  }
  std::cout << "\n########### Done! ###########" << std::endl;

  posterior["irf"] = irf;
  return posterior;
}

// Estimate historical decomposition
//
// [[Rcpp:interface(cpp)]]
arma::mat historical_decomposition_cpp(arma::cube irf,
                                       arma::mat E,
                                       int var_i,
                                       int start,
                                       int end) {
  int periods = end - start + 1;
  int N = E.n_cols;

  mat hd = zeros(periods, N);
  for(int shock_j = 0; shock_j < N; shock_j++) {
    for(int h = 0; h < periods; h++) {
      for(int l = 0; l <= h; l++) {
        hd(h, shock_j) +=
            irf(var_i - 1, shock_j, l) * E(start + h - l - 1, shock_j);
      }
    }
  }

  return hd;
}

// Estimate historical decomposition
//
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
Rcpp::List historical_decomposition_cpp(Rcpp::List posterior,
                                        int var_i,
                                        int start,
                                        int end) {
  cube A = posterior["A"];
  cube B = posterior["B"];
  cube U = posterior["U"];

  int T = U.n_rows;
  int N = A.n_cols;
  int S = A.n_slices;
  int periods = end - start + 1;

  std::cout << "#############################" << std::endl;
  std::cout << "## Computing hist. decomp. ##" << std::endl;
  int progress = floor(S / 10);

  cube E = zeros(T, N, S);
  cube hd = zeros(periods, N, S);
  for(int s = 0; s < S; s++) {
    cube irf = irf_cpp(A.slice(s), B.slice(s), periods);
    mat E_s = solve(B.slice(s), U.slice(s).t()).t();
    E.slice(s) = E_s;
    hd.slice(s) = historical_decomposition_cpp(irf, E_s, var_i, start, end);

    if(progress == 0 || (s + 1) % progress == 0) {
      std::cout << "## Progress: [" << s + 1 << "/" << S << "]\r";
    }
  }
  std::cout << "\n########### Done! ###########" << std::endl;

  posterior["E"] = E;
  posterior["hd"] = hd;
  return posterior;
}
