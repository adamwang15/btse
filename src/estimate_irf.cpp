#include <RcppArmadillo.h>

#include "utils.h"

using namespace arma;

//' Estimate impulse response function
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export(.estimate_irf_cpp)]]
Rcpp::List estimate_irf_cpp(Rcpp::List posterior,
                            const int& periods,
                            const bool& structural = true) {
  cube A = posterior["A"];
  cube B = posterior["B"];

  int N = A.n_cols;
  int K = A.n_rows;
  int S = A.n_slices;
  int k = (K - 1) / N;

  std::cout << "#############################" << std::endl;
  std::cout << "## Computing IRF for draws ##" << std::endl;
  int progress = S / 10;

  cube irf = zeros(N, N * periods, S);
  for(int s = 0; s < S; s++) {
    cube irf_s = zeros(N, N, periods);

    mat B_s;
    if(structural) {
      B_s = B.slice(s);
    } else {
      B_s = eye(N, N);
    }

    irf_s.slice(0) = B_s;
    for(int t = 2; t <= periods; t++) {
      for(int j = 1; j <= std::min(k, t - 1); j++) {
        mat A_j = A.slice(s).rows(1 + (j - 1) * N, j * N).t();
        irf_s.slice(t - 1) += A_j * irf_s.slice(t - j - 1);
      }
    }
    irf_s = reshape(irf_s, N, N * periods, 1);
    irf.slice(s) = irf_s.slice(0);

    if((s + 1) % progress == 0) {
      std::cout << "## Progress: [" << s + 1 << "/" << S << "]\r";
    }
  }
  std::cout << "\n########### Done! ###########" << std::endl;

  posterior["irf"] = irf;
  return posterior;
}
