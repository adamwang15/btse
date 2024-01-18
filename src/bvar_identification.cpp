#include <RcppArmadillo.h>

#include "bvar_analysis.h"
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

//' If matches traditional sign restriction
//'
// [[Rcpp:interface(cpp)]]
bool traditional_sign_cpp(const arma::mat& B, const arma::mat sign) {
  return accu(((B % sign) > 0)) == accu(abs(sign));
}

//' If matches narrative sign restriction
//'
// [[Rcpp:interface(cpp)]]
bool narrative_sign_cpp(const arma::cube irf,
                        const arma::mat E,
                        const arma::mat sign) {
  bool match = true;

  for(int i = 0; i < sign.n_rows; i++) {
    int type = sign(i, 0);
    int var_i = sign(i, 1);
    int shock_j = sign(i, 2);
    int start = sign(i, 3);
    int end = sign(i, 4);
    int is_greater = sign(i, 5);
    int periods = end - start + 1;

    if(type == 1) {
      // type A
      // check if shock_j is greater (less) than 0
      mat E_j = E.submat(start - 1, shock_j - 1, end - 1, shock_j - 1);
      std::cout << (accu(E_j) < 0) << std::endl;
      if(is_greater) {
        match = accu(E_j > 0) == periods;
      } else {
        match = accu(E_j < 0) == periods;
      }
    } else {
      // type B
      // check if historical decomposition is the largest (smallest)
      mat hd = historical_decomposition_cpp(irf, E, var_i, start, end);
      if(is_greater) {
        match = sum(index_max(abs(hd), 1) == (shock_j - 1)) == periods;
      } else {
        match = sum(index_min(abs(hd), 1) == (shock_j - 1)) == periods;
      }
    }

    if(match == false) {
      break;
    }
  }
  // return true;
  return match;
}

//' Sign restriction identification
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

//' Narrative sign restriction identification
//'
// [[Rcpp:interface(cpp)]]
// [[Rcpp::export]]
Rcpp::List identify_narrative_sign(Rcpp::List posterior,
                                   const arma::mat& traditional,
                                   const arma::mat& narrative) {
  cube Sigma = posterior["Sigma"];
  cube A = posterior["A"];
  cube B = cube(size(Sigma), fill::zeros);
  cube U = posterior["U"];

  int S = Sigma.n_slices;
  int N = Sigma.n_rows;
  int target = accu(abs(traditional));
  int max_periods = max(narrative.col(4) - narrative.col(3)) + 1;

  std::cout << "#############################" << std::endl;
  std::cout << "## Computing nar. sign res.##" << std::endl;
  int progress = floor(S / 10);

  for(int s = 0; s < S; s++) {
    mat P = chol(Sigma.slice(s)).t();  // Choleski decomposition of Sigma

    // draw Q till satisfied
    mat B_s;
    int trial = 1;
    while(true) {
      mat Q = qr_sign_cpp(mat(N, N, fill::randn));
      B_s = P * Q;
      cube irf = irf_cpp(A.slice(s), B_s, max_periods);
      mat E = solve(B_s, U.slice(s).t()).t();

      if(traditional_sign_cpp(B_s, traditional) &&
         narrative_sign_cpp(irf, E, narrative)) {
        break;
      }

      std::cout << "Q: " << trial << "\r";
      trial++;
    }
    B.slice(s) = B_s;

    if(progress == 0 || (s + 1) % progress == 0) {
      std::cout << "\n## Progress: [" << s + 1 << "/" << S << "]\r";
    }
  }

  std::cout << "\n########### Done! ###########" << std::endl;

  posterior["B"] = B;
  return posterior;
}
