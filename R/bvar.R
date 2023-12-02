#' @title
#' Bayesian Vector Autoregression
#'
#' @description
#' Estimate a homoskedastic Bayesian vector autoregression model with the
#' Normal-Inverse-Wishart conjugate prior.
#'
#' @details
#' Estimate a Bayesian Structural vector autoregression, base on a VAR(p) specification:
#' \deqn{y_t = \mu + A_1 y_{t-1} + \cdots + A_k y_{t-k} + u_t}
#' where \eqn{y_t} and \eqn{u_t} are \eqn{N \times 1}, \eqn{A_j} is \eqn{N \times N}.
#' We can stack the \eqn{A_j}'s vertically
#' such that \eqn{A=(\mu,A_1,\dots,A_k)'} which is \eqn{K \times N} where \eqn{K=1+Nk},
#' and \eqn{y_{t-j}}'s vertically such that \eqn{x_t=(y_{t-1}',
#' \cdots,y_{t-k}')'} which is \eqn{K \times 1}. We have
#' \deqn{y_t = A' x_t + u_t, \quad Y = XA + U}
#' where \eqn{Y/U} are vertically stacked \eqn{y_t'/u_t'} and \eqn{T \times N},
#' \eqn{X} is vertically stacked \eqn{x_t'} and \eqn{T \times K}. \cr
#' Now the homoskedastic conjugate Normal-Inverse-Wishart prior is
#' \eqn{u_t|\Sigma\sim MN(0,I,\Sigma)} and
#' \deqn{Y|A,\Sigma,X\sim MN(XA,I,\Sigma),\quad
#' A|\Sigma\sim MN(\underline A,\underline V,\Sigma),\quad
#' \Sigma\sim IW(\underline S,\underline v)}
#' For SVAR, the forecast errors \eqn{u_t} can be decomposed to structural shocks \eqn{e_t}
#' such that
#' \deqn{u_t = B e_t = P Q e_t}
#' where \eqn{B,P,Q} are \eqn{N \times N}, \eqn{P=\text{chol}(\Sigma)} and
#' \eqn{Q} is some orthogonal matrix (usually just the identity matrix \eqn{I}) since
#' \deqn{\text{var}(PQe_t)=PQQ'P'=PP'=\Sigma=\text{var}(u_t) \Rightarrow P=\text{chol}(\Sigma)}
#'
#' @param specification list of output from \code{\link{minnesota_prior}}.
#' @param S total draws.
#' @param burn number of burnin draws.
#' @param thin thinning parameter for storing MCMC draws.
#'
#' @examples
#' data("fomc")
#' sign_restrictions <- matrix(c(1, -1, 1, 1), nrow = 2)
#' fomc |>
#'   minnesota_prior(k = 12) |>
#'   bvar(S = 1000) |>
#'   identify_sign(sign_restrictions) |>
#'   irf(periods = 36, shock_names = c("mp", "cbi"))
#'
#' @export
bvar <- function(specification, S = 200, burn = 50, thin = 10) {
  posterior <- .bvar_cpp(
    Y = as.matrix(specification$Y), k = specification$k, model = specification$model,
    S = S, burn = burn, thin = thin, prior = specification$prior
  )
  posterior$data$names <- colnames(specification$Y)
  posterior
}
