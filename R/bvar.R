#' Bayesian vector autoregression
#'
#' Bayesian vector autoregression
#'
#' @param Y dependent variables, n by m
#' @param p number of lags
#' @param S number of draws
#' @param prior list of priors
#'
#' @export
bvar <- function(Y, p, S, prior) {
  # Minnesota-like Prior
  if (is.null(prior)) {
    m <- ncol(Y)
    mp <- m * p

    A_0 <- matrix(0, mp, m)
    A_0[1:m, 1:m] <- diag(m)
    V_0 <- diag(rep((1:p)^-2, each = m))
    Lambda_0 <- diag(m)
    nu_0 <- 1
    prior <- list(A_0 = A_0, V_0 = V_0, Lambda_0 = Lambda_0, nu_0 = nu_0)
  }

  .bvar_cpp(as.matrix(Y), p, S, prior)
}


#' Bayesian Structural vector autoregression
#'
#' Estimate a Bayesian Structural vector autoregression, base on the VAR(p) specification:
#' \deqn{y_t = A_1 y_{t-1} + \cdots + A_p y_{t-p} + u_t}
#' where \eqn{y_t} and \eqn{u_t} are \eqn{m \times 1}, \eqn{A_j} is \eqn{m \times m}.
#' We can stack the \eqn{A_j}'s vertically
#' such that \eqn{A=(A_1,\dots,A_p)'} which is \eqn{mp \times m}, and \eqn{y_{t-j}}'s vertically
#' such that \eqn{x_t=(y_{t-1}',\cdots,y_{t-p}')'} which is \eqn{mp \times 1}. We have
#' \deqn{y_t = A' x_t + u_t, \quad Y = XA + U}
#' where \eqn{Y/U} are vertically stacked \eqn{y_t'/u_t'} and \eqn{n \times m},
#' \eqn{X} is vertically stacked \eqn{x_t'} and \eqn{n \times mp}. \cr
#' Now the conjugate homoskedastic Normal-Inverse-Wishart prior is
#' \eqn{u_t|\Sigma\sim MN(0,I,\Sigma)} and
#' \deqn{Y|A,\Sigma,X\sim MN(XA,I,\Sigma),\quad
#' A|\Sigma,X\sim MN(A_0,V_0,\Sigma),\quad
#' \Sigma|X\sim IW(\Lambda_0,\nu_0)}
#' For SVAR, the forecast errors \eqn{u_t} can be decomposed to structural shocks \eqn{e_t}
#' such that
#' \deqn{u_t = P Q e_t}
#' where \eqn{P,Q} are \eqn{m \times m}, \eqn{P=\text{chol}(\Sigma)} and
#' \eqn{Q} is some orthogonal matrix (usually just the identity matrix \eqn{I}) since
#' \deqn{\text{var}(PQe_t)=PQQ'P'=PP'=\Sigma=\text{var}(u_t) \Rightarrow P=\text{chol}(\Sigma)}
#' Then, the model can be compactly written as
#' \deqn{y_t = A' x_t + P Q e_t}
#'
#' @param Y variables, n by m, n observations/periods, m variables
#' @param p number of lags
#' @param S number of draws
#' @param prior list of priors \eqn{(A_0,V_0,\Lambda_0,\nu_0)}
#'
#' @export
bsvar <- function(Y, p, S, prior, identification, sign = NULL) {
  posterior <- bvar(Y, p, S, prior)

  if (identification == "short") {
    posterior <- posterior |> .identify_shortrun_cpp()
  } else if (identification == "long") {
    posterior <- posterior |> .identify_longrun_cpp()
  } else {
    posterior <- posterior |> .identify_sign_cpp(sign)
  }

  posterior |> .estimate_irf_cpp(24)
}
