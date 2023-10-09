#' Bayesian vector autoregression
#'
#' Bayesian vector autoregression
#'
#' @param Y dependent variables, n by m
#' @param p number of lags
#' @param S number of draws
#' @param prior list of priors \eqn{(B_0,V_0,\Lambda_0,\nu_0)}
#'
#' @export
bvar <- function(Y, p, S, prior) {
  # Minnesota-like Prior
  if (is.null(prior)) {
    m <- ncol(Y)
    mp <- m * p

    B_0 <- matrix(0, mp, m)
    B_0[1:m, 1:m] <- diag(m)
    V_0 <- diag(rep((1:p)^-2, each = m))
    Lambda_0 <- diag(m)
    nu_0 <- 1
    prior <- list(B_0 = B_0, V_0 = V_0, Lambda_0 = Lambda_0, nu_0 = nu_0)
  }

  .bvar_cpp(as.matrix(Y), p, S, prior)
}

#' Bayesian Structural vector autoregression
#'
#' Bayesian Structural vector autoregression
#'
#' @param Y dependent variables, n by m
#' @param p number of lags
#' @param S number of draws
#' @param prior list of priors \eqn{(B_0,V_0,\Lambda_0,\nu_0)}
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
