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
    V_0 <- diag((1:mp)^-2) * 10
    Lambda_0 <- diag(m)
    nu_0 <- 1
    prior <- list(B_0 = B_0, V_0 = V_0, Lambda_0 = Lambda_0, nu_0 = nu_0)
  }

  posterior <- .bvar_cpp(as.matrix(Y), p, S, prior) |>
    .identify_sign_cpp(matrix(c(1, -1, 1, 1), nrow = 2)) |>
    # .identify_longrun_cpp() |>
    .estimate_irf_cpp(24)

  posterior
}
