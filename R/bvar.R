#' Bayesian vector autoregression
#'
#' Bayesian vector autoregression
#'
#' @param Y dependent variables, n by m
#' @param k number of lags
#' @param S number of draws
#' @param prior list of priors \eqn{(B_0,V_0,\Lambda_0,\nu_0)}
#'
#' @export
bvar <- function(Y, k, S, prior) {
  draws <- bvar_cpp(as.matrix(Y), k, S, prior)
  result <- apply(draws$B, c(1, 2), mean) |> round(2)
  # rownames(result) <- colnames(X)
  # colnames(result) <- colnames(Y)
  result
}
