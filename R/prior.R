#' @title
#' Minnesota prior
#'
#' @description
#' Given data and lag order, construct Minnesota prior for Bayesian VAR.
#'
#' @param Y data matrix
#' @param p lag order
#'
#' @export
minnesota_prior <- function(Y, p, non_stationary, lambda = 5) {
  m <- ncol(Y)
  mp <- m * p

  ar_sigma2 <- vector(length = m)
  for (i in 1:m) {
    ar_sigma2[i] <- ar(Y[, i], order.max = p)$resid |> var()
  }
  ar_sigma2[is.na(ar_sigma2)] <- 1

  A_0 <- matrix(0, mp, m)
  A_0[1:m, 1:m] <- diag(non_stationary)
  V_0 <- lambda^-1 * diag(rep((1:p)^-1, each = m) * rep(ar_sigma2^-1, p))

  Lambda_0 <- diag(ar_sigma2)
  Lambda_0 <- diag(m)
  nu_0 <- m + 2

  list(A_0 = A_0, V_0 = V_0, Lambda_0 = Lambda_0, nu_0 = nu_0)
}
