#' @title
#' Minnesota prior
#'
#' @description
#' Given data and lag order, construct Minnesota prior for Bayesian VAR.
#'
#' @param Y data matrix
#' @param k lag order
#' @param lambda0 prior variance of intercept
#' @param lambda1 prior variance of coefficients
#' @param lambda2 prior decay rate of coefficients
#' @param prior list of priors, might contain "non_stationary" for custom shrinkage
#'
#' @export
minnesota_prior <- function(Y,
                            k,
                            lambda0 = 1 / .Machine$double.eps,
                            lambda1 = 0.2,
                            lambda2 = 1,
                            prior = NULL) {
  N <- ncol(Y)
  K <- 1 + N * k

  ar_sigma2 <- vector(length = N)
  for (n in 1:N) {
    ar_sigma2[n] <- ar(Y[, n], order.max = k)$resid |> var()
  }
  ar_sigma2[is.na(ar_sigma2)] <- 1

  if ("non_stationary" %in% names(prior)) {
    non_stationary <- prior$non_stationary
  } else {
    non_stationary <- rep(1, N)
  }

  A <- matrix(0, K, N)
  A[2:(N + 1), 1:N] <- diag(non_stationary)

  V <- matrix(0, K, K)
  V[1, 1] <- lambda0
  V[2:K, 2:K] <- diag(lambda1^2 * rep((1:k)^-(2 * lambda2), each = N) * rep(ar_sigma2^-1, k))

  S <- diag(ar_sigma2)
  v <- N + 2

  list(A_l = A, V_l = V, S_l = S, v_l = v)
}
