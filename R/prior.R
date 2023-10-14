#' @title
#' Minnesota prior
#'
#' @description
#' Given data and lag order, construct Minnesota prior for Bayesian VAR.
#'
#' @param Y data matrix
#' @param k lag order
#' @param model type of model, either "natural" or "independent"
#' @param lambda_1 prior variance of intercept
#' @param lambda_2 prior variance of coefficients
#' @param lambda_3 prior decay rate of coefficients
#' @param prior list of priors, might contain "non_stationary" for custom shrinkage
#'
#' @export
minnesota_prior <- function(Y,
                            k,
                            non_stationary = NULL,
                            model = "conjugate",
                            lambda_1 = 1 / .Machine$double.eps,
                            lambda_2 = 0.2,
                            lambda_3 = 1) {
  T <- nrow(Y)
  N <- ncol(Y)
  K <- 1 + N * k

  ar_s2 <- vector(length = N)
  for (n in 1:N) {
    resid <- ar(Y[, n], order.max = k)$resid |> na.omit()
    ar_s2[n] <- t(resid) %*% resid / (T - k - 1)
  }

  if (is.null(non_stationary)) {
    non_stationary <- vector(length = N)
    for (n in 1:N) {
      non_stationary[n] <- suppressWarnings(as.numeric(tseries::adf.test(Y[, n])$p.value > 0.05))
    }
  }

  A <- matrix(0, K, N)
  A[2:(N + 1), 1:N] <- diag(non_stationary)

  V <- matrix(0, K, K)
  V[1, 1] <- lambda_1
  V[2:K, 2:K] <- diag(lambda_2^2 * rep((1:k)^-(2 * lambda_3), each = N) * rep(ar_s2^-1, k))
  if (model != "conjugate") {
    V <- kronecker(diag(ar_s2), V)
    for (n in 1:N) {
      idx <- (n - 1) * K + 1
      V[idx, idx] <- lambda_1
    }
  }

  S <- diag(ar_s2)
  v <- N + 2

  list(
    Y = Y, k = k, model = model,
    prior = list(A_l = A, V_l = V, S_l = S, v_l = v)
  )
}
