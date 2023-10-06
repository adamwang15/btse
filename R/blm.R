#' Bayesian linear model
#'
#' Bayesian linear regression with normal-inverse-Wishart conjugate prior
#' \deqn{
#' Y|B,\Sigma,X \sim MN(XB,I,\Sigma),\
#' B|\Sigma \sim MN(B_0, V_0, \Sigma),\
#' \Sigma \sim IW(\Lambda_0, \nu_0)
#' }
#' where \eqn{MN} is the matrix normal distribution.
#'
#' @param Y dependent variables, n by m
#' @param X independent variables, n by k
#' @param S number of draws
#' @param prior list of priors \eqn{(B_0,V_0,\Lambda_0,\nu_0)}
#'
#' @export
blm <- function(Y, X, S, prior) {
  if (is.null(prior)) {
    m <- ncol(Y)
    k <- ncol(X)

    B_0 <- matrix(0, k, m)
    V_0 <- diag(k)
    Lambda_0 <- diag(m)
    nu_0 <- 10
    prior <- list(B_0 = B_0, V_0 = V_0, Lambda_0 = Lambda_0, nu_0 = nu_0)
  }

  draws <- blm_cpp(as.matrix(Y), as.matrix(X), S, prior)
  result <- apply(draws$B, c(1, 2), mean) |> round(2)
  rownames(result) <- colnames(X)
  colnames(result) <- colnames(Y)
  result
}
