#' @title
#' Structural shocks identification
#'
#' @description
#' Needed for IRF estimation.
#'
#' @param posterior draws from bvar
#' @param strategy identification strategy, short/long/sign
#' @param sign sign matrix, only used when strategy is sign
#'
#' @export
identify <- function(posterior, strategy, sign = NULL) {
  if (strategy == "short") {
    posterior <- .identify_shortrun_cpp(posterior)
  } else if (strategy == "long") {
    posterior <- .identify_longrun_cpp(posterior)
  } else if (strategy == "sign") {
    posterior <- .identify_sign_cpp(posterior, sign)
  } else {
    stop("strategy must be short/long/sign")
  }
  posterior
}
