#' @title
#' Impulse response function
#'
#' @description
#' Given posterior draws from bsvar model, plot the impulse response function.
#'
#' @param posterior posterior draws with impulse response function
#' @param shock_names custom shock names in order
#'
#' @export
irf <- function(posterior, shock_names = NULL) {
  irf <- posterior$irf
  N <- nrow(irf)
  periods <- ncol(irf) / N

  # naming by order if no input
  if (is.null(shock_names)) {
    shock_names <- sapply(1:N, \(x) paste("shock", x))
  }
  variable_index <- shock_index <- 1:N

  # median and credible intervals 68% / 95%
  irf <- apply(irf, c(1, 2), \(x) {
    quantile(x, probs = c(0.500, 0.025, 0.975, 0.160, 0.840))
  })

  # some data manipulation for desired ggplot2 format
  temp <- c()
  for (i in 1:dim(irf)[3]) {
    temp <- rbind(temp, as.matrix(irf[, , i]))
  }
  irf <- dplyr::as_tibble(temp)
  colnames(irf) <- variable_index
  irf["CI"] <- rep(
    c("median", "CI_95_lower", "CI_95_upper", "CI_68_lower", "CI_68_upper"),
    N * periods
  )
  irf["shock"] <- rep(rep(shock_index, each = 5), periods)
  irf["period"] <- rep(0:(periods - 1), each = 5 * N)
  irf <- irf |>
    tidyr::pivot_longer(variable_index, names_to = "variable") |>
    tidyr::pivot_wider(names_from = "CI", values_from = "value")

  # labeller functions, otherwise the ordering will be messy
  labeller_variable <- labeller_shock <- c()
  for (i in 1:N) {
    labeller_variable[as.character(i)] <- posterior$names[i]
    labeller_shock[as.character(i)] <- shock_names[i]
  }

  # plot
  irf |> ggplot2::ggplot(ggplot2::aes(x = period)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_line(ggplot2::aes(y = median)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = CI_95_lower, ymax = CI_95_upper),
      color = NA, fill = "dodgerblue", alpha = 0.4
    ) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = CI_68_lower, ymax = CI_68_upper),
      color = NA, fill = "dodgerblue", alpha = 0.5
    ) +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(variable ~ shock,
      switch = "y",
      scales = "free_y",
      labeller = ggplot2::labeller(
        variable = labeller_variable,
        shock = labeller_shock
      )
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )
}
