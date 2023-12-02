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
irf <- function(posterior, periods, shock_names = NULL) {
  posterior <- .irf_cpp(posterior, periods)
  irf <- posterior$irf
  N <- nrow(irf)
  periods <- ncol(irf) / N

  # naming by order if no input
  if (is.null(shock_names)) {
    shock_names <- posterior$data$names
  }
  variable_index <- shock_index <- 1:N

  # median and credible intervals 68% / 90%
  irf <- apply(irf, c(1, 2), \(x) {
    quantile(x, probs = c(0.500, 0.05, 0.95, 0.160, 0.840))
  })

  # some data manipulation for desired ggplot2 format
  temp <- c()
  for (i in 1:dim(irf)[3]) {
    temp <- rbind(temp, as.matrix(irf[, , i]))
  }
  irf <- data.frame(temp)
  colnames(irf) <- variable_index
  irf["CI"] <- rep(
    c("median", "ci_wide_lower", "ci_wide_upper", "ci_narrow_lower", "ci_narrow_upper"),
    N * periods
  )
  irf["period"] <- rep(0:(periods - 1), each = 5 * N)
  irf["shock"] <- rep(rep(shock_index, each = 5), periods)

  irf <- irf |>
    dplyr::filter(shock <= length(shock_names)) |>
    tidyr::pivot_longer(variable_index, names_to = "variable") |>
    tidyr::pivot_wider(names_from = "CI", values_from = "value") |>
    dplyr::arrange(period, shock, variable)

  # if (cummulative) {
  #   irf <- irf |>
  #     dplyr::group_by(shock, variable) |>
  #     dplyr::mutate(
  #       ci_wide_lower = cumsum(ci_wide_lower),
  #       ci_wide_upper = cumsum(ci_wide_upper),
  #       ci_narrow_lower = cumsum(ci_narrow_lower),
  #       ci_narrow_upper = cumsum(ci_narrow_upper),
  #       median = cumsum(median)
  #     ) |>
  #     dplyr::ungroup()
  # }

  # labeller functions, otherwise the ordering will be messy
  labeller_variable <- labeller_shock <- c()
  for (i in 1:N) {
    labeller_variable[as.character(i)] <- posterior$data$names[i]
    labeller_shock[as.character(i)] <- shock_names[i]
  }

  fill_color <- "#b9d3ee"
  # plot
  irf |>
    ggplot2::ggplot(ggplot2::aes(x = period)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_wide_lower, ymax = ci_wide_upper),
      color = NA, fill = fill_color, alpha = 0.7
    ) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_narrow_lower, ymax = ci_narrow_upper),
      color = NA, fill = fill_color, alpha = 1
    ) +
    ggplot2::geom_line(ggplot2::aes(y = median), linetype = 7, color = "#021b8b", alpha = 0.8) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0))) +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(variable ~ shock,
      scales = "free_y",
      labeller = ggplot2::labeller(
        variable = labeller_variable,
        shock = labeller_shock
      )
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
    )
}
