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
irf <- function(posterior, horizons, shock_names) {
  if (missing(shock_names)) {
    shock_names <- 1:length(posterior$data$names)
  }

  posterior <- irf_cpp(posterior, horizons)

  # change armadillo field to array
  S <- length(posterior$irf)
  irf <- array(NA, c(dim(posterior$irf[[1]]), S))
  for (s in 1:S) {
    irf[, , , s] <- posterior$irf[[s]]
  }
  posterior$irf <- irf

  print(plot_irf(irf, posterior$data$names, shock_names))
  posterior
}

#' @title
#' Plot impulse response function
#'
#' @description
#' Given draws of irf, plot the impulse response function.
#'
#' @param irf draws of impulse response function
#' @param variable_names variable names
#' @param shock_names structural shock names
#'
#' @export
plot_irf <- function(irf, variable_names, shock_names) {
  N <- dim(irf)[1]
  h <- dim(irf)[3]
  S <- dim(irf)[4]
  indices <- 1:N

  stacked_irf <- array(dim = c(N, N * h, S))
  for (s in 1:S) {
    stacked_irf[, , s] <- matrix(irf[, , , s], N, N * h)
  }
  irf <- stacked_irf

  # credible intervals 68%, 50%, 90% and some manipulation for ploting
  irf <- apply(irf, c(1, 2), \(x) {
    quantile(x, probs = c(0.05, 0.16, 0.50, 0.84, 0.95))
  }) |>
    aperm(c(1, 3, 2)) |>
    matrix(nrow = 5 * N * h, ncol = N) |>
    as.data.frame() |>
    dplyr::rename_with(~ as.character(indices)) |>
    dplyr::mutate(ci = rep(c(
      "ci_05", "ci_16", "ci_50", "ci_84", "ci_95"
    ), N * h)) |>
    dplyr::mutate(horizon = rep(0:(h - 1), each = 5 * N)) |>
    dplyr::mutate(shock = rep(rep(indices, each = 5), h)) |>
    dplyr::filter(shock <= length(shock_names)) |>
    tidyr::pivot_longer(indices, names_to = "variable") |>
    tidyr::pivot_wider(names_from = "ci", values_from = "value") |>
    dplyr::arrange("horizon", "shock", "variable")

  # label variables and shocks
  labeller_variable <- labeller_shock <- c()
  for (i in 1:N) {
    labeller_variable[as.character(i)] <- variable_names[i]
    labeller_shock[as.character(i)] <- shock_names[i]
  }

  # plot
  fill_color <- "#b9d3ee"
  irf |>
    ggplot2::ggplot(ggplot2::aes(x = horizon)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_05, ymax = ci_95),
      color = NA,
      fill = fill_color,
      alpha = 0.7
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_16, ymax = ci_84),
      color = NA,
      fill = fill_color,
      alpha = 1
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = ci_50),
      linetype = 7,
      color = "#021b8b",
      alpha = 0.8
    ) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0))) +
    ggplot2::theme_bw() +
    ggplot2::facet_grid(
      variable ~ shock,
      scales = "free_y",
      labeller = ggplot2::labeller(variable = labeller_variable, shock = labeller_shock)
    ) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
    )
}
