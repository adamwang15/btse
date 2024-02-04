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
  plot_irf(posterior$irf, posterior$data$names, shock_names) |> print()
  posterior
}

plot_irf <- function(irf, variable_names, shock_names) {
  dim_irf <- dim(irf[[1]])
  N <- dim_irf[1]
  h <- dim_irf[3]
  S <- length(irf)
  indices <- 1:N

  temp <- array(dim = c(N, N * h, S))
  for (s in 1:S) {
    temp[, , s] <- matrix(irf[[s]], N, N * h)
  }
  irf <- temp

  # credible intervals 68%, 50%, 90%
  irf <- apply(irf, c(1, 2), \(x) {
    quantile(x, probs = c(0.05, 0.16, 0.50, 0.84, 0.95))
  })

  # some data manipulation for desired ggplot2 format
  temp <- c()
  for (i in 1:dim(irf)[3]) {
    temp <- rbind(temp, as.matrix(irf[, , i]))
  }
  irf <- data.frame(temp)

  colnames(irf) <- indices
  irf["ci"] <- rep(c("ci_05", "ci_16", "ci_50", "ci_84", "ci_95"),
                   N * h)
  irf["horizon"] <- rep(0:(h - 1), each = 5 * N)
  irf["shock"] <- rep(rep(indices, each = 5), h)

  irf <- irf |>
    dplyr::filter(shock <= length(shock_names)) |>
    tidyr::pivot_longer(indices, names_to = "variable") |>
    tidyr::pivot_wider(names_from = "ci", values_from = "value") |>
    dplyr::arrange("horizon", "shock", "variable")

  # label functions, otherwise the ordering will be messy
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