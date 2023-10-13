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
plot_irf <- function(posterior, shock_names = NULL) {
  IRF <- posterior$IRF
  m <- nrow(IRF)
  periods <- ncol(IRF) / m

  # naming by order if no input
  if (is.null(shock_names)) {
    shock_names <- sapply(1:m, \(x) paste("shock", x))
  }
  variable_index <- shock_index <- 1:m

  # median and credible intervals 68% / 95%
  IRF <- apply(IRF, c(1, 2), \(x) {
    quantile(x, probs = c(0.500, 0.025, 0.975, 0.160, 0.840))
  })

  # some data manipulation for desired ggplot2 format
  temp <- c()
  for (i in 1:dim(IRF)[3]) {
    temp <- rbind(temp, as.matrix(IRF[, , i]))
  }
  IRF <- dplyr::as_tibble(temp)
  colnames(IRF) <- variable_index
  IRF["CI"] <- rep(
    c("median", "CI_95_lower", "CI_95_upper", "CI_68_lower", "CI_68_upper"),
    m * periods
  )
  IRF["shock"] <- rep(rep(shock_index, each = 5), periods)
  IRF["period"] <- rep(0:(periods - 1), each = 5 * m)
  IRF <- IRF |>
    tidyr::pivot_longer(variable_index, names_to = "variable") |>
    tidyr::pivot_wider(names_from = "CI", values_from = "value")

  # labeller functions, otherwise the ordering will be messy
  labeller_variable <- labeller_shock <- c()
  for (i in 1:m) {
    labeller_variable[as.character(i)] <- posterior$names[i]
    labeller_shock[as.character(i)] <- shock_names[i]
  }

  # plot
  IRF |> ggplot2::ggplot(ggplot2::aes(x = period)) +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_line(ggplot2::aes(y = median)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = CI_95_lower, ymax = CI_95_upper),
      color = NA, fill = "dodgerblue", alpha = 0.2
    ) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = CI_68_lower, ymax = CI_68_upper),
      color = NA, fill = "dodgerblue", alpha = 0.5
    ) +
    ggplot2::theme_bw() +
    ggplot2::facet_wrap(~ variable + shock,
      scales = "free_y",
      labeller = ggplot2::labeller(
        variable = labeller_variable,
        shock = labeller_shock
      )
    )
}
