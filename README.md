# Greetings

Hello, World! This is an R package about Bayesian Econometrics.

# Installation

``` r
devtools::install_git("https://github.com/adamwang15/btse.git")
```

# Example

A replication of Jarociński and Karadi (2020)

``` r
data("fomc")
sign_restrictions <- matrix(c(1, -1, 1, 1), nrow = 2)
fomc |>
  minnesota_prior(k = 12) |>
  bvar(S = 1000) |>
  identify_sign(sign_restrictions) |>
  irf(periods = 36, shock_names = c("mp", "cbi"))
```
