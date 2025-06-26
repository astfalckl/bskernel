#' Inverse Fourier Transform of a B-spline Basis Function
#'
#' Computes the inverse Fourier transform of a B-spline basis function of degree \code{k},
#' indexed by \code{i}, using its truncated power expansion representation. This is used
#' to construct spectral representations of B-spline basis functions.
#'
#' @param tau A numeric vector of frequencies at which to evaluate the inverse Fourier transform.
#' @inheritParams bspline_alpha_rule
#'
#' @return A complex-valued vector of the same length as \code{tau}, giving the inverse Fourier transform
#' of the basis function evaluated at each frequency.
#'
#' @examples
#' knots <- c(0, 1, 2, 3, 4, 5)
#' tau <- seq(-5, 5, length.out = 200)
#' rho <- inverse_fourier_truncated_power(knots, i = 1, k = 1, tau = tau)
#' plot(tau, Re(rho), type = "l", main = "Real part of Fourier transform")
#'
#' @export
inverse_fourier_truncated_power <- function(knots, i, k, tau) {
  alpha <- bspline_alpha_rule(knots, i, k)
  alpha <- alpha * bspline_normalisation_constant(knots, i, k)

  stopifnot(length(alpha) == k + 2)
  stopifnot(i + k + 2 <= length(knots))

  kappa_w <- knots[i + k + 2]  # upper limit of support
  two_pi_i_tau <- 2 * pi * 1i * tau
  denom <- (two_pi_i_tau)^(k + 1)

  rho <- rep(0 + 0i, length(tau))

  for (j in 0:(k + 1)) {
    kappa_j <- knots[i + j + 1]
    alpha_j <- alpha[j + 1]

    exp1 <- exp(two_pi_i_tau * kappa_j)
    exp2 <- exp(two_pi_i_tau * (kappa_w - kappa_j))

    l_sum <- rep(0 + 0i, length(tau))
    for (l in 0:k) {
      l_sum <- l_sum + ((-two_pi_i_tau * (kappa_w - kappa_j))^l) / factorial(l)
    }

    rho <- rho + alpha_j * exp1 * (exp2 * l_sum - 1)
  }

  prefactor <- (-1)^k * factorial(k) / denom
  return(prefactor * rho)
}

#' Reconstruct the Autocovariance Function (ACF) from Spectral Coefficients
#'
#' Given a set of spectral coefficients and a B-spline basis, this function reconstructs
#' the autocovariance function via inverse Fourier transform.
#'
#' @param c A numeric vector of spectral coefficients (length equals number of basis functions).
#' @inheritParams bspline_alpha_rule
#' @inheritParams inverse_fourier_truncated_power
#'
#' @return A real-valued numeric vector of the same length as \code{tau}, representing the reconstructed autocovariance.
#'
#' @details
#' The function computes a weighted sum of inverse Fourier transforms of the B-spline basis functions.
#' The zero-lag autocovariance (at \code{tau[1]}) is manually set to the total spectral mass \code{sum(c)},
#' ensuring numerical consistency.
#'
#' @examples
#' knots <- c(0, 1, 2, 3, 4, 5)
#' c <- c(1, 0.5, 0.2)
#' tau <- seq(0, 10, length.out = 200)
#' acf <- reconstruct_acf(c, knots, k = 1, tau = tau)
#' plot(tau, acf, type = "l")
#'
#' @export
reconstruct_acf <- function(c, knots, k, tau) {
  n_basis <- length(c)
  acf <- rep(0 + 0i, length(tau))

  for (i in seq_len(n_basis) - 1) {
    rho_i <- inverse_fourier_truncated_power(knots, i, k, tau)
    acf <- acf + c[i + 1] * rho_i
  }

  acf[1] <- sum(c)  # enforce correct variance at lag zero
  return(2 * Re(acf))   # autocovariance is real-valued
}
