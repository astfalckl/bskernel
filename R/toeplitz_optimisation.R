#' Log-Likelihood under Toeplitz Gaussian Model
#'
#' Computes the exact log-likelihood of a univariate process \code{y} under a
#' Gaussian process model with autocovariance structure derived from a B-spline
#' spectral basis. The autocovariance function is reconstructed and used to define
#' a Toeplitz covariance matrix, enabling efficient log-likelihood computation.
#'
#' @param y A numeric vector representing the observed time series.
#' @param c A numeric vector of coefficients for the ACF basis functions.
#' @param knots A numeric vector of knot locations defining the B-spline basis.
#' @param k The degree of the B-spline basis functions (here, 0, 1, or 2).
#'
#' @return A numeric scalar: the log-likelihood of the data under the implied Toeplitz covariance model.
#'
#' @details
#' The autocovariance function is reconstructed using \code{\link{reconstruct_acf}} and
#' passed to the \code{SuperGauss::NormalToeplitz} class for efficient evaluation of the
#' Gaussian likelihood with Toeplitz structure. This method avoids explicit matrix inversion
#' and is well-suited for large univariate processes.
#'
#' @seealso \code{\link{reconstruct_acf}}, \code{SuperGauss::NormalToeplitz}
#'
#' @examples
#' knots <- seq(0, 0.5, length.out = 6)
#' k <- 1
#' c <- rep(1, length(knots) - k - 1)
#' y <- arima.sim(n = 100, model = list(ar = 0.8))
#' compute_toeplitz_loglik(y, c, knots, k)
#'
#' @export
compute_toeplitz_loglik <- function(y, c, knots, k) {
  n <- length(y)
  tau <- 0:(n - 1)
  acf <- Re(reconstruct_acf(c, knots, k, tau))
  nt <- SuperGauss::NormalToeplitz$new(N = n)
  ll <- nt$logdens(z = y, acf = acf)
  return(ll)
}


#' Gradient of Toeplitz Log-Likelihood
#'
#' Computes the gradient of the Gaussian log-likelihood of a univariate process
#' under a model with autocovariance structure defined via B-spline basis functions,
#' either with respect to the coefficients or their logarithms.
#'
#' @param y A numeric vector representing the observed time series.
#' @param c A numeric vector of spectral coefficients.
#' @param knots A numeric vector of knot locations defining the B-spline basis.
#' @param k The degree of the B-spline basis functions.
#' @param log_coef Logical; if \code{TRUE} (default), returns the gradient with respect to \code{log(c)}.
#'
#' @return A numeric vector: the gradient of the log-likelihood with respect to \code{c} or \code{log(c)}.
#'
#' @details
#' When \code{log_coef = TRUE}, this function returns the gradient with respect to the log-transformed
#' spectral coefficients \eqn{\theta = \log(c)}, applying the chain rule. This is useful for optimization
#' in an unconstrained parameter space. When \code{log_coef = FALSE}, it returns the gradient with respect
#' to the original coefficients \code{c}.
#'
#' @seealso \code{\link{compute_toeplitz_loglik}}, \code{\link{reconstruct_acf}}, \code{SuperGauss::NormalToeplitz}
#'
#' @export
compute_toeplitz_loglik_grad <- function(
  y, c, knots, k, log_coef = TRUE
  ) {

  n <- length(y)
  tau <- 0:(n - 1)
  n_basis <- length(c)

  acf <- Re(reconstruct_acf(c, knots, k, tau))

  dacf <- sapply(seq_len(n_basis) - 1, function(i) {
    Re(inverse_fourier_truncated_power(knots, i, k, tau))
  })

  # Correct first row: acf[1] = sum(c), so ∂acf[1]/∂c_i = 1
  dacf[1, ] <- 1

  nt <- SuperGauss::NormalToeplitz$new(N = n)
  grad <- nt$grad(z = y, dz = matrix(0, n, n_basis), acf = acf, dacf = dacf)

  if (log_coef) {
    return(grad * c)
  } else {
    return(grad)
  }
}


#' Maximum Likelihood Estimation via Toeplitz Gaussian Likelihood
#'
#' Estimates the spectral basis coefficients that maximize the Gaussian log-likelihood
#' of a univariate time series with autocovariance structure modeled by B-spline basis functions.
#'
#' @param c_init A numeric vector of positive initial values for the coefficients.
#' @param knots A numeric vector of knot locations defining the B-spline basis.
#' @param k The degree of the B-spline basis functions.
#' @param y A numeric vector representing the observed time series.
#'
#' @return An object of class \code{\link[stats]{optim}}, containing the MLE estimates
#' (on the log scale), convergence diagnostics, and gradient information.
#'
#' @details
#' This function uses \code{\link[stats]{optim}} with the BFGS method to maximize
#' the log-likelihood under a Gaussian process model with a Toeplitz covariance matrix,
#' using \code{\link{compute_toeplitz_loglik}} and \code{\link{compute_toeplitz_loglik_grad}}.
#' The optimization is unconstrained by reparameterizing \code{c} on the log scale.
#'
#' The final estimates (in log-scale) are returned in \code{$par}, and can be transformed
#' back via \code{exp(result$par)}.
#'
#' @seealso \code{\link{compute_toeplitz_loglik}}, \code{\link{compute_toeplitz_loglik_grad}}, \code{\link[stats]{optim}}
#'
#'
#' @export
optim_toeplitz_mle <- function(c_init, knots, k, y) {
  stats::optim(
    par = log(c_init),
    fn = function(log_c) {
      -compute_toeplitz_loglik(y, exp(log_c), knots, k)
    },
    gr = function(log_c) {
      -compute_toeplitz_loglik_grad(y, exp(log_c), knots, k, log_coef = TRUE)
    },
    method = "BFGS"
  )
}