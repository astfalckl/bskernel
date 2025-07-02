#' Whittle Log-Likelihood Objective
#'
#' Computes the negative Whittle log-likelihood for a given coefficient vector.
#'
#' @param c Coefficient vector for the spectral basis expansion.
#' @param B Design matrix (rows = frequencies, cols = basis functions). Easiest to use the output from `build_bspline_design_matrix`
#' @param periodogram Vector of periodogram values.
#'
#' @return A scalar value: the negative Whittle log-likelihood.
#' @export
whittle_loglik <- function(c, B, periodogram) {
  f <- pmax(as.numeric(B %*% c), 1e-8)  # ensure strictly positive
  if (any(f <= 0)) return(Inf)
  -sum(log(f) + periodogram / f)
}


#' Gradient of Whittle Log-Likelihood
#'
#' Computes the gradient of the negative Whittle log-likelihood with respect to the coefficients.
#'
#' @param c Coefficient vector for the spectral basis expansion.
#' @param B Design matrix (rows = frequencies, cols = basis functions).
#' @param periodogram Vector of periodogram values.
#'
#' @return A numeric vector: the gradient with respect to \code{c}.
#' @export
whittle_loglik_grad <- function(c, B, periodogram) {
  f <- pmax(as.numeric(B %*% c), 1e-8)
  df <- (1 / f) - (periodogram / f^2)
  -as.numeric(t(B) %*% df)
}

#' Optimize Whittle Likelihood
#'
#' Fits spline coefficients by minimizing the negative Whittle log-likelihood.
#'
#' @param omega Vector of frequencies.
#' @param periodogram Vector of periodogram values.
#' @param B Design matrix (rows = frequencies, cols = basis functions).
#' @param c0 Optional starting values for the coefficients. Defaults to 0.1 for each basis.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{c}}{The estimated coefficient vector.}
#'   \item{\code{fitted}}{The fitted spectral density values at each frequency.}
#'   \item{\code{optim}}{The full output from \code{\link[stats]{optim}}.}
#' }
#' @export
optim_whittle <- function(omega, periodogram, B, c0 = NULL) {
  n_basis <- ncol(B)
  if (is.null(c0)) c0 <- rep(0.1, n_basis)

  fit <- stats::optim(
    par = c0,
    fn = function(c) -whittle_loglik(c, B, periodogram),
    gr = function(c) -whittle_loglik_grad(c, B, periodogram),
    method = "L-BFGS-B",
    lower = rep(0, n_basis)
  )

  list(
    c = fit$par,
    fitted = as.numeric(B %*% fit$par),
    optim = fit
  )
}
