#' Truncated Power Function
#'
#' Computes the truncated power basis function \eqn{(x - \kappa)^k_+}, which is zero for values of \code{x < kappa}
#' and equal to \eqn{(x - \kappa)^k_+} otherwise. Commonly used in spline basis construction.
#'
#' @param x A numeric vector of input values.
#' @param kappa A numeric scalar specifying the truncation point.
#' @param k A non-negative integer specifying the exponent.
#'
#' @return A numeric vector of the same length as \code{x}, with each element transformed by the truncated power function.
#' 
#' @examples
#' x <- seq(0, 2, by = 0.1)
#' truncated_power(x, kappa = 1, k = 2)
#'
#' @export
truncated_power <- function(x, kappa, k) {
  out <- (x - kappa)^k
  out[x < kappa] <- 0
  return(out)
}

#' Offset Logarithm Transform
#'
#' Computes a log-transformation of \code{x} with a small positive bias \code{b},
#' and an optional logarithmic \code{base}.
#'
#' @param x A numeric vector or matrix of values to transform.
#' @param b A small positive offset added to \code{x} before applying the logarithm. Defaults to \code{1e-6}.
#' @param base The base of the logarithm. Defaults to \code{10}.
#'
#' @return A numeric vector or matrix with the log-transformed values.
#'
#' @examples
#' logbp(1:10)
#' logbp(0, b = 1e-4)
#' logbp(c(0.1, 1, 10), base = exp(1))
#'
#' @export
logbp <- function(x, b = 1e-6, base = 10^1) log(x + b, base)
