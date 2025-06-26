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
