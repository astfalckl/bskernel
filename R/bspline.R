#' Alpha Coefficients for Truncated Power Expansion of B-spline Basis
#'
#' Computes the coefficients \eqn{\alpha_j} in the truncated power representation of a B-spline basis function
#' of degree \code{k}, valid for \code{k = 0, 1, 2}. These coefficients allow the basis function to be expressed as a
#' linear combination of truncated power functions over the relevant knot interval.
#'
#' @param knots A numeric vector of knots.
#' @param i An integer index indicating the position in the knot sequence for which to compute the coefficients.
#' @param k The degree of the B-spline basis function (must be 0, 1, or 2).
#'
#' @return A numeric vector of coefficients \eqn{\alpha_j} for the truncated power expansion.
#' 
#' @details
#' For:
#' \itemize{
#'   \item \code{k = 0}, returns \eqn{(1, -1)} for the characteristic function over a single interval.
#'   \item \code{k = 1}, returns a 3-vector corresponding to linear B-spline basis functions.
#'   \item \code{k = 2}, returns a 4-vector for quadratic B-spline basis functions.
#' }
#' The function throws an error if any required knot difference is zero (i.e. repeated knots), or if \code{k > 2}.
#'
#' @examples
#' knots <- c(0, 1, 2, 3, 4, 5)
#' bspline_alpha_rule(knots, i = 1, k = 1)
#'
#' @export
bspline_alpha_rule <- function(knots, i, k) {

  if (!is.numeric(knots) || anyNA(knots)) stop("knots must be a numeric vector with no NAs")

  if (!is.numeric(i) || length(i) != 1) stop("i must be a single integer index")

  if (!is.numeric(k) || length(k) != 1 || !(k %in% 0:2)) stop("k must be 0, 1, or 2")

  if (k == 0) {
    return(c(1, -1))
  }

  if (k == 1) {
    d1 <- knots[i + 3] - knots[i + 2]
    d0 <- knots[i + 2] - knots[i + 1]
    if (d0 == 0 || d1 == 0) stop("Zero-length knot interval for k = 1")
    return(c(1 / d0, - (1 / d0 + 1 / d1), 1 / d1))
  }

  if (k == 2) {
    ks <- knots[(i + 1):(i + 4)]
    d01 <- ks[2] - ks[1]
    d02 <- ks[3] - ks[1]
    d13 <- ks[4] - ks[3]
    d12 <- ks[3] - ks[2]
    d14 <- ks[4] - ks[1]
    d24 <- ks[4] - ks[2]
    d34 <- ks[4] - ks[3]

    if (any(c(d01, d02, d13, d12, d14, d24, d34) == 0)) stop("Zero-length knot interval for k = 2")

    c1 <- 1 / d01 / d02
    c2 <- ((d13)^2 / d24 / d34 - c1 * (ks[3] - ks[1])^2) / d12^2
    c3 <- (-c1 * d14^2 - c2 * d24^2) / d34^2
    c4 <- -(c1 + c2 + c3)

    return(c(c1, c2, c3, c4))
  }
}

#' Normalisation Constant for Truncated Power Expansion of B-spline
#'
#' Computes the scaling factor required to normalise a B-spline basis function,
#' using its truncated power expansion, such that the basis integrates to 1 over its support.
#'
#' @inheritParams bspline_alpha_rule
#'
#' @return A numeric scalar giving the inverse of the integral under the truncated power representation - i.e., a
#' multiplicative constant \eqn{Z^{-1}} that makes the basis function integrate to 1.
#'
#' @details
#' This function uses \code{\link{bspline_alpha_rule}} to compute the coefficients
#' \eqn{\alpha_j} for the truncated power expansion, and integrates over the support.
#'
#' @examples
#' knots <- c(0, 1, 2, 3, 4, 5)
#' bspline_normalisation_constant(knots, i = 1, k = 1)
#'
#' @export
bspline_normalisation_constant <- function(knots, i, k) {
  alpha <- bspline_alpha_rule(knots, i, k)

  support_upper <- knots[i + k + 2]
  kappas <- knots[(i + 1):(i + k + 2)]

  deltas <- support_upper - kappas
  deltas[deltas < 0] <- 0  # ensure finite support only

  integral_terms <- alpha * deltas^(k + 1) / (k + 1)
  Z <- sum(integral_terms)
  if (Z <= 0) stop("Invalid integral - check knots and alpha")

  return(1 / Z)  # scaling factor to multiply alpha by
}

#' Evaluate a Normalised B-spline Basis Function
#'
#' Evaluates the B-spline basis function of degree \code{k}, indexed by \code{i},
#' at a vector of input locations \code{omega}, using its truncated power expansion.
#' The result is scaled so that the basis integrates to 1 over its support.
#'
#' @param omega A numeric vector of input locations at which to evaluate the basis function.
#' @inheritParams bspline_alpha_rule
#'
#' @return A numeric vector of the same length as \code{omega}, giving the evaluated basis function values.
#'
#' @examples
#' knots <- c(0, 1, 2, 3, 4, 5)
#' omega <- seq(0, 5, length.out = 100)
#' y <- evaluate_bspline_basis(omega, knots, i = 1, k = 1)
#' plot(omega, y, type = "l")
#'
#' @export
evaluate_bspline_basis <- function(omega, knots, i, k) {
  alpha <- bspline_alpha_rule(knots, i, k)
  knot_segment <- knots[(i + 1):(i + k + 2)]

  B_val <- numeric(length(omega))
  for (j in seq_along(alpha)) {
    B_val <- B_val + alpha[j] * truncated_power(omega, knot_segment[j], k)
  }

  norm <- bspline_normalisation_constant(knots, i, k)
  return(norm * B_val)
}

#' Construct a B-spline Design Matrix
#'
#' Builds the B-spline basis matrix by evaluating each normalised truncated power B-spline
#' basis function at the input locations \code{omega}. Each column of the matrix corresponds
#' to one basis function.
#'
#' @inheritParams evaluate_bspline_basis
#' @inheritParams bspline_alpha_rule
#'
#' @return A numeric matrix with \code{length(omega)} rows and \code{length(knots) - k - 1} columns,
#' where each column is a normalised B-spline basis function evaluated at \code{omega}.
#'
#' @details
#' The number of basis functions is given by \code{length(knots) - k - 1}. Each column of the matrix
#' is computed using \code{\link{evaluate_bspline_basis}}.
#'
#' @examples
#' knots <- c(0, 1, 2, 3, 4, 5)
#' omega <- seq(0, 5, length.out = 100)
#' B <- build_bspline_design_matrix(omega, knots, k = 1)
#' matplot(omega, B, type = "l", lty = 1, col = 1:ncol(B))
#'
#' @export
build_bspline_design_matrix <- function(omega, knots, k) {
  n_basis <- length(knots) - (k + 1)
  stopifnot(n_basis > 0)

  B_mat <- matrix(0, nrow = length(omega), ncol = n_basis)

  for (i in seq_len(n_basis) - 1) {
    B_mat[, i + 1] <- evaluate_bspline_basis(omega, knots, i, k)
  }

  return(B_mat)
}

#' Plot B-spline Basis Functions Using ggplot2
#'
#' Plots a B-spline design matrix using \code{ggplot2}, optionally applying log scales
#' to the x or y axes. Each column of the matrix is treated as a separate basis function.
#'
#' @param omega A numeric vector of input values at which the basis is evaluated.
#' @param design_matrix A numeric matrix (rows = length of \code{omega}, columns = basis functions).
#' @param k (Optional) The degree of the B-spline basis. Used only for labelling.
#' @param logx Logical; if \code{TRUE}, plot the x-axis on a log10 scale.
#' @param logy Logical; if \code{TRUE}, plot the y-axis on a log10 scale.
#'
#' @return A \code{ggplot} object displaying each basis function as a line.
#'
#'
#' @importFrom ggplot2 ggplot aes geom_line labs scale_x_continuous scale_y_continuous theme_minimal theme
#' @importFrom tidyr pivot_longer
#' @export
plot_bspline_basis_ggplot <- function(omega,design_matrix, k = NULL, logx = FALSE, logy = FALSE) {
  stopifnot(length(omega) == nrow(design_matrix))

  basis_df <- as.data.frame(design_matrix)
  basis_df$omega <- omega

  basis_long <- pivot_longer(
    basis_df,
    cols = -omega,
    names_to = "basis",
    values_to = "value"
  )

  plt <- ggplot(basis_long, aes(x = omega, y = value, group = basis, linetype = basis)) +
    geom_line() +
    labs(
      x = expression(omega),
      y = "Basis value",
      title = if (!is.null(k)) paste("B-spline basis functions ( k =", k, ")") else "B-spline basis functions"
    ) +
    theme_minimal() +
    theme(legend.position = "none")

  if (logx) plt <- plt + scale_x_continuous(trans = "log10")
  if (logy) plt <- plt + scale_y_continuous(trans = "log10")

  plt
}
