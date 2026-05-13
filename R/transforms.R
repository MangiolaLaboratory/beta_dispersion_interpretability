# =============================================================================
# R/transforms.R - Basic transforms and clamping helpers
# =============================================================================
#
# Source this via the top-level functions.R, which sources all R/*.R modules
# in the correct dependency order.
#
# Pure base-R utilities; no external package dependencies.
#
# =============================================================================

#' Softmax: convert log-linear predictors to compositional probabilities
#'
#' Exponentiates and normalizes row-wise so each row sums to 1. Used to map
#' design-matrix linear predictors (e.g. \code{design_matrix \%*\% coeffs}) to
#' taxon proportions. This is the standard softmax (not inverse softmax /
#' centered log-ratio).
#'
#' Numerical note: no max-shift is applied here. If your inputs can be very
#' large or very negative, prefer \code{\link{inv_softmax}} which subtracts the
#' row max for numerical stability.
#'
#' @param log_linear_predictors Numeric matrix of shape
#'   \code{n_samples x n_taxa} containing log-space values.
#'
#' @return Numeric matrix of probabilities with the same shape as
#'   \code{log_linear_predictors}; each row sums to 1.
#'
#' @seealso \code{\link{inv_softmax}} for the numerically stable variant that
#'   also accepts vectors.
#' @keywords internal
softmax <- function(log_linear_predictors) {
  exp_log <- exp(log_linear_predictors)
  row_sums <- rowSums(exp_log)
  probabilities <- exp_log / row_sums
  return(probabilities)
}

#' Clamp values into the open unit interval \eqn{[\epsilon, 1-\epsilon]}
#'
#' Internal helper used by \code{\link{arcsin_sqrt}}, \code{\link{inv_softmax}}
#' and (when present) variance / posterior-mean helpers, to keep proportions
#' away from the exact 0/1 boundary where downstream transforms (e.g.
#' arcsin-sqrt, logit) are ill-conditioned.
#'
#' @param x Numeric vector or matrix.
#' @param eps Non-negative numeric scalar; the boundary half-width. \code{0}
#'   disables clamping.
#'
#' @return Numeric object the same shape as \code{x}, with values clamped to
#'   \code{[eps, 1 - eps]}.
#' @noRd
.clamp01 <- function(x, eps = 0) {
  x <- pmax(x, 0 + eps)
  x <- pmin(x, 1 - eps)
  x
}

#' Arcsin-sqrt variance-stabilising transform for proportions
#'
#' Computes \eqn{g(p) = \arcsin(\sqrt{p})}. For binomial proportions this
#' transform approximately stabilises variance to \eqn{1/(4N)}.
#'
#' @param p Numeric vector or matrix of proportions in \eqn{[0, 1]}.
#' @param eps Non-negative numeric scalar. If \code{> 0}, \code{p} is clamped
#'   to \code{[eps, 1 - eps]} before the transform via \code{.clamp01}.
#'
#' @return Numeric object of the same shape as \code{p} containing
#'   \code{asin(sqrt(p))}.
#'
#' @seealso \code{\link{inv_softmax}}, the unit-interval helpers in this file.
#' @keywords internal
arcsin_sqrt <- function(p, eps = 0) {
  if (eps > 0) p <- .clamp01(p, eps = eps)
  asin(sqrt(p))
}

#' Numerically stable softmax for unconstrained log-linear predictors
#'
#' In this repo's naming, \code{"mu_inv_softmax"} refers to an *unconstrained*
#' vector (log-linear predictors). Applying \code{inv_softmax} converts it to
#' proportions. Mathematically this is the standard softmax with a row-wise
#' max-shift for numerical stability.
#'
#' Accepts both vectors and matrices:
#' \itemize{
#'   \item Vector input -> vector output that sums to 1.
#'   \item Matrix input -> row-wise softmax; each row sums to 1 and the input
#'     dimnames are preserved.
#' }
#'
#' Throws an informative error if any entry is \code{NA}; the caller is
#' expected to handle missing values upstream.
#'
#' @param x Numeric vector or matrix of unconstrained (real-valued) values.
#' @param eps Non-negative numeric scalar. If \code{> 0}, the returned
#'   probabilities are clamped to \code{[eps, 1 - eps]}.
#'
#' @return Numeric object of the same shape as \code{x} with probabilities in
#'   \eqn{(0, 1)} that sum to 1 along each row (or in total, for vector input).
#'
#' @seealso \code{\link{softmax}} for the non-stabilised matrix-only variant.
#' @keywords internal
inv_softmax <- function(x, eps = 0) {
  if (!is.numeric(x)) stop("x must be numeric.")

  softmax_vec <- function(v) {
    if (anyNA(v)) stop("x contains NA; please handle missing values before softmax.")
    v <- v - max(v) # numerical stability
    ex <- exp(v)
    out <- ex / sum(ex)
    if (eps > 0) out <- .clamp01(out, eps = eps)
    out
  }

  if (is.matrix(x)) {
    out <- t(apply(x, 1, softmax_vec))
    dimnames(out) <- dimnames(x)
    return(out)
  }

  if (is.vector(x)) {
    return(softmax_vec(as.numeric(x)))
  }

  stop("x must be a numeric vector or matrix.")
}
