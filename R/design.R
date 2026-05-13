# =============================================================================
# R/design.R - Design matrix and group-factor helpers
# =============================================================================
#
# Source this via the top-level functions.R, which sources all R/*.R modules
# in the correct dependency order.
#
# Pure base-R utilities; no external package dependencies.
#
# =============================================================================

#' Build an intercept + binary-group design matrix
#'
#' Returns a two-column design matrix \code{cbind(Intercept = 1, Group)} where
#' \code{Group} is a 0/1 indicator. Used to feed
#' \code{\link{simulate_compositional_bb}} when one wants the simplest possible
#' two-cohort design.
#'
#' Two input modes are supported:
#' \describe{
#'   \item{\strong{Balanced}}{\code{n_samples_per_group} of length 1 -> the
#'     given count is used for both groups (total \code{2n} rows; first \code{n}
#'     are Group 0, next \code{n} are Group 1).}
#'   \item{\strong{Unbalanced}}{\code{n_samples_per_group} of length 2 -> the
#'     two entries give the per-group sample sizes; first block is Group 0,
#'     second is Group 1.}
#' }
#'
#' @param n_samples_per_group Positive integer scalar (balanced) or length-2
#'   integer vector (per-group sample counts).
#'
#' @return Numeric matrix with two columns named \code{"Intercept"} and
#'   \code{"Group"}, and \code{sum(n_samples_per_group)} (or \code{2*n}) rows.
#'
#' @seealso \code{\link{simulate_compositional_bb}},
#'   \code{\link{synthetic_sim_add_two_group_factor}}.
#' @keywords internal
design_matrix_from_groups <- function(n_samples_per_group) {
  if (length(n_samples_per_group) == 1L) {
    Group <- rep(c(0L, 1L), each = n_samples_per_group)
  } else if (length(n_samples_per_group) == 2L) {
    Group <- rep(c(0L, 1L), times = n_samples_per_group)
  } else {
    stop("n_samples_per_group must be length 1 (balanced) or length 2 (per-group counts).")
  }
  cbind(
    Intercept = 1,
    Group = Group
  )
}

#' Annotate a simulation result with a two-level study-specific group factor
#'
#' Takes the integer \code{Group} column produced by
#' \code{\link{design_matrix_from_groups}} (values \code{0L}/\code{1L}) and
#' adds a human-readable \code{group} factor with caller-supplied labels to
#' both \code{sim_r$sample_metadata} and \code{sim_r$count_long}. Used in the
#' per-dataset Quarto reports so plots and PERMANOVA output show e.g.
#' \code{"non-IBD"}/\code{"IBD"} or \code{"healthy"}/\code{"case"} instead of
#' \code{0}/\code{1}.
#'
#' The function mutates two known list elements in place (returning a new list);
#' it does not modify any other fields.
#'
#' @param sim_r List returned by \code{\link{simulate_compositional_bb}}. Must
#'   contain a \code{sample_metadata} and a \code{count_long} data frame, each
#'   with an integer \code{Group} column.
#' @param group_levels Length-2 character vector. \code{group_levels[[1]]}
#'   labels \code{Group == 0L}; \code{group_levels[[2]]} labels
#'   \code{Group == 1L}. The factor levels are set in this order.
#'
#' @return The same list \code{sim_r} with \code{sample_metadata$group} and
#'   \code{count_long$group} added/overwritten as a factor.
#'
#' @seealso \code{\link{design_matrix_from_groups}},
#'   \code{\link{simulate_compositional_bb}}.
#' @keywords internal
synthetic_sim_add_two_group_factor <- function(sim_r, group_levels) {
  if (length(group_levels) != 2) {
    stop("group_levels must have length 2 (first = Group 0, second = Group 1).")
  }
  to_factor <- function(g) factor(group_levels[g + 1L], levels = group_levels)
  sim_r$sample_metadata$group <- to_factor(sim_r$sample_metadata$Group)
  sim_r$count_long$group      <- to_factor(sim_r$count_long$Group)
  sim_r
}
