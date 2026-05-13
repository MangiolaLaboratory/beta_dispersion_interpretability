# =============================================================================
# R/betadisper.R - betadisper input builders and permutest wrappers
# =============================================================================
#
# Source this via the top-level functions.R, which sources all R/*.R modules
# in the correct dependency order.
#
# Dependencies (see @importFrom on each function):
#   vegan   - permutest, vegdist, decostand
#   stats   - dist
#   dplyr   - select
#   tidyr   - pivot_wider
#   tibble  - column_to_rownames
#
# =============================================================================

#' Fixed-budget permutation test for a betadisper object
#'
#' Thin convenience wrapper around \code{vegan::permutest} that runs a fixed
#' number of permutations on a \code{vegan::betadisper} object and returns
#' only the headline group-effect p-value plus the permutation count.
#'
#' @param bd A \code{betadisper} object (the output of
#'   \code{vegan::betadisper}).
#' @param permutations Positive integer; number of permutations to run.
#'   Default \code{999}.
#'
#' @return Named list with:
#' \describe{
#'   \item{\code{p_value}}{Numeric in \eqn{[0, 1]}; the first \code{Pr(>F)}
#'     entry from \code{pt$tab}, i.e. the overall test of homogeneity of
#'     multivariate dispersions.}
#'   \item{\code{n_perm}}{Integer; echoes \code{permutations}.}
#' }
#'
#' @seealso \code{\link{permutest_betadisper_stable}} for the adaptive
#'   variant that keeps adding permutations until a target SE is reached.
#' @keywords internal
#' @importFrom vegan permutest
permutest_betadisper_fixed <- function(bd, permutations = 999) {
  pt <- vegan::permutest(bd, pairwise = TRUE, permutations = permutations)
  list(
    p_value = unname(pt$tab$`Pr(>F)`[1]),
    n_perm = permutations
  )
}

#' Adaptive permutation test that stops when the p-value SE is small enough
#'
#' Repeatedly calls \code{vegan::permutest} in batches, accumulating
#' permutations until either the Monte Carlo standard error of the headline
#' p-value falls below \code{target_se} or the budget \code{max_perm} is
#' exhausted. The SE estimate is \eqn{\sqrt{p(1-p)/m}} for the running
#' permutation count \eqn{m}.
#'
#' Note: each batch is an independent permutest call, so the final p-value is
#' the one from the \emph{last} batch (not a pooled estimate). This matches
#' the historical behaviour of the function and is fine when the only goal is
#' a stable headline p-value.
#'
#' @param bd A \code{betadisper} object.
#' @param max_perm Positive integer >= 999; the maximum number of
#'   permutations to attempt. Default \code{99999}.
#' @param batch Positive integer >= 100; permutations added per iteration.
#'   Default \code{5000}.
#' @param target_se Positive numeric; stop once the running SE drops below
#'   this. Default \code{0.0025}.
#'
#' @return Named list with:
#' \describe{
#'   \item{\code{p_value}}{Numeric; the headline \code{Pr(>F)} from the last
#'     batch.}
#'   \item{\code{p_se}}{Numeric; the running SE estimate at termination.}
#'   \item{\code{n_perm}}{Integer; total permutations run.}
#' }
#'
#' @seealso \code{\link{permutest_betadisper_fixed}}.
#' @keywords internal
#' @importFrom vegan permutest
permutest_betadisper_stable <- function(bd, max_perm = 99999, batch = 5000, target_se = 0.0025) {
  if (max_perm  < 999) stop("max_perm must be >= 999")
  if (batch     < 100) stop("batch must be >= 100")
  if (target_se <= 0)  stop("target_se must be > 0")

  m  <- 0L
  p  <- NA_real_
  se <- Inf
  while (m < max_perm && se > target_se) {
    m  <- as.integer(min(max_perm, m + batch))
    pt <- vegan::permutest(bd, pairwise = TRUE, permutations = m)
    p  <- unname(pt$tab$`Pr(>F)`[1])
    se <- if (is.finite(p)) sqrt(p * (1 - p) / m) else Inf
  }
  list(p_value = p, p_se = se, n_perm = m)
}

#' Build a standard distance + group factor for vegan::betadisper
#'
#' Pivots a tall \code{count_long} table into a sample-by-taxon count matrix,
#' applies a distance-method-specific pre-processing step, computes the
#' distance, and aligns a group factor so the result can be fed straight to
#' \code{vegan::betadisper(dist, group)}.
#'
#' Distance methods supported:
#' \describe{
#'   \item{\code{"aitchison"}, \code{"robust.aitchison"}}{Computed by
#'     \code{vegan::vegdist}. Zero/negative counts are replaced by
#'     \code{clr_pseudocount} before the call to keep the log finite.}
#'   \item{\code{"hellinger"}}{Hellinger transform (\code{vegan::decostand})
#'     followed by Euclidean \code{stats::dist}.}
#'   \item{\code{"jaccard"}}{Total-sum normalisation
#'     (\code{vegan::decostand(..., method = "total")}) then
#'     \code{vegan::vegdist(..., method = "jaccard", binary = FALSE)}
#'     (i.e. abundance-Jaccard rather than presence-Jaccard).}
#'   \item{\emph{Anything else}}{Forwarded directly to
#'     \code{vegan::vegdist(matrix, method = distance_method)}; supports
#'     bray, kulczynski, gower, ...}
#' }
#'
#' Samples with zero total counts are silently dropped before distance
#' computation. The returned group factor uses the order of
#' \code{group_levels} so downstream plots/tests have a stable reference.
#'
#' @param count_long Data frame with at least \code{sample_id}, \code{taxon_id}
#'   and \code{count} columns.
#' @param sample_metadata Data frame with at least \code{sample_id} and
#'   \code{group} columns (one row per sample).
#' @param distance_method Character scalar selecting the distance; see the
#'   list above.
#' @param group_levels Character vector defining the factor levels for the
#'   returned group factor (order matters for downstream contrasts).
#' @param clr_pseudocount Positive numeric scalar; pseudocount used for the
#'   Aitchison family. Ignored for other distances. Default \code{1}.
#'
#' @return Named list with:
#' \describe{
#'   \item{\code{dist}}{A \code{"dist"} object suitable for
#'     \code{vegan::betadisper}.}
#'   \item{\code{group}}{Factor aligned to the rows of the count matrix that
#'     produced \code{dist}, with levels \code{group_levels}.}
#' }
#'
#' @seealso \code{\link{run_analysis_job}} for the main consumer.
#' @keywords internal
#' @importFrom dplyr select all_of
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom vegan vegdist decostand
#' @importFrom stats dist
build_standard_betadisper_inputs_paper <- function(
  count_long,
  sample_metadata,
  distance_method,
  group_levels,
  clr_pseudocount = 1
) {
  missing_cols <- setdiff(c("sample_id", "taxon_id", "count"), colnames(count_long))
  if (length(missing_cols) > 0L) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!all(c("sample_id", "group") %in% colnames(sample_metadata))) {
    stop("sample_metadata must contain columns: sample_id, group")
  }

  count_matrix <- count_long %>%
    dplyr::select(dplyr::all_of(c("sample_id", "taxon_id", "count"))) %>%
    tidyr::pivot_wider(names_from = "taxon_id", values_from = "count") %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()
  count_matrix <- count_matrix[rowSums(count_matrix) > 0, , drop = FALSE]

  dist_obj <- switch(
    distance_method,
    "aitchison" = ,
    "robust.aitchison" = {
      if (clr_pseudocount <= 0) stop("clr_pseudocount must be positive.")
      vegan::vegdist(pmax(count_matrix, clr_pseudocount), method = distance_method)
    },
    "hellinger" = stats::dist(
      vegan::decostand(count_matrix, method = "hellinger", MARGIN = 1),
      method = "euclidean"
    ),
    "jaccard" = vegan::vegdist(
      vegan::decostand(count_matrix, method = "total", MARGIN = 1),
      method = "jaccard", binary = FALSE
    ),
    vegan::vegdist(count_matrix, method = distance_method)  # default branch
  )

  if (any(!is.finite(dist_obj))) {
    stop("Distance matrix contains NA/Inf values.")
  }

  group_factor <- factor(
    sample_metadata$group[match(rownames(count_matrix), sample_metadata$sample_id)],
    levels = group_levels
  )
  list(dist = dist_obj, group = group_factor)
}
