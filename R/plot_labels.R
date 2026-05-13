# =============================================================================
# R/plot_labels.R - Pretty method labels for plots
# =============================================================================
#
# Source this via the top-level functions.R, which sources all R/*.R modules
# in the correct dependency order.
#
# Pure base-R utilities; no external package dependencies.
#
# =============================================================================

#' Map an internal method key to a human-readable label
#'
#' Translates the per-method identifier strings produced by the analysis
#' pipeline (e.g. \code{"permdisp_bray"}, \code{"hellinger"}) into
#' typographically clean display labels used on plots (e.g.
#' \code{"Bray-Curtis"}, \code{"Hellinger"}). Unknown keys are returned
#' unchanged so the function is safe to use as a generic
#' \code{ggplot2::scale_*(labels = ...)} hook.
#'
#' Supported keys (both the raw distance name and the \code{permdisp_}-prefixed
#' variant are recognised):
#' \itemize{
#'   \item \code{bray} / \code{permdisp_bray} -> "Bray-Curtis"
#'   \item \code{jaccard} / \code{permdisp_jaccard} -> "Jaccard"
#'   \item \code{hellinger} / \code{permdisp_hellinger} -> "Hellinger"
#'   \item \code{aitchison} / \code{permdisp_aitchison} -> "Aitchison"
#'   \item \code{robust.aitchison} / \code{permdisp_robust.aitchison} -> "robust Aitchison"
#'   \item \code{permdisp_euclidean} -> "Euclidean"
#' }
#'
#' @param method_key Character scalar; the internal method identifier.
#'
#' @return Character scalar; the display label, or \code{method_key} if no
#'   match was found.
#'
#' @seealso \code{\link{pretty_method_label_roc}},
#'   \code{\link{pretty_method_label_box}}.
#' @keywords internal
pretty_method_label <- function(method_key) {
  labels <- c(
    bray                      = "Bray\u2013Curtis",
    permdisp_bray             = "Bray\u2013Curtis",
    jaccard                   = "Jaccard",
    permdisp_jaccard          = "Jaccard",
    hellinger                 = "Hellinger",
    permdisp_hellinger        = "Hellinger",
    aitchison                 = "Aitchison",
    permdisp_aitchison        = "Aitchison",
    robust.aitchison          = "robust Aitchison",
    permdisp_robust.aitchison = "robust Aitchison",
    permdisp_euclidean        = "Euclidean"
  )
  if (method_key %in% names(labels)) unname(labels[method_key]) else method_key
}

#' Pretty method label for ROC-curve legends (currently an alias)
#'
#' Thin wrapper around \code{\link{pretty_method_label}} kept as a named
#' indirection so the ROC plots can later diverge from the box-plot labels
#' (e.g. with annotated AUC) without touching call sites.
#'
#' @param method_key Character scalar; internal method identifier.
#' @return Character scalar; same as \code{\link{pretty_method_label}}.
#'
#' @seealso \code{\link{pretty_method_label}}.
#' @keywords internal
pretty_method_label_roc <- function(method_key) {
  pretty_method_label(method_key)
}

#' Pretty method label for box-plot legends (currently an alias)
#'
#' Thin wrapper around \code{\link{pretty_method_label}} kept as a named
#' indirection (see \code{\link{pretty_method_label_roc}}).
#'
#' @param method_key Character scalar; internal method identifier.
#' @return Character scalar; same as \code{\link{pretty_method_label}}.
#'
#' @seealso \code{\link{pretty_method_label}}.
#' @keywords internal
pretty_method_label_box <- function(method_key) {
  pretty_method_label(method_key)
}
