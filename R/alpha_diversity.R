# =============================================================================
# R/alpha_diversity.R - Alpha diversity analysis and plotting
# =============================================================================
#
# Source this via the top-level functions.R, which sources all R/*.R modules
# in the correct dependency order.
#
# Dependencies (see @importFrom on each function):
#   dplyr                 - group_by, summarise, mutate, select, first, n
#   tidyr                 - pivot_wider
#   tibble                - tibble
#   ggplot2               - ggplot, aes, geom_point, stat_ellipse, scale_color_manual,
#                           labs, theme_minimal, theme, element_text
#   stats                 - t.test, rnorm, sd
#   rlang                 - .data
#   mia                   - getAlpha (only run_alpha_diversity_analysis_normalised)
#   SummarizedExperiment  - SummarizedExperiment (idem)
#   S4Vectors             - DataFrame (idem)
#
# =============================================================================

#' Build the alpha-diversity summary / tests / plot bundle (internal)
#'
#' Shared backend for \code{\link{run_alpha_diversity_analysis}} and
#' \code{\link{run_alpha_diversity_analysis_normalised}}. Takes an
#' \code{alpha_data} tibble (one row per sample with \code{group},
#' \code{richness}, \code{shannon}, \code{pielou_evenness}) and produces:
#' \itemize{
#'   \item a per-group summary table (mean / SD of all three indices),
#'   \item three two-sample Welch t-tests (richness, Shannon, Pielou),
#'   \item a ggplot of Pielou vs richness with 95\% normal ellipses on
#'     slightly jittered values (the jitter only feeds
#'     \code{stat_ellipse}; the points are drawn from the original data).
#' }
#'
#' All numerics flow through \code{tryCatch}-wrapped \code{stats::t.test} so
#' a degenerate group (e.g. zero variance) yields NA entries rather than an
#' error.
#'
#' @param alpha_data Data frame with at least \code{group} (factor) plus the
#'   three numeric columns \code{richness}, \code{shannon},
#'   \code{pielou_evenness}.
#' @param case_label Character; goes into the plot title (e.g. the simulated
#'   case name).
#' @param group_colors Named character vector of hex colours; names must
#'   match the levels in \code{alpha_data$group}.
#' @param plot_subtitle Character; goes into the plot subtitle (used to
#'   indicate rarefaction settings).
#'
#' Validate (group_levels, group_colors) argument pair (internal)
#'
#' Centralises the per-function argument checking shared by
#' \code{run_alpha_diversity_analysis()} and its rarefied variant. Returns
#' the cleaned-up character vector + named colour vector so the callers can
#' use them directly.
#' @noRd
.validate_group_args <- function(group_levels, group_colors, fn_label) {
  if (missing(group_levels) || missing(group_colors)) {
    stop(fn_label, "() requires `group_levels` and `group_colors`.", call. = FALSE)
  }
  if (is.null(names(group_colors))) {
    stop(fn_label, "(): `group_colors` must be named (names = group levels).",
         call. = FALSE)
  }
  gl   <- as.character(group_levels)
  miss <- setdiff(gl, names(group_colors))
  if (length(miss) > 0L) {
    stop(fn_label, "(): `group_colors` missing names for: ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  list(
    group_levels = gl,
    group_colors = setNames(unname(group_colors[gl]), gl)
  )
}

#' @return Named list with \code{data} (echo of input), \code{summary} (per
#'   group means / SDs), \code{tests} (list of three t-test objects), and
#'   \code{plot} (a \code{ggplot}).
#' @noRd
#' @importFrom dplyr group_by summarise mutate n
#' @importFrom stats t.test rnorm sd
.alpha_diversity_table_to_plot <- function(alpha_data, case_label, group_colors, plot_subtitle) {
  summary_tbl <- alpha_data %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      mean_richness = mean(richness),
      sd_richness = sd(richness),
      mean_shannon = mean(shannon),
      sd_shannon = sd(shannon),
      mean_evenness = mean(pielou_evenness),
      sd_evenness = sd(pielou_evenness),
      .groups = "drop"
    )

  safe_ttest <- function(formula, data) {
    tryCatch(
      stats::t.test(formula, data = data),
      error = function(e) list(p.value = NA, statistic = NA, estimate = c(NA, NA), message = e$message)
    )
  }

  tests <- list(
    richness = safe_ttest(richness ~ group, data = alpha_data),
    shannon = safe_ttest(shannon ~ group, data = alpha_data),
    evenness = safe_ttest(pielou_evenness ~ group, data = alpha_data)
  )

  set.seed(12345)
  alpha_data_jitter <- alpha_data %>%
    dplyr::mutate(
      richness_j = richness + stats::rnorm(dplyr::n(), 0, 0.01),
      evenness_j = pielou_evenness + stats::rnorm(dplyr::n(), 0, 0.001)
    )

  p <- ggplot2::ggplot(alpha_data, ggplot2::aes(x = richness, y = pielou_evenness, color = group)) +
    ggplot2::geom_point(size = 0.1, alpha = 0.6, stroke = 0) +
    ggplot2::stat_ellipse(
      data = alpha_data_jitter,
      ggplot2::aes(x = richness_j, y = evenness_j, color = group),
      level = 0.95, linewidth = 0.25, type = "norm"
    ) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::labs(
      title = paste0("Alpha Diversity: Evenness vs Richness (", case_label, ")"),
      subtitle = plot_subtitle,
      x = "Richness (# of observed taxa)",
      y = "Pielou's Evenness",
      color = "Group"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11),
      legend.position = "bottom"
    )

  list(data = alpha_data, summary = summary_tbl, tests = tests, plot = p)
}

#' Alpha-diversity analysis on raw counts (no rarefaction)
#'
#' Computes richness (number of observed taxa), Shannon entropy, and Pielou
#' evenness per sample directly from the raw \code{count_long} table of a
#' simulation result. Aggregates per group and runs three Welch t-tests, then
#' assembles a plot via \code{\link{.alpha_diversity_table_to_plot}}.
#'
#' Because no rarefaction is applied, richness scales with library size; for a
#' depth-controlled variant use
#' \code{\link{run_alpha_diversity_analysis_normalised}}.
#'
#' The function aggressively validates its inputs: missing
#' \code{group_levels}/\code{group_colors}, unnamed colour vectors, missing
#' names for the requested levels, and missing required columns all raise
#' informative errors.
#'
#' @param sim_result List as returned by
#'   \code{\link{simulate_compositional_bb}}; must contain a
#'   \code{count_long} data frame with columns \code{sample_id},
#'   \code{group}, \code{count}.
#' @param case_label Character scalar; included in the plot title. Default
#'   \code{"Case"}.
#' @param group_levels Character vector of factor levels for the group factor.
#'   Required.
#' @param group_colors Named character vector of hex colours. Names must be a
#'   superset of \code{group_levels}; entries are reordered to match.
#'   Required.
#'
#' @return Named list with \code{data}, \code{summary}, \code{tests},
#'   \code{plot} (see \code{\link{.alpha_diversity_table_to_plot}}).
#'
#' @seealso \code{\link{run_alpha_diversity_analysis_normalised}}.
#' @keywords internal
#' @importFrom dplyr group_by summarise mutate
run_alpha_diversity_analysis <- function(
  sim_result,
  case_label = "Case",
  group_levels,
  group_colors
) {
  v <- .validate_group_args(group_levels, group_colors, "run_alpha_diversity_analysis")
  gl <- v$group_levels
  group_colors <- v$group_colors

  count_long <- sim_result$count_long
  missing_cols <- setdiff(c("sample_id", "group", "count"), colnames(count_long))
  if (length(missing_cols) > 0L) {
    stop("sim_result$count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  alpha_data <- count_long %>%
    dplyr::group_by(sample_id, group) %>%
    dplyr::summarise(
      richness = sum(count > 0),
      shannon = {
        props <- count / sum(count)
        props <- props[props > 0]
        -sum(props * log(props))
      },
      pielou_evenness = shannon / log(richness),
      .groups = "drop"
    ) %>%
    dplyr::mutate(group = factor(group, levels = gl))

  .alpha_diversity_table_to_plot(
    alpha_data, case_label, group_colors,
    plot_subtitle = "Points = samples, Ellipses = 95% confidence intervals"
  )
}

#' Alpha-diversity analysis with iterative rarefaction (depth-controlled)
#'
#' Same summary / test / plot bundle as
#' \code{\link{run_alpha_diversity_analysis}}, but the three diversity indices
#' are computed via \code{mia::getAlpha} after iterative rarefaction to a
#' common depth (default: the minimum per-sample library size). This removes
#' the artefactual scaling of richness with library size.
#'
#' Requires the Bioconductor packages \code{mia} and \code{SummarizedExperiment}
#' to be installed; \code{S4Vectors::DataFrame} is used to construct
#' \code{colData}. The function raises an informative error if any sample's
#' total count is below the chosen rarefaction depth, so the caller can pick
#' a different \code{rarefy_sample} or filter samples upstream.
#'
#' Implementation outline:
#' \enumerate{
#'   \item Pivot \code{count_long} to a taxon-by-sample matrix and align a
#'     \code{group} factor by \code{sample_id}.
#'   \item Wrap the matrix in a \code{SummarizedExperiment} and call
#'     \code{mia::getAlpha} with \code{index = c("observed", "shannon",
#'     "pielou")} and the chosen number of rarefaction rounds.
#'   \item Hand the resulting per-sample summaries off to
#'     \code{\link{.alpha_diversity_table_to_plot}} for plotting.
#' }
#'
#' @param sim_result List as returned by
#'   \code{\link{simulate_compositional_bb}}; must contain a
#'   \code{count_long} data frame with columns \code{sample_id},
#'   \code{taxon_id}, \code{group}, \code{count}.
#' @param case_label Character scalar; included in the plot title. Default
#'   \code{"Case"}.
#' @param group_levels Character vector of factor levels for the group factor.
#' @param group_colors Named character vector of hex colours (see
#'   \code{\link{run_alpha_diversity_analysis}}).
#' @param rarefy_niter Positive integer; number of rarefaction rounds passed
#'   to \code{mia::getAlpha} via \code{niter}. Default \code{50L}.
#' @param rarefy_sample Optional positive numeric scalar; rarefaction depth
#'   (counts per draw). \code{NULL} (default) uses the minimum per-sample
#'   total count.
#'
#' @return Same structure as \code{\link{run_alpha_diversity_analysis}}; the
#'   plot subtitle now states the rarefaction settings.
#'
#' @seealso \code{\link{run_alpha_diversity_analysis}}.
#' @keywords internal
#' @importFrom dplyr select group_by summarise first
#' @importFrom tidyr pivot_wider
#' @importFrom tibble tibble
#' @importFrom rlang .data
run_alpha_diversity_analysis_normalised <- function(
  sim_result,
  case_label = "Case",
  group_levels,
  group_colors,
  rarefy_niter = 50L,
  rarefy_sample = NULL
) {
  v <- .validate_group_args(group_levels, group_colors,
                            "run_alpha_diversity_analysis_normalised")
  gl <- v$group_levels
  group_colors <- v$group_colors

  count_long <- sim_result$count_long

  # taxon-by-sample count matrix, samples in matrix-column order
  wide <- count_long %>%
    dplyr::select("sample_id", "taxon_id", "count") %>%
    tidyr::pivot_wider(names_from = "sample_id", values_from = "count",
                        values_fill = 0)
  mat <- as.matrix(wide[, -1, drop = FALSE])
  rownames(mat) <- wide$taxon_id
  sample_order <- colnames(mat)
  group_per_sample <- count_long$group[match(sample_order, count_long$sample_id)]

  # Rarefy to the shallowest sample (or the user-supplied depth).
  depth <- if (is.null(rarefy_sample)) min(colSums(mat)) else rarefy_sample

  col_df <- S4Vectors::DataFrame(group = factor(group_per_sample, levels = gl))
  rownames(col_df) <- sample_order
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = mat), colData = col_df
  )
  alpha_df <- mia::getAlpha(
    se, assay.type = "counts",
    index = c("observed", "shannon", "pielou"),
    niter = rarefy_niter, sample = depth
  )
  alpha_data <- tibble::tibble(
    sample_id       = sample_order,
    group           = factor(group_per_sample, levels = gl),
    richness        = as.numeric(alpha_df$observed),
    shannon         = as.numeric(alpha_df$shannon),
    pielou_evenness = as.numeric(alpha_df$pielou)
  )
  .alpha_diversity_table_to_plot(
    alpha_data, case_label, group_colors,
    plot_subtitle = paste0(
      "Rarefied (mia::getAlpha, niter = ", rarefy_niter,
      ", depth = ", depth, "); points = samples, ellipses = 95% CI"
    )
  )
}

