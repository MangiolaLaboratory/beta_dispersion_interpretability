# =============================================================================
# R/pipeline.R - Simulation + analysis pipeline (per-job runners and summaries)
# =============================================================================
#
# Dependencies (see @importFrom on each function):
#   dplyr   - group_by, summarise, filter, pull, bind_rows, n
#   purrr   - set_names, imap, imap_dfr
#   tibble  - tibble
#   vegan   - betadisper
#   rlang   - .data
#   internal:
#     simulate_compositional_bb, design_matrix_from_groups,
#     synthetic_sim_add_two_group_factor,
#     build_standard_betadisper_inputs_paper, permutest_betadisper_stable,
#     pretty_method_label_roc.
#
# =============================================================================

#' Run a single two-group simulation with taxon-level dispersion intercepts
#'
#' Thin wrapper around \code{\link{simulate_compositional_bb}} that fixes the
#' design to "intercept + two-group" (built by
#' \code{\link{design_matrix_from_groups}}) and lets the caller supply
#' \emph{either} a scalar/length-\code{n_taxa} dispersion intercept
#' (\code{intercept_dispersion}) \emph{or} a 2-by-\code{n_taxa} matrix
#' (\code{intercept_dispersion_matrix}) -- but not both.
#'
#' Useful when you want to vary the dispersion intercept by group (e.g. for
#' differential-dispersion null/alternative cases) without rebuilding the
#' design matrix every time.
#'
#' @param slope_vector,mu_inv_softmax,log_dispersion_assoc,sd_log_overdispersion
#'   Forwarded to \code{\link{simulate_compositional_bb}}; see its
#'   documentation.
#' @param intercept_dispersion,intercept_dispersion_matrix Exactly one must be
#'   non-\code{NULL}. The vector form is length 1 or \code{n_taxa}; the
#'   matrix form is \code{2 x n_taxa} (rows = groups).
#' @param seed Integer; forwarded as the seed to
#'   \code{simulate_compositional_bb}.
#' @param n_taxa Positive integer.
#' @param n_samples_per_group Length-1 (balanced) or length-2 integer vector
#'   forwarded to \code{\link{design_matrix_from_groups}}; the resulting
#'   total sample count is forwarded as \code{n_samples}.
#' @param library_size_mean,library_size_sd Forwarded to
#'   \code{simulate_compositional_bb}.
#'
#' @return The simulation result list returned by
#'   \code{\link{simulate_compositional_bb}}.
#'
#' @seealso \code{\link{simulate_compositional_bb}},
#'   \code{\link{design_matrix_from_groups}},
#'   \code{\link{run_simulation_job_brito}}.
#' @export
simulate_case_realistic_taxa <- function(
  slope_vector,
  mu_inv_softmax,
  log_dispersion_assoc,
  intercept_dispersion = NULL,
  intercept_dispersion_matrix = NULL,
  sd_log_overdispersion,
  seed,
  n_taxa,
  n_samples_per_group,
  library_size_mean,
  library_size_sd
) {
  if (is.null(intercept_dispersion) == is.null(intercept_dispersion_matrix)) {
    stop("Provide exactly one of intercept_dispersion (length n_taxa) or intercept_dispersion_matrix (2 x n_taxa).")
  }
  n_samples <- if (length(n_samples_per_group) == 1L) 2L * n_samples_per_group
               else sum(n_samples_per_group)

  simulate_compositional_bb(
    slope_vector                = slope_vector,
    mu_inv_softmax              = mu_inv_softmax,
    log_dispersion_assoc        = log_dispersion_assoc,
    n_taxa                      = n_taxa,
    n_samples                   = n_samples,
    sd_log_overdispersion       = sd_log_overdispersion,
    intercept_dispersion        = intercept_dispersion,
    intercept_dispersion_matrix = intercept_dispersion_matrix,
    library_size_mean           = library_size_mean,
    library_size_sd             = library_size_sd,
    design_matrix               = design_matrix_from_groups(n_samples_per_group),
    seed                        = seed
  )
}

#' Filter \code{count_long} by mean proportional abundance
#'
#' Drops taxa whose mean proportional abundance (\code{count / library_size})
#' across samples is below \code{threshold}. Returns the filtered tall table.
#'
#' If \code{threshold} is \code{NULL}, the function falls back to a
#' free variable named \code{abundance_threshold} expected to exist in the
#' calling environment (this is the convention used by the per-dataset
#' Quarto reports). Prefer passing \code{threshold} explicitly in new code.
#'
#' @param count_long Tall data frame with at least \code{taxon_id},
#'   \code{count} and \code{library_size} columns.
#' @param threshold Numeric scalar in \eqn{[0, 1]}; minimum mean proportional
#'   abundance to keep a taxon. \code{NULL} (default) looks up
#'   \code{abundance_threshold} in the calling scope.
#'
#' @return Filtered tall data frame with the same columns as \code{count_long}
#'   but possibly fewer taxa.
#'
#' @seealso \code{\link{run_permanova_job}} which uses this to pre-filter
#'   before computing distances.
#' @export
filter_by_abundance_cases <- function(count_long, threshold = NULL) {
  if (is.null(threshold)) threshold <- abundance_threshold
  # Per-row proportion, averaged within taxon via stats::ave(), then keep
  # rows whose taxon mean is >= threshold.
  taxon_mean_prop <- stats::ave(
    count_long$count / count_long$library_size,
    count_long$taxon_id,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  count_long[taxon_mean_prop >= threshold, , drop = FALSE]
}

#' Extract \code{seed} and \code{n_da_taxa} from a job-grid row
#'
#' Pulls the two identifying fields out of a single-row data frame / tibble
#' (or list) coming out of the simulation jobs grid. Coerces to data.frame
#' for uniform indexing.
#'
#' @param job_row A single-row data frame / tibble / list with \code{seed}
#'   and \code{n_da_taxa} fields.
#'
#' @return Named list with \code{seed} (numeric) and \code{n_da_taxa}
#'   (integer).
#'
#' @seealso \code{\link{run_simulation_job_brito}}.
#' @export
extract_job_params <- function(job_row) {
  if (is.null(job_row)) stop("job_row is NULL")
  job_row <- as.data.frame(job_row)
  if (nrow(job_row) == 0L) stop("job_row is an empty data.frame")

  pick <- function(field) {
    if (!field %in% names(job_row)) {
      stop("job_row missing '", field, "' column. Available columns: ",
           paste(names(job_row), collapse = ", "))
    }
    val <- job_row[[field]][1]
    if (is.na(val)) stop("Could not extract ", field, " from job_row: ", deparse(val))
    val
  }

  list(
    seed      = as.numeric(pick("seed")),
    n_da_taxa = as.integer(pick("n_da_taxa"))
  )
}

#' Run one simulation job with taxon-level dispersion intercepts (Brito variant)
#'
#' Single-job entry point used by the per-dataset \code{targets} pipelines.
#' Given a job-grid row and the empirical (composition, dispersion) profile
#' produced by \code{\link{build_simulation_params_brito}}, it:
#' \enumerate{
#'   \item Resolves \code{seed} and \code{n_da_taxa} via
#'     \code{\link{extract_job_params}}.
#'   \item Picks \code{n_da_taxa} random taxa and assigns each a
#'     non-zero slope sampled from \code{slope_effects_distribution}
#'     (mean-centred to preserve the compositional constraint).
#'   \item Forwards everything to \code{\link{simulate_compositional_bb}}
#'     with the two-group design matrix from
#'     \code{\link{design_matrix_from_groups}}.
#'   \item Records the job's \code{seed} / \code{n_da_taxa} back on the
#'     simulation list, and optionally annotates a human-readable
#'     \code{group} factor via
#'     \code{\link{synthetic_sim_add_two_group_factor}}.
#' }
#'
#' The "brito" qualifier distinguishes it from a previous variant that took a
#' scalar dispersion intercept; here every taxon may have its own (so the
#' empirical mean-dispersion association is preserved).
#'
#' Reproducibility note: a deterministic offset is added to the input
#' \code{seed} when picking DA taxa (\code{seed + 10000}) so that the DA
#' indices are stable per \code{(seed, n_da_taxa)} pair, independent of the
#' downstream count-sampling seed.
#'
#' @param job_row Single row of the job grid; see
#'   \code{\link{extract_job_params}} for accepted shapes.
#' @param composition_intercept_by_taxon Numeric vector of length
#'   \code{n_taxa}; baseline logits used as \code{mu_inv_softmax} (must sum
#'   to 0; not re-centred here).
#' @param dispersion_intercept_by_taxon Numeric scalar or length-\code{n_taxa}
#'   vector; per-taxon dispersion intercept \eqn{c} (\eqn{\log\sigma} at
#'   \eqn{\mu = 0}). Scalars are recycled.
#' @param mean_dispersion_assoc_slope Numeric scalar; the slope \eqn{k} in
#'   \eqn{\log\sigma = -k\,\log\mu + c}.
#' @param slope_effects_distribution Numeric vector of empirical slope effect
#'   sizes; sampled with replacement to populate DA taxa.
#' @param sd_log_overdispersion Non-negative numeric scalar; SD of additional
#'   log-sigma noise.
#' @param n_taxa,n_samples,n_samples_per_group,library_size_mean,library_size_sd
#'   Standard simulator inputs; \code{n_samples_per_group} also feeds
#'   \code{\link{design_matrix_from_groups}}.
#' @param group_levels Optional length-2 character vector. If supplied, the
#'   simulation result is post-processed with
#'   \code{\link{synthetic_sim_add_two_group_factor}}.
#'
#' @return The simulation result list from
#'   \code{\link{simulate_compositional_bb}} with two extra fields
#'   \code{seed} and \code{n_da_taxa} appended, and (when \code{group_levels}
#'   is supplied) the \code{group} factor populated.
#'
#' @seealso \code{\link{extract_job_params}},
#'   \code{\link{simulate_compositional_bb}},
#'   \code{\link{run_analysis_job}},
#'   \code{\link{build_simulation_params_brito}}.
#' @export
run_simulation_job_brito <- function(
  job_row,
  composition_intercept_by_taxon,
  dispersion_intercept_by_taxon,
  mean_dispersion_assoc_slope,
  slope_effects_distribution,
  sd_log_overdispersion,
  n_taxa,
  n_samples,
  n_samples_per_group,
  library_size_mean,
  library_size_sd,
  group_levels = NULL
) {
  if (length(composition_intercept_by_taxon) != n_taxa) {
    stop("composition_intercept_by_taxon must have length n_taxa.")
  }
  if (length(slope_effects_distribution) == 0L) {
    stop("slope_effects_distribution is empty.")
  }
  alpha_intercept_by_taxon <- rep_len(dispersion_intercept_by_taxon, n_taxa)

  params <- extract_job_params(job_row)
  seed <- params$seed
  n_da <- params$n_da_taxa

  # Build slope vector: n_da random taxa get sampled empirical effects.
  # Seed offset keeps DA picks stable per (seed, n_da) independent of count sampling.
  slope_vec <- rep(0, n_taxa)
  if (n_da > 0) {
    set.seed(seed + 10000)
    da_idx <- sample.int(n_taxa, size = min(n_da, n_taxa))
    slope_vec[da_idx] <- sample(slope_effects_distribution, length(da_idx),
                                replace = TRUE)
    slope_vec <- slope_vec - mean(slope_vec)
  }

  sim_r <- simulate_compositional_bb(
    slope_vector = slope_vec,
    mu_inv_softmax = composition_intercept_by_taxon,
    log_dispersion_assoc = mean_dispersion_assoc_slope,
    n_taxa = n_taxa,
    n_samples = n_samples,
    sd_log_overdispersion = sd_log_overdispersion,
    intercept_dispersion = alpha_intercept_by_taxon,
    library_size_mean = library_size_mean,
    library_size_sd = library_size_sd,
    design_matrix = design_matrix_from_groups(n_samples_per_group),
    seed = seed
  )
  sim_r$seed <- seed
  sim_r$n_da_taxa <- n_da

  if (!is.null(group_levels)) {
    sim_r <- synthetic_sim_add_two_group_factor(sim_r, group_levels)
  }
  sim_r
}

#' Run betadisper across several standard distances on one simulation
#'
#' Single source of truth for the per-simulation betadisper analysis used by
#' all per-dataset benchmarks. For each requested distance method it:
#' \enumerate{
#'   \item Builds a \code{(dist, group)} pair via
#'     \code{\link{build_standard_betadisper_inputs_paper}}.
#'   \item Fits \code{vegan::betadisper}.
#'   \item Runs an adaptive permutation test via
#'     \code{\link{permutest_betadisper_stable}} (so the headline p-value's
#'     Monte Carlo SE is bounded).
#'   \item Captures the per-group mean distances-to-centroid (\code{group1},
#'     \code{group2}, \code{delta}) and the per-sample distances.
#' }
#'
#' Two output modes via \code{return_full}:
#' \describe{
#'   \item{\code{FALSE} (default)}{Returns a long tibble with one row per
#'     method (\code{seed}, \code{n_da_taxa}, \code{method}, \code{p_value}).
#'     This is what the FPR / TPR summaries consume.}
#'   \item{\code{TRUE}}{Returns the full \code{method_results} list plus
#'     \code{seed} / \code{n_da_taxa}, retaining the fitted \code{bd}
#'     objects and per-sample distances. Used for plotting.}
#' }
#'
#' The function requires that the caller has populated
#' \code{sample_metadata$group} as a factor (typically via
#' \code{\link{synthetic_sim_add_two_group_factor}} run from the Quarto
#' document with the appropriate per-dataset labels).
#'
#' @param sim_r Simulation result list as returned by
#'   \code{\link{run_simulation_job_brito}} (or compatible), with
#'   \code{seed}, \code{n_da_taxa} and \code{sample_metadata$group} fields.
#' @param max_perm,batch,target_se Forwarded to
#'   \code{\link{permutest_betadisper_stable}}. Defaults are tuned for the
#'   benchmark sweeps.
#' @param return_full Logical; if \code{TRUE} returns the rich list instead
#'   of the FPR-summary tibble. Default \code{FALSE}.
#' @param standard_distance_methods Character vector of distance methods
#'   passed to \code{\link{build_standard_betadisper_inputs_paper}}. Default
#'   covers Bray-Curtis, Aitchison family, Jaccard and Hellinger.
#'
#' @return Either a tibble (\code{return_full = FALSE}) or a list
#'   (\code{return_full = TRUE}); see Description.
#'
#' @seealso \code{\link{lighten_sweep_analysis_result_for_fpr}},
#'   \code{\link{summarise_fpr}}.
#' @export
#' @importFrom tibble tibble
#' @importFrom vegan betadisper
run_analysis_job <- function(
  sim_r,
  max_perm = 19999,
  batch = 2000,
  target_se = 0.005,
  return_full = FALSE,
  standard_distance_methods = c("bray", "aitchison", "robust.aitchison", "jaccard", "hellinger")
) {
  if (!("group" %in% names(sim_r$sample_metadata))) {
    stop("sim_r$sample_metadata must include factor column `group`. ",
         "Add it in the Quarto document after simulation (see benchmark_*_paper_analyses*.qmd).")
  }
  if (!is.factor(sim_r$sample_metadata$group)) {
    stop("sim_r$sample_metadata$group must be a factor (set in Quarto after simulation).")
  }
  group_levels <- levels(sim_r$sample_metadata$group)
  method_keys  <- paste0("permdisp_", standard_distance_methods)

  method_results <- setNames(
    lapply(standard_distance_methods, function(dm) {
      inp <- build_standard_betadisper_inputs_paper(
        count_long = sim_r$count_long,
        sample_metadata = sim_r$sample_metadata,
        distance_method = dm,
        group_levels = group_levels
      )
      bd   <- vegan::betadisper(inp$dist, inp$group)
      stab <- permutest_betadisper_stable(bd, max_perm = max_perm, batch = batch,
                                          target_se = target_se)
      list(
        method_key = paste0("permdisp_", dm),
        bd         = bd,
        group1    = unname(bd$group.distances[1]),
        group2    = unname(bd$group.distances[2]),
        delta     = unname(bd$group.distances[2] - bd$group.distances[1]),
        p_value   = unname(stab$p_value),
        p_se      = unname(stab$p_se),
        n_perm    = unname(stab$n_perm),
        distances = unname(bd$distances),
        group     = bd$group
      )
    }),
    method_keys
  )

  if (return_full) {
    list(method_results = method_results,
         seed = sim_r$seed, n_da_taxa = sim_r$n_da_taxa)
  } else {
    tibble::tibble(
      seed      = sim_r$seed,
      n_da_taxa = sim_r$n_da_taxa,
      method    = vapply(method_keys, pretty_method_label_roc, character(1L), USE.NAMES = FALSE),
      p_value   = vapply(method_results, `[[`, numeric(1L), "p_value", USE.NAMES = FALSE)
    )
  }
}

#' Shrink a \code{run_analysis_job(return_full = TRUE)} result for FPR storage
#'
#' Keeps only the fields that \code{\link{summarise_fpr}} actually consumes
#' (\code{seed}, \code{n_da_taxa}, optional \code{threshold}, and a
#' \code{method_results} list reduced to a per-method \code{list(p_value =
#' ...)} stub). Drops the fitted \code{betadisper} objects and per-sample
#' distances, which would otherwise dominate the on-disk size of a
#' \code{targets} \code{pattern = map()} sweep.
#'
#' If the input is not a \code{run_analysis_job} list (i.e. lacks
#' \code{method_results}) it is returned unchanged so the function is safe
#' to map blindly over heterogeneous worker outputs.
#'
#' @param x A list returned by \code{\link{run_analysis_job}} with
#'   \code{return_full = TRUE}, optionally annotated with \code{threshold}.
#'
#' @return A trimmed list with the same top-level structure but per-method
#'   \code{p_value}-only entries.
#'
#' @seealso \code{\link{run_analysis_job}}, \code{\link{summarise_fpr}}.
#' @export
#' @importFrom purrr imap
lighten_sweep_analysis_result_for_fpr <- function(x) {
  if (!is.list(x) || !("method_results" %in% names(x))) {
    return(x)
  }
  method_results_light <- purrr::imap(x$method_results, ~ list(p_value = .x$p_value))
  out <- list(
    method_results = method_results_light,
    seed = x$seed,
    n_da_taxa = x$n_da_taxa
  )
  if ("threshold" %in% names(x)) {
    out$threshold <- x$threshold
  }
  out
}

#' Summarise sweep results into a false-positive rate (FPR) table
#'
#' Aggregates the outputs of a null-setting sweep (\code{n_da_taxa == 0}) or
#' a mixed sweep into a per-\code{(n_da_taxa, method, threshold)} table with:
#' \itemize{
#'   \item \code{n_sims}: count of replicates contributing,
#'   \item \code{n_fp}: count with \code{p < 0.05},
#'   \item \code{FPR}: \code{n_fp / n_sims}.
#' }
#'
#' Two input shapes are supported:
#' \describe{
#'   \item{Auto-bound tibble}{A single tibble with columns \code{seed},
#'     \code{n_da_taxa}, \code{method}, \code{p_value} (this is what
#'     \code{targets} produces when worker outputs are tibbles).}
#'   \item{List of per-job lists}{Each element is a
#'     \code{run_analysis_job(return_full = TRUE)} result (optionally
#'     trimmed by \code{\link{lighten_sweep_analysis_result_for_fpr}})
#'     containing \code{method_results}, \code{seed}, \code{n_da_taxa} and
#'     optional \code{threshold}. Per-method rows are pivoted out of
#'     \code{method_results} before aggregation.}
#' }
#'
#' Raises informative errors for empty input or unrecognised element types so
#' a misconfigured sweep fails fast rather than silently returning an empty
#' table.
#'
#' @param sweep_results Either a tibble of the shape described above, or a
#'   list of per-job results (see Description).
#'
#' @return Tibble grouped by \code{(n_da_taxa, method, threshold)} with the
#'   columns \code{n_sims}, \code{n_fp}, \code{FPR}.
#'
#' @seealso \code{\link{run_analysis_job}},
#'   \code{\link{lighten_sweep_analysis_result_for_fpr}}.
#' @export
#' @importFrom dplyr group_by summarise bind_rows n
#' @importFrom purrr imap_dfr
#' @importFrom tibble tibble
#' @importFrom rlang .data
summarise_fpr <- function(sweep_results) {
  if (is.data.frame(sweep_results) &&
      all(c("seed", "n_da_taxa", "method", "p_value") %in% names(sweep_results))) {
    fpr_data <- sweep_results
  } else {
    if (length(sweep_results) == 0L) stop("sweep_results is empty")
    fpr_data <- dplyr::bind_rows(lapply(sweep_results, function(x) {
      purrr::imap_dfr(x$method_results, function(mr, method_key) {
        tibble::tibble(
          seed      = x$seed,
          n_da_taxa = x$n_da_taxa,
          threshold = if ("threshold" %in% names(x)) x$threshold else NA_real_,
          method    = method_key,
          p_value   = mr$p_value
        )
      })
    }))
  }

  fpr_data %>%
    dplyr::group_by(.data$n_da_taxa, .data$method, .data$threshold) %>%
    dplyr::summarise(
      n_sims = dplyr::n(),
      n_fp   = sum(.data$p_value < 0.05, na.rm = TRUE),
      FPR    = n_fp / n_sims,
      .groups = "drop"
    )
}

