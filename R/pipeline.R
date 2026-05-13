# =============================================================================
# R/pipeline.R - Simulation + analysis pipeline (per-job runners and summaries)
# =============================================================================
#
# Source this via the top-level functions.R, which sources all R/*.R modules
# in the correct dependency order.
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
#' @keywords internal
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
  n_groups <- 2
  n_samples <- if (length(n_samples_per_group) == 1L) {
    n_samples_per_group * n_groups
  } else {
    sum(n_samples_per_group)
  }
  has_vec <- !is.null(intercept_dispersion)
  has_mat <- !is.null(intercept_dispersion_matrix)
  if (has_vec == has_mat) {
    stop("Provide exactly one of intercept_dispersion (length n_taxa) or intercept_dispersion_matrix (2 x n_taxa).")
  }
  simulate_compositional_bb( # nolint
    slope_vector = slope_vector,
    mu_inv_softmax = mu_inv_softmax,
    log_dispersion_assoc = log_dispersion_assoc,
    n_taxa = n_taxa,
    n_samples = n_samples,
    sd_log_overdispersion = sd_log_overdispersion,
    intercept_dispersion = intercept_dispersion,
    intercept_dispersion_matrix = intercept_dispersion_matrix,
    library_size_mean = library_size_mean,
    library_size_sd = library_size_sd,
    design_matrix = design_matrix_from_groups(n_samples_per_group),
    seed = seed
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
#' @keywords internal
#' @importFrom dplyr group_by summarise filter pull
#' @importFrom rlang .data
filter_by_abundance_cases <- function(count_long, threshold = NULL) {
  if (is.null(threshold)) {
    threshold <- abundance_threshold
  }
  mean_abundance <- count_long %>%
    dplyr::group_by(.data$taxon_id) %>%
    dplyr::summarise(
      mean_prop = mean(.data$count / .data$library_size, na.rm = TRUE),
      .groups = "drop"
    )
  keep_taxa <- mean_abundance %>%
    dplyr::filter(.data$mean_prop >= threshold) %>%
    dplyr::pull(.data$taxon_id)
  count_long %>% dplyr::filter(.data$taxon_id %in% keep_taxa)
}

#' Robustly extract \code{seed} and \code{n_da_taxa} from a job-grid row
#'
#' Shared helper used by the simulation runners. Accepts data frames,
#' tibbles, and lists (the latter arising from \code{targets} workers that
#' round-trip the row through serialisation), and copes with column-name
#' mangling that sometimes happens when nested tibbles are coerced
#' (\code{X1.seed}, \code{X2.seed}, ...).
#'
#' Failure modes are all wrapped in informative \code{stop()} calls:
#' \itemize{
#'   \item \code{NULL} or empty row,
#'   \item neither a \code{seed} (or \code{X\\d+\\.seed}) column nor an
#'     \code{n_da_taxa} (or \code{X\\d+\\.n_da_taxa}) column,
#'   \item missing/\code{NA} values after extraction.
#' }
#'
#' @param job_row A single-row data frame / tibble / list with at minimum
#'   \code{seed} and \code{n_da_taxa} fields.
#'
#' @return Named list with \code{seed} (numeric) and \code{n_da_taxa}
#'   (integer).
#'
#' @seealso \code{\link{run_simulation_job_brito}}.
#' @keywords internal
extract_job_params <- function(job_row) {
  # Normalize job_row to data.frame to handle serialization issues
  if (is.null(job_row)) {
    stop("job_row is NULL")
  }
  
  # Convert to data.frame if needed (handles tibbles, lists, etc.)
  if (!is.data.frame(job_row)) {
    tryCatch({
      job_row <- as.data.frame(job_row)
    }, error = function(e) {
      stop("Could not convert job_row to data.frame. Class: ", paste(class(job_row), collapse=", "), ". Error: ", e$message)
    })
  }
  
  if (is.data.frame(job_row)) {
    if (nrow(job_row) == 0) {
      stop("job_row is an empty data.frame")
    }
    # Handle nested column names (e.g., "X2.seed" instead of "seed")
    seed_col <- NULL
    n_da_col <- NULL
    
    if ("seed" %in% names(job_row)) {
      seed_col <- "seed"
    } else {
      # Look for nested column names like "X2.seed", "X1.seed", etc.
      seed_cols <- grep("^X\\d+\\.seed$|seed$", names(job_row), value = TRUE)
      if (length(seed_cols) > 0) {
        seed_col <- seed_cols[1]
      }
    }
    
    if ("n_da_taxa" %in% names(job_row)) {
      n_da_col <- "n_da_taxa"
    } else {
      # Look for nested column names
      n_da_cols <- grep("^X\\d+\\.n_da_taxa$|n_da_taxa$", names(job_row), value = TRUE)
      if (length(n_da_cols) > 0) {
        n_da_col <- n_da_cols[1]
      }
    }
    
    if (is.null(seed_col)) {
      stop("job_row missing 'seed' column. Available columns: ", paste(names(job_row), collapse=", "))
    }
    if (is.null(n_da_col)) {
      stop("job_row missing 'n_da_taxa' column. Available columns: ", paste(names(job_row), collapse=", "))
    }
    
    seed_val <- job_row[[seed_col]][1]
    n_da_val <- job_row[[n_da_col]][1]
    if (is.null(seed_val) || length(seed_val) == 0 || is.na(seed_val)) {
      stop("Could not extract seed from job_row: ", deparse(seed_val))
    }
    if (is.null(n_da_val) || length(n_da_val) == 0 || is.na(n_da_val)) {
      stop("Could not extract n_da_taxa from job_row: ", deparse(n_da_val))
    }
    seed <- as.numeric(seed_val)
    n_da <- as.integer(n_da_val)
  } else if (is.list(job_row)) {
    # Handle case where tibble was converted to list during serialization
    # Check if it has tibble-like structure (seed and n_da_taxa as vectors)
    if ("seed" %in% names(job_row) && "n_da_taxa" %in% names(job_row)) {
      seed_val <- job_row$seed
      n_da_val <- job_row$n_da_taxa
      # If they're vectors, take first element
      if (length(seed_val) > 0) seed_val <- seed_val[1]
      if (length(n_da_val) > 0) n_da_val <- n_da_val[1]
      seed <- as.numeric(seed_val)
      n_da <- as.integer(n_da_val)
    } else {
      # Try to convert back to data.frame if possible
      tryCatch({
        job_df <- as.data.frame(job_row)
        if ("seed" %in% names(job_df) && "n_da_taxa" %in% names(job_df)) {
          seed <- as.numeric(job_df$seed[1])
          n_da <- as.integer(job_df$n_da_taxa[1])
        } else {
          stop("job_row list missing 'seed' or 'n_da_taxa'. Available names: ", paste(names(job_row), collapse=", "))
        }
      }, error = function(e) {
        stop("job_row list missing 'seed' or 'n_da_taxa'. Available names: ", paste(names(job_row), collapse=", "), ". Error: ", e$message)
      })
    }
  } else {
    stop("job_row must be a data.frame or list, got: ", paste(class(job_row), collapse=", "))
  }
  
  # Final validation
  if (is.na(seed) || length(seed) == 0) {
    stop("seed is invalid after extraction: ", deparse(seed))
  }
  if (is.na(n_da) || length(n_da) == 0) {
    stop("n_da_taxa is invalid after extraction: ", deparse(n_da))
  }
  
  list(seed = seed, n_da_taxa = n_da)
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
#'   vector; per-taxon dispersion intercept \eqn{c} (\code{log\sigma} at
#'   \code{mu = 0}). Scalars are recycled.
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
#' @keywords internal
run_simulation_job_brito <- function(
  job_row,                      # One row from jobs grid (must include seed, n_da_taxa)
  composition_intercept_by_taxon,  # Composition intercept by taxon (baseline logits; length n_taxa)
  dispersion_intercept_by_taxon,   # Dispersion intercept by taxon (alpha/log-sigma; length n_taxa; scalar allowed)
  mean_dispersion_assoc_slope,     # Mean-variance association slope k in log(sigma) = -k*log(mu) + alpha
  slope_effects_distribution,   # Empirical DA effect-size pool to sample non-zero taxa slopes
  sd_log_overdispersion,        # Additional random log-sigma noise SD
  n_taxa,                       # Number of taxa (must match composition/dispersion intercept vector lengths)
  n_samples,                    # Total simulated samples
  n_samples_per_group,          # Samples per group (used to build design matrix)
  library_size_mean,            # Mean sequencing depth
  library_size_sd,              # SD sequencing depth
  group_levels = NULL            # If non-NULL, length-2 vector maps Group 0/1 to factor levels (caller supplies labels)
) {
  # Aliases with explicit intent for readability.
  mu_inv_softmax_base <- composition_intercept_by_taxon
  log_dispersion_assoc <- mean_dispersion_assoc_slope
  alpha_intercept_by_taxon <- dispersion_intercept_by_taxon

  # Composition and dispersion intercepts are taxon-length vectors.
  if (length(mu_inv_softmax_base) != n_taxa) {
    stop(
      sprintf(
        "composition_intercept_by_taxon must be length n_taxa=%d (got %d).",
        n_taxa, length(mu_inv_softmax_base)
      )
    )
  }
  if (length(alpha_intercept_by_taxon) == 1) {
    alpha_intercept_by_taxon <- rep(alpha_intercept_by_taxon, n_taxa)
  }
  if (length(alpha_intercept_by_taxon) != n_taxa) {
    stop(
      sprintf(
        paste0(
          "dispersion_intercept_by_taxon must be length n_taxa=%d ",
          "(or scalar, which is expanded); got %d."
        ),
        n_taxa, length(alpha_intercept_by_taxon)
      )
    )
  }

  # Extract seed and n_da_taxa
  params <- extract_job_params(job_row)
  seed <- params$seed
  n_da <- params$n_da_taxa

  # Ensure slope_effects_distribution is not empty
  if (length(slope_effects_distribution) == 0) {
    stop("slope_effects_distribution is empty - cannot sample from empty vector")
  }

  # Create slope vector
  slope_vec <- rep(0, n_taxa)
  if (n_da > 0) {
    set.seed(seed + 10000)
    da_indices <- sample(1:n_taxa, size = min(n_da, n_taxa), replace = FALSE)
    slope_vec[da_indices] <- sample(slope_effects_distribution, length(da_indices), replace = TRUE)
    slope_vec <- slope_vec - mean(slope_vec)
  }

  # Run simulation
  sim_r <- simulate_compositional_bb( # nolint
    slope_vector = slope_vec,
    mu_inv_softmax = mu_inv_softmax_base,
    log_dispersion_assoc = log_dispersion_assoc,
    n_taxa = n_taxa,
    n_samples = n_samples,
    sd_log_overdispersion = sd_log_overdispersion,
    intercept_dispersion = alpha_intercept_by_taxon,
    library_size_mean = library_size_mean,
    library_size_sd = library_size_sd,
    design_matrix = design_matrix_from_groups(n_samples_per_group),
    seed = seed
  )

  # Add seed and n_da_taxa to simulation for tracking
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
#' @keywords internal
#' @importFrom purrr set_names imap imap_dfr
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
    stop(
      "sim_r$sample_metadata must include factor column `group`. ",
      "Add it in the Quarto document after simulation (see benchmark_*_paper_analyses*.qmd)."
    )
  }
  g <- sim_r$sample_metadata$group
  if (!is.factor(g)) {
    stop("sim_r$sample_metadata$group must be a factor (set in Quarto after simulation).")
  }
  group_levels <- levels(g)
  
  seed <- sim_r$seed
  n_da <- sim_r$n_da_taxa
  
  inputs_std <- purrr::set_names(
    lapply(standard_distance_methods, \(dm) {
      build_standard_betadisper_inputs_paper(
        count_long = sim_r$count_long,
        sample_metadata = sim_r$sample_metadata,
        distance_method = dm,
        group_levels = group_levels
      )
    }),
    nm = paste0("permdisp_", standard_distance_methods)
  )

  # Run analysis for each method using adaptive permutations
  method_results <- purrr::imap(inputs_std, \(inp, method_key) {
    bd <- vegan::betadisper(inp$dist, inp$group)
    
    # Always use adaptive permutations for consistent methodology
    stab <- permutest_betadisper_stable(bd, max_perm = max_perm, batch = batch, target_se = target_se)
    p_value <- unname(stab$p_value)
    p_se <- unname(stab$p_se)
    n_perm <- unname(stab$n_perm)
    
    list(
      method_key = method_key,
      bd = bd,
      group1 = unname(bd$group.distances[1]),
      group2 = unname(bd$group.distances[2]),
      delta = unname(bd$group.distances[2] - bd$group.distances[1]),
      p_value = p_value,
      p_se = p_se,
      n_perm = n_perm,
      distances = unname(bd$distances),
      group = bd$group
    )
  })
  
  # Return format depends on return_full flag
  if (return_full) {
    # Return full results for plotting
    list(
      method_results = method_results,
      seed = seed,
      n_da_taxa = n_da
    )
  } else {
    # Return minimal results for FPR calculations (sweeps)
    purrr::imap_dfr(method_results, \(mr, method_key) {
      tibble::tibble(
        seed = seed,
        n_da_taxa = n_da,
        method = pretty_method_label_roc(method_key),
        p_value = mr$p_value
      )
    })
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
#' @keywords internal
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
#' @keywords internal
#' @importFrom dplyr group_by summarise bind_rows n
#' @importFrom purrr imap_dfr
#' @importFrom tibble tibble
#' @importFrom rlang .data
summarise_fpr <- function(sweep_results) {
  # Handle case where targets auto-bound tibbles into a single tibble
  if (is.data.frame(sweep_results) && all(c("seed", "n_da_taxa", "method", "p_value") %in% names(sweep_results))) {
    fpr_data <- sweep_results
  } else {
    # Extract fpr from list structure
    if (length(sweep_results) == 0) {
      stop("sweep_results is empty")
    }
    
    # Extract FPR data, handling mixed structures
    fpr_list <- lapply(sweep_results, \(x) {
      if (!is.list(x) || !("method_results" %in% names(x))) {
        stop("summarise_fpr expects full analysis objects (return_full = TRUE).")
      }
      purrr::imap_dfr(x$method_results, \(mr, method_key) {
        tibble::tibble(
          seed = x$seed,
          n_da_taxa = x$n_da_taxa,
          threshold = if ("threshold" %in% names(x)) x$threshold else NA_real_,
          method = method_key,
          p_value = mr$p_value
        )
      })
    })
    
    # Remove NULLs and bind
    fpr_list <- fpr_list[!sapply(fpr_list, is.null)]
    
    if (length(fpr_list) == 0) {
      stop("No FPR data found in sweep_results")
    }
    
    fpr_data <- dplyr::bind_rows(fpr_list)
  }
  
  if (nrow(fpr_data) == 0) {
    stop("No FPR data extracted from sweep_results")
  }
  
  fpr_data %>%
    dplyr::group_by(.data$n_da_taxa, .data$method, .data$threshold) %>%
    dplyr::summarise(
      n_sims = dplyr::n(),
      n_fp = sum(.data$p_value < 0.05, na.rm = TRUE),
      FPR = sum(.data$p_value < 0.05, na.rm = TRUE) / dplyr::n(),
      .groups = "drop"
    )
}

