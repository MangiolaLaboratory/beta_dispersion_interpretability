# =============================================================================
# R/permanova.R - PERMANOVA job runner and TPR summariser
# =============================================================================
#
# Source this via the top-level functions.R, which sources all R/*.R modules
# in the correct dependency order.
#
# Dependencies (see @importFrom on each function):
#   vegan   - adonis2
#   purrr   - set_names, imap_dfr
#   tibble  - tibble
#   dplyr   - filter, group_by, summarise, if_else
#   rlang   - .data
#   internal - filter_by_abundance_cases (R/pipeline.R),
#              build_standard_betadisper_inputs_paper (R/betadisper.R)
#
# =============================================================================

#' Run a per-simulation PERMANOVA across several standard distances
#'
#' For one simulated dataset, this function:
#' \enumerate{
#'   \item Optionally drops low-prevalence/abundance taxa via
#'     \code{\link{filter_by_abundance_cases}}.
#'   \item Builds a \code{betadisper}-style \code{(dist, group)} pair for
#'     each requested distance via
#'     \code{\link{build_standard_betadisper_inputs_paper}}.
#'   \item Runs \code{vegan::adonis2(dist ~ group, ...)} on each and returns
#'     the headline \code{Pr(>F)} p-value alongside the (seed, n_da_taxa,
#'     threshold, method) identifiers.
#' }
#'
#' Edge cases handled:
#' \itemize{
#'   \item If the filtered table has fewer than 3 distinct taxa, an empty
#'     tibble with the correct columns is returned (so downstream
#'     \code{bind_rows} still works).
#'   \item \code{adonis2} errors are caught and surface as \code{NA} p-values
#'     for that method (rather than failing the whole job).
#' }
#'
#' @param sim_r Simulation result list (as returned by
#'   \code{\link{simulate_compositional_bb}}, optionally annotated by
#'   \code{\link{synthetic_sim_add_two_group_factor}} so
#'   \code{sample_metadata$group} is a factor).
#' @param threshold Abundance / prevalence threshold passed to
#'   \code{\link{filter_by_abundance_cases}}. Use \code{NULL} or \code{0} for
#'   no filtering.
#' @param seed Integer; the seed used to generate \code{sim_r}, recorded for
#'   provenance in the returned tibble.
#' @param n_da_taxa Integer; the number of differentially-abundant taxa
#'   injected when building \code{sim_r}; also for provenance.
#' @param permutations Positive integer; permutation count for
#'   \code{vegan::adonis2}. Default \code{199}.
#' @param standard_distance_methods Character vector of distance method names
#'   recognised by \code{\link{build_standard_betadisper_inputs_paper}}.
#'   Default covers Bray-Curtis, Aitchison family, Jaccard and Hellinger.
#'
#' @return Tibble with one row per distance method:
#'   \code{seed}, \code{n_da_taxa}, \code{threshold}, \code{method}
#'   (e.g. \code{"permanova_bray"}), \code{p_value}.
#'
#' @seealso \code{\link{summarise_permanova_tpr}}.
#' @keywords internal
#' @importFrom purrr set_names imap_dfr
#' @importFrom tibble tibble
#' @importFrom vegan adonis2
run_permanova_job <- function(
  sim_r,
  threshold,
  seed,
  n_da_taxa,
  permutations = 199,
  standard_distance_methods = c("bray", "aitchison", "robust.aitchison", "jaccard", "hellinger")
) {
  sim_filtered <- sim_r
  sim_filtered$count_long <- filter_by_abundance_cases(sim_r$count_long, threshold = threshold)

  if (length(unique(sim_filtered$count_long$taxon_id)) < 3) {
    return(tibble::tibble(
      seed = integer(0),
      n_da_taxa = integer(0),
      threshold = numeric(0),
      method = character(0),
      p_value = numeric(0)
    ))
  }

  g <- sim_filtered$sample_metadata$group
  if (!is.factor(g)) {
    stop("sample_metadata$group must be a factor (set in Quarto after simulation).")
  }
  group_levels <- levels(g)

  inputs_std <- purrr::set_names(
    lapply(standard_distance_methods, \(dm) {
      build_standard_betadisper_inputs_paper(
        count_long = sim_filtered$count_long,
        sample_metadata = sim_filtered$sample_metadata,
        distance_method = dm,
        group_levels = group_levels
      )
    }),
    nm = paste0("permanova_", standard_distance_methods)
  )

  purrr::imap_dfr(inputs_std, \(inp, method_key) {
    p_val <- tryCatch({
      ad <- vegan::adonis2(inp$dist ~ group, data = data.frame(group = inp$group), permutations = permutations)
      unname(ad$`Pr(>F)`[1])
    }, error = function(e) {
      NA_real_
    })

    tibble::tibble(
      seed = seed,
      n_da_taxa = n_da_taxa,
      threshold = threshold,
      method = method_key,
      p_value = p_val
    )
  })
}

#' Summarise PERMANOVA results into a true-positive rate (TPR)
#'
#' Aggregates a long table of per-simulation PERMANOVA results (as produced
#' by \code{\link{run_permanova_job}}) into per-\code{(n_da_taxa, method,
#' threshold)} groups, computing:
#' \itemize{
#'   \item \code{n_sims}: number of simulations contributing a non-\code{NA}
#'     p-value,
#'   \item \code{n_tp}: number of those with \code{p < 0.05},
#'   \item \code{TPR}: \code{n_tp / n_sims} (NA if \code{n_sims == 0}).
#' }
#'
#' Rows with \code{n_da_taxa == 0} (the null setting) are excluded -- those
#' belong in the FPR summary, not the TPR one.
#'
#' @param permanova_results Tibble with columns \code{n_da_taxa},
#'   \code{method}, \code{threshold}, \code{p_value} (typically the result of
#'   \code{bind_rows()}-ing many \code{\link{run_permanova_job}} outputs).
#'
#' @return Tibble grouped by \code{(n_da_taxa, method, threshold)} with
#'   \code{n_sims}, \code{n_tp}, \code{TPR}.
#'
#' @seealso \code{\link{run_permanova_job}}.
#' @keywords internal
#' @importFrom dplyr filter group_by summarise if_else
#' @importFrom rlang .data
summarise_permanova_tpr <- function(permanova_results) {
  permanova_results %>%
    dplyr::filter(.data$n_da_taxa > 0) %>%
    dplyr::group_by(.data$n_da_taxa, .data$method, .data$threshold) %>%
    dplyr::summarise(
      n_sims = sum(!is.na(.data$p_value)),
      n_tp = sum(.data$p_value < 0.05, na.rm = TRUE),
      TPR = dplyr::if_else(.data$n_sims > 0, .data$n_tp / .data$n_sims, NA_real_),
      .groups = "drop"
    )
}
