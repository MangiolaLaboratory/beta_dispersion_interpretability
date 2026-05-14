# =============================================================================
# R/sccomp.R - sccomp parameter extraction and simulation-param builders
# =============================================================================
#
# Dependencies (see @importFrom on each function):
#   stats  - sd
#   Other inputs come from sccomp::sccomp_estimate() results (the `fit` object
#   uses CmdStanR-style `$metadata()` and `$summary()` methods).
#
# =============================================================================

#' Extract precision-trend posterior summaries from an sccomp Stan fit
#'
#' Pulls the mean and SD of the precision-trend intercept/slope variables out
#' of a CmdStanR fit produced by \code{sccomp}. Tries three name conventions
#' in order, returning as soon as one set is found:
#' \enumerate{
#'   \item Indexed names \code{prec_intercept_1} / \code{prec_slope_1}
#'     (modern sccomp, "first mixture component").
#'   \item Indexed names \code{prec_intercept_2} / \code{prec_slope_2}
#'     (bimodal models).
#'   \item Older scalar \code{prec_intercept} + any name starting with
#'     \code{prec_slope}.
#' }
#'
#' If none of these are present the function falls through to the legacy
#' branch and will raise an error from \code{[[1]]} indexing of an empty
#' vector -- this is the historical behaviour and is preserved here.
#'
#' @param fit A CmdStanR-style fit object exposing \code{$metadata()} and
#'   \code{$summary()} (e.g. as attached to an \code{sccomp} result).
#'
#' @return Named list with:
#' \describe{
#'   \item{\code{prec_intercept}}{Posterior mean of the precision intercept.}
#'   \item{\code{prec_slope}}{Posterior mean of the precision slope.}
#'   \item{\code{k_sd_real}}{Posterior SD of the slope (used to calibrate
#'     downstream draws).}
#'   \item{\code{slope_parameter}}{Character; the actual Stan variable name
#'     resolved to (useful for traceability).}
#' }
#'
#' @seealso \code{\link{extract_sccomp_params_brito}} which can optionally
#'   consume the values returned here.
#' @export
extract_precision_trend_params <- function(fit) {
  fit_vars <- fit$metadata()$stan_variables

  # Try modern indexed names first (j = 1, then 2), fall back to legacy
  # unindexed `prec_intercept` + first `prec_slope*` match.
  for (j in c(1L, 2L)) {
    int_name <- paste0("prec_intercept_", j)
    slo_name <- paste0("prec_slope_",     j)
    if (int_name %in% fit_vars && slo_name %in% fit_vars) {
      int_s <- fit$summary(int_name)
      slo_s <- fit$summary(slo_name)
      return(list(
        prec_intercept  = int_s$mean[1],
        prec_slope      = slo_s$mean[1],
        k_sd_real       = slo_s$sd[1],
        slope_parameter = slo_name
      ))
    }
  }
  slope_name <- sort(grep("^prec_slope", fit_vars, value = TRUE))[[1]]
  slo_s <- fit$summary(slope_name)
  list(
    prec_intercept  = fit$summary("prec_intercept")$mean[1],
    prec_slope      = slo_s$mean[1],
    k_sd_real       = slo_s$sd[1],
    slope_parameter = slope_name
  )
}

#' Summarise per-sample library sizes from an sccomp model input
#'
#' Reads the library-size vector out of the \code{exposure} field of an
#' \code{sccomp} \code{model_input} object and returns a mean / SD pair. Used
#' to seed library-size sampling in downstream simulators so the synthetic
#' counts have a realistic depth distribution.
#'
#' Robustness rules:
#' \itemize{
#'   \item Missing/NULL inputs or empty/zero exposures -> \code{NA_real_} pair
#'     (callers fall back to defaults via
#'     \code{\link{.sim_library_size_or_default}}).
#'   \item If the empirical SD is zero or non-finite, it is replaced by
#'     \code{max(0.05 * mean, 1)} to avoid degenerate distributions.
#'   \item Otherwise SD is floored at 5\% of the mean to avoid pathologically
#'     tight library-size distributions.
#' }
#'
#' @param model_input A list with at least an \code{exposure} numeric vector
#'   (typically \code{attr(sccomp_result, "model_input")}).
#'
#' @return Named list with elements \code{library_size_mean} and
#'   \code{library_size_sd} (both numeric scalars, possibly \code{NA_real_}).
#'
#' @seealso \code{\link{.sim_library_size_or_default}},
#'   \code{\link{extract_sccomp_params_brito}}.
#' @export
#' @importFrom stats sd
library_size_from_sccomp_model_input <- function(model_input) {
  if (is.null(model_input) || is.null(model_input$exposure)) {
    return(list(library_size_mean = NA_real_, library_size_sd = NA_real_))
  }
  ex <- as.numeric(model_input$exposure)
  ex <- ex[is.finite(ex) & ex > 0]
  if (length(ex) == 0L) {
    return(list(library_size_mean = NA_real_, library_size_sd = NA_real_))
  }
  m <- mean(ex)
  s <- stats::sd(ex)
  # Floor SD at 5% of the mean to avoid pathologically tight or degenerate
  # library-size distributions when extracting from real cohorts.
  s_floor <- max(0.05 * m, 1)
  list(library_size_mean = m,
       library_size_sd   = if (is.finite(s) && s > 0) max(s, s_floor) else s_floor)
}

#' Pick simulation library-size parameters with a safe default
#'
#' Internal helper that returns the library-size mean/SD recorded inside an
#' \code{sccomp_params} list, falling back to a fixed default
#' (\code{mean = 15125}, \code{sd = 5000}) when the mean is missing or
#' non-positive. If only the SD is bad, it is replaced by
#' \code{max(0.05 * mean, 1)}.
#'
#' @param sccomp_params A list with (optional) numeric scalars
#'   \code{library_size_mean} and \code{library_size_sd}, typically produced
#'   by \code{\link{extract_sccomp_params_brito}}.
#'
#' @return Named list with elements \code{library_size_mean} and
#'   \code{library_size_sd}, both finite positive numeric scalars.
#' @noRd
.sim_library_size_or_default <- function(sccomp_params) {
  m <- sccomp_params$library_size_mean
  s <- sccomp_params$library_size_sd
  if (is.null(m) || !is.finite(m) || m <= 0) {
    return(list(library_size_mean = 15125, library_size_sd = 5000))
  }
  if (is.null(s) || !is.finite(s) || s <= 0) {
    s <- max(0.05 * m, 1)
  }
  list(library_size_mean = m, library_size_sd = s)
}

#' Extract simulation-ready parameters from an sccomp result (Brito variant)
#'
#' Pulls the composition and variability effects out of an \code{sccomp}
#' result tibble (one row per (\code{cell_group}, \code{parameter})) so they
#' can seed \code{\link{simulate_compositional_bb}} via
#' \code{\link{build_simulation_params_brito}}.
#'
#' This is the "Brito-style" extractor used by every per-dataset benchmark
#' analysis. Key differences from any future non-Brito extractor:
#' \itemize{
#'   \item It supports models whose non-intercept coefficient name is not
#'     fixed (e.g. \code{body_site}-style terms) via \code{target_parameter}.
#'   \item It does NOT consume the precision-trend slope/intercept from the
#'     Stan fit. \code{k_real} is always \code{0}, and
#'     \code{c_real}/\code{prec_sd_real} are \code{NA}. Dispersion intercepts
#'     come from the per-taxon \code{v_effect} column via
#'     \code{alpha_intercept_effects <- v_intercept} in
#'     \code{\link{build_simulation_params_brito}}.
#' }
#'
#' Variability semantics: sccomp parameterises precision as
#' \eqn{\phi = \exp(-v_{\text{effect}})}, while
#' \code{simulate_compositional_bb} uses dispersion
#' \eqn{\sigma} with concentration \eqn{\kappa = 1/\sigma}. To preserve
#' "higher precision -> lower dispersion" we map
#' \eqn{\log\sigma = v_{\text{effect}}}.
#'
#' @param result Tibble returned by \code{sccomp::sccomp_estimate()} (or a
#'   downstream test); must have columns \code{parameter}, \code{c_effect},
#'   \code{c_lower}, \code{c_upper}, \code{v_effect}, and (for
#'   \code{n_taxa_real} fallback) \code{cell_group}. Optionally carries a
#'   \code{model_input} attribute used to read \code{N}, \code{M} and
#'   \code{exposure}.
#' @param use_precision_trend_from_fit Kept only for backward compatibility;
#'   currently ignored. The Brito path always returns \code{k_real = 0}.
#' @param target_parameter Optional character scalar naming which
#'   non-intercept parameter to treat as the group effect. If \code{NULL}
#'   (default), the first non-intercept parameter is used. Must be one of
#'   the parameters present in \code{result}.
#'
#' @return A list with these scalar fields:
#'   \code{n_samples_real}, \code{n_taxa_real}, \code{k_real} (always 0),
#'   \code{k_sd_real} (NA), \code{c_real} (NA), \code{prec_sd_real} (NA),
#'   \code{group_parameter} (resolved name), \code{library_size_mean},
#'   \code{library_size_sd}; and these per-taxon vectors:
#'   \code{intercept_effects}, \code{slope_effects},
#'   \code{slope_effects_significant} (only taxa with significant
#'   \code{c_effect}), \code{v_intercept}, \code{v_slope},
#'   \code{alpha_intercept_effects} (= \code{v_intercept}; see semantics above).
#'
#' @seealso \code{\link{build_simulation_params_brito}} which consumes this
#'   output; \code{\link{library_size_from_sccomp_model_input}}.
#' @export
extract_sccomp_params_brito <- function(
    result,
    use_precision_trend_from_fit = FALSE,  # kept for back-compat; ignored
    target_parameter = NULL) {
  model_input <- attr(result, "model_input")
  if (!is.null(model_input) && is.list(model_input)) {
    n_samples_real <- model_input$N
    n_taxa_real    <- model_input$M
  } else {
    n_samples_real <- 178
    n_taxa_real    <- length(unique(result$cell_group))
  }

  non_intercept <- unique(result$parameter[result$parameter != "Intercept"])
  if (length(non_intercept) == 0L) {
    stop("No non-intercept parameter found in sccomp result.")
  }
  group_parameter <- if (is.null(target_parameter)) non_intercept[1] else target_parameter
  if (!group_parameter %in% result$parameter) {
    stop("group_parameter '", group_parameter, "' not present in sccomp result.")
  }

  int_rows <- result[result$parameter == "Intercept",       , drop = FALSE]
  slp_rows <- result[result$parameter == group_parameter,   , drop = FALSE]
  sig_mask <- slp_rows$c_lower > 0 | slp_rows$c_upper < 0

  lib_ls <- library_size_from_sccomp_model_input(model_input)

  # sccomp variability: precision phi = exp(-v_effect). Simulator's dispersion
  # sigma has concentration kappa = 1/sigma; mapping log(sigma) = v_effect
  # preserves "higher precision -> lower dispersion".
  list(
    n_samples_real            = n_samples_real,
    n_taxa_real               = n_taxa_real,
    k_real                    = 0,
    k_sd_real                 = NA_real_,
    c_real                    = NA_real_,
    prec_sd_real              = NA_real_,
    intercept_effects         = int_rows$c_effect,
    slope_effects             = slp_rows$c_effect,
    slope_effects_significant = slp_rows$c_effect[sig_mask],
    v_intercept               = int_rows$v_effect,
    v_slope                   = slp_rows$v_effect,
    alpha_intercept_effects   = int_rows$v_effect,
    group_parameter           = group_parameter,
    library_size_mean         = lib_ls$library_size_mean,
    library_size_sd           = lib_ls$library_size_sd
  )
}

#' Build simulation parameters from sccomp-derived inputs (Brito variant)
#'
#' Takes the output of \code{\link{extract_sccomp_params_brito}} and turns it
#' into a self-contained list of inputs ready to feed
#' \code{\link{simulate_compositional_bb}} (intercepts, slopes, library-size
#' distribution, dispersion intercept(s), etc.) plus job-replication settings
#' (\code{n_reps_auc}, \code{perm_reps_auc}, \code{rep_seeds_auc}) consumed
#' by the per-dataset pipelines.
#'
#' Highlights of the construction:
#' \itemize{
#'   \item Resampling is \strong{paired}: for each synthetic taxon a row of
#'     \code{(intercept_effects, alpha_intercept_effects)} is drawn together,
#'     so the empirical mean-variability dependence between composition and
#'     variability intercepts is preserved.
#'   \item The mean-centred sampled intercepts feed
#'     \code{mu_inv_softmax_base_realistic}.
#'   \item The paired \code{alpha_intercept_effects} drive the per-taxon log
#'     sigma intercepts; their mean / SD become
#'     \code{intercept_disp_realistic} and
#'     \code{sd_log_overdispersion_realistic} respectively.
#'   \item Slopes are independently resampled (with mean centring) so the
#'     two-group effect remains compositional.
#'   \item Two \code{set.seed(123)} calls fix the resampling so two reports
#'     that share the same \code{sccomp_params} simulate the same parameter
#'     set; per-rep variability comes from the \code{rep_seeds_auc} block.
#' }
#'
#' @param sccomp_params List returned by
#'   \code{\link{extract_sccomp_params_brito}}.
#' @param n_reps_auc Positive integer; number of simulation replicates that
#'   the downstream AUC sweep should run. Default \code{20}.
#' @param group_levels Length-2 character vector; the per-group labels used
#'   for downstream factor conversion via
#'   \code{\link{synthetic_sim_add_two_group_factor}}.
#' @param n_samples_per_group Optional. \code{NULL} (default) sets a balanced
#'   design with \code{floor(n_samples_real / 2)} per group. Length-1 gives a
#'   balanced design at the supplied size; length-2 gives an unbalanced design.
#'
#' @return Named list of simulation inputs and replication settings:
#'   \code{n_taxa}, \code{n_samples_per_group}, \code{n_groups},
#'   \code{n_samples}, \code{group_levels}, \code{library_size_mean},
#'   \code{library_size_sd}, \code{mu_inv_softmax_base_realistic},
#'   \code{intercept_disp_realistic}, \code{intercept_disp_realistic_taxon},
#'   \code{sigma_realistic}, \code{slope_realistic},
#'   \code{sd_log_overdispersion_realistic}, \code{alpha_intercept_sampled},
#'   \code{n_reps_auc}, \code{perm_reps_auc}, \code{rep_seeds_auc}.
#'
#' @seealso \code{\link{extract_sccomp_params_brito}},
#'   \code{\link{simulate_compositional_bb}},
#'   \code{\link{synthetic_sim_add_two_group_factor}}.
#' @export
#' @importFrom stats sd
build_simulation_params_brito <- function(
  sccomp_params,
  n_reps_auc = 20,
  group_levels,
  n_samples_per_group = NULL
) {
  if (length(group_levels) != 2L) {
    stop("group_levels must have length 2 (Group==0 then Group==1 in design_matrix_from_groups).")
  }
  if (is.null(n_samples_per_group)) {
    n_samples_per_group <- floor(sccomp_params$n_samples_real / 2)
  }
  if (!length(n_samples_per_group) %in% c(1L, 2L)) {
    stop("n_samples_per_group must be length 1 (balanced) or length 2 (per-group counts).")
  }
  if (length(sccomp_params$intercept_effects) != length(sccomp_params$alpha_intercept_effects)) {
    stop("intercept_effects and alpha_intercept_effects must have the same length ",
         "for paired resampling.")
  }

  n_taxa    <- sccomp_params$n_taxa_real
  n_samples <- if (length(n_samples_per_group) == 1L) 2L * n_samples_per_group
               else sum(n_samples_per_group)
  lib       <- .sim_library_size_or_default(sccomp_params)

  # Paired resampling: keep empirical (composition, variability) dependence.
  # Two `set.seed(123)` calls (one per sample()) are deliberate: they make the
  # paired-intercept draw and the slope draw each start from the same fixed
  # state, so the two are independently reproducible.
  set.seed(123)
  paired_idx         <- sample.int(length(sccomp_params$intercept_effects),
                                   size = n_taxa, replace = TRUE)
  sampled_intercepts <- sccomp_params$intercept_effects[paired_idx]
  sampled_alpha      <- sccomp_params$alpha_intercept_effects[paired_idx]
  set.seed(123)
  slope_realistic    <- sample(sccomp_params$slope_effects, n_taxa, replace = TRUE)

  list(
    n_taxa                          = n_taxa,
    n_samples_per_group             = n_samples_per_group,
    n_groups                        = 2L,
    n_samples                       = n_samples,
    group_levels                    = group_levels,
    library_size_mean               = lib$library_size_mean,
    library_size_sd                 = lib$library_size_sd,
    mu_inv_softmax_base_realistic   = sampled_intercepts - mean(sampled_intercepts),
    intercept_disp_realistic        = mean(sampled_alpha),
    intercept_disp_realistic_taxon  = sampled_alpha,
    sigma_realistic                 = exp(mean(sampled_alpha)),
    slope_realistic                 = slope_realistic - mean(slope_realistic),
    sd_log_overdispersion_realistic = stats::sd(sampled_alpha),
    alpha_intercept_sampled         = sampled_alpha,
    n_reps_auc                      = n_reps_auc,
    perm_reps_auc                   = 199,
    rep_seeds_auc                   = 12000 + seq_len(n_reps_auc)
  )
}
