# =============================================================================
# R/simulation.R - Beta-binomial compositional count simulation
# =============================================================================
#
# Dependencies (see @importFrom on each function):
#   VGAM      - rbetabinom.ab
#   dplyr     - left_join, select, any_of
#   purrr     - pmap_dbl
#   stats     - rnorm
#   internal  - softmax (R/transforms.R)
#
# =============================================================================

#' Draw a beta-binomial sample
#'
#' Draws \code{n_sim} counts from a beta-binomial distribution parameterised by
#' mean \code{mu} and dispersion \code{sigma}. The concentration is
#' \eqn{\kappa = 1/\sigma}, with shape parameters
#' \eqn{\alpha = \mu \kappa},
#' \eqn{\beta = (1-\mu) \kappa}, so higher \code{sigma} produces more
#' overdispersion relative to a binomial.
#'
#' Calls \code{VGAM::rbetabinom.ab} for the actual draw.
#'
#' @param n Positive integer scalar (or vector recycled by VGAM); number of
#'   trials i.e. the library size.
#' @param mu Numeric scalar in \eqn{(0, 1)} giving the mean proportion.
#' @param sigma Positive numeric scalar; the dispersion parameter
#'   (\eqn{\kappa = 1/\sigma}).
#' @param n_sim Positive integer; number of independent draws to produce
#'   (default \code{1}).
#'
#' @return Integer vector of length \code{n_sim} containing simulated counts
#'   in \eqn{\{0, 1, \dots, n\}}.
#'
#' @seealso \code{\link{simulate_compositional_bb}} which uses this routine
#'   per (sample, taxon) cell.
#' @export
simulate_beta_binomial <- function(n, mu, sigma, n_sim = 1) {
  kappa <- 1 / sigma
  VGAM::rbetabinom.ab(n = n_sim, size = n, shape1 = mu * kappa, shape2 = (1 - mu) * kappa)
}

#' Simulate compositional beta-binomial microbiome counts
#'
#' Generates microbiome-like count data using a three-level generative model:
#' \enumerate{
#'   \item Compositional means \eqn{\mu_{gj}} come from
#'     \code{softmax(design_matrix \%*\% coeffs)}, so each sample's taxon
#'     proportions sum to 1.
#'   \item Per-taxon log-dispersion follows a mean-dispersion association
#'     \eqn{\log(\sigma_{gj}) = -k\,\mathrm{logit}(\mu_{gj}) + c + \epsilon_j},
#'     with \eqn{\epsilon_j \sim N(0, \tau^2)} fixed across cohorts.
#'   \item Counts are drawn taxon-wise from BetaBinom\eqn{(N_i, \mu_{gj},
#'     \sigma_{gj})} via \code{\link{simulate_beta_binomial}}, with library
#'     sizes \eqn{N_i} sampled from \eqn{N(\text{library\_size\_mean},
#'     \text{library\_size\_sd}^2)} (clipped at 100).
#' }
#'
#' Parameters are computed once per unique cohort (unique row of
#' \code{design_matrix}) and then expanded to the sample level via the
#' cohort -> sample mapping.
#'
#' Compositional constraints (slope and baseline both sum to 0) are validated
#' up-front; the function raises an informative error and suggests a centring
#' fix if they are violated.
#'
#' @param slope_vector Numeric vector of length \code{n_taxa}; per-taxon
#'   group/covariate effect on the linear predictor scale. Must sum to 0.
#' @param mu_inv_softmax Numeric vector of length \code{n_taxa}; baseline
#'   log-linear predictors. Must sum to 0 (compositional constraint).
#' @param log_dispersion_assoc Numeric scalar \eqn{k} in the mean-dispersion
#'   relationship \eqn{\log(\sigma) = -k\,\mathrm{logit}(\mu) + c}. Larger
#'   values yield a stronger inverse mean-dispersion relationship.
#' @param n_taxa Positive integer; number of taxa to simulate.
#' @param n_samples Positive integer; total number of samples (allocated
#'   across cohorts by the design matrix).
#' @param sd_log_overdispersion Non-negative numeric scalar; SD of the
#'   taxon-level \eqn{\log(\sigma)} noise \eqn{\tau}.
#' @param intercept_dispersion Optional length-\code{n_taxa} (or length-1
#'   recycling) numeric vector for the dispersion intercept \eqn{c}. Defaults
#'   to \code{rep(2.0, n_taxa)} when \code{NULL}. Mutually exclusive with
#'   \code{intercept_dispersion_matrix}.
#' @param intercept_dispersion_matrix Optional \code{n_cohorts x n_taxa}
#'   numeric matrix overriding \code{intercept_dispersion}; useful for
#'   differential overdispersion between groups.
#' @param library_size_mean,library_size_sd Numeric scalars; Gaussian
#'   library-size distribution parameters. Library sizes are rounded and
#'   floored at 100.
#' @param design_matrix Numeric matrix of shape \code{n_samples x 2} (the
#'   intercept + binary-group design built by
#'   \code{\link{design_matrix_from_groups}}). Required.
#' @param seed Optional integer; if non-\code{NULL} sets the random seed once
#'   before any sampling occurs.
#'
#' @return Named list with elements:
#' \describe{
#'   \item{\code{ground_truth_params}}{Data frame with one row per
#'     \code{(cohort, taxon)}: baseline intercept, slope, log-linear
#'     predictor, \eqn{\mu}, \eqn{\log\sigma}, \eqn{\sigma}.}
#'   \item{\code{count_long}}{Tall data frame with one row per
#'     \code{(sample, taxon)}: simulated \code{count}, \code{library_size},
#'     \code{mu}, \code{sigma}, \code{log_sigma},
#'     \code{unconstrained_log_mu}, plus all design columns.}
#'   \item{\code{cohort_log_linear_predictors}, \code{cohort_log_sigma},
#'     \code{cohort_sigma}, \code{cohort_mu}}{Cohort-by-taxon matrices used
#'     to build the per-sample data.}
#'   \item{\code{unique_design}}{The unique rows of \code{design_matrix}.}
#'   \item{\code{sample_to_cohort}}{Integer vector mapping each sample to its
#'     cohort row.}
#'   \item{\code{n_cohorts}}{Integer; number of unique cohorts.}
#'   \item{\code{sample_metadata}, \code{taxon_metadata}, \code{design_matrix},
#'     \code{mu_inv_softmax}, \code{slope_vector}}{Inputs preserved for
#'     downstream joining.}
#'   \item{\code{parameters}}{List echoing all numeric inputs for
#'     reproducibility.}
#' }
#'
#' @seealso \code{\link{simulate_beta_binomial}} for the per-cell draw;
#'   \code{\link{softmax}} for the link used in step 1; the alpha-diversity
#'   and betadisper helpers consume the returned \code{count_long}.
#' @export
#' @importFrom dplyr left_join select mutate rename any_of
#' @importFrom purrr pmap_dbl
#' @importFrom stats rnorm
simulate_compositional_bb <- function(
  slope_vector,
  mu_inv_softmax,
  log_dispersion_assoc,
  n_taxa,
  n_samples,
  sd_log_overdispersion,
  intercept_dispersion = NULL,
  intercept_dispersion_matrix = NULL,
  library_size_mean = 10000,
  library_size_sd = 2000,
  design_matrix,
  seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  # Compositional / shape constraints (must remain: callers rely on them).
  if (length(slope_vector)   != n_taxa) stop("slope_vector length must equal n_taxa")
  if (length(mu_inv_softmax) != n_taxa) stop("mu_inv_softmax length must equal n_taxa")
  if (abs(sum(mu_inv_softmax)) > 1e-10) stop("mu_inv_softmax must sum to 0 (compositional constraint).")
  if (abs(sum(slope_vector))   > 1e-10) stop("slope_vector must sum to 0 (compositional constraint).")
  if (!is.null(intercept_dispersion_matrix) && !is.null(intercept_dispersion)) {
    stop("Provide only one of intercept_dispersion or intercept_dispersion_matrix.")
  }

  # ---- cohort linear predictors -------------------------------------------
  # Design is always [Intercept, Group] (see design_matrix_from_groups()), so
  # coefficient matrix has exactly two rows: baseline logits and slope.
  coeff_matrix   <- rbind(mu_inv_softmax, slope_vector)
  unique_design  <- unique(design_matrix)
  n_cohorts      <- nrow(unique_design)
  cohort_log_linear_predictors <- unique_design %*% coeff_matrix
  colnames(cohort_log_linear_predictors) <- paste0("Taxon_", seq_len(n_taxa))

  # ---- dispersion intercept (vector or matrix path) -----------------------
  # `intercept_dispersion_matrix` lets the qmd "DA off" null sims set per-group
  # dispersion intercepts (n_cohorts x n_taxa); the vector path is the usual
  # paired-resampling case driven by build_simulation_params_brito().
  if (!is.null(intercept_dispersion_matrix)) {
    intercept_dispersion_matrix <- as.matrix(intercept_dispersion_matrix)
    if (!identical(dim(intercept_dispersion_matrix), c(n_cohorts, n_taxa))) {
      stop(sprintf("intercept_dispersion_matrix must be n_cohorts x n_taxa (%d x %d); got %d x %d.",
                   n_cohorts, n_taxa,
                   nrow(intercept_dispersion_matrix), ncol(intercept_dispersion_matrix)))
    }
    intercept_dispersion_param <- NULL
  } else {
    if (is.null(intercept_dispersion)) intercept_dispersion <- 2.0
    if (!length(intercept_dispersion) %in% c(1L, n_taxa)) {
      stop("intercept_dispersion must have length n_taxa (", n_taxa,
           "), or length 1 to recycle; got ", length(intercept_dispersion), ".")
    }
    intercept_dispersion        <- rep_len(intercept_dispersion, n_taxa)
    intercept_dispersion_matrix <- matrix(intercept_dispersion, n_cohorts, n_taxa, byrow = TRUE)
    intercept_dispersion_param  <- intercept_dispersion
  }

  # Per-taxon log(sigma) noise (shared across cohorts, so it shifts the
  # mean-dispersion intercept but not its cohort contrast).
  cohort_log_sigma_errors <- matrix(
    rnorm(n_taxa, sd = sd_log_overdispersion),
    nrow = n_cohorts, ncol = n_taxa, byrow = TRUE
  )

  # log(sigma_gj) = -k * log_mu_gj + c_j + epsilon_j
  cohort_log_sigma <- -log_dispersion_assoc * cohort_log_linear_predictors +
                       intercept_dispersion_matrix + cohort_log_sigma_errors
  cohort_sigma     <- exp(cohort_log_sigma)
  cohort_mu        <- softmax(cohort_log_linear_predictors)

  # ---- sample <-> cohort map ----------------------------------------------
  sample_to_cohort <- match(
    apply(design_matrix, 1, paste, collapse = "\r"),
    apply(unique_design, 1, paste, collapse = "\r")
  )
  library_sizes <- pmax(round(rnorm(n_samples, library_size_mean, library_size_sd)), 100)

  taxon_ids <- paste0("Taxon_", seq_len(n_taxa))
  sample_metadata <- data.frame(
    sample_id    = paste0("Sample_", seq_len(n_samples)),
    library_size = library_sizes,
    cohort_idx   = sample_to_cohort,
    design_matrix,
    check.names  = FALSE
  )

  # ground_truth_params: one row per (cohort, taxon), built by row-replicating
  # the n_cohorts x n_taxa parameter matrices in cohort-major order.
  ground_truth_params <- data.frame(
    taxon_id                     = rep(taxon_ids,           times = n_cohorts),
    cohort_idx                   = rep(seq_len(n_cohorts),  each  = n_taxa),
    baseline_intercept           = rep(mu_inv_softmax,      times = n_cohorts),
    slope                        = rep(slope_vector,        times = n_cohorts),
    cohort_log_linear_predictors = as.vector(t(cohort_log_linear_predictors)),
    cohort_mu                    = as.vector(t(cohort_mu)),
    cohort_log_sigma             = as.vector(t(cohort_log_sigma)),
    cohort_sigma                 = as.vector(t(cohort_sigma)),
    stringsAsFactors             = FALSE
  )

  # count_long: cross-join of (sample, taxon) with matrix-indexed cohort
  # lookups, then per-cell BB draw.
  count_long <- data.frame(
    sample_id    = rep(sample_metadata$sample_id, times = n_taxa),
    taxon_id     = rep(taxon_ids,                 each  = n_samples),
    library_size = rep(library_sizes,             times = n_taxa),
    cohort_idx   = rep(sample_to_cohort,          times = n_taxa),
    stringsAsFactors = FALSE
  )
  ti  <- match(count_long$taxon_id, taxon_ids)
  ci  <- count_long$cohort_idx
  count_long$mu                   <- cohort_mu[cbind(ci, ti)]
  count_long$sigma                <- cohort_sigma[cbind(ci, ti)]
  count_long$unconstrained_log_mu <- cohort_log_linear_predictors[cbind(ci, ti)]
  count_long$log_sigma            <- cohort_log_sigma[cbind(ci, ti)]
  count_long$count <- purrr::pmap_dbl(
    list(count_long$library_size, count_long$mu, count_long$sigma),
    function(n, mu, sigma) simulate_beta_binomial(n = n, mu = mu, sigma = sigma, n_sim = 1L)
  )
  count_long <- count_long %>%
    dplyr::left_join(sample_metadata,
                     by = c("sample_id", "library_size", "cohort_idx")) %>%
    dplyr::select(-dplyr::any_of("cohort_idx"))

  list(
    ground_truth_params          = ground_truth_params,
    count_long                   = count_long,
    cohort_log_linear_predictors = cohort_log_linear_predictors,
    cohort_log_sigma             = cohort_log_sigma,
    cohort_sigma                 = cohort_sigma,
    cohort_mu                    = cohort_mu,
    unique_design                = unique_design,
    sample_to_cohort             = sample_to_cohort,
    n_cohorts                    = n_cohorts,
    sample_metadata              = sample_metadata,
    design_matrix                = design_matrix,
    mu_inv_softmax               = mu_inv_softmax,
    slope_vector                 = slope_vector,
    parameters = list(
      slope_vector                = slope_vector,
      mu_inv_softmax              = mu_inv_softmax,
      baseline_intercepts         = mu_inv_softmax,
      log_dispersion_assoc        = log_dispersion_assoc,
      n_taxa                      = n_taxa,
      n_samples                   = n_samples,
      n_cohorts                   = n_cohorts,
      sd_log_overdispersion       = sd_log_overdispersion,
      intercept_dispersion        = intercept_dispersion_param,
      intercept_dispersion_matrix = intercept_dispersion_matrix,
      library_size_mean           = library_size_mean,
      library_size_sd             = library_size_sd
    )
  )
}

