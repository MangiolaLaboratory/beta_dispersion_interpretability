# =============================================================================
# R/simulation.R - Beta-binomial compositional count simulation
# =============================================================================
#
# Source this via the top-level functions.R, which sources all R/*.R modules
# in the correct dependency order.
#
# Dependencies (see @importFrom on each function):
#   VGAM      - rbetabinom.ab
#   rlang     - check_installed
#   dplyr     - left_join, select, mutate, rename, any_of
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
#' Calls \code{VGAM::rbetabinom.ab} for the actual draw; a friendly
#' \code{rlang::check_installed} guard is used to surface a clear error if
#' VGAM is missing.
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
#' @keywords internal
#' @importFrom rlang check_installed
simulate_beta_binomial <- function(n, mu, sigma, n_sim = 1) {
  # Concentration kappa = 1/sigma; alpha = mu*kappa, beta = (1-mu)*kappa
  # Higher dispersion → Lower concentration → More overdispersion
  concentration <- 1 / sigma
  
  # Convert to alpha, beta parameters
  alpha <- mu * concentration
  beta <- (1 - mu) * concentration
  
  # Check if VGAM is available (required for both simulation and fitting)
  rlang::check_installed("VGAM", reason = "for beta-binomial simulation and fitting")
  
  # Use VGAM's direct beta-binomial simulator
  # This is more efficient and accurate than the two-step rbeta + rbinom approach
  counts <- VGAM::rbetabinom.ab(n = n_sim, size = n, shape1 = alpha, shape2 = beta)
  
  return(counts)
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
#' @param design_matrix Optional numeric matrix of shape
#'   \code{n_samples x n_covariates}. If \code{NULL} (default), an intercept
#'   plus a single standard-normal covariate is used.
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
#' @keywords internal
#' @importFrom dplyr left_join select mutate rename any_of
#' @importFrom purrr pmap_dbl
#' @importFrom stats rnorm
simulate_compositional_bb <- function(
  slope_vector,           # Vector of slopes for each taxon (length = n_taxa)
                          # Determines how each taxon responds to covariates
  mu_inv_softmax,         # Vector of baseline log-linear predictors (length = n_taxa)
                          # These are in log space and sum to 0 (compositional constraint)
                          # Will be transformed to probabilities via softmax
  log_dispersion_assoc,   # Association parameter: log(σ) = -k·logit(μ) + c
                          # This is the 'k' parameter (slope of inverse relationship)
                          # Higher k → stronger inverse relationship
  n_taxa,                 # Number of taxa/species
  n_samples,               # Number of samples
  sd_log_overdispersion,   # SD of log(dispersion) around regression line
                          # Controls variability: log(σ) = -k·logit(μ) + c + N(0, sd^2)
  intercept_dispersion = NULL,  # Length n_taxa (or 1 to recycle); see roxygen
  intercept_dispersion_matrix = NULL,  # Optional n_cohorts x n_taxa; if set, overrides vector
  library_size_mean = 10000,   # Mean library size per sample
  library_size_sd = 2000,     # SD of library size
  design_matrix = NULL,        # Optional design matrix (n_samples x n_covariates)
                                # If NULL, creates intercept + one random covariate
  seed = NULL                   # Optional random seed
) {
  
  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Validate inputs
  if (length(slope_vector) != n_taxa) {
    stop("slope_vector length must equal n_taxa")
  }
  if (length(mu_inv_softmax) != n_taxa) {
    stop("mu_inv_softmax length must equal n_taxa")
  }
  
  # Check compositional constraints (sum to 0)
  tolerance <- 1e-10
  mu_sum <- sum(mu_inv_softmax)
  if (abs(mu_sum) > tolerance) {
    stop(sprintf(
      "mu_inv_softmax must sum to 0 (compositional constraint). Current sum = %.10f\n  Use: mu_inv_softmax <- mu_inv_softmax - mean(mu_inv_softmax)",
      mu_sum
    ))
  }
  
  slope_sum <- sum(slope_vector)
  if (abs(slope_sum) > tolerance) {
    stop(sprintf(
      "slope_vector must sum to 0 (compositional constraint). Current sum = %.10f\n  Use: slope_vector <- slope_vector - mean(slope_vector)",
      slope_sum
    ))
  }
  
  # ========================================================================
  # Step 1: Generate compositional means using log-linear predictors + softmax
  # ========================================================================
  # Note on terminology:
  # - mu_inv_softmax = log-linear predictors (in log space, sum to 0)
  # - These represent the "inverse" of softmax (i.e., what you'd get from centered log-ratio transform)
  # - We apply softmax() to convert them to probabilities that sum to 1
  
  # Create design matrix if not provided
  if (is.null(design_matrix)) {
    # Default: intercept + one covariate
    design_matrix <- cbind(1, rnorm(n_samples))
    colnames(design_matrix) <- c("Intercept", "Covariate1")
  }
  
  n_covariates <- ncol(design_matrix)
  
  # Use mu_inv_softmax as baseline intercepts for each taxon
  baseline_intercepts <- mu_inv_softmax
  
  # Build coefficient matrix for linear model
  # Coefficient matrix: (n_covariates x n_taxa)
  # Row 1: baseline_intercepts (for intercept term)
  # Row 2: slope_vector (for first covariate/group effect)
  # Additional rows: random effects for additional covariates (if any)
  
  # Start with baseline intercepts and slopes
  coeff_matrix <- rbind(
    baseline_intercepts,
    if (n_covariates > 1) slope_vector else NULL
  )
  
  # Add random effects for additional covariates (if any)
  if (n_covariates > 2) {
    additional_effects <- matrix(
      rnorm((n_covariates - 2) * n_taxa, mean = 0, sd = 0.5),
      nrow = n_covariates - 2,
      ncol = n_taxa
    )
    coeff_matrix <- rbind(coeff_matrix, additional_effects)
  }
  
  # ========================================================================
  # NEW APPROACH: Generate parameters at COHORT level, then sample
  # ========================================================================
  
  # Step 1.1: Get unique cohorts from design matrix
  unique_design <- unique(design_matrix)
  n_cohorts <- nrow(unique_design)
  
  cat("Generating parameters for", n_cohorts, "unique cohorts\n")
  
  # Step 1.2: Calculate log-linear predictors at COHORT level
  # Matrix multiplication: (n_cohorts x n_covariates) %*% (n_covariates x n_taxa)
  # Result: (n_cohorts x n_taxa)
  cohort_log_linear_predictors <- unique_design %*% coeff_matrix
  colnames(cohort_log_linear_predictors) <- paste0("Taxon_", 1:n_taxa)
  
  # ========================================================================
  # Step 2: Simulate log(sigma) at COHORT level
  # ========================================================================
  
  # Random errors for log(sigma): taxon-specific, shared across cohorts so
  # sd_log_overdispersion does not by itself induce cohort-varying mean dispersion.
  taxon_log_sigma_errors <- rnorm(n_taxa, mean = 0, sd = sd_log_overdispersion)
  cohort_log_sigma_errors <- matrix(
    taxon_log_sigma_errors,
    nrow = n_cohorts, ncol = n_taxa,
    byrow = TRUE
  )
  
  if (!is.null(intercept_dispersion_matrix) && !is.null(intercept_dispersion)) {
    stop("Provide only one of intercept_dispersion or intercept_dispersion_matrix.")
  }
  
  if (!is.null(intercept_dispersion_matrix)) {
    intercept_dispersion_matrix <- as.matrix(intercept_dispersion_matrix)
    if (nrow(intercept_dispersion_matrix) != n_cohorts ||
        ncol(intercept_dispersion_matrix) != n_taxa) {
      stop(
        sprintf(
          "intercept_dispersion_matrix must be n_cohorts x n_taxa (%d x %d); got %d x %d.",
          n_cohorts, n_taxa,
          nrow(intercept_dispersion_matrix), ncol(intercept_dispersion_matrix)
        )
      )
    }
    intercept_dispersion_param <- NULL
  } else {
    if (is.null(intercept_dispersion)) {
      intercept_dispersion <- rep(2.0, n_taxa)
    } else if (length(intercept_dispersion) == 1L) {
      intercept_dispersion <- rep(intercept_dispersion, n_taxa)
    }
    if (length(intercept_dispersion) != n_taxa) {
      stop(
        sprintf(
          "intercept_dispersion must have length n_taxa (%d), or length 1 to recycle; got %d.",
          n_taxa, length(intercept_dispersion)
        )
      )
    }
    intercept_dispersion_matrix <- matrix(
      intercept_dispersion,
      nrow = n_cohorts,
      ncol = n_taxa,
      byrow = TRUE
    )
    intercept_dispersion_param <- intercept_dispersion
  }

  # Calculate log(sigma) at cohort level
  # log(sigma) = -k * cohort_log_mu + c + error
  cohort_log_sigma <- -log_dispersion_assoc * cohort_log_linear_predictors + 
                       intercept_dispersion_matrix + cohort_log_sigma_errors
  
  # Convert to sigma (vectorized)
  cohort_sigma <- exp(cohort_log_sigma)
  
  # Calculate compositional probabilities (mu) at cohort level using softmax
  cohort_mu <- softmax(cohort_log_linear_predictors)
  
  # ========================================================================
  # Step 3: Map each sample to its cohort and expand parameters
  # ========================================================================
  
  # For each sample, find which cohort (row of unique_design) it belongs to
  # This will be used to join samples with ground_truth_params
  sample_to_cohort <- apply(design_matrix, 1, function(row) {
    which(apply(unique_design, 1, function(urow) all(row == urow)))[1]
  })
  
  # ========================================================================
  # Step 3: Generate library sizes
  # ========================================================================
  
  library_sizes <- round(rnorm(n_samples, mean = library_size_mean, 
                                sd = library_size_sd))
  library_sizes <- pmax(library_sizes, 100)  # Ensure minimum library size
  
  # ========================================================================
  # Step 5: Create output data structures - SIMPLIFIED APPROACH
  # ========================================================================
  
  # Create sample metadata by binding design matrix directly
  sample_metadata <- data.frame(
    sample_id = paste0("Sample_", 1:n_samples),
    library_size = library_sizes,
    cohort_idx = sample_to_cohort,  # Add cohort index for joining
    design_matrix,  # Add all design matrix columns at once
    check.names = FALSE  # Preserve column names
  )
  
  # ========================================================================
  # Create Ground Truth Parameters Table from COHORT-level parameters
  # ========================================================================
  # This table shows the cohort-level parameter structure:
  # - baseline_intercept and slope: user inputs (sum-to-0)
  # - cohort_log_linear_predictors: baseline + slope * design (per cohort-taxon)
  # - cohort_mu: compositional probabilities (per cohort-taxon)
  # - cohort_log_sigma: calculated from log-linear predictors + noise
  
  # Create one row per cohort per taxon using the cohort-level matrices
  ground_truth_params_list <- vector("list", n_cohorts)
  
  for (cohort_idx in 1:n_cohorts) {
    # Extract cohort-level parameters directly from matrices
    cohort_log_linear_pred <- cohort_log_linear_predictors[cohort_idx, ]
    cohort_log_sig <- cohort_log_sigma[cohort_idx, ]
    cohort_sigma_vals <- cohort_sigma[cohort_idx, ]
    cohort_mu_vals <- cohort_mu[cohort_idx, ]
    
    # Create simplified data frame for this cohort
    ground_truth_params_list[[cohort_idx]] <- data.frame(
      taxon_id = paste0("Taxon_", 1:n_taxa),
      cohort_idx = cohort_idx,
      baseline_intercept = baseline_intercepts,
      slope = slope_vector,
      cohort_log_linear_predictors = cohort_log_linear_pred,
      cohort_mu = cohort_mu_vals,
      cohort_log_sigma = cohort_log_sig,
      cohort_sigma = cohort_sigma_vals,  # Add sigma for easier simulation
      stringsAsFactors = FALSE
    )
  }
  
  ground_truth_params <- do.call(rbind, ground_truth_params_list)
  
  # ========================================================================
  # Create count_long by expanding ground_truth_params to sample level
  # ========================================================================
  # Simple approach: For each sample, join with ground_truth_params based on cohort
  
  # Create a cross join of samples and taxa
  count_long <- expand.grid(
    sample_id = sample_metadata$sample_id,
    taxon_id = paste0("Taxon_", 1:n_taxa),
    stringsAsFactors = FALSE
  )
  
  # Join with sample metadata to get cohort_idx and library_size
  count_long <- count_long %>%
    dplyr::left_join(
      sample_metadata %>% dplyr::select("sample_id", "library_size", "cohort_idx"),
      by = "sample_id"
    )

  # Join with ground_truth_params to get cohort-level parameters
  count_long <- count_long %>%
    dplyr::left_join(
      ground_truth_params %>% dplyr::select(
        taxon_id, cohort_idx,
        cohort_mu, cohort_sigma,
        cohort_log_linear_predictors, cohort_log_sigma
      ),
      by = c("taxon_id", "cohort_idx")
    )

  # Rename cohort_ columns to sample-level columns (they're the same within a cohort).
  # Use dplyr::rename explicitly: S4Vectors (loaded by sccomp) also provides rename().
  count_long <- count_long %>%
    dplyr::rename(
      mu = cohort_mu,
      sigma = cohort_sigma,
      unconstrained_log_mu = cohort_log_linear_predictors,
      log_sigma = cohort_log_sigma
    )

  # Simulate counts using pmap
  count_long <- count_long %>%
    dplyr::mutate(
      count = pmap_dbl(
        list(library_size, mu, sigma),
        function(n, mu, sigma) simulate_beta_binomial(n = n, mu = mu, sigma = sigma, n_sim = 1)
      )
    )

  # Add remaining sample metadata (design matrix columns)
  count_long <- count_long %>%
    dplyr::left_join(sample_metadata, by = c("sample_id", "library_size", "cohort_idx"))

  count_long <- count_long %>% dplyr::select(-dplyr::any_of("cohort_idx"))

  # Factor `group` for betadisper is *not* set here — add it in the Quarto report after
  # simulate_compositional_bb() so labels match each benchmark dataset.

  # Also create legacy taxon_metadata for backward compatibility
  taxon_metadata <- ground_truth_params
  
  # ========================================================================
  # Return results
  # ========================================================================
  
  return(list(
    # *** Ground Truth Parameters Table ***
    # This is the key table showing the full parameter generation pipeline
    # Use this for all theoretical plots in reports
    # Structure: one row per cohort per taxon
    ground_truth_params = ground_truth_params,
    
    # *** Sample-level data ***
    # Expanded from ground_truth_params: one row per sample per taxon
    # Contains all parameters and simulated counts
    count_long = count_long,
    
    # Cohort-level parameters (matrices for advanced users)
    cohort_log_linear_predictors = cohort_log_linear_predictors,  # (n_cohorts x n_taxa)
    cohort_log_sigma = cohort_log_sigma,                           # (n_cohorts x n_taxa)
    cohort_sigma = cohort_sigma,                                   # (n_cohorts x n_taxa)
    cohort_mu = cohort_mu,                                         # (n_cohorts x n_taxa)
    unique_design = unique_design,                                 # (n_cohorts x n_covariates)
    sample_to_cohort = sample_to_cohort,                           # (n_samples) - cohort index per sample
    n_cohorts = n_cohorts,
    
    # Metadata
    sample_metadata = sample_metadata,
    taxon_metadata = taxon_metadata,  # Legacy name, same as ground_truth_params
    design_matrix = design_matrix,
    
    # Input parameters (for reference)
    mu_inv_softmax = mu_inv_softmax,  # Input baseline intercepts (vector, length = n_taxa)
    slope_vector = slope_vector,       # Input slopes (vector, length = n_taxa)
    parameters = list(
      slope_vector = slope_vector,
      mu_inv_softmax = mu_inv_softmax,
      baseline_intercepts = baseline_intercepts,
      log_dispersion_assoc = log_dispersion_assoc,
      n_taxa = n_taxa,
      n_samples = n_samples,
      n_cohorts = n_cohorts,
      sd_log_overdispersion = sd_log_overdispersion,
      intercept_dispersion = intercept_dispersion_param,
      intercept_dispersion_matrix = intercept_dispersion_matrix,
      library_size_mean = library_size_mean,
      library_size_sd = library_size_sd
    )
  ))
}

