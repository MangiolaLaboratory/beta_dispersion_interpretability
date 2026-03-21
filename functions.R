# =============================================================================
# functions.R - Simulation and Beta-Dispersion Analysis for Compositional Data
# =============================================================================
#
# This file provides functions for simulating compositional count data (e.g.,
# microbiome) and building association-adjusted distance matrices for beta-
# dispersion analysis (vegan::betadisper). It supports the paper:
# "Beta dispersion lacks interpretability for differential stochastic dispersion
# analyses" (beta-dispersion interpretability project).
#
# =============================================================================

library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)

# ============================================================================
# Helper Functions
# ============================================================================

#' Softmax: convert log-linear predictors to compositional probabilities
#'
#' Exponentiates and normalizes so each row sums to 1. Used to map design-
#' matrix effects to taxon proportions. This is the standard softmax (not
#' inverse softmax / centered log-ratio).
#'
#' @param log_linear_predictors Matrix (n_samples x n_taxa) of log-space values.
#' @return Matrix of probabilities (n_samples x n_taxa), rows sum to 1.
softmax <- function(log_linear_predictors) {
  exp_log <- exp(log_linear_predictors)
  row_sums <- rowSums(exp_log)
  probabilities <- exp_log / row_sums
  return(probabilities)
}

#' Simulate beta-binomial counts
#'
#' Draws counts from BetaBinom(size=n, mu, sigma) using VGAM. Higher sigma
#' implies more overdispersion relative to binomial.
#'
#' @param n Number of trials (library size).
#' @param mu Mean probability in (0, 1).
#' @param sigma Dispersion parameter (higher = more overdispersion).
#' @param n_sim Number of draws (default 1).
#' @return Numeric vector of counts (length n_sim).
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

#' Standard deviation of beta-binomial proportion from overdispersion
#'
#' Computes SD of p_hat = X/n under BetaBinom(n, mu, overdispersion). Supports
#' sigma (dispersion), factor (n+kappa)/(kappa+1), or ratio (Var/Var_binomial).
#'
#' @param mu Mean probability in (0, 1).
#' @param overdispersion Overdispersion parameter (interpretation by type).
#' @param n Number of trials.
#' @param overdispersion_type One of "sigma", "factor", "ratio".
#' @return Numeric scalar SD.
calculate_sd_from_overdispersion <- function(mu, overdispersion, n, 
                                            overdispersion_type = c("sigma", "factor", "ratio")) {
  overdispersion_type <- match.arg(overdispersion_type)
  
  # Calculate concentration (κ) based on overdispersion type
  if (overdispersion_type == "sigma") {
    # sigma represents dispersion: concentration = 1/sigma
    sigma <- overdispersion
    concentration <- 1 / sigma
  } else if (overdispersion_type == "factor") {
    # Overdispersion factor: (n + κ) / (κ + 1) = overdispersion
    # where κ is concentration (κ = 1/σ for dispersion σ)
    # Solve for κ: overdispersion * (κ + 1) = n + κ
    # overdispersion * κ + overdispersion = n + κ
    # overdispersion * κ - κ = n - overdispersion
    # κ * (overdispersion - 1) = n - overdispersion
    # κ = (n - overdispersion) / (overdispersion - 1)
    if (overdispersion <= 1) {
      stop("Overdispersion factor must be > 1")
    }
    concentration <- (n - overdispersion) / (overdispersion - 1)
    # Ensure concentration is positive
    if (concentration <= 0) {
      stop("Invalid overdispersion factor: results in non-positive concentration")
    }
  } else if (overdispersion_type == "ratio") {
    # Overdispersion ratio: Var / Var_binomial = overdispersion
    # Var = n * μ * (1-μ) * (n + κ) / (κ + 1)
    # Var_binomial = n * μ * (1-μ)
    # Ratio = (n + κ) / (κ + 1) = overdispersion
    # This is the same as factor, so use same calculation
    if (overdispersion <= 1) {
      stop("Overdispersion ratio must be > 1")
    }
    concentration <- (n - overdispersion) / (overdispersion - 1)
    if (concentration <= 0) {
      stop("Invalid overdispersion ratio: results in non-positive concentration")
    }
  }
  
  # Calculate variance: Var = n * μ * (1-μ) * (n + κ) / (κ + 1)
  # where κ is concentration
  variance <- n * mu * (1 - mu) * (n + concentration) / (concentration + 1)
  
  # Calculate standard deviation
  sd <- sqrt(variance)
  
  return(sd)
}

#' @describeIn calculate_sd_from_overdispersion Wrapper using sigma (dispersion).
calculate_sd_from_sigma <- function(mu, sigma, n) {
  return(calculate_sd_from_overdispersion(mu, sigma, n, overdispersion_type = "sigma"))
}

#' @describeIn calculate_sd_from_overdispersion Wrapper using overdispersion factor.
calculate_sd_from_overdispersion_factor <- function(mu, overdispersion_factor, n) {
  return(calculate_sd_from_overdispersion(mu, overdispersion_factor, n, 
                                          overdispersion_type = "factor"))
}

#' @describeIn calculate_sd_from_overdispersion Wrapper using overdispersion ratio.
calculate_sd_from_overdispersion_ratio <- function(mu, overdispersion_ratio, n) {
  return(calculate_sd_from_overdispersion(mu, overdispersion_ratio, n, 
                                          overdispersion_type = "ratio"))
}

# ============================================================================
# Main Simulation Function
# ============================================================================

#' Simulate compositional beta-binomial count data
#'
#' Generates microbiome-like count data where: (1) compositional means come
#' from softmax of log-linear predictors (design_matrix %*% coeffs); (2) counts
#' are beta-binomial with mean-dispersion association log(sigma) = -k*logit(mu)
#' + c + N(0, sd^2); (3) library sizes vary per sample.
#'
#' @param slope_vector Length-n_taxa slopes (must sum to 0). Group/covariate effect per taxon.
#' @param mu_inv_softmax Length-n_taxa baseline log-predictors (must sum to 0).
#' @param log_dispersion_assoc Scalar k in log(sigma) = -k*logit(mu) + c.
#' @param n_taxa Number of taxa.
#' @param n_samples Total samples (split by design).
#' @param sd_log_overdispersion SD of log(sigma) noise (taxon-level heterogeneity).
#' @param intercept_dispersion Scalar or per-cohort intercept c. Higher = more overdispersion.
#' @param library_size_mean,library_size_sd Library size distribution.
#' @param design_matrix Optional; default intercept + one group covariate.
#' @param seed Optional random seed.
#' @return List with count_long, sample_metadata, ground_truth_params, cohort matrices, etc.
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
  intercept_dispersion = 2.0,  # Intercept 'c' in log(σ) = -k·logit(μ) + c
                          # NOTE: σ = dispersion (higher σ = MORE overdispersion)
                          # Negative values like -1 give LOW overdispersion
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
  
  # Generate random errors for log(sigma) that are TAXON-specific but shared across cohorts.
  # This ensures that `sd_log_overdispersion > 0` creates baseline OD heterogeneity across taxa
  # without introducing differential OD between conditions when `intercept_dispersion` is scalar.
  taxon_log_sigma_errors <- rnorm(n_taxa, mean = 0, sd = sd_log_overdispersion)
  cohort_log_sigma_errors <- matrix(
    taxon_log_sigma_errors,
    nrow = n_cohorts, ncol = n_taxa,
    byrow = TRUE
  )
  
  # Allow cohort/group-specific dispersion intercepts
  # - scalar: same intercept for all cohorts
  # - vector length n_cohorts: one intercept per cohort (e.g., Group1 vs Group2)
  intercept_dispersion_vec <- intercept_dispersion
  if (length(intercept_dispersion_vec) == 1) {
    intercept_dispersion_vec <- rep(intercept_dispersion_vec, n_cohorts)
  }
  if (length(intercept_dispersion_vec) != n_cohorts) {
    stop(
      sprintf(
        "intercept_dispersion must be length 1 or length n_cohorts=%d (got %d).",
        n_cohorts, length(intercept_dispersion_vec)
      )
    )
  }
  intercept_dispersion_matrix <- matrix(intercept_dispersion_vec, nrow = n_cohorts, ncol = n_taxa, byrow = FALSE)

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
  
  # Add cohort label to sample_metadata
  if ("Group" %in% colnames(unique_design)) {
    sample_metadata$cohort_label <- ifelse(sample_metadata$Group == 1, "IBD", "non-IBD")
    # Also add lowercase 'group' column for convenience
    sample_metadata$group <- sample_metadata$cohort_label
  } else {
    sample_metadata$cohort_label <- paste0("Cohort_", sample_metadata$cohort_idx)
  }
  
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
    # Get design matrix values for this cohort
    cohort_design <- unique_design[cohort_idx, , drop = TRUE]
    
    # Create a name/label for this cohort
    if ("Group" %in% colnames(unique_design)) {
      # Design matrix: non-IBD has Group=0, IBD has Group=1
      cohort_label <- ifelse(cohort_design["Group"] == 1, "IBD", "non-IBD")
    } else {
      cohort_label <- paste0("Cohort_", cohort_idx)
    }
    
    # Extract cohort-level parameters directly from matrices
    cohort_log_linear_pred <- cohort_log_linear_predictors[cohort_idx, ]
    cohort_log_sig <- cohort_log_sigma[cohort_idx, ]
    cohort_sigma_vals <- cohort_sigma[cohort_idx, ]
    cohort_mu_vals <- cohort_mu[cohort_idx, ]
    
    # Create simplified data frame for this cohort
    ground_truth_params_list[[cohort_idx]] <- data.frame(
      taxon_id = paste0("Taxon_", 1:n_taxa),
      cohort_label = cohort_label,
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
  
  # Rename 'cohort_label' to 'group' if it's a group-based design for backward compatibility
  if ("Group" %in% colnames(unique_design)) {
    colnames(ground_truth_params)[colnames(ground_truth_params) == "cohort_label"] <- "group"
  }
  
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
    left_join(sample_metadata %>% select(sample_id, library_size, cohort_idx), 
              by = "sample_id")
  
  # Join with ground_truth_params to get cohort-level parameters
  count_long <- count_long %>%
    left_join(ground_truth_params %>% select(taxon_id, cohort_idx, 
                                              cohort_mu, cohort_sigma, 
                                              cohort_log_linear_predictors, cohort_log_sigma),
              by = c("taxon_id", "cohort_idx"))
  
  # Rename cohort_ columns to sample-level columns (they're the same within a cohort)
  count_long <- count_long %>%
    rename(
      mu = cohort_mu,
      sigma = cohort_sigma,
      unconstrained_log_mu = cohort_log_linear_predictors,
      log_sigma = cohort_log_sigma
    )
  
  # Simulate counts using pmap
  count_long <- count_long %>%
    mutate(
      count = pmap_dbl(
        list(library_size, mu, sigma),
        function(n, mu, sigma) simulate_beta_binomial(n = n, mu = mu, sigma = sigma, n_sim = 1)
      )
    )
  
  # Add remaining sample metadata (design matrix columns)
  count_long <- count_long %>%
    left_join(sample_metadata, by = c("sample_id", "library_size", "cohort_idx"))
  
  # Remove cohort_idx and cohort_label as they're internal
  count_long <- count_long %>%
    select(-cohort_idx, -cohort_label)
  
  # Add lowercase 'group' column for backward compatibility (if Group column exists)
  if ("Group" %in% colnames(count_long)) {
    count_long <- count_long %>%
      mutate(group = ifelse(Group == 1, "IBD", "non-IBD"))
  }
  
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
      intercept_dispersion = intercept_dispersion_vec,
      library_size_mean = library_size_mean,
      library_size_sd = library_size_sd
    )
  ))
}

# ============================================================================
# Visualization Functions
# ============================================================================

#' Plot mean-dispersion relationship from simulation output
#'
#' Scatter of log(sigma) vs logit(mu) with fitted regression line. Visualizes
#' the inverse mean-dispersion association and residual heterogeneity.
#'
#' @param sim_result Output from simulate_compositional_bb.
#' @return ggplot object.
plot_mean_dispersion_relationship <- function(sim_result) {
  df <- sim_result$count_long %>%
    mutate(
      logit_mu = log(mu / (1 - mu)),
      log_sigma = log(sigma)
    )
  
  p <- ggplot(df, aes(x = logit_mu, y = log_sigma)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1.2) +
    labs(
      title = "Mean-Dispersion Relationship",
      subtitle = "log(σ) vs logit(μ) with regression line",
      x = "logit(μ)",
      y = "log(σ)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(p)
}

#' Plot stacked compositional abundances (bar chart)
#'
#' Shows proportion of each taxon per sample for a subset of samples.
#'
#' @param sim_result Output from simulate_compositional_bb.
#' @param n_samples_plot Number of samples to show (default 20).
#' @return ggplot object.
plot_compositional_abundances <- function(sim_result, n_samples_plot = 20) {
  count_long <- sim_result$count_long %>%
    group_by(sample_id) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup()
  
  # Select subset of samples
  sample_subset <- unique(count_long$sample_id)[1:min(n_samples_plot, 
                                                       length(unique(count_long$sample_id)))]
  count_long_subset <- count_long %>%
    filter(sample_id %in% sample_subset)
  
  p <- ggplot(count_long_subset, aes(x = sample_id, y = proportion, fill = taxon_id)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
      title = "Compositional Abundances",
      subtitle = paste("First", length(sample_subset), "samples"),
      x = "Sample",
      y = "Proportion"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  return(p)
}

#' Plot mean sigma (dispersion) by taxon with error bars
#'
#' Summarizes dispersion per taxon across samples. Useful for checking
#' taxon-level heterogeneity in overdispersion.
#'
#' @param sim_result Output from simulate_compositional_bb.
#' @return ggplot object.
plot_sigma_by_taxon <- function(sim_result) {
  df <- sim_result$count_long %>%
    group_by(taxon_id) %>%
    summarise(
      mean_sigma = mean(sigma),
      median_sigma = median(sigma),
      sd_sigma = sd(sigma),
      mean_mu = mean(mu)
    )
  
  p <- ggplot(df, aes(x = taxon_id, y = mean_sigma)) +
    geom_point(size = 2, color = "steelblue") +
    geom_errorbar(aes(ymin = mean_sigma - sd_sigma, 
                      ymax = mean_sigma + sd_sigma), 
                  width = 0.2, color = "steelblue") +
    labs(
      title = "Mean Sigma by Taxon",
      subtitle = "With standard deviation error bars",
      x = "Taxon",
      y = "Mean Sigma (σ)"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# ============================================================================
# Association-Adjusted Distance Builders (for vegan::betadisper)
# ============================================================================
#
# All build_*_betadisper_inputs functions return list(dist, group) suitable for
# vegan::betadisper(dist, group). They compute distance matrices that account
# for the mean-dispersion association (rare taxa have higher variance). Without
# adjustment, standard distances conflate differential abundance with differential
# dispersion.
#
# Variants:
# - build_assoc_adjusted_*: Arcsin-sqrt residuals scaled by sigma ratio (slope-only or pi-aware)
# - build_assoc_adjusted_piaware_*: Uses beta-binomial variance (delta method) for scaling
# - postmean: Regularizes p_hat with posterior mean under Beta prior (avoids 0/1 boundary)
# - postmeanT: Posterior mean on arcsin-sqrt scale (count-dependent shrinkage)
# - prevfilter / highprevfilter: Drops low-prevalence taxa before distance
# - Hellinger: Uses sqrt(p) transform instead of arcsin-sqrt
# - build_standard_*: Conventional distances (Bray, Aitchison, etc.) for comparison
#
# ============================================================================

#' Association-adjusted arcsin-sqrt distance (non-pi-aware)
#'
#' Core method: arcsin-sqrt residuals scaled by sigma ratio from mean-dispersion
#' association. Does not use beta-binomial variance (pi-aware); scaling uses
#' slope-only adjust_residual_assoc_arcsin_bb. Requires count_long with mu,
#' unconstrained_log_mu. Returns list(dist, group) for vegan::betadisper.
#'
#' Pipeline: arcsin-sqrt residuals -> scale by sigma ratio -> add centroid -> Euclidean dist.
build_assoc_adjusted_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  group_levels = c("non-IBD", "IBD")
) {
  # Expected columns
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  # Map intercept to each row (scalar or per-group)
  .intercept_by_group <- function(group_chr) {
    if (!is.numeric(log_sigma_intercept)) stop("log_sigma_intercept must be numeric.")
    if (length(log_sigma_intercept) == 1) return(rep(log_sigma_intercept, length(group_chr)))

    if (!is.null(names(log_sigma_intercept))) {
      out <- unname(log_sigma_intercept[group_chr])
      if (anyNA(out)) {
        stop("log_sigma_intercept has names but is missing entries for some groups.")
      }
      return(out)
    }

    if (length(log_sigma_intercept) == length(group_levels)) {
      named <- stats::setNames(log_sigma_intercept, group_levels)
      return(unname(named[group_chr]))
    }

    stop("log_sigma_intercept must be length 1, or named by group, or length(group_levels).")
  }

  df <- count_long %>%
    mutate(
      p_hat = count_to_proportion(count = count, library_size = library_size, eps = 0),
      obs_t = arcsin_sqrt(p_hat, eps = 0),
      exp_t = arcsin_sqrt(mu, eps = 0),
      residual = residual_subtract(obs_t, exp_t),
      # Association adjustment for arcsin-sqrt residuals should be done on the
      # *variance scale of the arcsin-sqrt residuals* (beta-binomial), not by
      # dividing by sqrt(sigma_rel).
      assoc_adj_residual = adjust_residual_assoc_arcsin_bb(
        residual = residual,
        mu_inv_softmax = unconstrained_log_mu,
        log_dispersion_assoc = log_dispersion_assoc,
        library_size = library_size,
        log_sigma_intercept = .intercept_by_group(as.character(group))
      ),
      centroid = arcsin_sqrt(mean(mu), eps = 0),
      .by = c(group, taxon_id)
    ) %>%
    mutate(
      assoc_adj_with_location = translate_residual_location(
        residual = assoc_adj_residual,
        expected_location = centroid
      )
    )

  mat <- df %>%
    select(sample_id, taxon_id, assoc_adj_with_location) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = assoc_adj_with_location) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  dist_obj <- stats::dist(mat, method = "euclidean")

  group_factor <- df %>%
    distinct(sample_id, group) %>%
    arrange(match(sample_id, rownames(mat))) %>%
    pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

#' Pi-aware association-adjusted arcsin-sqrt distance
#'
#' Like build_assoc_adjusted_betadisper_inputs but uses beta-binomial variance
#' (delta method) for residual scaling: sd0/sd_exp from bb_arcsin_sqrt_sd_piaware.
#' Better calibrated when library size and mean vary. Requires mu, unconstrained_log_mu.
#' Returns list(dist, group) for vegan::betadisper.
build_assoc_adjusted_piaware_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  .intercept_by_group <- function(group_chr) {
    if (!is.numeric(log_sigma_intercept)) stop("log_sigma_intercept must be numeric.")
    if (length(log_sigma_intercept) == 1) return(rep(log_sigma_intercept, length(group_chr)))
    if (!is.null(names(log_sigma_intercept))) {
      out <- unname(log_sigma_intercept[group_chr])
      if (anyNA(out)) stop("log_sigma_intercept (named) missing some groups.")
      return(out)
    }
    if (length(log_sigma_intercept) == length(group_levels)) {
      named <- stats::setNames(log_sigma_intercept, group_levels)
      return(unname(named[group_chr]))
    }
    stop("log_sigma_intercept must be length 1, or named by group, or length(group_levels).")
  }

  df <- count_long %>%
    dplyr::mutate(
      p_hat = count_to_proportion(count = count, library_size = library_size, eps = 0),
      obs_t = arcsin_sqrt(p_hat, eps = 0),
      exp_t = arcsin_sqrt(mu, eps = 0),
      residual = residual_subtract(obs_t, exp_t),

      c = .intercept_by_group(as.character(group)),
      sigma0 = exp(c),
      rho0 = sigma_to_rho(sigma0),

      sigma_exp = exp(-log_dispersion_assoc * unconstrained_log_mu + c),
      rho_exp = sigma_to_rho(sigma_exp),

      # pi-aware SD in arcsin-sqrt space (second-order delta method)
      sd0 = bb_arcsin_sqrt_sd_piaware(pi = mu, library_size = library_size, rho = rho0, eps = eps_pi),
      sd_exp = bb_arcsin_sqrt_sd_piaware(pi = mu, library_size = library_size, rho = rho_exp, eps = eps_pi),

      assoc_adj_residual = residual * (sd0 / sd_exp),
      centroid = arcsin_sqrt(mean(mu), eps = 0),
      .by = c(group, taxon_id)
    ) %>%
    dplyr::mutate(
      assoc_adj_with_location = translate_residual_location(
        residual = assoc_adj_residual,
        expected_location = centroid
      )
    )

  mat <- df %>%
    dplyr::select(sample_id, taxon_id, assoc_adj_with_location) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = assoc_adj_with_location) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  dist_obj <- stats::dist(mat, method = "euclidean")

  group_factor <- df %>%
    dplyr::distinct(sample_id, group) %>%
    dplyr::arrange(match(sample_id, rownames(mat))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

# Build pi-aware association-adjusted arcsin-sqrt distance using posterior-mean regularization.
#
# Key idea:
# - Replace p_hat = X/N with posterior mean p_tilde = E[p | X] under Beta prior:
#     p ~ Beta(alpha, beta),  X|p ~ Binomial(N, p)
#   where prior mean is pi (= mu) and prior strength is implied by rho(mu).
# - This ensures p_tilde ∈ (0,1) even when X ∈ {0, N}, avoiding boundary blow-ups.
# - Then compute residuals in arcsin-sqrt space and apply the same pi-aware SD-ratio scaling.
#
# Returns: list(dist, group)
build_assoc_adjusted_piaware_postmean_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  eps_rho = 1e-8,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  .intercept_by_group <- function(group_chr) {
    if (!is.numeric(log_sigma_intercept)) stop("log_sigma_intercept must be numeric.")
    if (length(log_sigma_intercept) == 1) return(rep(log_sigma_intercept, length(group_chr)))
    if (!is.null(names(log_sigma_intercept))) {
      out <- unname(log_sigma_intercept[group_chr])
      if (anyNA(out)) stop("log_sigma_intercept (named) missing some groups.")
      return(out)
    }
    if (length(log_sigma_intercept) == length(group_levels)) {
      named <- stats::setNames(log_sigma_intercept, group_levels)
      return(unname(named[group_chr]))
    }
    stop("log_sigma_intercept must be length 1, or named by group, or length(group_levels).")
  }

  df <- count_long %>%
    dplyr::mutate(
      c = .intercept_by_group(as.character(group)),
      sigma0 = exp(c),
      rho0 = sigma_to_rho(sigma0),

      sigma_exp = exp(-log_dispersion_assoc * unconstrained_log_mu + c),
      rho_exp = sigma_to_rho(sigma_exp),

      # Clamp rho_exp slightly to avoid rare boundary issues (rho must be in (0,1))
      rho_exp = pmax(pmin(rho_exp, 1 - eps_rho), eps_rho),

      # Posterior-mean regularization of p_hat under Beta prior implied by (pi=mu, rho=rho_exp)
      p_tilde = bb_posterior_mean_p(
        count = count,
        library_size = library_size,
        pi = mu,
        rho = rho_exp,
        eps = eps_pi
      ),

      obs_t = arcsin_sqrt(p_tilde, eps = 0),
      exp_t = arcsin_sqrt(mu, eps = 0),
      residual = residual_subtract(obs_t, exp_t),

      # pi-aware SD ratio (second-order delta method)
      sd0 = bb_arcsin_sqrt_sd_piaware(pi = mu, library_size = library_size, rho = rho0, eps = eps_pi),
      sd_exp = bb_arcsin_sqrt_sd_piaware(pi = mu, library_size = library_size, rho = rho_exp, eps = eps_pi),

      assoc_adj_residual = residual * (sd0 / sd_exp),
      centroid = arcsin_sqrt(mean(mu), eps = 0),
      .by = c(group, taxon_id)
    ) %>%
    dplyr::mutate(
      assoc_adj_with_location = translate_residual_location(
        residual = assoc_adj_residual,
        expected_location = centroid
      )
    )

  mat <- df %>%
    dplyr::select(sample_id, taxon_id, assoc_adj_with_location) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = assoc_adj_with_location) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  dist_obj <- stats::dist(mat, method = "euclidean")

  group_factor <- df %>%
    dplyr::distinct(sample_id, group) %>%
    dplyr::arrange(match(sample_id, rownames(mat))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

# Build pi-aware association-adjusted distance using posterior-mean regularization
# on the *arcsin-sqrt* scale (count-dependent EB-style shrinkage).
#
# This uses bb_posterior_mean_arcsin_sqrt(), which approximates E[asin(sqrt(p)) | X]
# under the Beta posterior implied by (pi=mu, rho=rho_exp). Unlike shrinking p then
# transforming, this introduces explicit dependence on X via the posterior.
#
# Returns: list(dist, group)
build_assoc_adjusted_piaware_postmeanT_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  eps_rho = 1e-8,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  .intercept_by_group <- function(group_chr) {
    if (!is.numeric(log_sigma_intercept)) stop("log_sigma_intercept must be numeric.")
    if (length(log_sigma_intercept) == 1) return(rep(log_sigma_intercept, length(group_chr)))
    if (!is.null(names(log_sigma_intercept))) {
      out <- unname(log_sigma_intercept[group_chr])
      if (anyNA(out)) stop("log_sigma_intercept (named) missing some groups.")
      return(out)
    }
    if (length(log_sigma_intercept) == length(group_levels)) {
      named <- stats::setNames(log_sigma_intercept, group_levels)
      return(unname(named[group_chr]))
    }
    stop("log_sigma_intercept must be length 1, or named by group, or length(group_levels).")
  }

  df <- count_long %>%
    dplyr::mutate(
      c = .intercept_by_group(as.character(group)),
      sigma0 = exp(c),
      rho0 = sigma_to_rho(sigma0),

      sigma_exp = exp(-log_dispersion_assoc * unconstrained_log_mu + c),
      rho_exp = sigma_to_rho(sigma_exp),
      rho_exp = pmax(pmin(rho_exp, 1 - eps_rho), eps_rho),

      # Posterior mean of transformed value (count-dependent)
      obs_t = bb_posterior_mean_arcsin_sqrt(
        count = count,
        library_size = library_size,
        pi = mu,
        rho = rho_exp,
        eps = eps_pi
      ),
      exp_t = arcsin_sqrt(mu, eps = 0),
      residual = residual_subtract(obs_t, exp_t),

      # pi-aware SD ratio (second-order delta method)
      sd0 = bb_arcsin_sqrt_sd_piaware(pi = mu, library_size = library_size, rho = rho0, eps = eps_pi),
      sd_exp = bb_arcsin_sqrt_sd_piaware(pi = mu, library_size = library_size, rho = rho_exp, eps = eps_pi),

      assoc_adj_residual = residual * (sd0 / sd_exp),
      centroid = arcsin_sqrt(mean(mu), eps = 0),
      .by = c(group, taxon_id)
    ) %>%
    dplyr::mutate(
      assoc_adj_with_location = translate_residual_location(
        residual = assoc_adj_residual,
        expected_location = centroid
      )
    )

  mat <- df %>%
    dplyr::select(sample_id, taxon_id, assoc_adj_with_location) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = assoc_adj_with_location) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  dist_obj <- stats::dist(mat, method = "euclidean")

  group_factor <- df %>%
    dplyr::distinct(sample_id, group) %>%
    dplyr::arrange(match(sample_id, rownames(mat))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

# Build pi-aware association-adjusted arcsin-sqrt distance using posterior-mean
# regularization of p_hat under the beta-binomial implied by (mu, rho_exp).
#
# Motivation: avoid exact p_hat in {0,1} and reduce boundary leverage without
# filtering/weighting, by using the conjugate Beta posterior mean:
#   p_tilde = E[p | X] where prior p ~ Beta(pi, rho_exp) (ICC parameterization).
#
# Steps:
# 1) Compute rho_exp from the mean–dispersion association at each mu_inv_softmax.
# 2) Compute posterior-mean proportion p_tilde = bb_posterior_mean_p(X, N, pi=mu, rho=rho_exp).
# 3) Compute arcsin-sqrt residuals using p_tilde instead of p_hat.
# 4) Apply the same pi-aware SD scaling (sd0/sd_exp) as Option A, then add centroid back.
#
# Returns: list(dist, group)
build_assoc_adjusted_piaware_postmean_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  eps_rho = 1e-8,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  .intercept_by_group <- function(group_chr) {
    if (!is.numeric(log_sigma_intercept)) stop("log_sigma_intercept must be numeric.")
    if (length(log_sigma_intercept) == 1) return(rep(log_sigma_intercept, length(group_chr)))
    if (!is.null(names(log_sigma_intercept))) {
      out <- unname(log_sigma_intercept[group_chr])
      if (anyNA(out)) stop("log_sigma_intercept (named) missing some groups.")
      return(out)
    }
    if (length(log_sigma_intercept) == length(group_levels)) {
      named <- stats::setNames(log_sigma_intercept, group_levels)
      return(unname(named[group_chr]))
    }
    stop("log_sigma_intercept must be length 1, or named by group, or length(group_levels).")
  }

  df <- count_long %>%
    dplyr::mutate(
      c = .intercept_by_group(as.character(group)),
      sigma0 = exp(c),
      rho0 = sigma_to_rho(sigma0),

      sigma_exp = exp(-log_dispersion_assoc * unconstrained_log_mu + c),
      rho_exp = sigma_to_rho(sigma_exp),
      rho_exp = pmax(pmin(rho_exp, 1 - eps_rho), eps_rho),

      # Regularized observed proportion via BB posterior mean under (mu, rho_exp)
      p_tilde = bb_posterior_mean_p(
        count = count,
        library_size = library_size,
        pi = mu,
        rho = rho_exp,
        eps = eps_pi
      ),

      obs_t = arcsin_sqrt(p_tilde, eps = eps_pi),
      exp_t = arcsin_sqrt(mu, eps = eps_pi),
      residual = residual_subtract(obs_t, exp_t),

      # pi-aware SD in arcsin-sqrt space (second-order delta method)
      sd0 = bb_arcsin_sqrt_sd_piaware(pi = mu, library_size = library_size, rho = rho0, eps = eps_pi),
      sd_exp = bb_arcsin_sqrt_sd_piaware(pi = mu, library_size = library_size, rho = rho_exp, eps = eps_pi),

      assoc_adj_residual = residual * (sd0 / sd_exp),
      centroid = arcsin_sqrt(mean(mu), eps = eps_pi),
      .by = c(group, taxon_id)
    ) %>%
    dplyr::mutate(
      assoc_adj_with_location = translate_residual_location(
        residual = assoc_adj_residual,
        expected_location = centroid
      )
    )

  mat <- df %>%
    dplyr::select(sample_id, taxon_id, assoc_adj_with_location) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = assoc_adj_with_location) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  dist_obj <- stats::dist(mat, method = "euclidean")

  group_factor <- df %>%
    dplyr::distinct(sample_id, group) %>%
    dplyr::arrange(match(sample_id, rownames(mat))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

# Build prevalence-filtered, pi-aware posterior-mean association-adjusted distance.
# Build covariance-whitened association-adjusted distance + aligned group factor.
#
# Option B: take Option A features and apply a pooled covariance whitening transform
# Build pooled-prevalence-weighted association-adjusted arcsin-sqrt distance (non-π-aware).
#
# Like `assoc_adj_prevweight`, but weights are computed from pooled prevalence across

# Build prevalence-filtered association-adjusted arcsin-sqrt distance (non-π-aware).
#
# Filter-only (no weighting): drops taxa with low prevalence, without altering the
# metric among retained taxa. This is often preferable for betadisper because
# feature weighting can itself induce/flip dispersion differences.
#
# Returns: list(dist, group)
build_assoc_adjusted_prevfilter_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  min_prevalence = 0.05,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(min_prevalence) || length(min_prevalence) != 1 || min_prevalence < 0 || min_prevalence > 1) {
    stop("min_prevalence must be a numeric scalar in [0,1].")
  }

  # prevalence by taxon (min across groups)
  prev_tbl <- count_long %>%
    dplyr::mutate(present = count > 0) %>%
    dplyr::group_by(group, taxon_id) %>%
    dplyr::summarise(prevalence = mean(present), .groups = "drop") %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::summarise(prevalence_min = min(prevalence), .groups = "drop")

  keep_taxa <- prev_tbl %>%
    dplyr::filter(prevalence_min >= min_prevalence) %>%
    dplyr::pull(taxon_id)

  if (length(keep_taxa) == 0) {
    stop("No taxa remain after prevalence filtering. Lower min_prevalence.")
  }

  # filter and delegate to assoc_adj builder (metric unchanged for retained taxa)
  count_long_f <- count_long %>% dplyr::filter(taxon_id %in% keep_taxa)

  build_assoc_adjusted_betadisper_inputs(
    count_long = count_long_f,
    log_dispersion_assoc = log_dispersion_assoc,
    log_sigma_intercept = log_sigma_intercept,
    group_levels = group_levels
  )
}

# Build prevalence-filtered pi-aware association-adjusted distance (Option A filter-only).
#
# Returns: list(dist, group)
build_assoc_adjusted_piaware_prevfilter_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  min_prevalence = 0.05,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(min_prevalence) || length(min_prevalence) != 1 || min_prevalence < 0 || min_prevalence > 1) {
    stop("min_prevalence must be a numeric scalar in [0,1].")
  }

  prev_tbl <- count_long %>%
    dplyr::mutate(present = count > 0) %>%
    dplyr::group_by(group, taxon_id) %>%
    dplyr::summarise(prevalence = mean(present), .groups = "drop") %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::summarise(prevalence_min = min(prevalence), .groups = "drop")

  keep_taxa <- prev_tbl %>%
    dplyr::filter(prevalence_min >= min_prevalence) %>%
    dplyr::pull(taxon_id)

  if (length(keep_taxa) == 0) {
    stop("No taxa remain after prevalence filtering. Lower min_prevalence.")
  }

  count_long_f <- count_long %>% dplyr::filter(taxon_id %in% keep_taxa)

  build_assoc_adjusted_piaware_betadisper_inputs(
    count_long = count_long_f,
    log_dispersion_assoc = log_dispersion_assoc,
    log_sigma_intercept = log_sigma_intercept,
    eps_pi = eps_pi,
    group_levels = group_levels
  )
}

# Build high-prevalence-filtered pi-aware association-adjusted distance
# Filters to only species present in most/all samples (e.g., >90% or 100%)
# This tests whether rare species drive false positives when richness differs
build_assoc_adjusted_piaware_highprevfilter_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  min_prevalence = 0.9,  # Default: present in at least 90% of samples
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(min_prevalence) || length(min_prevalence) != 1 || min_prevalence < 0 || min_prevalence > 1) {
    stop("min_prevalence must be a numeric scalar in [0,1].")
  }

  # Calculate prevalence across ALL samples (not per group)
  # This ensures we only keep species that are consistently present
  n_samples_total <- length(unique(count_long$sample_id))
  
  prev_tbl <- count_long %>%
    dplyr::mutate(present = count > 0) %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::summarise(
      prevalence_overall = sum(present) / n_samples_total,
      .groups = "drop"
    )

  keep_taxa <- prev_tbl %>%
    dplyr::filter(prevalence_overall >= min_prevalence) %>%
    dplyr::pull(taxon_id)

  if (length(keep_taxa) == 0) {
    stop("No taxa remain after high-prevalence filtering. Lower min_prevalence.")
  }

  count_long_f <- count_long %>% dplyr::filter(taxon_id %in% keep_taxa)

  build_assoc_adjusted_piaware_betadisper_inputs(
    count_long = count_long_f,
    log_dispersion_assoc = log_dispersion_assoc,
    log_sigma_intercept = log_sigma_intercept,
    eps_pi = eps_pi,
    group_levels = group_levels
  )
}

# Build pi-aware RAW (no residuals) association-adjusted distance
# Applies SD adjustment but NO residual calculation (no centroid subtraction)
# This is the arcsin-sqrt equivalent of build_assoc_adjusted_piaware_hellinger_raw
build_assoc_adjusted_piaware_raw_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  # Use arcsin-sqrt transformation but no residuals
  .intercept_by_group <- function(group_chr) {
    if (!is.numeric(log_sigma_intercept)) stop("log_sigma_intercept must be numeric.")
    if (length(log_sigma_intercept) == 1) return(rep(log_sigma_intercept, length(group_chr)))
    if (!is.null(names(log_sigma_intercept))) {
      out <- unname(log_sigma_intercept[group_chr])
      if (anyNA(out)) stop("log_sigma_intercept (named) missing some groups.")
      return(out)
    }
    if (length(log_sigma_intercept) == length(group_levels)) {
      named <- stats::setNames(log_sigma_intercept, group_levels)
      return(unname(named[group_chr]))
    }
    stop("log_sigma_intercept must be length 1, or named by group, or length(group_levels).")
  }

  df <- count_long %>%
    dplyr::mutate(
      p_hat = count_to_proportion(count = count, library_size = library_size, eps = 0),
      obs_t = arcsin_sqrt(p_hat, eps = 0),
      # NO residual calculation - use raw transformed value

      c = .intercept_by_group(as.character(group)),
      sigma0 = exp(c),
      rho0 = sigma_to_rho(sigma0),

      sigma_exp = exp(-log_dispersion_assoc * unconstrained_log_mu + c),
      rho_exp = sigma_to_rho(sigma_exp),

      # pi-aware SD in arcsin-sqrt space
      sd0 = bb_arcsin_sqrt_sd_piaware(pi = mu, library_size = library_size, rho = rho0, eps = eps_pi),
      sd_exp = bb_arcsin_sqrt_sd_piaware(pi = mu, library_size = library_size, rho = rho_exp, eps = eps_pi),

      # Apply SD adjustment directly to raw transformed value (no residual)
      assoc_adj_raw = obs_t * (sd0 / sd_exp),
      .by = c(group, taxon_id)
    )

  mat <- df %>%
    dplyr::select(sample_id, taxon_id, assoc_adj_raw) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = assoc_adj_raw) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  dist_obj <- stats::dist(mat, method = "euclidean")

  group_factor <- df %>%
    dplyr::distinct(sample_id, group) %>%
    dplyr::arrange(match(sample_id, rownames(mat))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

# Build high-prevalence-filtered pi-aware RAW (no residuals) association-adjusted distance
# Filters to species present in >90% of samples, applies SD adjustment but NO residual calculation
build_assoc_adjusted_piaware_highprevfilter_raw_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  min_prevalence = 0.9,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(min_prevalence) || length(min_prevalence) != 1 || min_prevalence < 0 || min_prevalence > 1) {
    stop("min_prevalence must be a numeric scalar in [0,1].")
  }

  n_samples_total <- length(unique(count_long$sample_id))
  
  prev_tbl <- count_long %>%
    dplyr::mutate(present = count > 0) %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::summarise(
      prevalence_overall = sum(present) / n_samples_total,
      .groups = "drop"
    )

  keep_taxa <- prev_tbl %>%
    dplyr::filter(prevalence_overall >= min_prevalence) %>%
    dplyr::pull(taxon_id)

  if (length(keep_taxa) == 0) {
    stop("No taxa remain after high-prevalence filtering. Lower min_prevalence.")
  }

  count_long_f <- count_long %>% dplyr::filter(taxon_id %in% keep_taxa)

  # Use arcsin-sqrt transformation (like standard IBD) but no residuals
  .intercept_by_group <- function(group_chr) {
    if (!is.numeric(log_sigma_intercept)) stop("log_sigma_intercept must be numeric.")
    if (length(log_sigma_intercept) == 1) return(rep(log_sigma_intercept, length(group_chr)))
    if (!is.null(names(log_sigma_intercept))) {
      out <- unname(log_sigma_intercept[group_chr])
      if (anyNA(out)) stop("log_sigma_intercept (named) missing some groups.")
      return(out)
    }
    if (length(log_sigma_intercept) == length(group_levels)) {
      named <- stats::setNames(log_sigma_intercept, group_levels)
      return(unname(named[group_chr]))
    }
    stop("log_sigma_intercept must be length 1, or named by group, or length(group_levels).")
  }

  df <- count_long_f %>%
    dplyr::mutate(
      p_hat = count_to_proportion(count = count, library_size = library_size, eps = 0),
      obs_t = arcsin_sqrt(p_hat, eps = 0),
      # NO residual calculation - use raw transformed value

      c = .intercept_by_group(as.character(group)),
      sigma0 = exp(c),
      rho0 = sigma_to_rho(sigma0),

      sigma_exp = exp(-log_dispersion_assoc * unconstrained_log_mu + c),
      rho_exp = sigma_to_rho(sigma_exp),

      # pi-aware SD in arcsin-sqrt space
      sd0 = bb_arcsin_sqrt_sd_piaware(pi = mu, library_size = library_size, rho = rho0, eps = eps_pi),
      sd_exp = bb_arcsin_sqrt_sd_piaware(pi = mu, library_size = library_size, rho = rho_exp, eps = eps_pi),

      # Apply SD adjustment directly to raw transformed value (no residual)
      assoc_adj_raw = obs_t * (sd0 / sd_exp),
      .by = c(group, taxon_id)
    )

  mat <- df %>%
    dplyr::select(sample_id, taxon_id, assoc_adj_raw) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = assoc_adj_raw) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  dist_obj <- stats::dist(mat, method = "euclidean")

  group_factor <- df %>%
    dplyr::distinct(sample_id, group) %>%
    dplyr::arrange(match(sample_id, rownames(mat))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

# Build high-prevalence-filtered RAW distance with NO adjustment (no residuals, no SD adjustment)
# Filters to species present in >90% of samples, uses raw arcsin-sqrt transformation only
build_assoc_adjusted_piaware_highprevfilter_raw_noadj_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  min_prevalence = 0.9,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(min_prevalence) || length(min_prevalence) != 1 || min_prevalence < 0 || min_prevalence > 1) {
    stop("min_prevalence must be a numeric scalar in [0,1].")
  }

  n_samples_total <- length(unique(count_long$sample_id))
  
  prev_tbl <- count_long %>%
    dplyr::mutate(present = count > 0) %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::summarise(
      prevalence_overall = sum(present) / n_samples_total,
      .groups = "drop"
    )

  keep_taxa <- prev_tbl %>%
    dplyr::filter(prevalence_overall >= min_prevalence) %>%
    dplyr::pull(taxon_id)

  if (length(keep_taxa) == 0) {
    stop("No taxa remain after high-prevalence filtering. Lower min_prevalence.")
  }

  count_long_f <- count_long %>% dplyr::filter(taxon_id %in% keep_taxa)

  # Just arcsin-sqrt transformation, no residuals, no SD adjustment
  df <- count_long_f %>%
    dplyr::mutate(
      p_hat = count_to_proportion(count = count, library_size = library_size, eps = 0),
      obs_t = arcsin_sqrt(p_hat, eps = 0),
      # NO residual calculation, NO SD adjustment - just raw transformed value
      .by = c(group, taxon_id)
    )

  mat <- df %>%
    dplyr::select(sample_id, taxon_id, obs_t) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = obs_t) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  dist_obj <- stats::dist(mat, method = "euclidean")

  group_factor <- df %>%
    dplyr::distinct(sample_id, group) %>%
    dplyr::arrange(match(sample_id, rownames(mat))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

# =============================================================================
# IBD Hellinger variant: uses sqrt(p) transformation instead of arcsin(sqrt(p))
# =============================================================================

# Hellinger transformation: sqrt(p)
hellinger_transform <- function(p, eps = 0) {
  p <- pmax(eps, pmin(1 - eps, p))
  sqrt(p)
}

# SD of Hellinger-transformed beta-binomial (delta method approximation)
# Var(sqrt(p)) ≈ (1/(4*p)) * Var(p) for p not too close to 0
# For beta-binomial: Var(p) = mu*(1-mu)*(1 + (n-1)*rho) / n
bb_hellinger_sd_piaware <- function(pi, library_size, rho, eps = 1e-8) {
  pi <- pmax(eps, pmin(1 - eps, pi))
  rho <- pmax(eps, pmin(1 - eps, rho))
  n <- library_size
  
  # Beta-binomial variance of proportion
  var_p <- pi * (1 - pi) * (1 + (n - 1) * rho) / n
  
  # Delta method: Var(sqrt(p)) ≈ Var(p) / (4*p)
  # But sqrt(p) has derivative 1/(2*sqrt(p)), so Var(sqrt(p)) = Var(p) / (4*p)
  var_hellinger <- var_p / (4 * pi)
  
  sqrt(pmax(var_hellinger, eps))
}

# Build pi-aware association-adjusted Hellinger distance
build_assoc_adjusted_piaware_hellinger_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  .intercept_by_group <- function(group_chr) {
    if (!is.numeric(log_sigma_intercept)) stop("log_sigma_intercept must be numeric.")
    if (length(log_sigma_intercept) == 1) return(rep(log_sigma_intercept, length(group_chr)))
    if (!is.null(names(log_sigma_intercept))) {
      out <- unname(log_sigma_intercept[group_chr])
      if (anyNA(out)) stop("log_sigma_intercept (named) missing some groups.")
      return(out)
    }
    if (length(log_sigma_intercept) == length(group_levels)) {
      named <- stats::setNames(log_sigma_intercept, group_levels)
      return(unname(named[group_chr]))
    }
    stop("log_sigma_intercept must be length 1, or named by group, or length(group_levels).")
  }

  df <- count_long %>%
    dplyr::mutate(
      p_hat = count_to_proportion(count = count, library_size = library_size, eps = 0),
      obs_t = hellinger_transform(p_hat, eps = 0),
      exp_t = hellinger_transform(mu, eps = 0),
      residual = residual_subtract(obs_t, exp_t),

      c = .intercept_by_group(as.character(group)),
      sigma0 = exp(c),
      rho0 = sigma_to_rho(sigma0),

      sigma_exp = exp(-log_dispersion_assoc * unconstrained_log_mu + c),
      rho_exp = sigma_to_rho(sigma_exp),

      # pi-aware SD in Hellinger space (delta method)
      sd0 = bb_hellinger_sd_piaware(pi = mu, library_size = library_size, rho = rho0, eps = eps_pi),
      sd_exp = bb_hellinger_sd_piaware(pi = mu, library_size = library_size, rho = rho_exp, eps = eps_pi),

      assoc_adj_residual = residual * (sd0 / sd_exp),
      centroid = hellinger_transform(mean(mu), eps = 0),
      .by = c(group, taxon_id)
    ) %>%
    dplyr::mutate(
      assoc_adj_with_location = translate_residual_location(
        residual = assoc_adj_residual,
        expected_location = centroid
      )
    )

  mat <- df %>%
    dplyr::select(sample_id, taxon_id, assoc_adj_with_location) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = assoc_adj_with_location) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  dist_obj <- stats::dist(mat, method = "euclidean")

  group_factor <- df %>%
    dplyr::distinct(sample_id, group) %>%
    dplyr::arrange(match(sample_id, rownames(mat))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

# Build pi-aware association-adjusted Hellinger distance (RAW - no residuals)
# This version applies the SD adjustment directly to transformed values,
# without subtracting expected values (no residual calculation).
build_assoc_adjusted_piaware_hellinger_raw_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  .intercept_by_group <- function(group_chr) {
    if (!is.numeric(log_sigma_intercept)) stop("log_sigma_intercept must be numeric.")
    if (length(log_sigma_intercept) == 1) return(rep(log_sigma_intercept, length(group_chr)))
    if (!is.null(names(log_sigma_intercept))) {
      out <- unname(log_sigma_intercept[group_chr])
      if (anyNA(out)) stop("log_sigma_intercept (named) missing some groups.")
      return(out)
    }
    if (length(log_sigma_intercept) == length(group_levels)) {
      named <- stats::setNames(log_sigma_intercept, group_levels)
      return(unname(named[group_chr]))
    }
    stop("log_sigma_intercept must be length 1, or named by group, or length(group_levels).")
  }

  df <- count_long %>%
    dplyr::mutate(
      p_hat = count_to_proportion(count = count, library_size = library_size, eps = 0),
      obs_t = hellinger_transform(p_hat, eps = 0),
      # NO residual calculation - use raw transformed value

      c = .intercept_by_group(as.character(group)),
      sigma0 = exp(c),
      rho0 = sigma_to_rho(sigma0),

      sigma_exp = exp(-log_dispersion_assoc * unconstrained_log_mu + c),
      rho_exp = sigma_to_rho(sigma_exp),

      # pi-aware SD in Hellinger space (delta method)
      sd0 = bb_hellinger_sd_piaware(pi = mu, library_size = library_size, rho = rho0, eps = eps_pi),
      sd_exp = bb_hellinger_sd_piaware(pi = mu, library_size = library_size, rho = rho_exp, eps = eps_pi),

      # Apply SD adjustment directly to raw transformed value (no residual)
      assoc_adj_raw = obs_t * (sd0 / sd_exp),
      .by = c(group, taxon_id)
    )

  mat <- df %>%
    dplyr::select(sample_id, taxon_id, assoc_adj_raw) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = assoc_adj_raw) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  dist_obj <- stats::dist(mat, method = "euclidean")

  group_factor <- df %>%
    dplyr::distinct(sample_id, group) %>%
    dplyr::arrange(match(sample_id, rownames(mat))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

# Build prevalence-filtered pi-aware Hellinger RAW association-adjusted distance
build_assoc_adjusted_piaware_hellinger_raw_prevfilter_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  min_prevalence = 0.05,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(min_prevalence) || length(min_prevalence) != 1 || min_prevalence < 0 || min_prevalence > 1) {
    stop("min_prevalence must be a numeric scalar in [0,1].")
  }

  prev_tbl <- count_long %>%
    dplyr::mutate(present = count > 0) %>%
    dplyr::group_by(group, taxon_id) %>%
    dplyr::summarise(prevalence = mean(present), .groups = "drop") %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::summarise(prevalence_min = min(prevalence), .groups = "drop")

  keep_taxa <- prev_tbl %>%
    dplyr::filter(prevalence_min >= min_prevalence) %>%
    dplyr::pull(taxon_id)

  if (length(keep_taxa) == 0) {
    stop("No taxa remain after prevalence filtering. Lower min_prevalence.")
  }

  count_long_f <- count_long %>% dplyr::filter(taxon_id %in% keep_taxa)

  build_assoc_adjusted_piaware_hellinger_raw_betadisper_inputs(
    count_long = count_long_f,
    log_dispersion_assoc = log_dispersion_assoc,
    log_sigma_intercept = log_sigma_intercept,
    eps_pi = eps_pi,
    group_levels = group_levels
  )
}

# Build prevalence-filtered pi-aware Hellinger association-adjusted distance
build_assoc_adjusted_piaware_hellinger_prevfilter_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  min_prevalence = 0.05,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(min_prevalence) || length(min_prevalence) != 1 || min_prevalence < 0 || min_prevalence > 1) {
    stop("min_prevalence must be a numeric scalar in [0,1].")
  }

  prev_tbl <- count_long %>%
    dplyr::mutate(present = count > 0) %>%
    dplyr::group_by(group, taxon_id) %>%
    dplyr::summarise(prevalence = mean(present), .groups = "drop") %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::summarise(prevalence_min = min(prevalence), .groups = "drop")

  keep_taxa <- prev_tbl %>%
    dplyr::filter(prevalence_min >= min_prevalence) %>%
    dplyr::pull(taxon_id)

  if (length(keep_taxa) == 0) {
    stop("No taxa remain after prevalence filtering. Lower min_prevalence.")
  }

  count_long_f <- count_long %>% dplyr::filter(taxon_id %in% keep_taxa)

  build_assoc_adjusted_piaware_hellinger_betadisper_inputs(
    count_long = count_long_f,
    log_dispersion_assoc = log_dispersion_assoc,
    log_sigma_intercept = log_sigma_intercept,
    eps_pi = eps_pi,
    group_levels = group_levels
  )
}

# Build high-prevalence-filtered pi-aware Hellinger association-adjusted distance
# Filters to only species present in most/all samples (e.g., >90% or 100%)
# Uses Hellinger transformation with residuals and SD adjustment
build_assoc_adjusted_piaware_hellinger_highprevfilter_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  min_prevalence = 0.9,  # Default: present in at least 90% of samples
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(min_prevalence) || length(min_prevalence) != 1 || min_prevalence < 0 || min_prevalence > 1) {
    stop("min_prevalence must be a numeric scalar in [0,1].")
  }

  # Calculate prevalence across ALL samples (not per group)
  n_samples_total <- length(unique(count_long$sample_id))
  
  prev_tbl <- count_long %>%
    dplyr::mutate(present = count > 0) %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::summarise(
      prevalence_overall = sum(present) / n_samples_total,
      .groups = "drop"
    )

  keep_taxa <- prev_tbl %>%
    dplyr::filter(prevalence_overall >= min_prevalence) %>%
    dplyr::pull(taxon_id)

  if (length(keep_taxa) == 0) {
    stop("No taxa remain after high-prevalence filtering. Lower min_prevalence.")
  }

  count_long_f <- count_long %>% dplyr::filter(taxon_id %in% keep_taxa)

  build_assoc_adjusted_piaware_hellinger_betadisper_inputs(
    count_long = count_long_f,
    log_dispersion_assoc = log_dispersion_assoc,
    log_sigma_intercept = log_sigma_intercept,
    eps_pi = eps_pi,
    group_levels = group_levels
  )
}

# Build high-prevalence-filtered pi-aware Hellinger RAW association-adjusted distance
# Combines: Hellinger transformation + high-prevalence filtering + no residuals + SD adjustment
# This hybrid approach aims to work well for both asymmetric and symmetric DA
build_assoc_adjusted_piaware_hellinger_highprevfilter_raw_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  min_prevalence = 0.9,  # Default: present in at least 90% of samples
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(min_prevalence) || length(min_prevalence) != 1 || min_prevalence < 0 || min_prevalence > 1) {
    stop("min_prevalence must be a numeric scalar in [0,1].")
  }

  # Calculate prevalence across ALL samples (not per group)
  n_samples_total <- length(unique(count_long$sample_id))
  
  prev_tbl <- count_long %>%
    dplyr::mutate(present = count > 0) %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::summarise(
      prevalence_overall = sum(present) / n_samples_total,
      .groups = "drop"
    )

  keep_taxa <- prev_tbl %>%
    dplyr::filter(prevalence_overall >= min_prevalence) %>%
    dplyr::pull(taxon_id)

  if (length(keep_taxa) == 0) {
    stop("No taxa remain after high-prevalence filtering. Lower min_prevalence.")
  }

  count_long_f <- count_long %>% dplyr::filter(taxon_id %in% keep_taxa)

  # Use Hellinger raw (no residuals, with SD adjustment) on high-prevalence filtered data
  build_assoc_adjusted_piaware_hellinger_raw_betadisper_inputs(
    count_long = count_long_f,
    log_dispersion_assoc = log_dispersion_assoc,
    log_sigma_intercept = log_sigma_intercept,
    eps_pi = eps_pi,
    group_levels = group_levels
  )
}

#' Standard distance matrix for betadisper (Bray-Curtis, Aitchison, etc.)
#'
#' Uses vegan::vegdist (or decostand + dist for Hellinger). No association
#' adjustment; useful as baseline comparison. Supports "bray", "aitchison",
#' "robust.aitchison", "jaccard", "hellinger".
#'
#' @param count_long Long-format counts (sample_id, taxon_id, count).
#' @param sample_metadata Must have sample_id, group.
#' @param distance_method One of bray, aitchison, robust.aitchison, jaccard, hellinger.
#' @param group_levels Factor levels for group.
#' @return list(dist, group) for vegan::betadisper.
build_standard_betadisper_inputs <- function(
  count_long,
  sample_metadata,
  distance_method = "bray",
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c("sample_id", "taxon_id", "count")
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!("sample_id" %in% colnames(sample_metadata)) || !("group" %in% colnames(sample_metadata))) {
    stop("sample_metadata must contain columns: sample_id, group")
  }

  count_matrix <- count_long %>%
    dplyr::select(sample_id, taxon_id, count) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = count) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  nonzero_samples <- rowSums(count_matrix) > 0
  count_matrix <- count_matrix[nonzero_samples, , drop = FALSE]

  dist_obj <- vegan::vegdist(count_matrix, method = distance_method)
  if (any(is.na(dist_obj)) || any(is.infinite(dist_obj))) {
    stop("Distance matrix contains NA/Inf values.")
  }

  group_factor <- sample_metadata %>%
    dplyr::filter(sample_id %in% rownames(count_matrix)) %>%
    dplyr::arrange(match(sample_id, rownames(count_matrix))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

#' Arcsin-sqrt residual distance with group centroid restored
#'
#' Computes r_ij = asin(sqrt(p_hat)) - asin(sqrt(mu)), then adds back the
#' group×taxon centroid so groups separate in PCoA. Preserves within-group
#' dispersion; used for residual-based beta-dispersion without association
#' adjustment. Requires count_long with mu (e.g. from simulation).
#'
#' @param count_long Long-format with sample_id, taxon_id, count, library_size, mu.
#' @param sample_metadata Must have sample_id, group.
#' @param group_levels Factor levels.
#' @return list(dist, group) for vegan::betadisper.
build_arcsin_residual_with_location_betadisper_inputs <- function(
  count_long,
  sample_metadata,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c("sample_id", "taxon_id", "count", "library_size", "mu")
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!("sample_id" %in% colnames(sample_metadata)) || !("group" %in% colnames(sample_metadata))) {
    stop("sample_metadata must contain columns: sample_id, group")
  }

  df <- count_long %>%
    dplyr::mutate(
      p_hat = count_to_proportion(count = count, library_size = library_size, eps = 0),
      obs_t = arcsin_sqrt(p_hat, eps = 0),
      exp_t = arcsin_sqrt(mu, eps = 0),
      resid = residual_subtract(obs_t, exp_t)
    )

  # Use group from count_long if present; otherwise join from sample_metadata.
  if ("group" %in% colnames(df)) {
    df <- df %>% dplyr::mutate(group = as.character(group))
  } else {
    df <- df %>%
      dplyr::left_join(
        sample_metadata %>% dplyr::select(sample_id, group),
        by = "sample_id"
      ) %>%
      dplyr::mutate(group = as.character(group))
  }

  df <- df %>%
    dplyr::mutate(
      centroid = mean(exp_t),
      .by = c(group, taxon_id)
    ) %>%
    dplyr::mutate(
      resid_with_location = translate_residual_location(
        residual = resid,
        expected_location = centroid
      )
    )

  mat <- df %>%
    dplyr::select(sample_id, taxon_id, resid_with_location) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = resid_with_location) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  dist_obj <- stats::dist(mat, method = "euclidean")

  group_factor <- df %>%
    dplyr::distinct(sample_id, group) %>%
    dplyr::arrange(match(sample_id, rownames(mat))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

#' Alpha diversity analysis (richness, Shannon, evenness) with group comparison
#'
#' Computes richness (observed taxa), Shannon entropy, and Pielou evenness per
#' sample. Summarizes by group and runs t-tests. Returns a list with summary
#' table, tests, plot (PCoA of alpha metrics with stat_ellipse), and raw alpha_data.
#'
#' @param sim_result Output from simulate_compositional_bb.
#' @param case_label Label for plot title.
#' @param group_levels Factor levels for group (default c("non-IBD", "IBD")).
#' @return List with summary_tbl, tests, plot, alpha_data.
run_alpha_diversity_analysis <- function(
  sim_result,
  case_label = "Case",
  group_levels = c("non-IBD", "IBD")
) {
  count_long <- sim_result$count_long

  required_cols <- c("sample_id", "group", "count")
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
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
    dplyr::mutate(group = factor(group, levels = group_levels))

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

  # Wrap t.tests in tryCatch to handle constant data (e.g., equal baseline)
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

  # Add tiny noise to help stat_ellipse when variance is extremely low
  # (without noise, covariance matrix can be singular and ellipse fails silently)
  set.seed(12345)  # Reproducible noise
  alpha_data_jitter <- alpha_data %>%
    dplyr::mutate(
      richness_j = richness + stats::rnorm(dplyr::n(), 0, 0.01),
      evenness_j = pielou_evenness + stats::rnorm(dplyr::n(), 0, 0.001)
    )
  
  p <- ggplot2::ggplot(alpha_data, ggplot2::aes(x = richness, y = pielou_evenness, color = group)) +
    ggplot2::geom_point(size = 0.1, alpha = 0.6, stroke = 0) +
    # Use jittered data for ellipse calculation only
    ggplot2::stat_ellipse(
      data = alpha_data_jitter,
      ggplot2::aes(x = richness_j, y = evenness_j, color = group),
      level = 0.95, linewidth = 0.25, type = "norm"
    ) +
    ggplot2::scale_color_manual(values = c("non-IBD" = "#E74C3C", "IBD" = "#3498DB")) +
    ggplot2::labs(
      title = paste0("Alpha Diversity: Evenness vs Richness (", case_label, ")"),
      subtitle = "Points = samples, Ellipses = 95% confidence intervals",
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

# =============================================================================
# IBD Hellinger Dampened: reduces sensitivity to high k values
# Uses (sd0/sd_exp)^damping with damping < 1 to reduce adjustment aggressiveness
# =============================================================================

build_assoc_adjusted_piaware_hellinger_highprev_dampened_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  min_prevalence = 0.9,
  damping = 0.5,  # Damping factor: 1 = full adjustment, 0 = no adjustment
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "mu", "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  # High-prevalence filtering
  n_samples_total <- length(unique(count_long$sample_id))
  prev_tbl <- count_long %>%
    dplyr::mutate(present = count > 0) %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::summarise(prevalence_overall = sum(present) / n_samples_total, .groups = "drop")

  keep_taxa <- prev_tbl %>%
    dplyr::filter(prevalence_overall >= min_prevalence) %>%
    dplyr::pull(taxon_id)

  if (length(keep_taxa) == 0) {
    stop("No taxa remain after high-prevalence filtering.")
  }

  count_long_f <- count_long %>% dplyr::filter(taxon_id %in% keep_taxa)

  .intercept_by_group <- function(group_chr) {
    if (!is.numeric(log_sigma_intercept)) stop("log_sigma_intercept must be numeric.")
    if (length(log_sigma_intercept) == 1) return(rep(log_sigma_intercept, length(group_chr)))
    if (!is.null(names(log_sigma_intercept))) {
      out <- unname(log_sigma_intercept[group_chr])
      if (anyNA(out)) stop("log_sigma_intercept (named) missing some groups.")
      return(out)
    }
    if (length(log_sigma_intercept) == length(group_levels)) {
      named <- stats::setNames(log_sigma_intercept, group_levels)
      return(unname(named[group_chr]))
    }
    stop("log_sigma_intercept must be length 1, or named by group, or length(group_levels).")
  }

  df <- count_long_f %>%
    dplyr::mutate(
      p_hat = count_to_proportion(count = count, library_size = library_size, eps = 0),
      obs_t = hellinger_transform(p_hat, eps = 0),
      # No residual calculation - use raw transformed value

      c = .intercept_by_group(as.character(group)),
      sigma0 = exp(c),
      rho0 = sigma_to_rho(sigma0),

      sigma_exp = exp(-log_dispersion_assoc * unconstrained_log_mu + c),
      rho_exp = sigma_to_rho(sigma_exp),

      sd0 = bb_hellinger_sd_piaware(pi = mu, library_size = library_size, rho = rho0, eps = eps_pi),
      sd_exp = bb_hellinger_sd_piaware(pi = mu, library_size = library_size, rho = rho_exp, eps = eps_pi),

      # DAMPENED adjustment: (sd0/sd_exp)^damping instead of (sd0/sd_exp)
      assoc_adj_raw = obs_t * ((sd0 / sd_exp)^damping),
      .by = c(group, taxon_id)
    )

  mat <- df %>%
    dplyr::select(sample_id, taxon_id, assoc_adj_raw) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = assoc_adj_raw) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  dist_obj <- stats::dist(mat, method = "euclidean")

  group_factor <- df %>%
    dplyr::distinct(sample_id, group) %>%
    dplyr::arrange(match(sample_id, rownames(mat))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}
