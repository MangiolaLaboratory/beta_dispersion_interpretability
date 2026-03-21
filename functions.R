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


# Functions used in targets pipeline
extract_sccomp_params <- function(result) {
  fit <- attr(result, "fit")
  model_input <- attr(result, "model_input")

  if (!is.null(model_input) && is.list(model_input)) {
    n_samples_real <- model_input$N
    n_taxa_real <- model_input$M
  } else {
    n_taxa_real <- length(unique(result$cell_group))
    n_samples_real <- 178
  }

  prec_coeff_summary <- fit$summary("prec_coeff")
  prec_intercept <- prec_coeff_summary$mean[1]
  prec_slope <- prec_coeff_summary$mean[2]
  k_real <- prec_slope
  k_sd_real <- prec_coeff_summary$sd[2]
  c_real <- -prec_intercept

  prec_sd_summary <- fit$summary("prec_sd")
  prec_sd_real <- prec_sd_summary$mean[1]

  intercept_effects <- result %>%
    filter(.data$parameter == "(Intercept)") %>%
    pull(.data$c_effect)

  slope_effects <- result %>%
    filter(.data$parameter == "groupIBD") %>%
    pull(.data$c_effect)

  slope_data <- result %>%
    filter(.data$parameter == "groupIBD")

  slope_effects_significant <- slope_data %>%
    filter(.data$c_lower > 0 | .data$c_upper < 0) %>%
    pull(.data$c_effect)

  v_intercept <- result %>%
    filter(.data$parameter == "(Intercept)") %>%
    pull(.data$v_effect)

  v_slope <- result %>%
    filter(.data$parameter == "groupIBD") %>%
    pull(.data$v_effect)

  list(
    n_samples_real = n_samples_real,
    n_taxa_real = n_taxa_real,
    k_real = k_real,
    k_sd_real = k_sd_real,
    c_real = c_real,
    prec_sd_real = prec_sd_real,
    intercept_effects = intercept_effects,
    slope_effects = slope_effects,
    slope_effects_significant = slope_effects_significant,
    v_intercept = v_intercept,
    v_slope = v_slope
  )
}

build_simulation_params <- function(sccomp_params, n_reps_auc = 20) {
  set.seed(123)
  n_taxa <- min(sccomp_params$n_taxa_real, 200)
  n_samples_per_group <- floor(sccomp_params$n_samples_real / 2)
  n_groups <- 2
  n_samples <- n_samples_per_group * n_groups
  group_levels <- c("non-IBD", "IBD")

  library_size_mean <- 15125
  library_size_sd <- 5000

  sampled_intercepts <- sample(sccomp_params$intercept_effects, n_taxa, replace = TRUE)
  mu_inv_softmax_base_realistic <- sampled_intercepts - mean(sampled_intercepts)

  intercept_disp_realistic <- sccomp_params$c_real
  sigma_realistic <- exp(intercept_disp_realistic)

  set.seed(123)
  slope_realistic <- sample(sccomp_params$slope_effects, n_taxa, replace = TRUE)
  slope_realistic <- slope_realistic - mean(slope_realistic)

  perm_reps_auc <- 199
  rep_seeds_auc <- 12000 + seq_len(n_reps_auc)

  list(
    n_taxa = n_taxa,
    n_samples_per_group = n_samples_per_group,
    n_groups = n_groups,
    n_samples = n_samples,
    group_levels = group_levels,
    library_size_mean = library_size_mean,
    library_size_sd = library_size_sd,
    mu_inv_softmax_base_realistic = mu_inv_softmax_base_realistic,
    intercept_disp_realistic = intercept_disp_realistic,
    sigma_realistic = sigma_realistic,
    slope_realistic = slope_realistic,
    n_reps_auc = n_reps_auc,
    perm_reps_auc = perm_reps_auc,
    rep_seeds_auc = rep_seeds_auc
  )
}

design_matrix_from_groups <- function(n_samples_per_group) {
  group <- rep(c("non-IBD", "IBD"), each = n_samples_per_group)
  cbind(
    Intercept = 1,
    Group = ifelse(group == "IBD", 1, 0)
  )
}

detect_n_cores <- function(max_cores = NULL) {
  available <- parallelly::availableCores()
  if (is.null(max_cores)) {
    return(available)
  } else {
    return(min(available, max_cores))
  }
}



pretty_method_label <- function(method_key) {
  if (method_key == "permdisp_bray" || method_key == "bray") return("Bray–Curtis")
  if (method_key == "permdisp_jaccard" || method_key == "jaccard") return("Jaccard")
  if (method_key == "permdisp_hellinger" || method_key == "hellinger") return("Hellinger")
  if (method_key == "permdisp_aitchison" || method_key == "aitchison") return("Aitchison")
  if (method_key == "permdisp_robust.aitchison" || method_key == "robust.aitchison") return("robust Aitchison")
  if (method_key == "permdisp_euclidean") return("Euclidean")
  return(method_key)
}

pretty_method_label_roc <- function(method_key) {
  pretty_method_label(method_key)
}

pretty_method_label_box <- function(method_key) {
  pretty_method_label(method_key)
}

permutest_betadisper_fixed <- function(bd, permutations = 999) {
  pt <- vegan::permutest(bd, pairwise = TRUE, permutations = permutations)
  list(
    p_value = unname(pt$tab$`Pr(>F)`[1]),
    n_perm = permutations
  )
}

permutest_betadisper_stable <- function(bd, max_perm = 99999, batch = 5000, target_se = 0.0025) {
  if (!is.numeric(max_perm) || length(max_perm) != 1 || max_perm < 999) stop("max_perm must be >= 999")
  if (!is.numeric(batch) || length(batch) != 1 || batch < 100) stop("batch must be >= 100")
  if (!is.numeric(target_se) || length(target_se) != 1 || target_se <= 0) stop("target_se must be > 0")

  m <- 0L
  p <- NA_real_
  se <- Inf

  while (m < max_perm && se > target_se) {
    m <- as.integer(min(max_perm, m + batch))
    pt <- vegan::permutest(bd, pairwise = TRUE, permutations = m)
    p <- unname(pt$tab$`Pr(>F)`[1])
    if (!is.finite(p)) {
      se <- Inf
    } else {
      se <- sqrt(p * (1 - p) / m)
    }
  }

  list(p_value = p, p_se = se, n_perm = m)
}

build_standard_betadisper_inputs_paper <- function(
  count_long,
  sample_metadata,
  distance_method,
  group_levels = c("non-IBD", "IBD"),
  clr_pseudocount = 1
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
    dplyr::select(.data$sample_id, .data$taxon_id, .data$count) %>%
    tidyr::pivot_wider(names_from = "taxon_id", values_from = "count") %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  nonzero_samples <- rowSums(count_matrix) > 0
  count_matrix <- count_matrix[nonzero_samples, , drop = FALSE]

  if (distance_method %in% c("aitchison", "robust.aitchison")) {
    if (!is.numeric(clr_pseudocount) || length(clr_pseudocount) != 1 || !is.finite(clr_pseudocount) || clr_pseudocount <= 0) {
      stop("clr_pseudocount must be a positive finite numeric scalar (or NULL).")
    }
    if (any(count_matrix <= 0, na.rm = TRUE)) {
      count_matrix[count_matrix <= 0] <- clr_pseudocount
    }
  }

  if (distance_method == "jaccard") {
    count_matrix <- vegan::decostand(count_matrix, method = "total", MARGIN = 1)
  }

  dist_obj <- if (distance_method == "hellinger") {
    hel <- vegan::decostand(count_matrix, method = "hellinger", MARGIN = 1)
    stats::dist(hel, method = "euclidean")
  } else if (distance_method == "jaccard") {
    vegan::vegdist(count_matrix, method = "jaccard", binary = FALSE)
  } else {
    vegan::vegdist(count_matrix, method = distance_method)
  }

  if (any(is.na(dist_obj)) || any(is.infinite(dist_obj))) {
    stop("Distance matrix contains NA/Inf values.")
  }

  group_factor <- sample_metadata %>%
    dplyr::filter(.data$sample_id %in% rownames(count_matrix)) %>%
    dplyr::arrange(match(.data$sample_id, rownames(count_matrix))) %>%
    dplyr::pull(.data$group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

simulate_case_realistic_taxa <- function(
  slope_vector,
  mu_inv_softmax,
  log_dispersion_assoc,
  intercept_dispersion,
  sd_log_overdispersion,
  seed,
  n_taxa,
  n_samples_per_group,
  library_size_mean,
  library_size_sd
) {
  n_groups <- 2
  n_samples <- n_samples_per_group * n_groups
  simulate_compositional_bb( # nolint
    slope_vector = slope_vector,
    mu_inv_softmax = mu_inv_softmax,
    log_dispersion_assoc = log_dispersion_assoc,
    n_taxa = n_taxa,
    n_samples = n_samples,
    sd_log_overdispersion = sd_log_overdispersion,
    intercept_dispersion = intercept_dispersion,
    library_size_mean = library_size_mean,
    library_size_sd = library_size_sd,
    design_matrix = design_matrix_from_groups(n_samples_per_group),
    seed = seed
  )
}

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

run_da_taxa_sweep <- function(
  intercept_dispersion,
  mu_base,
  k_fixed,
  n_da_taxa_grid,
  slope_effects_distribution,
  setting_name = "DA Taxa Sweep",
  representative_n_da = 5,
  representative_seed = 12001,
  sd_log_overdispersion = 0,
  n_taxa,
  n_samples,
  n_samples_per_group,
  library_size_mean,
  library_size_sd,
  rep_seeds_auc,
  perm_reps_auc
) {
  standard_distance_methods <- c("bray", "aitchison", "robust.aitchison", "jaccard", "hellinger")
  group_levels <- c("non-IBD", "IBD")

  simulate_case_local <- function(slope_vec, mu_inv_softmax, log_dispersion_assoc, seed) {
    simulate_compositional_bb( # nolint
      slope_vector = slope_vec,
      mu_inv_softmax = mu_inv_softmax,
      log_dispersion_assoc = log_dispersion_assoc,
      n_taxa = n_taxa,
      n_samples = n_samples,
      sd_log_overdispersion = sd_log_overdispersion,
      intercept_dispersion = intercept_dispersion,
      library_size_mean = library_size_mean,
      library_size_sd = library_size_sd,
      design_matrix = design_matrix_from_groups(n_samples_per_group),
      seed = seed
    )
  }

  jobs <- tidyr::crossing(
    seed = rep_seeds_auc,
    n_da_taxa = n_da_taxa_grid
  )

  run_job <- function(job_row) {
    seed <- job_row$seed
    n_da <- job_row$n_da_taxa
    slope_vec <- rep(0, n_taxa)
    if (n_da > 0) {
      set.seed(seed + 10000)
      da_indices <- sample(1:n_taxa, size = min(n_da, n_taxa), replace = FALSE)
      slope_vec[da_indices] <- sample(slope_effects_distribution, length(da_indices), replace = TRUE)
      slope_vec <- slope_vec - mean(slope_vec)
    }

    sim_r <- simulate_case_local(slope_vec, mu_base, k_fixed, seed)

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

    purrr::imap_dfr(inputs_std, \(inp, method_key) {
      bd <- vegan::betadisper(inp$dist, inp$group)
      p <- permutest_betadisper_fixed(bd, permutations = perm_reps_auc)$p_value
      tibble::tibble(
        seed = seed,
        n_da_taxa = n_da,
        method = pretty_method_label_roc(method_key),
        p_value = p
      )
    })
  }

  sweep_results <- dplyr::bind_rows(lapply(split(jobs, seq_len(nrow(jobs))), run_job))

  fpr_data <- sweep_results %>%
    dplyr::group_by(.data$n_da_taxa, .data$method) %>%
    dplyr::summarise(
      n_sims = dplyr::n(),
      n_fp = sum(.data$p_value < 0.05, na.rm = TRUE),
      FPR = sum(.data$p_value < 0.05, na.rm = TRUE) / dplyr::n(),
      .groups = "drop"
    )

  set.seed(representative_seed)
  slope_vec_rep <- rep(0, n_taxa)
  if (representative_n_da > 0) {
    da_indices_rep <- sample(1:n_taxa, size = min(representative_n_da, n_taxa), replace = FALSE)
    slope_vec_rep[da_indices_rep] <- sample(slope_effects_distribution, length(da_indices_rep), replace = TRUE)
    slope_vec_rep <- slope_vec_rep - mean(slope_vec_rep)
  }

  sim_rep <- simulate_case_local(slope_vec_rep, mu_base, k_fixed, representative_seed)

  list(
    fpr_data = fpr_data,
    representative = sim_rep
  )
}

run_da_taxa_sweep_abundance_filtered <- function(
  intercept_dispersion,
  mu_base,
  k_fixed,
  n_da_taxa_grid,
  slope_effects_distribution,
  setting_name = "DA Taxa Sweep (Filtered)",
  representative_n_da = 5,
  representative_seed = 12001,
  sd_log_overdispersion = 0,
  n_taxa,
  n_samples,
  n_samples_per_group,
  library_size_mean,
  library_size_sd,
  rep_seeds_auc,
  perm_reps_auc,
  abundance_threshold = abundance_threshold
) {
  standard_distance_methods <- c("bray", "aitchison", "robust.aitchison", "jaccard", "hellinger")
  group_levels <- c("non-IBD", "IBD")

  simulate_case_local <- function(slope_vec, mu_inv_softmax, log_dispersion_assoc, seed) {
    simulate_compositional_bb( # nolint
      slope_vector = slope_vec,
      mu_inv_softmax = mu_inv_softmax,
      log_dispersion_assoc = log_dispersion_assoc,
      n_taxa = n_taxa,
      n_samples = n_samples,
      sd_log_overdispersion = sd_log_overdispersion,
      intercept_dispersion = intercept_dispersion,
      library_size_mean = library_size_mean,
      library_size_sd = library_size_sd,
      design_matrix = design_matrix_from_groups(n_samples_per_group),
      seed = seed
    )
  }

  jobs <- tidyr::crossing(
    seed = rep_seeds_auc,
    n_da_taxa = n_da_taxa_grid
  )

  run_job <- function(job_row) {
    seed <- job_row$seed
    n_da <- job_row$n_da_taxa
    slope_vec <- rep(0, n_taxa)
    if (n_da > 0) {
      set.seed(seed + 10000)
      da_indices <- sample(1:n_taxa, size = min(n_da, n_taxa), replace = FALSE)
      slope_vec[da_indices] <- sample(slope_effects_distribution, length(da_indices), replace = TRUE)
      slope_vec <- slope_vec - mean(slope_vec)
    }

    sim_r <- simulate_case_local(slope_vec, mu_base, k_fixed, seed)
    sim_r$count_long <- filter_by_abundance_cases(sim_r$count_long, threshold = abundance_threshold)

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

    purrr::imap_dfr(inputs_std, \(inp, method_key) {
      bd <- vegan::betadisper(inp$dist, inp$group)
      p <- permutest_betadisper_fixed(bd, permutations = perm_reps_auc)$p_value
      tibble::tibble(
        seed = seed,
        n_da_taxa = n_da,
        method = pretty_method_label_roc(method_key),
        p_value = p
      )
    })
  }

  sweep_results <- dplyr::bind_rows(lapply(split(jobs, seq_len(nrow(jobs))), run_job))

  fpr_data <- sweep_results %>%
    dplyr::group_by(.data$n_da_taxa, .data$method) %>%
    dplyr::summarise(
      n_sims = dplyr::n(),
      n_fp = sum(.data$p_value < 0.05, na.rm = TRUE),
      FPR = sum(.data$p_value < 0.05, na.rm = TRUE) / dplyr::n(),
      .groups = "drop"
    )

  set.seed(representative_seed)
  slope_vec_rep <- rep(0, n_taxa)
  if (representative_n_da > 0) {
    da_indices_rep <- sample(1:n_taxa, size = min(representative_n_da, n_taxa), replace = FALSE)
    slope_vec_rep[da_indices_rep] <- sample(slope_effects_distribution, length(da_indices_rep), replace = TRUE)
    slope_vec_rep <- slope_vec_rep - mean(slope_vec_rep)
  }

  sim_rep <- simulate_case_local(slope_vec_rep, mu_base, k_fixed, representative_seed)
  sim_rep$count_long <- filter_by_abundance_cases(sim_rep$count_long, threshold = abundance_threshold)

  list(
    fpr_data = fpr_data,
    representative = sim_rep
  )
}

# Extract seed and n_da_taxa from job_row (shared helper)
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

# Run simulation only (no analysis)
run_simulation_job <- function(
  job_row,
  intercept_dispersion,
  mu_base,
  k_fixed,
  slope_effects_distribution,
  sd_log_overdispersion,
  n_taxa,
  n_samples,
  n_samples_per_group,
  library_size_mean,
  library_size_sd
) {
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
    mu_inv_softmax = mu_base,
    log_dispersion_assoc = k_fixed,
    n_taxa = n_taxa,
    n_samples = n_samples,
    sd_log_overdispersion = sd_log_overdispersion,
    intercept_dispersion = intercept_dispersion,
    library_size_mean = library_size_mean,
    library_size_sd = library_size_sd,
    design_matrix = design_matrix_from_groups(n_samples_per_group),
    seed = seed
  )
  
  # Add seed and n_da_taxa to simulation for tracking
  sim_r$seed <- seed
  sim_r$n_da_taxa <- n_da
  
  sim_r
}

# Run analysis on a simulation (betadisper) - SINGLE SOURCE OF TRUTH for analysis
# This is the only function that performs betadisper analysis
# Uses adaptive permutations for consistent, accurate p-values across all analyses
run_analysis_job <- function(
  sim_r,
  max_perm = 19999,
  batch = 2000,
  target_se = 0.005,
  return_full = FALSE,
  standard_distance_methods = c("bray", "aitchison", "robust.aitchison", "jaccard", "hellinger")
) {
  group_levels <- c("non-IBD", "IBD")
  
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

  inputs_std <- purrr::set_names(
    lapply(standard_distance_methods, \(dm) {
      build_standard_betadisper_inputs_paper(
        count_long = sim_filtered$count_long,
        sample_metadata = sim_filtered$sample_metadata,
        distance_method = dm,
        group_levels = c("non-IBD", "IBD")
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


# Basic utility functions
## Utilities: proportion + variance-stabilizing transforms for compositional counts
##
## These helpers are intentionally dependency-free (base R only) so they can be
## sourced from scripts and Quarto documents.
##
## Conventions used in this repo:
## - X_ij (count) with N_i (library_size)
## - p_hat = X_ij / N_i
## - arcsin-sqrt transform: asin(sqrt(p))

.clamp01 <- function(x, eps = 0) {
  x <- pmax(x, 0 + eps)
  x <- pmin(x, 1 - eps)
  x
}

#' Convert counts to proportions p_hat = count / library_size
#'
#' @param count numeric vector (>= 0)
#' @param library_size numeric vector (> 0), recycled to length(count) if needed
#' @param eps numeric scalar. If > 0, proportions are clamped to [eps, 1-eps].
#'
#' @return numeric vector of proportions in [0, 1] (or [eps, 1-eps] if eps > 0)
count_to_proportion <- function(count, library_size, eps = 0) {
  if (any(library_size <= 0, na.rm = TRUE)) {
    stop("library_size must be > 0 for all entries.")
  }
  p_hat <- count / library_size
  if (eps > 0) p_hat <- .clamp01(p_hat, eps = eps)
  p_hat
}

#' Arcsin-sqrt variance-stabilizing transform for proportions
#'
#' @param p numeric vector of proportions in [0, 1]
#' @param eps numeric scalar. If > 0, p is clamped to [eps, 1-eps] before transform.
#'
#' @return numeric vector asin(sqrt(p))
arcsin_sqrt <- function(p, eps = 0) {
  if (eps > 0) p <- .clamp01(p, eps = eps)
  asin(sqrt(p))
}

#' Residual (generic): subtraction on any scale
#'
#' This is intentionally transformation-agnostic. Whatever transformation you want
#' (e.g., proportions, arcsin-sqrt, CLR, etc.) should be applied BEFORE calling this.
#'
#' @param observed numeric vector (already on the chosen scale)
#' @param expected numeric vector (already on the same scale as observed)
#'
#' @return observed - expected
residual_subtract <- function(observed, expected) {
  observed - expected
}

#' Translate residuals back to a group/taxon expected location
#'
#' This is the "add the centroid back" step:
#' \deqn{y = r + c}
#'
#' where:
#' - \eqn{r} is a residual on some scale (e.g. arcsin–sqrt),
#' - \eqn{c} is the expected abundance/location on the *same scale*.
#'
#' Note: the expected abundances should be transformed externally with the same
#' transform used for the observed proportions (e.g. `arcsin_sqrt()`), so this
#' function remains transformation-agnostic.
#'
#' @param residual numeric vector/matrix of residuals.
#' @param expected_location numeric vector/matrix of expected locations on the same scale.
#'   Must be conformable with `residual` (same length or same dimensions).
#'
#' @return translated values `residual + expected_location`.
translate_residual_location <- function(residual, expected_location) {
  if (!is.numeric(residual)) stop("residual must be numeric.")
  if (!is.numeric(expected_location)) stop("expected_location must be numeric.")
  residual + expected_location
}

#' Convert dispersion (sigma) to beta-binomial ICC parameter rho
#'
#' In this repo, `sigma` denotes dispersion (higher = more overdispersion).
#' A convenient reparameterization is:
#' \deqn{\rho = \sigma / (1 + \sigma)}
#'
#' @param sigma numeric vector/matrix of dispersion values (> 0)
#' @return rho in (0, 1)
sigma_to_rho <- function(sigma) {
  if (!is.numeric(sigma)) stop("sigma must be numeric.")
  if (any(sigma <= 0, na.rm = TRUE)) stop("sigma must be > 0.")
  # Numerically stable for very large sigma (including Inf):
  # sigma/(1+sigma) can become Inf/Inf -> NaN when sigma is Inf.
  rho <- sigma / (1 + sigma)
  rho[is.infinite(sigma)] <- 1

  # Clamp away from exact boundaries to keep downstream variance formulas well-defined.
  eps <- 1e-12
  rho <- pmax(pmin(rho, 1 - eps), eps)
  rho
}

#' Beta-binomial variance of the sample proportion p_hat = X/N
#'
#' Using rho parameterization:
#' \deqn{Var(p_hat) = pi(1-pi) * (1 + (N-1)rho) / N}
#'
#' @param pi numeric vector in [0,1]
#' @param library_size numeric vector of N (>0)
#' @param rho numeric vector in (0,1)
#' @param eps clamp for pi away from 0/1
#'
#' @return numeric vector variance of p_hat
bb_var_p_hat <- function(pi, library_size, rho, eps = 1e-12) {
  if (!is.numeric(pi)) stop("pi must be numeric.")
  if (!is.numeric(library_size)) stop("library_size must be numeric.")
  if (!is.numeric(rho)) stop("rho must be numeric.")
  if (any(library_size <= 0, na.rm = TRUE)) stop("library_size must be > 0.")
  if (any(rho <= 0 | rho >= 1, na.rm = TRUE)) stop("rho must be in (0, 1).")

  pi <- .clamp01(pi, eps = eps)
  pi * (1 - pi) * (1 + (library_size - 1) * rho) / library_size
}

#' Posterior mean of the underlying proportion p under a beta-binomial model
#'
#' Model:
#' - p ~ Beta(alpha, beta)
#' - X | p ~ Binomial(N, p)
#'
#' Using the common beta-binomial ICC parameterization:
#' \deqn{\rho = 1 / (\alpha + \beta + 1)}
#'
#' Given mean \eqn{\pi} and \eqn{\rho} we have:
#' \deqn{\phi = \alpha + \beta = 1/\rho - 1,\ \alpha = \pi\phi,\ \beta = (1-\pi)\phi}
#'
#' Posterior:
#' \deqn{p | X \sim Beta(\alpha + X, \beta + N - X)}
#' Posterior mean:
#' \deqn{E[p|X] = (\alpha + X) / (\alpha + \beta + N)}
#'
#' This provides principled shrinkage of extreme observed proportions (0/1)
#' toward \eqn{\pi}, with strength implied by \eqn{\rho}.
#'
#' @param count numeric vector of counts X (0..N)
#' @param library_size numeric vector of N (>0)
#' @param pi numeric vector of prior mean in [0,1]
#' @param rho numeric vector of ICC in (0,1)
#' @param eps clamp for pi and rho away from boundaries
#'
#' @return numeric vector posterior mean in (0,1)
bb_posterior_mean_p <- function(count, library_size, pi, rho, eps = 1e-8) {
  if (!is.numeric(count)) stop("count must be numeric.")
  if (!is.numeric(library_size)) stop("library_size must be numeric.")
  if (!is.numeric(pi)) stop("pi must be numeric.")
  if (!is.numeric(rho)) stop("rho must be numeric.")
  if (any(library_size <= 0, na.rm = TRUE)) stop("library_size must be > 0.")

  # Clamp away from exact boundaries for numerical stability
  rho <- pmax(pmin(rho, 1 - eps), eps)
  pi <- .clamp01(pi, eps = eps)

  # Guard counts within [0, N]
  count <- pmax(0, pmin(count, library_size))

  phi <- (1 / rho) - 1
  alpha <- pi * phi
  beta <- (1 - pi) * phi

  (alpha + count) / (alpha + beta + library_size)
}

#' Posterior mean of arcsin-sqrt transformed proportion under Beta posterior (delta approx)
#'
#' Instead of shrinking on the proportion scale (p), we approximate shrinkage directly on
#' the arcsin-sqrt scale by approximating \(E[g(p)\mid X]\) where \(g(p)=\arcsin\sqrt{p}\).
#' This introduces explicit **count-dependence** (via the Beta posterior) so that
#' low-information observations (e.g. very small X, or boundary outcomes X∈{0,N})
#' are regularized more strongly.
#'
#' Approximation (second-order delta method around posterior mean m):
#' \deqn{E[g(p)\mid X] \approx g(m) + \tfrac{1}{2} g''(m)\,Var(p\mid X)}
#'
#' @param count numeric vector of counts X (0..N)
#' @param library_size numeric vector of N (>0)
#' @param pi numeric vector prior mean in [0,1]
#' @param rho numeric vector ICC in (0,1)
#' @param eps clamp for numerical stability
#'
#' @return numeric vector approx \(E[\arcsin(\sqrt{p})\mid X]\)
bb_posterior_mean_arcsin_sqrt <- function(count, library_size, pi, rho, eps = 1e-8) {
  if (!is.numeric(count)) stop("count must be numeric.")
  if (!is.numeric(library_size)) stop("library_size must be numeric.")
  if (!is.numeric(pi)) stop("pi must be numeric.")
  if (!is.numeric(rho)) stop("rho must be numeric.")
  if (any(library_size <= 0, na.rm = TRUE)) stop("library_size must be > 0.")

  rho <- pmax(pmin(rho, 1 - eps), eps)
  pi <- .clamp01(pi, eps = eps)
  count <- pmax(0, pmin(count, library_size))

  phi <- (1 / rho) - 1
  alpha <- pi * phi
  beta <- (1 - pi) * phi

  a_post <- alpha + count
  b_post <- beta + (library_size - count)
  ab_post <- a_post + b_post

  m <- a_post / ab_post
  m <- .clamp01(m, eps = eps)

  # Var(p | X) for Beta(a,b): ab / ((a+b)^2 (a+b+1))
  v <- (a_post * b_post) / ((ab_post^2) * (ab_post + 1))

  g <- asin(sqrt(m))
  # g''(p) = (2p-1) / (4*(p(1-p))^(3/2))
  gpp <- (2 * m - 1) / (4 * (m * (1 - m))^(3/2))

  g + 0.5 * gpp * v
}

#' Approximate SD on arcsin-sqrt scale (pi-aware; second-order delta method)
#'
#' First-order delta method cancels pi dependence:
#' \deqn{Var(asin(sqrt(p_hat))) \approx (1 + (N-1)rho) / (4N)}
#'
#' Near boundaries (pi close to 0/1), first-order can be inaccurate. This
#' pi-aware approximation adds a second-order term:
#' \deqn{Var(g(p_hat)) \approx g'(pi)^2 Var(p_hat) + 0.5 g''(pi)^2 Var(p_hat)^2}
#'
#' where g(p) = asin(sqrt(p)).
#'
#' @param pi numeric vector in [0,1]
#' @param library_size numeric vector of N (>0)
#' @param rho numeric vector in (0,1)
#' @param eps clamp for pi away from 0/1
#'
#' @return numeric vector SD approximation
bb_arcsin_sqrt_sd_piaware <- function(pi, library_size, rho, eps = 1e-8) {
  pi <- .clamp01(pi, eps = eps)
  v_p <- bb_var_p_hat(pi = pi, library_size = library_size, rho = rho, eps = eps)

  # g'(p) = 1/(2*sqrt(p(1-p)))
  gp <- 1 / (2 * sqrt(pi * (1 - pi)))
  # g''(p) = (2p-1) / (4*(p(1-p))^(3/2))
  gpp <- (2 * pi - 1) / (4 * (pi * (1 - pi))^(3/2))

  v1 <- (gp^2) * v_p
  v2 <- 0.5 * (gpp^2) * (v_p^2)
  sqrt(v1 + v2)
}

#' Association-adjust residuals in arcsin-sqrt space using beta-binomial SD
#'
#' This is the "correct" variance-scale adjustment when residuals are computed as:
#' \deqn{r_{ij} = \arcsin\sqrt{\hat p_{ij}} - \arcsin\sqrt{\pi_{gj}}}
#'
#' We assume a mean–dispersion association in unconstrained space:
#' \deqn{\log(\sigma_{ij}) = -k \cdot \mu_{gj} + c}
#'
#' and adjust residuals by the *relative* arcsin-sqrt SD implied by the association:
#' \deqn{r^*_{ij} = r_{ij} \cdot \frac{\mathrm{SD}_0(N_i)}{\mathrm{SD}(N_i, \rho_{ij})}}
#'
#' where \eqn{\mathrm{SD}_0(N_i)} is the SD at \eqn{\mu = 0} (i.e. \eqn{\sigma_0 = \exp(c)}),
#' so this is "rotation around the origin" (no global inflation/deflation).
#'
#' @param residual numeric residuals on arcsin-sqrt scale
#' @param mu_inv_softmax numeric unconstrained predictors (same shape as residual)
#' @param log_dispersion_assoc numeric scalar k
#' @param library_size numeric vector of N (same shape as residual, or recycled)
#' @param log_sigma_intercept optional numeric scalar c. If NULL, it will be estimated from
#'   `log_sigma_observed` as mean(log_sigma_observed + k * mu_inv_softmax).
#' @param log_sigma_observed optional numeric vector of observed log_sigma (same shape as residual).
#'   Used only if log_sigma_intercept is NULL.
#'
#' @return numeric adjusted residuals
adjust_residual_assoc_arcsin_bb <- function(
  residual,
  mu_inv_softmax,
  log_dispersion_assoc,
  library_size,
  log_sigma_intercept
) {
  if (!is.numeric(residual)) stop("residual must be numeric.")
  if (!is.numeric(mu_inv_softmax)) stop("mu_inv_softmax must be numeric.")
  if (!is.numeric(log_dispersion_assoc) || length(log_dispersion_assoc) != 1) {
    stop("log_dispersion_assoc must be a numeric scalar.")
  }
  if (!is.numeric(library_size)) stop("library_size must be numeric.")
  if (!is.numeric(log_sigma_intercept)) stop("log_sigma_intercept must be numeric.")

  # Model: log_sigma = -k * mu + c
  # Here we assume `log_sigma_intercept` (c) is KNOWN (e.g. from simulation inputs).
  # This keeps the adjustment fully algebraic (no estimation).

  # Baseline at mu = 0
  sigma0 <- exp(log_sigma_intercept)
  rho0 <- sigma_to_rho(sigma0)

  # Expected sigma at each mu under the model
  log_sigma_exp <- -log_dispersion_assoc * mu_inv_softmax + log_sigma_intercept
  sigma_exp <- exp(log_sigma_exp)
  rho_exp <- sigma_to_rho(sigma_exp)

  var0 <- (1 + (library_size - 1) * rho0) / (4 * library_size)
  var_exp <- (1 + (library_size - 1) * rho_exp) / (4 * library_size)

  residual * sqrt(var0 / var_exp)
}

#' "inv_softmax" transform: map unconstrained log-linear predictors to proportions
#'
#' In this repo's naming, "mu_inv_softmax" refers to an *unconstrained* vector
#' (log-linear predictors). Applying "inv_softmax" converts it to proportions.
#' Mathematically, this is the SOFTMAX transform.
#'
#' @param x numeric vector or matrix of unconstrained values (real-valued).
#'   - If a vector: returns a vector of the same length that sums to 1.
#'   - If a matrix: applies row-wise; each row sums to 1.
#' @param eps numeric scalar. If > 0, probabilities are clamped to [eps, 1-eps].
#'
#' @return numeric vector or matrix of probabilities in (0,1), same shape as x.
inv_softmax <- function(x, eps = 0) {
  if (!is.numeric(x)) stop("x must be numeric.")

  softmax_vec <- function(v) {
    if (anyNA(v)) stop("x contains NA; please handle missing values before softmax.")
    v <- v - max(v) # numerical stability
    ex <- exp(v)
    out <- ex / sum(ex)
    if (eps > 0) out <- .clamp01(out, eps = eps)
    out
  }

  if (is.matrix(x)) {
    out <- t(apply(x, 1, softmax_vec))
    dimnames(out) <- dimnames(x)
    return(out)
  }

  if (is.vector(x)) {
    return(softmax_vec(as.numeric(x)))
  }

  stop("x must be a numeric vector or matrix.")
}