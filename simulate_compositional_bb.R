# Simulate Compositional Beta-Binomial Dataset
# 
# This script simulates compositional count data (e.g., microbiome) where:
# - Counts follow beta-binomial distribution
# - Means are compositional (inverse softmax of linear predictors)
# - Log-dispersion (sigma) is inversely proportional to logit-mean (sccomp-style)
# - There's variability around the mean-dispersion relationship
#
# Parameters:
# - slope_vector: Vector of slopes for each taxon (length = n_taxa)
# - inv_softmax_mean: Intercept for inverse softmax mean (scalar)
# - log_dispersion_assoc: Association parameter k in log(σ) = -k·logit(μ) + c
# - n_taxa: Number of taxa/species
# - n_samples: Number of samples
# - sd_log_overdispersion: SD of log(sigma) around the regression line
#
# Dependencies:
# - ggplot2 (for visualization)
# - dplyr (for data manipulation)
# - purrr (for functional programming)
# - tibble (for tibble data structures)
# - VGAM (for beta-binomial simulation and fitting)
#   Alternative packages for beta-binomial simulation:
#   - extraDistr::rbbinom()
#   - emdbook::rbetabinom()

library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)

# ============================================================================
# Helper Functions
# ============================================================================

# Softmax: converts log-linear predictors to probabilities
# Given log-linear predictors, exponentiates and normalizes to sum to 1
# Note: This is the standard softmax operation, NOT inverse softmax
# Inverse softmax would be: probabilities -> log-linear predictors (centered log-ratio)
softmax <- function(log_linear_predictors) {
  # log_linear_predictors: matrix of size (n_samples x n_taxa) in log space
  # Returns: matrix of probabilities (n_samples x n_taxa) that sum to 1 per row
  exp_log <- exp(log_linear_predictors)  # Exponentiate (back to regular space)
  row_sums <- rowSums(exp_log)
  probabilities <- exp_log / row_sums    # Normalize so each row sums to 1
  return(probabilities)
}

# Calculate sigma from mu under constraint: log(sigma) = -k·logit(mu) + c + error
calculate_sigma_constrained <- function(mu, k, c, error = 0) {
  # mu: probability (0 < mu < 1)
  # k: slope parameter (log-dispersion association)
  # c: intercept
  # error: random error term
  logit_mu <- log(mu / (1 - mu))
  log_sigma <- -k * logit_mu + c + error
  sigma <- exp(log_sigma)
  return(sigma)
}

# Simulate beta-binomial counts
simulate_beta_binomial <- function(n, mu, sigma, n_sim = 1) {
  # n: number of trials (library size)
  # mu: mean probability
  # sigma: dispersion parameter (higher = more overdispersion)
  # n_sim: number of simulations
  
  # Convert dispersion to concentration: concentration = 1/sigma
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

# Calculate standard deviation from unconstrained mu and overdispersion
calculate_sd_from_overdispersion <- function(mu, overdispersion, n, 
                                            overdispersion_type = c("sigma", "factor", "ratio")) {
  # mu: unconstrained mean probability (0 < mu < 1)
  # overdispersion: overdispersion parameter (interpretation depends on type)
  # n: number of trials
  # overdispersion_type: type of overdispersion parameter
  #   - "sigma": dispersion parameter σ (higher = more overdispersion)
  #   - "factor": overdispersion factor (n + κ) / (κ + 1) where κ = 1/σ
  #   - "ratio": overdispersion ratio Var / Var_binomial
  
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

# Convenience wrapper functions
calculate_sd_from_sigma <- function(mu, sigma, n) {
  # Wrapper for sigma (dispersion) parameter
  # sigma: dispersion parameter (higher = more overdispersion)
  return(calculate_sd_from_overdispersion(mu, sigma, n, overdispersion_type = "sigma"))
}

calculate_sd_from_overdispersion_factor <- function(mu, overdispersion_factor, n) {
  # Wrapper for overdispersion factor
  return(calculate_sd_from_overdispersion(mu, overdispersion_factor, n, 
                                          overdispersion_type = "factor"))
}

calculate_sd_from_overdispersion_ratio <- function(mu, overdispersion_ratio, n) {
  # Wrapper for overdispersion ratio
  return(calculate_sd_from_overdispersion(mu, overdispersion_ratio, n, 
                                          overdispersion_type = "ratio"))
}

# ============================================================================
# Beta-Binomial Fitting Function
# ============================================================================

fit_beta_binomial_per_taxon <- function(count_long, library_size_col = "library_size", group_col = "group") {
  # Fit beta-binomial distribution per taxon (and per group if group_col exists)
  # and estimate mu_inv_softmax and log_sigma (dispersion)
  #
  # Parameters:
  #   count_long: data frame with columns: taxon_id, sample_id, count, and library_size
  #   library_size_col: name of the column containing library sizes
  #   group_col: name of the column containing group information (optional)
  #
  # Returns:
  #   data frame with columns: taxon_id, group (if applicable), mu_inv_softmax, log_sigma, mu, sigma
  #   where sigma represents dispersion (higher = more overdispersion)
  #
  # Dependencies:
  #   - VGAM package for beta-binomial fitting
  rlang::check_installed("VGAM", reason = "for beta-binomial fitting")
  
  # Check if group column exists
  has_group <- group_col %in% colnames(count_long)
  
  # Get unique taxa
  taxa <- unique(count_long$taxon_id)
  n_taxa <- length(taxa)
  
  # Initialize results list
  results_list <- list()
  
  if (has_group) {
    # Fit per taxon per group
    groups <- unique(count_long[[group_col]])
    groups <- groups[!is.na(groups)]
    
    for (taxon in taxa) {
      for (group_val in groups) {
        # Get data for this taxon-group combination
        taxon_group_data <- count_long %>%
          filter(taxon_id == taxon, .data[[group_col]] == group_val) %>%
          mutate(
            library_size = .data[[library_size_col]],
            success = count,
            failure = library_size - count
          ) %>%
          # Ensure counts don't exceed library_size (can happen due to rounding)
          mutate(
            success = pmin(success, library_size),
            failure = pmax(0, library_size - success)
          ) %>%
          # Filter out invalid rows
          filter(library_size > 0 & success >= 0 & failure >= 0)
        
        # Skip if insufficient data
        if (nrow(taxon_group_data) < 3) {
          warning(paste("Insufficient data for taxon", taxon, "group", group_val, "- skipping"))
          next
        }
        
        # Fit beta-binomial model using VGAM
        tryCatch({
          # Create response matrix: cbind(success, failure)
          # VGAM automatically uses rowSums(y_matrix) = success + failure = library_size
          # as the total number of trials for each sample
          # This handles varying library sizes automatically
          y_matrix <- cbind(taxon_group_data$success, taxon_group_data$failure)
          
          # Additional validation: ensure no negative values
          if (any(y_matrix < 0)) {
            stop("Negative values in response matrix")
          }
          
          # Set initial values based on prior and data
          # Prior: log_sigma = -2 → sigma ≈ 0.135
          # Since sigma = rho / (1 - rho), we get: rho = sigma / (1 + sigma)
          prior_sigma <- exp(-2)  # sigma ≈ 0.135
          prior_rho <- prior_sigma / (1 + prior_sigma)  # rho ≈ 0.119
          
          # For mu, use empirical proportion
          empirical_mu <- mean(taxon_group_data$success / taxon_group_data$library_size)
          empirical_mu <- pmax(0.001, pmin(0.999, empirical_mu))  # Bound away from 0 and 1
          
          # VGAM betabinomial uses logit link for both parameters
          init_mu_logit <- qlogis(empirical_mu)
          init_rho_logit <- qlogis(prior_rho)
          
          # Fit with initial values and increased iterations
          fit <- VGAM::vglm(
            y_matrix ~ 1,
            family = VGAM::betabinomial(zero = NULL),
            coefstart = c(init_mu_logit, init_rho_logit),  # Initial values
            control = VGAM::vglm.control(
              maxit = 100,      # Increase from default (usually 10-20)
              epsilon = 1e-7    # Tighter convergence criterion
            )
          )
          
          coefs_link <- VGAM::coef(fit)
          
          if (length(coefs_link) < 2) {
            stop(paste("Expected 2 coefficients for taxon", taxon, "group", group_val, "but got", length(coefs_link)))
          }
          
          mu_logit <- coefs_link[1]
          mu <- plogis(mu_logit)
          
          rho_logit <- coefs_link[2]
          rho <- plogis(rho_logit)
          
          # Convert VGAM's rho to our dispersion parameter
          # VGAM's beta-binomial uses: rho = 1 / (alpha + beta + 1)
          # where alpha = mu * κ, beta = (1-mu) * κ, and κ is concentration
          # So: alpha + beta = κ (concentration)
          # Therefore: rho = 1 / (κ + 1)  =>  κ = (1 - rho) / rho
          # In our notation, dispersion σ = 1/κ (higher σ = more overdispersion)
          # So: σ = 1/κ = rho / (1 - rho)
          
          # Apply soft prior regularization on dispersion estimates
          # Instead of hard cutoffs, we use shrinkage toward a prior distribution
          # 
          # IMPORTANT: sigma = DISPERSION (not concentration!)
          #   - Higher sigma → More overdispersion → More variability
          #   - concentration = 1/sigma (inverse relationship)
          # 
          # Prior: log_sigma ~ N(prior_mean, prior_sd)
          # - prior_mean = -2 
          #   → sigma ≈ 0.14 (LOW dispersion, counts close to means)
          #   → concentration ≈ 7.4 (HIGH concentration)
          # - prior_sd = 3 (allows wide range: log_sigma roughly in [-8, 4])
          #
          # Shrinkage formula: log_sigma_regularized = w * log_sigma_mle + (1-w) * prior_mean
          # where w = 1 / (1 + lambda), and lambda increases with distance from prior
          
          if (is.na(rho) || rho <= 0 || rho >= 1 || is.nan(rho)) {
            warning(paste("Invalid rho for taxon", taxon, "group", group_val, ":", rho, "- using prior mean"))
            sigma <- exp(-2)  # Prior mean: log_sigma = -2 → sigma = 0.14 (LOW dispersion)
          } else {
            # Wide bounds on rho to avoid numerical issues: [0.0001, 0.9999]
            rho_bounded <- pmin(pmax(rho, 0.0001), 0.9999)
            
            # Convert to concentration: κ = (1 - rho) / rho
            concentration <- (1 - rho_bounded) / rho_bounded
            
            # Convert to dispersion: σ = 1/κ
            sigma_mle <- 1 / concentration
            
            # Compute log_sigma from MLE
            log_sigma_mle <- log(sigma_mle)
            
            # Apply soft prior regularization with shrinkage
            prior_mean <- -2    # Moderate low overdispersion
            prior_sd <- 3       # Wide prior (allows roughly [-8, 4] range)
            
            # Shrinkage weight: stronger shrinkage for extreme values
            # Distance from prior in standard deviations
            z_score <- abs(log_sigma_mle - prior_mean) / prior_sd
            
            # Weight decreases as z_score increases (more shrinkage for extreme values)
            # Using exponential decay: w = exp(-lambda * z_score^2)
            # lambda = 0.1 gives gentle shrinkage
            lambda <- 0.1
            weight <- exp(-lambda * z_score^2)
            
            # Regularized log_sigma: weighted average of MLE and prior
            log_sigma_regularized <- weight * log_sigma_mle + (1 - weight) * prior_mean
            
            # Convert back to sigma
            sigma <- exp(log_sigma_regularized)
            
            # Wide safety bounds to catch numerical issues only
            if (sigma <= exp(-15) || sigma >= exp(15) || is.infinite(sigma) || is.nan(sigma)) {
              warning(paste("Extreme dispersion for taxon", taxon, "group", group_val, ":", sigma, "- using prior mean"))
              sigma <- exp(-2)
            }
          }
          
          results_list[[length(results_list) + 1]] <- data.frame(
            taxon_id = taxon,
            group = as.character(group_val),
            mu_logit = mu_logit,
            mu = mu,
            log_sigma = log(sigma),
            sigma = sigma,
            stringsAsFactors = FALSE
          )
          
        }, error = function(e) {
          warning(paste("Error fitting taxon", taxon, "group", group_val, ":", e$message))
        })
      }
    }
    
    # Combine results
    if (length(results_list) == 0) {
      return(data.frame(
        taxon_id = character(0),
        group = character(0),
        mu_inv_softmax = numeric(0),
        log_sigma = numeric(0),
        mu = numeric(0),
        sigma = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    
    results <- do.call(rbind, results_list)
    
    # Convert estimated proportions to centered log-linear predictors per group
    # For each group, center the log-proportions so they sum to 0
    for (group_val in groups) {
      group_idx <- results$group == group_val
      if (sum(group_idx) > 0) {
        valid_idx <- group_idx & !is.na(results$mu) & results$mu > 0 & results$mu < 1
        if (sum(valid_idx) > 0) {
          log_props <- log(results$mu[valid_idx])
          mean_log_prop <- mean(log_props)
          results$mu_inv_softmax[valid_idx] <- log_props - mean_log_prop
        }
      }
    }
    
    # Remove temporary mu_logit column
    results$mu_logit <- NULL
    
  } else {
    # Fit per taxon only (no groups)
    results <- data.frame(
      taxon_id = character(n_taxa),
      mu_inv_softmax = numeric(n_taxa),
      log_sigma = numeric(n_taxa),
      mu = numeric(n_taxa),
      sigma = numeric(n_taxa),
      stringsAsFactors = FALSE
    )
    
    for (i in 1:n_taxa) {
      taxon <- taxa[i]
      
      taxon_data <- count_long %>%
        filter(taxon_id == taxon) %>%
        mutate(
          library_size = .data[[library_size_col]],
          success = count,
          failure = library_size - count
        ) %>%
        # Ensure counts don't exceed library_size
        mutate(
          success = pmin(success, library_size),
          failure = pmax(0, library_size - success)
        ) %>%
        # Filter out invalid rows
        filter(library_size > 0 & success >= 0 & failure >= 0)
      
      if (nrow(taxon_data) < 3) {
        warning(paste("Insufficient data for taxon", taxon, "- skipping"))
        results$taxon_id[i] <- taxon
        results$mu_inv_softmax[i] <- NA
        results$log_sigma[i] <- NA
        results$mu[i] <- NA
        results$sigma[i] <- NA
        next
      }
      
      tryCatch({
        # Create response matrix: cbind(success, failure)
        # VGAM automatically uses rowSums(y_matrix) = success + failure = library_size
        # as the total number of trials for each sample
        # This handles varying library sizes automatically
        y_matrix <- cbind(taxon_data$success, taxon_data$failure)
        
        # Additional validation: ensure no negative values
        if (any(y_matrix < 0)) {
          stop("Negative values in response matrix")
        }
        
        fit <- VGAM::vglm(
          y_matrix ~ 1,
          family = VGAM::betabinomial(zero = NULL)
        )
        
        coefs_link <- VGAM::coef(fit)
        
        if (length(coefs_link) < 2) {
          stop(paste("Expected 2 coefficients for taxon", taxon, "but got", length(coefs_link)))
        }
        
        mu_logit <- coefs_link[1]
        mu <- plogis(mu_logit)
        
        rho_logit <- coefs_link[2]
        rho <- plogis(rho_logit)
        
        # Convert VGAM's rho to our dispersion parameter
        # VGAM's beta-binomial uses: rho = 1 / (alpha + beta + 1)
        # where alpha = mu * κ, beta = (1-mu) * κ, and κ is concentration
        # So: alpha + beta = κ (concentration)
        # Therefore: rho = 1 / (κ + 1)  =>  κ = (1 - rho) / rho
        # In our notation, dispersion σ = 1/κ (higher σ = more overdispersion)
        # So: σ = 1/κ = rho / (1 - rho)
        
        # Apply reasonable bounds to rho to avoid extreme dispersion estimates
        # See grouped case above for detailed explanation of bounds
        # Keeping log_sigma in [-5, 5] range for interpretability
        
        if (is.na(rho) || rho <= 0 || rho >= 1 || is.nan(rho)) {
          warning(paste("Invalid rho for taxon", taxon, ":", rho, "- setting dispersion to default"))
          sigma <- 0.1  # Default dispersion
        } else {
          # Bound rho to reasonable range: [0.005, 0.995]
          rho_bounded <- pmin(pmax(rho, 0.005), 0.995)
          
          # Convert to concentration first: κ = (1 - rho) / rho
          concentration <- (1 - rho_bounded) / rho_bounded
          
          # Then to dispersion: σ = 1/κ
          sigma <- 1 / concentration
          
          # Apply strict bounds on sigma to ensure log_sigma stays in [-5, 5]
          sigma <- pmin(pmax(sigma, exp(-5)), exp(5))
          
          if (sigma <= 0 || is.infinite(sigma) || is.nan(sigma)) {
            warning(paste("Invalid dispersion for taxon", taxon, ":", sigma, "- setting to default"))
            sigma <- 0.1
          }
        }
        
        results$taxon_id[i] <- taxon
        results$mu_inv_softmax[i] <- mu_logit  # Temporary: will be centered later
        results$log_sigma[i] <- log(sigma)
        results$mu[i] <- mu
        results$sigma[i] <- sigma
        
      }, error = function(e) {
        warning(paste("Error fitting taxon", taxon, ":", e$message))
        results$taxon_id[i] <- taxon
        results$mu_inv_softmax[i] <- NA
        results$log_sigma[i] <- NA
        results$mu[i] <- NA
        results$sigma[i] <- NA
      })
    }
    
    # Convert estimated proportions to centered log-linear predictors
    valid_idx <- !is.na(results$mu) & results$mu > 0 & results$mu < 1
    if (sum(valid_idx) > 0) {
      log_props <- log(results$mu[valid_idx])
      mean_log_prop <- mean(log_props)
      results$mu_inv_softmax[valid_idx] <- log_props - mean_log_prop
    }
  }
  
  return(results)
}

# ============================================================================
# Main Simulation Function
# ============================================================================

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

# Plot mean-dispersion relationship
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

# Plot compositional abundances
plot_compositional_abundances <- function(sim_result, n_samples_plot = 20) {
  # Calculate proportions
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

# Plot sigma distribution by taxon
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
# Invariant beta-dispersion helpers (report utilities)
# ============================================================================

# Build the association-adjusted distance matrix and aligned group factor.
#
# This encapsulates the report pipeline:
# 1) compute transformed residuals (arcsin-sqrt by default),
# 2) scale residuals by the mean–dispersion association (slope-only),
# 3) translate residuals by group×taxon expected location (centroid),
# 4) return (distance, group_factor) for vegan::betadisper().
#
# Returns a list of length 2: $dist, $group
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

# Build pi-aware association-adjusted arcsin-sqrt distance + aligned group factor.
#
# Option A: pi-aware diagonal scaling (second-order delta method).
# Residuals are in arcsin-sqrt space; scaling uses both sigma(mu) and pi via
# a pi-aware SD approximation. Centroid is then added back (same as assoc_adj).
#
# Returns: list(dist, group)
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

# Build prevalence-filtered version of postmeanT (drops ultra-rare taxa).
#
# Returns: list(dist, group)
build_assoc_adjusted_piaware_postmeanT_prevfilter_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  eps_rho = 1e-8,
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

  build_assoc_adjusted_piaware_postmeanT_betadisper_inputs(
    count_long = count_long_f,
    log_dispersion_assoc = log_dispersion_assoc,
    log_sigma_intercept = log_sigma_intercept,
    eps_pi = eps_pi,
    eps_rho = eps_rho,
    group_levels = group_levels
  )
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
#
# Filter-only: drops taxa with low prevalence, without altering geometry among retained taxa.
# Returns: list(dist, group)
build_assoc_adjusted_piaware_postmean_prevfilter_betadisper_inputs <- function(
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

  build_assoc_adjusted_piaware_postmean_betadisper_inputs(
    count_long = count_long_f,
    log_dispersion_assoc = log_dispersion_assoc,
    log_sigma_intercept = log_sigma_intercept,
    eps_pi = eps_pi,
    group_levels = group_levels
  )
}

# Build covariance-whitened association-adjusted distance + aligned group factor.
#
# Option B: take Option A features and apply a pooled covariance whitening transform
# (Mahalanobis-like), then compute Euclidean distances in the whitened space.
#
# Returns: list(dist, group)
build_assoc_adjusted_whiten_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  ridge = 1e-6,
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

  # Pooled covariance (centered)
  mat_centered <- scale(mat, center = TRUE, scale = FALSE)
  Sigma <- stats::cov(mat_centered)
  if (any(!is.finite(Sigma))) stop("Covariance matrix has non-finite entries.")

  # Ridge for numerical stability
  lambda <- ridge * mean(diag(Sigma))
  Sigma_reg <- Sigma + diag(lambda, nrow(Sigma))

  # Whitening transform via Cholesky: Sigma_reg = R^T R
  R <- chol(Sigma_reg)
  W <- backsolve(R, diag(nrow(Sigma_reg))) # W = R^{-1}
  mat_white <- mat_centered %*% W

  dist_obj <- stats::dist(mat_white, method = "euclidean")

  group_factor <- df %>%
    dplyr::distinct(sample_id, group) %>%
    dplyr::arrange(match(sample_id, rownames(mat))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

# Build pooled-prevalence-weighted association-adjusted arcsin-sqrt distance (non-π-aware).
#
# Like `assoc_adj_prevweight`, but weights are computed from pooled prevalence across
# all samples (ignoring group), to reduce interaction with asymmetric composition.
#
# Returns: list(dist, group)
build_assoc_adjusted_prevweight_pooled_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  min_prevalence = 0.05,
  weight_power = 0.5,
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
  if (!is.numeric(weight_power) || length(weight_power) != 1 || weight_power <= 0) {
    stop("weight_power must be a positive numeric scalar.")
  }

  # Build base assoc_adj matrix (already anchored/scaled, with centroid translation)
  inputs <- build_assoc_adjusted_betadisper_inputs(
    count_long = count_long,
    log_dispersion_assoc = log_dispersion_assoc,
    log_sigma_intercept = log_sigma_intercept,
    group_levels = group_levels
  )

  # Rebuild the underlying matrix in the same way as assoc_adj (dist doesn't keep it)
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

  # Pooled prevalence by taxon (ignore group)
  prev_tbl <- count_long %>%
    dplyr::mutate(present = count > 0) %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::summarise(prevalence = mean(present), .groups = "drop")

  prev_vec <- prev_tbl$prevalence
  names(prev_vec) <- prev_tbl$taxon_id

  taxa <- colnames(mat)
  if (any(!(taxa %in% names(prev_vec)))) {
    missing <- taxa[!(taxa %in% names(prev_vec))]
    stop("Prevalence missing for taxa: ", paste(missing, collapse = ", "))
  }
  prev_vec <- prev_vec[taxa]

  keep <- prev_vec >= min_prevalence
  if (!any(keep)) stop("No taxa remain after prevalence filtering. Lower min_prevalence.")

  weights <- (pmax(prev_vec, 0))^weight_power
  mat_w <- mat[, keep, drop = FALSE]
  w_keep <- weights[keep]
  mat_w <- sweep(mat_w, 2, w_keep, `*`)

  dist_obj <- stats::dist(mat_w, method = "euclidean")

  group_factor <- df %>%
    dplyr::distinct(sample_id, group) %>%
    dplyr::arrange(match(sample_id, rownames(mat_w))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

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

# Build prevalence-filtered, posterior-mean-regularized pi-aware association-adjusted distance.
#
# Hybrid:
# - Regularize p_hat via posterior mean to avoid boundary blow-ups, AND
# - Drop ultra-rare taxa (min prevalence across groups) to reduce rare-taxa leverage.
#
# Returns: list(dist, group)
build_assoc_adjusted_piaware_postmean_prevfilter_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  eps_rho = 1e-8,
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

  build_assoc_adjusted_piaware_postmean_betadisper_inputs(
    count_long = count_long_f,
    log_dispersion_assoc = log_dispersion_assoc,
    log_sigma_intercept = log_sigma_intercept,
    eps_pi = eps_pi,
    eps_rho = eps_rho,
    group_levels = group_levels
  )
}

# Build prevalence-weighted, pi-aware association-adjusted arcsin-sqrt distance.
#
# Motivation: heavy zero-inflation / ultra-rare taxa can dominate Euclidean geometry
# after transforms via boundary effects (p_hat=0), even if mean–variance adjustment
# is algebraically correct. This flavor downweights low-prevalence taxa.
#
# Steps:
# 1) Build pi-aware adjusted matrix (Option A)
# 2) Compute per-taxon prevalence (min across groups)
# 3) Apply weights to columns (and optionally drop very-low-prevalence taxa)
#
# Returns: list(dist, group)
build_assoc_adjusted_piaware_prevweight_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  min_prevalence = 0.05,
  weight_power = 0.5,
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
  if (!is.numeric(weight_power) || length(weight_power) != 1 || weight_power <= 0) {
    stop("weight_power must be a positive numeric scalar.")
  }

  # 1) Base matrix from pi-aware method (already contains centroid translation)
  base <- build_assoc_adjusted_piaware_betadisper_inputs(
    count_long = count_long,
    log_dispersion_assoc = log_dispersion_assoc,
    log_sigma_intercept = log_sigma_intercept,
    eps_pi = eps_pi,
    group_levels = group_levels
  )

  # Rebuild the underlying matrix (stats::dist doesn't keep it)
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

  # 2) Prevalence by taxon (min across groups)
  prev_tbl <- count_long %>%
    dplyr::mutate(present = count > 0) %>%
    dplyr::group_by(group, taxon_id) %>%
    dplyr::summarise(prevalence = mean(present), .groups = "drop") %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::summarise(prevalence_min = min(prevalence), .groups = "drop")

  prev_vec <- prev_tbl$prevalence_min
  names(prev_vec) <- prev_tbl$taxon_id

  # Keep only taxa present in the matrix (should match, but be defensive)
  taxa <- colnames(mat)
  if (any(!(taxa %in% names(prev_vec)))) {
    missing <- taxa[!(taxa %in% names(prev_vec))]
    stop("Prevalence missing for taxa: ", paste(missing, collapse = ", "))
  }

  prev_vec <- prev_vec[taxa]

  # 3) Weights + optional filtering
  keep <- prev_vec >= min_prevalence
  if (!any(keep)) {
    stop("No taxa remain after prevalence filtering. Lower min_prevalence.")
  }

  weights <- (pmax(prev_vec, 0))^weight_power
  weights[!keep] <- 0

  mat_w <- mat[, keep, drop = FALSE]
  w_keep <- weights[keep]

  # Column-wise scaling
  mat_w <- sweep(mat_w, 2, w_keep, `*`)

  dist_obj <- stats::dist(mat_w, method = "euclidean")

  group_factor <- df %>%
    dplyr::distinct(sample_id, group) %>%
    dplyr::arrange(match(sample_id, rownames(mat_w))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

# Build prevalence-weighted association-adjusted arcsin-sqrt distance (non-π-aware).
#
# This is the same as `assoc_adj` (anchored variance ratio in arcsin-sqrt space),
# but with a prevalence-based downweighting/filter to reduce the leverage of
# ultra-rare / zero-heavy taxa that can dominate distances via boundary effects.
#
# Returns: list(dist, group)
build_assoc_adjusted_prevweight_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  min_prevalence = 0.05,
  weight_power = 0.5,
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
  if (!is.numeric(weight_power) || length(weight_power) != 1 || weight_power <= 0) {
    stop("weight_power must be a positive numeric scalar.")
  }

  # Rebuild assoc_adj matrix (stats::dist doesn't retain it)
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

  # Prevalence by taxon (min across groups)
  prev_tbl <- count_long %>%
    dplyr::mutate(present = count > 0) %>%
    dplyr::group_by(group, taxon_id) %>%
    dplyr::summarise(prevalence = mean(present), .groups = "drop") %>%
    dplyr::group_by(taxon_id) %>%
    dplyr::summarise(prevalence_min = min(prevalence), .groups = "drop")

  prev_vec <- prev_tbl$prevalence_min
  names(prev_vec) <- prev_tbl$taxon_id

  taxa <- colnames(mat)
  if (any(!(taxa %in% names(prev_vec)))) {
    missing <- taxa[!(taxa %in% names(prev_vec))]
    stop("Prevalence missing for taxa: ", paste(missing, collapse = ", "))
  }
  prev_vec <- prev_vec[taxa]

  keep <- prev_vec >= min_prevalence
  if (!any(keep)) stop("No taxa remain after prevalence filtering. Lower min_prevalence.")

  weights <- (pmax(prev_vec, 0))^weight_power
  mat_w <- mat[, keep, drop = FALSE]
  w_keep <- weights[keep]
  mat_w <- sweep(mat_w, 2, w_keep, `*`)

  dist_obj <- stats::dist(mat_w, method = "euclidean")

  group_factor <- df %>%
    dplyr::distinct(sample_id, group) %>%
    dplyr::arrange(match(sample_id, rownames(mat_w))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

# Build "inv_softmax-adjusted" distance + aligned group factor for betadisper.
#
# This is an alternative adjusted method where residuals are computed on the
# unconstrained (log-linear predictor) scale by applying an inverse-softmax-like
# transform to observed proportions:
#
#   p_hat_ij = X_ij / N_i
#   mu_hat_ij = log(p_hat_ij) - mean_j log(p_hat_ij)   (center per sample)
#   residual_ij = mu_hat_ij - mu_exp_ij
#
# Then apply the association adjustment (slope-only) and translate back to the
# expected location on the same scale:
#
#   residual_adj = residual / sqrt(exp(-k * mu_exp))
#   y_ij = residual_adj + centroid_gj   where centroid_gj = mean(mu_exp_ij) over group×taxon
#
# Returns a list of length 2: $dist, $group
build_inv_softmax_adjusted_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  eps = 1e-12,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size",
    "unconstrained_log_mu"
  )
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(eps) || length(eps) != 1 || eps <= 0) {
    stop("eps must be a positive numeric scalar.")
  }

  df <- count_long %>%
    dplyr::mutate(
      p_hat = count_to_proportion(count = count, library_size = library_size, eps = 0),
      p_hat_clamped = pmax(p_hat, eps),
      log_p = log(p_hat_clamped),
      mu_hat = log_p - mean(log_p),
      .by = sample_id
    ) %>%
    dplyr::mutate(
      residual = residual_subtract(mu_hat, unconstrained_log_mu),
      assoc_adj_residual = adjust_residual_assoc(
        residual = residual,
        mu_inv_softmax = unconstrained_log_mu,
        log_dispersion_assoc = log_dispersion_assoc
      ),
      centroid = mean(unconstrained_log_mu),
      .by = c(group, taxon_id)
    ) %>%
    dplyr::mutate(
      inv_softmax_adj_with_location = translate_residual_location(
        residual = assoc_adj_residual,
        expected_location = centroid
      )
    )

  mat <- df %>%
    dplyr::select(sample_id, taxon_id, inv_softmax_adj_with_location) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = inv_softmax_adj_with_location) %>%
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

# Build standard Bray–Curtis distance + aligned group factor for betadisper.
#
# Returns a list of length 2: $dist, $group
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

# Build arcsin-sqrt residual distance + aligned group factor for betadisper.
#
# This matches the report's "Approach 1: Arcsin(sqrt) residuals":
#   r_ij = asin(sqrt(p_hat_ij)) - asin(sqrt(mu_ij))
# and then uses Euclidean distance on the residual matrix.
#
# Returns a list of length 2: $dist, $group
build_arcsin_residual_betadisper_inputs <- function(
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

  mat <- df %>%
    dplyr::select(sample_id, taxon_id, resid) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = resid) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  dist_obj <- stats::dist(mat, method = "euclidean")

  group_factor <- sample_metadata %>%
    dplyr::filter(sample_id %in% rownames(mat)) %>%
    dplyr::arrange(match(sample_id, rownames(mat))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

# Build arcsin-sqrt residuals, then add back group×taxon centroid location.
#
# Motivation (paper figures):
# - Pure residuals (obs - expected) can visually overlap across groups in PCoA,
#   because the group-level location is removed by construction.
# - Adding back a constant centroid term per group×taxon shifts group centroids
#   without changing within-group dispersion in Euclidean space.
#
# Returns a list of length 2: $dist, $group
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

# Build Pearson residual distance + aligned group factor for betadisper.
#
# This matches the report's Pearson residual definition:
#   r_ij = (X_ij - N_i * mu_ij) / sqrt( N_i * mu_ij * (1-mu_ij) * (1 + (N_i-1)*rho_ij) )
# with rho_ij = sigma_ij / (1 + sigma_ij)
#
# Returns a list of length 2: $dist, $group
build_pearson_residual_betadisper_inputs <- function(
  count_long,
  sample_metadata,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c("sample_id", "taxon_id", "count", "library_size", "mu", "sigma")
  missing_cols <- setdiff(required_cols, colnames(count_long))
  if (length(missing_cols) > 0) {
    stop("count_long is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!("sample_id" %in% colnames(sample_metadata)) || !("group" %in% colnames(sample_metadata))) {
    stop("sample_metadata must contain columns: sample_id, group")
  }

  df <- count_long %>%
    dplyr::mutate(
      expected_count = library_size * mu,
      rho = sigma / (1 + sigma),
      expected_var = library_size * mu * (1 - mu) * (1 + (library_size - 1) * rho),
      resid = (count - expected_count) / sqrt(expected_var)
    )

  mat <- df %>%
    dplyr::select(sample_id, taxon_id, resid) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = resid) %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  dist_obj <- stats::dist(mat, method = "euclidean")

  group_factor <- sample_metadata %>%
    dplyr::filter(sample_id %in% rownames(mat)) %>%
    dplyr::arrange(match(sample_id, rownames(mat))) %>%
    dplyr::pull(group) %>%
    factor(levels = group_levels)

  list(dist = dist_obj, group = group_factor)
}

# Alpha diversity helper: richness, Shannon, Pielou evenness + plot
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

# =============================================================================
# IBD Asinh variant: uses asinh(p * scale) transformation
# Variance-stabilizing, handles full [0,1] range smoothly
# =============================================================================

# Asinh transformation: asinh(p * scale) where scale controls sensitivity
asinh_transform <- function(p, scale = 10, eps = 0) {
  p <- pmax(eps, pmin(1 - eps, p))
  asinh(p * scale)
}

# SD of asinh-transformed beta-binomial (delta method approximation)
# Var(asinh(p*s)) ≈ (d/dp asinh(p*s))^2 * Var(p)
# d/dp asinh(p*s) = s / sqrt(1 + (p*s)^2)
bb_asinh_sd_piaware <- function(pi, library_size, rho, scale = 10, eps = 1e-8) {
  pi <- pmax(eps, pmin(1 - eps, pi))
  rho <- pmax(eps, pmin(1 - eps, rho))
  n <- library_size
  
  # Beta-binomial variance of proportion
  var_p <- pi * (1 - pi) * (1 + (n - 1) * rho) / n
  
  # Derivative of asinh(p*s) w.r.t. p
  deriv <- scale / sqrt(1 + (pi * scale)^2)
  
  # Delta method: Var(f(p)) ≈ (f'(p))^2 * Var(p)
  var_asinh <- deriv^2 * var_p
  
  sqrt(pmax(var_asinh, eps))
}

# Build pi-aware association-adjusted Asinh distance with high-prevalence filtering
build_assoc_adjusted_piaware_asinh_highprev_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  scale = 10,
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
      obs_t = asinh_transform(p_hat, scale = scale, eps = 0),
      # No residuals - use raw transformed value with SD adjustment

      c = .intercept_by_group(as.character(group)),
      sigma0 = exp(c),
      rho0 = sigma_to_rho(sigma0),

      sigma_exp = exp(-log_dispersion_assoc * unconstrained_log_mu + c),
      rho_exp = sigma_to_rho(sigma_exp),

      sd0 = bb_asinh_sd_piaware(pi = mu, library_size = library_size, rho = rho0, scale = scale, eps = eps_pi),
      sd_exp = bb_asinh_sd_piaware(pi = mu, library_size = library_size, rho = rho_exp, scale = scale, eps = eps_pi),

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

# =============================================================================
# IBD Power variant: uses p^alpha transformation (tunable power)
# alpha=0.5 is Hellinger, lower values compress rare taxa less
# =============================================================================

# Power transformation: p^alpha
power_transform <- function(p, alpha = 0.4, eps = 0) {
  p <- pmax(eps, pmin(1 - eps, p))
  p^alpha
}

# SD of power-transformed beta-binomial (delta method approximation)
# Var(p^alpha) ≈ (alpha * p^(alpha-1))^2 * Var(p)
bb_power_sd_piaware <- function(pi, library_size, rho, alpha = 0.4, eps = 1e-8) {
  pi <- pmax(eps, pmin(1 - eps, pi))
  rho <- pmax(eps, pmin(1 - eps, rho))
  n <- library_size
  
  # Beta-binomial variance of proportion
  var_p <- pi * (1 - pi) * (1 + (n - 1) * rho) / n
  
  # Derivative of p^alpha w.r.t. p
  deriv <- alpha * pi^(alpha - 1)
  
  # Delta method
  var_power <- deriv^2 * var_p
  
  sqrt(pmax(var_power, eps))
}

# Build pi-aware association-adjusted Power distance with high-prevalence filtering
build_assoc_adjusted_piaware_power_highprev_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,
  log_sigma_intercept,
  eps_pi = 1e-8,
  alpha = 0.4,
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
      obs_t = power_transform(p_hat, alpha = alpha, eps = 0),
      # No residuals - use raw transformed value with SD adjustment

      c = .intercept_by_group(as.character(group)),
      sigma0 = exp(c),
      rho0 = sigma_to_rho(sigma0),

      sigma_exp = exp(-log_dispersion_assoc * unconstrained_log_mu + c),
      rho_exp = sigma_to_rho(sigma_exp),

      sd0 = bb_power_sd_piaware(pi = mu, library_size = library_size, rho = rho0, alpha = alpha, eps = eps_pi),
      sd_exp = bb_power_sd_piaware(pi = mu, library_size = library_size, rho = rho_exp, alpha = alpha, eps = eps_pi),

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

# =============================================================================
# IBD Rank variant: uses rank-based transformation
# Robust to outliers and magnitude, completely non-parametric
# =============================================================================

# Build rank-based distance with high-prevalence filtering
build_assoc_adjusted_piaware_rank_highprev_betadisper_inputs <- function(
  count_long,
  log_dispersion_assoc,  # Not used for ranks, but kept for API consistency

  log_sigma_intercept,   # Not used for ranks
  eps_pi = 1e-8,         # Not used for ranks
  min_prevalence = 0.9,
  group_levels = c("non-IBD", "IBD")
) {
  required_cols <- c(
    "sample_id", "taxon_id", "group",
    "count", "library_size"
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

  # Convert proportions to ranks within each sample
  df <- count_long_f %>%
    dplyr::mutate(
      p_hat = count_to_proportion(count = count, library_size = library_size, eps = 0)
    ) %>%
    dplyr::group_by(sample_id) %>%
    dplyr::mutate(
      rank_val = rank(p_hat, ties.method = "average")
    ) %>%
    dplyr::ungroup()

  mat <- df %>%
    dplyr::select(sample_id, taxon_id, rank_val) %>%
    tidyr::pivot_wider(names_from = taxon_id, values_from = rank_val) %>%
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

# ============================================================================
# Example Usage: Standard Deviation Calculation
# ============================================================================

if (FALSE) {  # Set to TRUE to run example
  # Example: Calculate SD from mu and sigma
  mu <- 0.3
  sigma <- 10
  n <- 1000
  
  sd_from_sigma <- calculate_sd_from_sigma(mu, sigma, n)
  cat("SD from sigma:", sd_from_sigma, "\n")
  
  # Example: Calculate SD from mu and overdispersion factor
  overdispersion_factor <- 1.5
  sd_from_factor <- calculate_sd_from_overdispersion_factor(mu, overdispersion_factor, n)
  cat("SD from overdispersion factor:", sd_from_factor, "\n")
  
  # Example: Calculate SD from mu and overdispersion ratio
  overdispersion_ratio <- 1.5
  sd_from_ratio <- calculate_sd_from_overdispersion_ratio(mu, overdispersion_ratio, n)
  cat("SD from overdispersion ratio:", sd_from_ratio, "\n")
  
  # Verify they give the same result when factor = ratio
  cat("Factor and ratio give same result:", 
      abs(sd_from_factor - sd_from_ratio) < 1e-10, "\n")
  
  # Example: Compare SD across different mu values
  mu_values <- seq(0.1, 0.9, by = 0.1)
  sigma_fixed <- 10
  n_fixed <- 1000
  
  sd_values <- sapply(mu_values, function(mu) {
    calculate_sd_from_sigma(mu, sigma_fixed, n_fixed)
  })
  
  plot(mu_values, sd_values, type = "b", 
       xlab = "Mean Probability (μ)", 
       ylab = "Standard Deviation",
       main = "SD vs μ (fixed sigma = 10, n = 1000)")
}

# ============================================================================
# Example Usage: Full Simulation
# ============================================================================

if (FALSE) {  # Set to TRUE to run example
  # Example parameters
  n_taxa <- 10
  n_samples <- 50
  
  # Slope vector: different taxa respond differently to covariates
  slope_vector <- c(-2, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3)
  
  # Inverse softmax mean (baseline intercept)
  inv_softmax_mean <- 0
  
  # Log-dispersion association (k parameter)
  log_dispersion_assoc <- 1.0
  
  # SD of log overdispersion around regression line
  sd_log_overdispersion <- 0.3
  
  # Run simulation
  sim_result <- simulate_compositional_bb(
    slope_vector = slope_vector,
    inv_softmax_mean = inv_softmax_mean,
    log_dispersion_assoc = log_dispersion_assoc,
    n_taxa = n_taxa,
    n_samples = n_samples,
    sd_log_overdispersion = sd_log_overdispersion,
    seed = 123
  )
  
  # Visualize results
  plot_mean_dispersion_relationship(sim_result)
  plot_compositional_abundances(sim_result)
  plot_sigma_by_taxon(sim_result)
  
  # Check summary statistics
  summary(sim_result$count_long$count)
  summary(sim_result$count_long$mu)
  summary(sim_result$count_long$sigma)
}

