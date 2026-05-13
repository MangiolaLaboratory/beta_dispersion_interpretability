# =============================================================================
# tests/testthat/helper-setup.R - Shared setup and fixtures
# =============================================================================
#
# testthat auto-sources every file matching `^helper.*\\.R$` in this directory
# before running tests. We use it to:
#
#   1. Locate the project root (so tests can be run from any cwd).
#   2. Source the top-level loader so every R/*.R module is in scope.
#   3. Define small, fast fixtures that several tests reuse.
#
# Keep fixtures tiny: simulation calls are CPU-bound; n_taxa = 4-6 and
# n_samples = 20-30 are usually enough to verify shapes and constraints
# without making the suite slow.
# =============================================================================

# Project root (works whether tests are run from the project root, the tests
# directory, or via Rscript tests/testthat.R).
.proj_root <- tryCatch(
  here::here(),
  error = function(e) {
    # Fall back: walk up from the test file until we hit functions.R.
    p <- normalizePath(getwd(), mustWork = FALSE)
    while (nchar(p) > 1 && !file.exists(file.path(p, "functions.R"))) {
      p <- dirname(p)
    }
    p
  }
)

# Source the loader (which sources every R/*.R module).
suppressMessages(source(file.path(.proj_root, "functions.R")))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Many of the simulators print progress via `cat()`. Wrap them so the test
# reporter stays readable.
.quietly <- function(expr) {
  utils::capture.output(out <- eval.parent(substitute(expr)))
  out
}

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

#' Build a tiny, deterministic compositional simulation result.
#'
#' Defaults: 5 taxa, 10 samples per group, balanced two-group design, seed 1.
#' The result has `sample_metadata$group` set as a factor with levels
#' `c("ctrl", "case")` so downstream betadisper / PERMANOVA tests work.
.make_tiny_sim <- function(
  n_taxa = 5L,
  n_samples_per_group = 10L,
  seed = 1L,
  group_levels = c("ctrl", "case"),
  slope_strength = 0.6,
  log_dispersion_assoc = 0.4,
  sd_log_overdispersion = 0.2,
  intercept_dispersion = 1.5,
  library_size_mean = 5000,
  library_size_sd = 500
) {
  slope <- c(slope_strength, -slope_strength, rep(0, n_taxa - 2))
  slope <- slope - mean(slope)
  mu <- seq(-1, 1, length.out = n_taxa)
  mu <- mu - mean(mu)

  sim_r <- .quietly({
    simulate_compositional_bb(
      slope_vector = slope,
      mu_inv_softmax = mu,
      log_dispersion_assoc = log_dispersion_assoc,
      n_taxa = n_taxa,
      n_samples = 2L * n_samples_per_group,
      sd_log_overdispersion = sd_log_overdispersion,
      intercept_dispersion = intercept_dispersion,
      library_size_mean = library_size_mean,
      library_size_sd = library_size_sd,
      design_matrix = design_matrix_from_groups(n_samples_per_group),
      seed = seed
    )
  })

  sim_r <- synthetic_sim_add_two_group_factor(sim_r, group_levels)
  sim_r$seed <- as.integer(seed)
  sim_r$n_da_taxa <- 2L
  sim_r
}

#' Build a tiny sccomp-result-shaped tibble for sccomp.R tests.
#'
#' Mimics the relevant columns of `sccomp::sccomp_estimate()` output:
#'   cell_group, parameter, c_effect, c_lower, c_upper, v_effect.
#' Optionally attaches a fake `model_input` attribute with N/M/exposure.
.make_fake_sccomp_result <- function(
  n_taxa = 6L,
  group_parameter = "groupcase",
  attach_model_input = TRUE,
  exposure_n = 30L,
  exposure_mean = 4000,
  seed = 7L
) {
  taxa <- paste0("Taxon_", seq_len(n_taxa))
  set.seed(seed)
  c_int <- stats::rnorm(n_taxa, 0, 1)
  c_int <- c_int - mean(c_int)
  # Deterministic slopes: exactly two taxa are "significant" (CI excludes 0),
  # the rest have wide CIs straddling 0 so they are NOT flagged.
  c_slp <- rep(0, n_taxa)
  c_slp[1] <-  0.30
  c_slp[2] <- -0.30
  c_low <- c_slp - 0.50
  c_upp <- c_slp + 0.50
  c_low[1] <-  0.05; c_upp[1] <-  0.50
  c_low[2] <- -0.50; c_upp[2] <- -0.05
  v_int <- stats::rnorm(n_taxa, 1.5, 0.3)
  v_slp <- stats::rnorm(n_taxa, 0, 0.2)

  result <- dplyr::bind_rows(
    tibble::tibble(
      cell_group = taxa,
      parameter = "Intercept",
      c_effect = c_int,
      c_lower = c_int - 0.1,
      c_upper = c_int + 0.1,
      v_effect = v_int
    ),
    tibble::tibble(
      cell_group = taxa,
      parameter = group_parameter,
      c_effect = c_slp,
      c_lower = c_low,
      c_upper = c_upp,
      v_effect = v_slp
    )
  )

  if (attach_model_input) {
    set.seed(seed + 1L)
    exposure <- pmax(round(stats::rnorm(exposure_n, exposure_mean, exposure_mean * 0.2)), 100)
    attr(result, "model_input") <- list(
      N = as.integer(exposure_n),
      M = as.integer(n_taxa),
      exposure = exposure
    )
  }

  result
}
