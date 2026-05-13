# =============================================================================
# tests/testthat/test-simulation.R - Tests for R/simulation.R
# =============================================================================

test_that("simulate_beta_binomial(): output shape and range", {
  set.seed(123)
  out <- simulate_beta_binomial(n = 1000, mu = 0.3, sigma = 0.5, n_sim = 500)
  expect_length(out, 500L)
  expect_true(all(out >= 0 & out <= 1000))
})

test_that("simulate_beta_binomial(): empirical mean is near n * mu", {
  set.seed(7)
  out <- simulate_beta_binomial(n = 1000, mu = 0.2, sigma = 0.05, n_sim = 5000)
  # n=1000, mu=0.2, low overdispersion -> mean very close to 200.
  expect_equal(mean(out), 200, tolerance = 5)
})

test_that("simulate_beta_binomial(): seed gives reproducible draws", {
  set.seed(11); a <- simulate_beta_binomial(n = 500, mu = 0.4, sigma = 0.2, n_sim = 50)
  set.seed(11); b <- simulate_beta_binomial(n = 500, mu = 0.4, sigma = 0.2, n_sim = 50)
  expect_identical(a, b)
})

test_that("simulate_compositional_bb(): produces expected list structure", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 6L, seed = 1L)

  expected <- c(
    "ground_truth_params", "count_long",
    "cohort_log_linear_predictors", "cohort_log_sigma",
    "cohort_sigma", "cohort_mu", "unique_design", "sample_to_cohort",
    "n_cohorts", "sample_metadata", "taxon_metadata", "design_matrix",
    "mu_inv_softmax", "slope_vector", "parameters"
  )
  expect_true(all(expected %in% names(sim_r)))

  # Two-group design -> two cohorts.
  expect_equal(sim_r$n_cohorts, 2L)
  expect_equal(dim(sim_r$cohort_mu), c(2L, 4L))
  expect_equal(unname(rowSums(sim_r$cohort_mu)), c(1, 1), tolerance = 1e-10)

  # 12 samples * 4 taxa = 48 rows in count_long.
  expect_equal(nrow(sim_r$count_long), 12L * 4L)
  expected_cols <- c("sample_id", "taxon_id", "count", "library_size",
                     "mu", "sigma", "log_sigma", "unconstrained_log_mu")
  expect_true(all(expected_cols %in% colnames(sim_r$count_long)))

  # Library size is floored at 100.
  expect_true(all(sim_r$count_long$library_size >= 100))

  # Counts are non-negative integers (as doubles via pmap_dbl).
  expect_true(all(sim_r$count_long$count >= 0))
  expect_true(all(sim_r$count_long$count <= sim_r$count_long$library_size))
})

test_that("simulate_compositional_bb(): mu_inv_softmax must sum to 0", {
  expect_error(
    .quietly({
      simulate_compositional_bb(
        slope_vector = c(0.5, -0.5),
        mu_inv_softmax = c(1, 1),    # sums to 2, not 0
        log_dispersion_assoc = 0.4,
        n_taxa = 2L,
        n_samples = 8L,
        sd_log_overdispersion = 0.1,
        intercept_dispersion = 1.5,
        library_size_mean = 1000,
        library_size_sd = 100,
        design_matrix = design_matrix_from_groups(4L),
        seed = 1L
      )
    }),
    regexp = "mu_inv_softmax must sum to 0"
  )
})

test_that("simulate_compositional_bb(): slope_vector must sum to 0", {
  expect_error(
    .quietly({
      simulate_compositional_bb(
        slope_vector = c(0.5, 0.5),   # sums to 1, not 0
        mu_inv_softmax = c(0.5, -0.5),
        log_dispersion_assoc = 0.4,
        n_taxa = 2L,
        n_samples = 8L,
        sd_log_overdispersion = 0.1,
        intercept_dispersion = 1.5,
        library_size_mean = 1000,
        library_size_sd = 100,
        design_matrix = design_matrix_from_groups(4L),
        seed = 1L
      )
    }),
    regexp = "slope_vector must sum to 0"
  )
})

test_that("simulate_compositional_bb(): vector / matrix dispersion intercepts are mutually exclusive", {
  expect_error(
    .quietly({
      simulate_compositional_bb(
        slope_vector = c(0.5, -0.5),
        mu_inv_softmax = c(0.5, -0.5),
        log_dispersion_assoc = 0.4,
        n_taxa = 2L,
        n_samples = 8L,
        sd_log_overdispersion = 0.1,
        intercept_dispersion = 1.5,
        intercept_dispersion_matrix = matrix(1.5, nrow = 2L, ncol = 2L),
        library_size_mean = 1000,
        library_size_sd = 100,
        design_matrix = design_matrix_from_groups(4L),
        seed = 1L
      )
    }),
    regexp = "Provide only one"
  )
})

test_that("simulate_compositional_bb(): reproducible given the same seed", {
  a <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 5L, seed = 99L)
  b <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 5L, seed = 99L)
  expect_equal(a$count_long$count, b$count_long$count)
})

test_that("simulate_compositional_bb(): wrong-shape intercept_dispersion_matrix errors", {
  expect_error(
    .quietly({
      simulate_compositional_bb(
        slope_vector = c(0.5, -0.5),
        mu_inv_softmax = c(0.5, -0.5),
        log_dispersion_assoc = 0.4,
        n_taxa = 2L,
        n_samples = 8L,
        sd_log_overdispersion = 0.1,
        intercept_dispersion = NULL,
        intercept_dispersion_matrix = matrix(1, nrow = 3L, ncol = 2L),  # 3 != n_cohorts
        library_size_mean = 1000,
        library_size_sd = 100,
        design_matrix = design_matrix_from_groups(4L),
        seed = 1L
      )
    }),
    regexp = "intercept_dispersion_matrix must be n_cohorts x n_taxa"
  )
})
