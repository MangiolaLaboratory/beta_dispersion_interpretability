# =============================================================================
# tests/testthat/test-simulation.R - Tests for R/simulation.R
# =============================================================================
#
# Strategy: rely on the simulator's `seed` argument (and explicit set.seed())
# to make outputs deterministic, then pin exact counts / matrix values, list
# element NAMES AND ORDER, column names of count_long, etc.
# =============================================================================

# ---------------------------------------------------------------------------
# simulate_beta_binomial()
# ---------------------------------------------------------------------------

test_that("simulate_beta_binomial(): exact reproducible draw snapshot", {
  set.seed(11)
  out <- simulate_beta_binomial(n = 500, mu = 0.4, sigma = 0.2, n_sim = 5)
  # Pinned snapshot for set.seed(11) + given parameters.
  expect_equal(out, c(140, 198, 59, 351, 91))
})

test_that("simulate_beta_binomial(): output type and range", {
  set.seed(123)
  out <- simulate_beta_binomial(n = 1000, mu = 0.3, sigma = 0.5, n_sim = 500)
  expect_type(out, "double")
  expect_length(out, 500L)
  expect_true(all(out >= 0 & out <= 1000))
})

test_that("simulate_beta_binomial(): empirical mean near n*mu under low dispersion", {
  set.seed(7)
  out <- simulate_beta_binomial(n = 1000, mu = 0.2, sigma = 0.05, n_sim = 5000)
  expect_equal(mean(out), 200, tolerance = 5)
})

test_that("simulate_beta_binomial(): identical RNG state -> identical draws", {
  set.seed(11); a <- simulate_beta_binomial(n = 500, mu = 0.4, sigma = 0.2, n_sim = 50)
  set.seed(11); b <- simulate_beta_binomial(n = 500, mu = 0.4, sigma = 0.2, n_sim = 50)
  expect_identical(a, b)
})

# ---------------------------------------------------------------------------
# simulate_compositional_bb(): structure and exact snapshot
# ---------------------------------------------------------------------------

test_that("simulate_compositional_bb(): return-list names and order are pinned", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 6L, seed = 1L)
  expect_named(
    sim_r,
    c("ground_truth_params", "count_long",
      "cohort_log_linear_predictors", "cohort_log_sigma", "cohort_sigma",
      "cohort_mu", "unique_design", "sample_to_cohort", "n_cohorts",
      "sample_metadata", "taxon_metadata", "design_matrix",
      "mu_inv_softmax", "slope_vector", "parameters",
      "seed", "n_da_taxa"),  # extra fields added by .make_tiny_sim
    ignore.order = FALSE
  )
})

test_that("simulate_compositional_bb(): count_long has the expected columns (set)", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 6L, seed = 1L)
  expect_setequal(
    colnames(sim_r$count_long),
    c("sample_id", "taxon_id", "library_size", "mu", "sigma",
      "unconstrained_log_mu", "log_sigma", "count",
      "Intercept", "Group", "group")
  )
  expect_equal(nrow(sim_r$count_long), 12L * 4L)
})

test_that("simulate_compositional_bb(): cohort matrices have the right shape and rowSums == 1", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 6L, seed = 1L)
  expect_equal(sim_r$n_cohorts, 2L)
  expect_equal(dim(sim_r$cohort_mu), c(2L, 4L))
  expect_equal(unname(rowSums(sim_r$cohort_mu)), c(1, 1), tolerance = 1e-12)
  expect_equal(dim(sim_r$cohort_sigma), c(2L, 4L))
  expect_equal(dim(sim_r$cohort_log_sigma), c(2L, 4L))
  expect_equal(dim(sim_r$cohort_log_linear_predictors), c(2L, 4L))
  # log_sigma <-> sigma consistency.
  expect_equal(exp(sim_r$cohort_log_sigma), sim_r$cohort_sigma, tolerance = 1e-12)
})

test_that("simulate_compositional_bb(): cohort_mu values are pinned (seed = 1)", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 6L, seed = 1L)
  expect_equal(
    sim_r$cohort_mu[1, ],
    c(Taxon_1 = 0.07076911, Taxon_2 = 0.13783941,
      Taxon_3 = 0.26847450, Taxon_4 = 0.52291696),
    tolerance = 1e-7
  )
  expect_equal(
    sim_r$cohort_mu[2, ],
    c(Taxon_1 = 0.12946902, Taxon_2 = 0.07595251,
      Taxon_3 = 0.26955570, Taxon_4 = 0.52502277),
    tolerance = 1e-7
  )
})

test_that("simulate_compositional_bb(): library sizes are pinned and >= 100", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 6L, seed = 1L)
  expect_equal(
    sim_r$sample_metadata$library_size,
    c(5165, 4590, 5244, 5369, 5288, 4847,
      5756, 5195, 4689, 3893, 5562, 4978)
  )
  expect_true(all(sim_r$sample_metadata$library_size >= 100))
})

test_that("simulate_compositional_bb(): first counts are pinned (deterministic via seed)", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 6L, seed = 1L)
  expect_equal(head(sim_r$count_long$count, 5L),
               c(2033, 0, 0, 0, 0))
  expect_true(all(sim_r$count_long$count >= 0))
  expect_true(all(sim_r$count_long$count <= sim_r$count_long$library_size))
})

test_that("simulate_compositional_bb(): ground_truth_params layout", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 6L, seed = 1L)
  expect_equal(nrow(sim_r$ground_truth_params), 2L * 4L)
  expect_setequal(
    colnames(sim_r$ground_truth_params),
    c("taxon_id", "cohort_idx", "baseline_intercept", "slope",
      "cohort_log_linear_predictors", "cohort_mu",
      "cohort_log_sigma", "cohort_sigma")
  )
  # Two cohorts, four taxa each.
  expect_equal(as.integer(table(sim_r$ground_truth_params$cohort_idx)),
               c(4L, 4L))
})

test_that("simulate_compositional_bb(): same seed -> identical counts", {
  a <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 5L, seed = 99L)
  b <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 5L, seed = 99L)
  expect_equal(a$count_long$count, b$count_long$count)
  expect_equal(a$sample_metadata$library_size, b$sample_metadata$library_size)
  expect_equal(a$cohort_mu, b$cohort_mu)
})

# ---------------------------------------------------------------------------
# simulate_compositional_bb(): validation
# ---------------------------------------------------------------------------

test_that("simulate_compositional_bb(): mu_inv_softmax must sum to 0", {
  expect_error(
    .quietly({
      simulate_compositional_bb(
        slope_vector = c(0.5, -0.5),
        mu_inv_softmax = c(1, 1),
        log_dispersion_assoc = 0.4,
        n_taxa = 2L, n_samples = 8L,
        sd_log_overdispersion = 0.1,
        intercept_dispersion = 1.5,
        library_size_mean = 1000, library_size_sd = 100,
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
        slope_vector = c(0.5, 0.5),
        mu_inv_softmax = c(0.5, -0.5),
        log_dispersion_assoc = 0.4,
        n_taxa = 2L, n_samples = 8L,
        sd_log_overdispersion = 0.1,
        intercept_dispersion = 1.5,
        library_size_mean = 1000, library_size_sd = 100,
        design_matrix = design_matrix_from_groups(4L),
        seed = 1L
      )
    }),
    regexp = "slope_vector must sum to 0"
  )
})

test_that("simulate_compositional_bb(): wrong-length vectors error out", {
  expect_error(
    .quietly({
      simulate_compositional_bb(
        slope_vector = c(0.5, -0.5, 0),       # length 3 != n_taxa = 2
        mu_inv_softmax = c(0.5, -0.5),
        log_dispersion_assoc = 0.4,
        n_taxa = 2L, n_samples = 8L,
        sd_log_overdispersion = 0.1,
        intercept_dispersion = 1.5,
        library_size_mean = 1000, library_size_sd = 100,
        design_matrix = design_matrix_from_groups(4L),
        seed = 1L
      )
    }),
    regexp = "slope_vector length must equal n_taxa"
  )

  expect_error(
    .quietly({
      simulate_compositional_bb(
        slope_vector = c(0.5, -0.5),
        mu_inv_softmax = c(0.5, -0.5, 0),     # length 3 != n_taxa = 2
        log_dispersion_assoc = 0.4,
        n_taxa = 2L, n_samples = 8L,
        sd_log_overdispersion = 0.1,
        intercept_dispersion = 1.5,
        library_size_mean = 1000, library_size_sd = 100,
        design_matrix = design_matrix_from_groups(4L),
        seed = 1L
      )
    }),
    regexp = "mu_inv_softmax length must equal n_taxa"
  )
})

test_that("simulate_compositional_bb(): cannot pass both intercept_dispersion and intercept_dispersion_matrix", {
  expect_error(
    .quietly({
      simulate_compositional_bb(
        slope_vector = c(0.5, -0.5),
        mu_inv_softmax = c(0.5, -0.5),
        log_dispersion_assoc = 0.4,
        n_taxa = 2L, n_samples = 8L,
        sd_log_overdispersion = 0.1,
        intercept_dispersion = 1.5,
        intercept_dispersion_matrix = matrix(1.5, nrow = 2L, ncol = 2L),
        library_size_mean = 1000, library_size_sd = 100,
        design_matrix = design_matrix_from_groups(4L),
        seed = 1L
      )
    }),
    regexp = "Provide only one"
  )
})

test_that("simulate_compositional_bb(): wrong-shape intercept_dispersion_matrix errors", {
  expect_error(
    .quietly({
      simulate_compositional_bb(
        slope_vector = c(0.5, -0.5),
        mu_inv_softmax = c(0.5, -0.5),
        log_dispersion_assoc = 0.4,
        n_taxa = 2L, n_samples = 8L,
        sd_log_overdispersion = 0.1,
        intercept_dispersion = NULL,
        intercept_dispersion_matrix = matrix(1, nrow = 3L, ncol = 2L),
        library_size_mean = 1000, library_size_sd = 100,
        design_matrix = design_matrix_from_groups(4L),
        seed = 1L
      )
    }),
    regexp = "intercept_dispersion_matrix must be n_cohorts x n_taxa"
  )
})

test_that("simulate_compositional_bb(): intercept_dispersion of wrong length errors", {
  expect_error(
    .quietly({
      simulate_compositional_bb(
        slope_vector = c(0.5, -0.5),
        mu_inv_softmax = c(0.5, -0.5),
        log_dispersion_assoc = 0.4,
        n_taxa = 2L, n_samples = 8L,
        sd_log_overdispersion = 0.1,
        intercept_dispersion = c(1.5, 1.5, 1.5),  # length 3 != n_taxa = 2
        library_size_mean = 1000, library_size_sd = 100,
        design_matrix = design_matrix_from_groups(4L),
        seed = 1L
      )
    }),
    regexp = "intercept_dispersion must have length n_taxa"
  )
})
