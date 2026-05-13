# =============================================================================
# tests/testthat/test-sccomp.R - Tests for R/sccomp.R
# =============================================================================

# ---------------------------------------------------------------------------
# library_size_from_sccomp_model_input
# ---------------------------------------------------------------------------

test_that("library_size_from_sccomp_model_input(): NULL or empty -> NA pair", {
  expect_equal(
    library_size_from_sccomp_model_input(NULL),
    list(library_size_mean = NA_real_, library_size_sd = NA_real_)
  )
  expect_equal(
    library_size_from_sccomp_model_input(list(exposure = numeric(0))),
    list(library_size_mean = NA_real_, library_size_sd = NA_real_)
  )
  expect_equal(
    library_size_from_sccomp_model_input(list(exposure = c(0, -1, NA))),
    list(library_size_mean = NA_real_, library_size_sd = NA_real_)
  )
})

test_that("library_size_from_sccomp_model_input(): normal exposure -> mean and floored SD", {
  ex <- c(1000, 1200, 1100, 1300, 900)
  out <- library_size_from_sccomp_model_input(list(exposure = ex))
  expect_equal(out$library_size_mean, mean(ex))
  # SD is real, finite, and >= 5% of the mean.
  expect_true(is.finite(out$library_size_sd))
  expect_gte(out$library_size_sd, 0.05 * mean(ex))
})

test_that("library_size_from_sccomp_model_input(): constant exposure -> SD floor kicks in", {
  ex <- rep(2000, 10)
  out <- library_size_from_sccomp_model_input(list(exposure = ex))
  expect_equal(out$library_size_mean, 2000)
  # sd(ex) == 0 -> replaced with max(0.05 * mean, 1) = 100.
  expect_equal(out$library_size_sd, 100)
})

# ---------------------------------------------------------------------------
# .sim_library_size_or_default
# ---------------------------------------------------------------------------

test_that(".sim_library_size_or_default(): falls back to fixed defaults for bad mean", {
  fixed <- list(library_size_mean = 15125, library_size_sd = 5000)
  expect_equal(.sim_library_size_or_default(list()), fixed)
  expect_equal(.sim_library_size_or_default(list(library_size_mean = NA_real_)), fixed)
  expect_equal(.sim_library_size_or_default(list(library_size_mean = 0)), fixed)
  expect_equal(.sim_library_size_or_default(list(library_size_mean = -50)), fixed)
})

test_that(".sim_library_size_or_default(): repairs bad SD but keeps mean", {
  out <- .sim_library_size_or_default(list(library_size_mean = 4000, library_size_sd = NA_real_))
  expect_equal(out$library_size_mean, 4000)
  expect_equal(out$library_size_sd, 0.05 * 4000)

  out2 <- .sim_library_size_or_default(list(library_size_mean = 4000, library_size_sd = -1))
  expect_equal(out2$library_size_sd, 0.05 * 4000)
})

test_that(".sim_library_size_or_default(): passes through good values", {
  out <- .sim_library_size_or_default(list(library_size_mean = 4000, library_size_sd = 800))
  expect_equal(out, list(library_size_mean = 4000, library_size_sd = 800))
})

# ---------------------------------------------------------------------------
# extract_sccomp_params_brito
# ---------------------------------------------------------------------------

test_that("extract_sccomp_params_brito(): pulls intercept/slope and library-size", {
  result <- .make_fake_sccomp_result(n_taxa = 6L)
  params <- extract_sccomp_params_brito(result)

  # Basic shape.
  expect_named(params, c(
    "n_samples_real", "n_taxa_real", "k_real", "k_sd_real",
    "c_real", "prec_sd_real", "intercept_effects", "slope_effects",
    "slope_effects_significant", "v_intercept", "v_slope",
    "alpha_intercept_effects", "group_parameter",
    "library_size_mean", "library_size_sd"
  ))
  expect_equal(params$n_taxa_real, 6L)
  expect_equal(params$n_samples_real, 30L)

  # k_real always 0 in the Brito path; uncertainty fields are NA.
  expect_equal(params$k_real, 0)
  expect_true(is.na(params$k_sd_real))
  expect_true(is.na(params$c_real))
  expect_true(is.na(params$prec_sd_real))

  # Intercept / slope / alpha vectors have length n_taxa.
  expect_length(params$intercept_effects, 6L)
  expect_length(params$slope_effects, 6L)
  expect_length(params$v_intercept, 6L)
  expect_length(params$v_slope, 6L)
  expect_length(params$alpha_intercept_effects, 6L)
  expect_identical(params$alpha_intercept_effects, params$v_intercept)

  # Two of the slope effects are significant (see fixture).
  expect_length(params$slope_effects_significant, 2L)

  # group_parameter is auto-detected.
  expect_equal(params$group_parameter, "groupcase")

  # Library size pulled from the attached model_input.
  expect_true(is.finite(params$library_size_mean))
  expect_true(is.finite(params$library_size_sd))
})

test_that("extract_sccomp_params_brito(): falls back when model_input is absent", {
  result <- .make_fake_sccomp_result(n_taxa = 5L, attach_model_input = FALSE)
  params <- extract_sccomp_params_brito(result)
  # Default fallback: n_samples_real = 178, n_taxa_real = #unique cell_groups.
  expect_equal(params$n_samples_real, 178)
  expect_equal(params$n_taxa_real, 5L)
  # And no library-size info available.
  expect_true(is.na(params$library_size_mean))
  expect_true(is.na(params$library_size_sd))
})

test_that("extract_sccomp_params_brito(): honours target_parameter override", {
  # Two non-intercept parameters: 'groupcase' and 'sex'.
  base <- .make_fake_sccomp_result(n_taxa = 4L, group_parameter = "groupcase",
                                    attach_model_input = FALSE)
  sex_rows <- base[base$parameter == "groupcase", , drop = FALSE]
  sex_rows$parameter <- "sex"
  result <- dplyr::bind_rows(base, sex_rows)
  params <- extract_sccomp_params_brito(result, target_parameter = "sex")
  expect_equal(params$group_parameter, "sex")
})

test_that("extract_sccomp_params_brito(): rejects bad target_parameter", {
  result <- .make_fake_sccomp_result(n_taxa = 4L, attach_model_input = FALSE)
  expect_error(
    extract_sccomp_params_brito(result, target_parameter = "does_not_exist"),
    regexp = "group_parameter"
  )
})

# ---------------------------------------------------------------------------
# build_simulation_params_brito
# ---------------------------------------------------------------------------

test_that("build_simulation_params_brito(): builds compositional inputs", {
  result <- .make_fake_sccomp_result(n_taxa = 8L)
  params <- extract_sccomp_params_brito(result)
  sim_params <- build_simulation_params_brito(
    sccomp_params = params,
    n_reps_auc = 4L,
    group_levels = c("ctrl", "case")
  )

  # Required fields.
  expect_true(all(c(
    "n_taxa", "n_samples_per_group", "n_groups", "n_samples",
    "library_size_mean", "library_size_sd",
    "mu_inv_softmax_base_realistic", "slope_realistic",
    "intercept_disp_realistic", "intercept_disp_realistic_taxon",
    "sigma_realistic", "sd_log_overdispersion_realistic",
    "alpha_intercept_sampled", "n_reps_auc", "perm_reps_auc",
    "rep_seeds_auc", "group_levels"
  ) %in% names(sim_params)))

  expect_equal(sim_params$n_taxa, 8L)
  expect_equal(sim_params$n_groups, 2L)
  # n_samples_real == 30 -> floor(30/2) == 15 per group.
  expect_equal(sim_params$n_samples_per_group, 15)
  expect_equal(sim_params$n_samples, 30)

  # mu_inv_softmax must satisfy compositional constraint.
  expect_equal(sum(sim_params$mu_inv_softmax_base_realistic), 0, tolerance = 1e-12)
  # Slope must also be centred.
  expect_equal(sum(sim_params$slope_realistic), 0, tolerance = 1e-12)

  # Sigma realistic is positive.
  expect_true(sim_params$sigma_realistic > 0)

  # Rep seeds form a contiguous block of `n_reps_auc` values.
  expect_length(sim_params$rep_seeds_auc, 4L)
  expect_equal(sim_params$rep_seeds_auc, 12000 + 1:4)
  expect_equal(sim_params$perm_reps_auc, 199)
})

test_that("build_simulation_params_brito(): honours custom n_samples_per_group", {
  result <- .make_fake_sccomp_result(n_taxa = 4L)
  params <- extract_sccomp_params_brito(result)

  bp <- build_simulation_params_brito(params, group_levels = c("a", "b"),
                                       n_samples_per_group = 12L)
  expect_equal(bp$n_samples_per_group, 12L)
  expect_equal(bp$n_samples, 24L)

  bp2 <- build_simulation_params_brito(params, group_levels = c("a", "b"),
                                        n_samples_per_group = c(8L, 14L))
  expect_equal(bp2$n_samples_per_group, c(8L, 14L))
  expect_equal(bp2$n_samples, 22L)
})

test_that("build_simulation_params_brito(): rejects bad group_levels / sample counts", {
  result <- .make_fake_sccomp_result(n_taxa = 4L)
  params <- extract_sccomp_params_brito(result)
  expect_error(
    build_simulation_params_brito(params, group_levels = c("a", "b", "c")),
    regexp = "length 2"
  )
  expect_error(
    build_simulation_params_brito(params, group_levels = c("a", "b"),
                                  n_samples_per_group = c(1L, 2L, 3L)),
    regexp = "length 1.*length 2"
  )
})
