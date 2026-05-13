# =============================================================================
# tests/testthat/test-sccomp.R - Tests for R/sccomp.R
# =============================================================================
#
# Strategy: the fixture in helper-setup.R has hard-coded slope values (so the
# significance test is deterministic) and a seeded `model_input$exposure`.
# Tests below pin the exact extracted parameters, list names and order.
# =============================================================================

# ---------------------------------------------------------------------------
# library_size_from_sccomp_model_input()
# ---------------------------------------------------------------------------

test_that("library_size_from_sccomp_model_input(): NULL / empty -> NA pair (exact list)", {
  expect_identical(
    library_size_from_sccomp_model_input(NULL),
    list(library_size_mean = NA_real_, library_size_sd = NA_real_)
  )
  expect_identical(
    library_size_from_sccomp_model_input(list(exposure = numeric(0))),
    list(library_size_mean = NA_real_, library_size_sd = NA_real_)
  )
  expect_identical(
    library_size_from_sccomp_model_input(list(exposure = c(0, -1, NA))),
    list(library_size_mean = NA_real_, library_size_sd = NA_real_)
  )
})

test_that("library_size_from_sccomp_model_input(): exact snapshot for the helper fixture", {
  result <- .make_fake_sccomp_result(n_taxa = 6L)
  out <- library_size_from_sccomp_model_input(attr(result, "model_input"))
  # Seeded snapshot (depends on .make_fake_sccomp_result's seed = 7).
  expect_equal(out$library_size_mean, 3919.6, tolerance = 1e-6)
  expect_equal(out$library_size_sd,   886.1921, tolerance = 1e-4)
})

test_that("library_size_from_sccomp_model_input(): mean / SD with floor on small SD", {
  # SD = 5 -> below 5% of mean = 50; floored at 50.
  ex <- c(1000, 1000, 1010, 990, 1000)
  out <- library_size_from_sccomp_model_input(list(exposure = ex))
  expect_equal(out$library_size_mean, 1000)
  expect_equal(out$library_size_sd, max(stats::sd(ex), 0.05 * mean(ex)),
               tolerance = 1e-12)
})

test_that("library_size_from_sccomp_model_input(): constant exposure -> SD = max(0.05*mean, 1) = 100", {
  out <- library_size_from_sccomp_model_input(list(exposure = rep(2000, 10)))
  expect_equal(out, list(library_size_mean = 2000, library_size_sd = 100))
})

# ---------------------------------------------------------------------------
# .sim_library_size_or_default()
# ---------------------------------------------------------------------------

test_that(".sim_library_size_or_default(): fixed fallback (15125 / 5000) for bad mean", {
  fixed <- list(library_size_mean = 15125, library_size_sd = 5000)
  expect_identical(.sim_library_size_or_default(list()), fixed)
  expect_identical(.sim_library_size_or_default(list(library_size_mean = NA_real_)), fixed)
  expect_identical(.sim_library_size_or_default(list(library_size_mean = 0)), fixed)
  expect_identical(.sim_library_size_or_default(list(library_size_mean = -50)), fixed)
})

test_that(".sim_library_size_or_default(): repairs bad SD but keeps mean (exact values)", {
  expect_identical(
    .sim_library_size_or_default(list(library_size_mean = 4000, library_size_sd = NA_real_)),
    list(library_size_mean = 4000, library_size_sd = 200)  # 0.05 * 4000
  )
  expect_identical(
    .sim_library_size_or_default(list(library_size_mean = 4000, library_size_sd = -1)),
    list(library_size_mean = 4000, library_size_sd = 200)
  )
})

test_that(".sim_library_size_or_default(): passes through good values verbatim", {
  out <- .sim_library_size_or_default(list(library_size_mean = 4000, library_size_sd = 800))
  expect_identical(out, list(library_size_mean = 4000, library_size_sd = 800))
})

# ---------------------------------------------------------------------------
# extract_sccomp_params_brito()
# ---------------------------------------------------------------------------

test_that("extract_sccomp_params_brito(): return-list names and order are pinned", {
  result <- .make_fake_sccomp_result(n_taxa = 6L)
  params <- extract_sccomp_params_brito(result)
  expect_named(
    params,
    c("n_samples_real", "n_taxa_real", "k_real", "k_sd_real", "c_real",
      "prec_sd_real", "intercept_effects", "slope_effects",
      "slope_effects_significant", "v_intercept", "v_slope",
      "alpha_intercept_effects", "group_parameter",
      "library_size_mean", "library_size_sd"),
    ignore.order = FALSE
  )
})

test_that("extract_sccomp_params_brito(): exact scalar field values for fixture", {
  result <- .make_fake_sccomp_result(n_taxa = 6L)
  params <- extract_sccomp_params_brito(result)
  expect_equal(params$n_samples_real, 30L)
  expect_equal(params$n_taxa_real, 6L)
  expect_equal(params$k_real, 0)
  expect_true(is.na(params$k_sd_real))
  expect_true(is.na(params$c_real))
  expect_true(is.na(params$prec_sd_real))
  expect_identical(params$group_parameter, "groupcase")
  expect_equal(params$library_size_mean, 3919.6, tolerance = 1e-6)
  expect_equal(params$library_size_sd,   886.1921, tolerance = 1e-4)
})

test_that("extract_sccomp_params_brito(): vectors have correct length and values for fixture", {
  result <- .make_fake_sccomp_result(n_taxa = 6L)
  params <- extract_sccomp_params_brito(result)

  # All per-taxon vectors length 6 (n_taxa).
  expect_length(params$intercept_effects, 6L)
  expect_length(params$slope_effects, 6L)
  expect_length(params$v_intercept, 6L)
  expect_length(params$v_slope, 6L)
  expect_length(params$alpha_intercept_effects, 6L)

  # Snapshot values (seeded by .make_fake_sccomp_result(seed = 7)).
  expect_equal(
    params$intercept_effects,
    c(2.60959104, -0.87442780, -0.37194863, -0.08994907, -0.64832946, -0.62493607),
    tolerance = 1e-7
  )
  # Fixture sets slope deterministically: only taxa 1 / 2 are non-zero.
  expect_equal(params$slope_effects, c(0.3, -0.3, 0.0, 0.0, 0.0, 0.0))

  # Exactly two slopes are "significant" (CI excludes 0).
  expect_equal(params$slope_effects_significant, c(0.3, -0.3))

  # alpha_intercept_effects == v_intercept (per the function's semantics).
  expect_identical(params$alpha_intercept_effects, params$v_intercept)
  expect_equal(
    params$v_intercept,
    c(1.724442, 1.464913, 1.545797, 2.156993, 1.607096, 2.315026),
    tolerance = 1e-6
  )
})

test_that("extract_sccomp_params_brito(): fallback when model_input is absent", {
  result <- .make_fake_sccomp_result(n_taxa = 5L, attach_model_input = FALSE)
  params <- extract_sccomp_params_brito(result)
  expect_equal(params$n_samples_real, 178)      # fixed legacy default
  expect_equal(params$n_taxa_real, 5L)
  expect_true(is.na(params$library_size_mean))
  expect_true(is.na(params$library_size_sd))
})

test_that("extract_sccomp_params_brito(): target_parameter override is honoured", {
  base <- .make_fake_sccomp_result(n_taxa = 4L, group_parameter = "groupcase",
                                    attach_model_input = FALSE)
  sex_rows <- base[base$parameter == "groupcase", , drop = FALSE]
  sex_rows$parameter <- "sex"
  result <- dplyr::bind_rows(base, sex_rows)
  params <- extract_sccomp_params_brito(result, target_parameter = "sex")
  expect_equal(params$group_parameter, "sex")
  # And slope_effects come from the SEX rows, not groupcase.
  expect_equal(params$slope_effects, params$slope_effects)  # tautology placeholder
  expect_equal(
    params$slope_effects,
    sex_rows$c_effect
  )
})

test_that("extract_sccomp_params_brito(): rejects unknown target_parameter", {
  result <- .make_fake_sccomp_result(n_taxa = 4L, attach_model_input = FALSE)
  expect_error(
    extract_sccomp_params_brito(result, target_parameter = "does_not_exist"),
    regexp = "group_parameter"
  )
})

# ---------------------------------------------------------------------------
# build_simulation_params_brito()
# ---------------------------------------------------------------------------

test_that("build_simulation_params_brito(): return-list names and order are pinned", {
  result <- .make_fake_sccomp_result(n_taxa = 6L)
  params <- extract_sccomp_params_brito(result)
  sim_params <- build_simulation_params_brito(
    sccomp_params = params,
    n_reps_auc = 4L,
    group_levels = c("ctrl", "case")
  )
  expect_named(
    sim_params,
    c("n_taxa", "n_samples_per_group", "n_groups", "n_samples",
      "group_levels", "library_size_mean", "library_size_sd",
      "mu_inv_softmax_base_realistic",
      "intercept_disp_realistic", "intercept_disp_realistic_taxon",
      "sigma_realistic", "slope_realistic",
      "sd_log_overdispersion_realistic", "alpha_intercept_sampled",
      "n_reps_auc", "perm_reps_auc", "rep_seeds_auc"),
    ignore.order = FALSE
  )
})

test_that("build_simulation_params_brito(): exact snapshot values for fixture", {
  result <- .make_fake_sccomp_result(n_taxa = 6L)
  params <- extract_sccomp_params_brito(result)
  sim_params <- build_simulation_params_brito(
    sccomp_params = params,
    n_reps_auc = 4L,
    group_levels = c("ctrl", "case")
  )

  # Sample-size fields.
  expect_equal(sim_params$n_taxa, 6L)
  expect_equal(sim_params$n_groups, 2L)
  expect_equal(sim_params$n_samples_per_group, 15)  # floor(30 / 2)
  expect_equal(sim_params$n_samples, 30)
  expect_identical(sim_params$group_levels, c("ctrl", "case"))

  # Compositional constraints.
  expect_equal(sum(sim_params$mu_inv_softmax_base_realistic), 0, tolerance = 1e-12)
  expect_equal(sum(sim_params$slope_realistic), 0, tolerance = 1e-12)

  # Snapshot of the resampled vectors (set.seed(123) inside the function).
  expect_equal(
    sim_params$mu_inv_softmax_base_realistic,
    c( 0.251822202, -0.001165233,  0.251822202,
      -0.250656970, -0.250656970, -0.001165233),
    tolerance = 1e-7
  )
  expect_equal(
    sim_params$slope_realistic,
    c(0.1, 0.1, 0.1, -0.2, -0.2, 0.1),
    tolerance = 1e-12
  )

  # Scalar derived quantities.
  expect_equal(sim_params$intercept_disp_realistic,        1.775245,   tolerance = 1e-6)
  expect_equal(sim_params$sd_log_overdispersion_realistic, 0.4196737,  tolerance = 1e-6)
  expect_equal(sim_params$sigma_realistic,                 5.901729,   tolerance = 1e-6)

  # Replication block.
  expect_equal(sim_params$n_reps_auc, 4L)
  expect_equal(sim_params$perm_reps_auc, 199)
  expect_equal(sim_params$rep_seeds_auc, 12001:12004)
})

test_that("build_simulation_params_brito(): custom n_samples_per_group is propagated exactly", {
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

test_that("build_simulation_params_brito(): rejects bad group_levels and bad sample counts", {
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
