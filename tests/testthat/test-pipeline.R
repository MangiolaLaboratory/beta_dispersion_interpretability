# =============================================================================
# tests/testthat/test-pipeline.R - Tests for R/pipeline.R
# =============================================================================
#
# Strategy: every output (tibbles, lists, named vectors) is pinned by NAME,
# ORDER, TYPE and VALUE wherever possible. Stochastic outputs (PERMANOVA-style
# p-values inside run_analysis_job) are made reproducible with set.seed() and
# their values snapshotted.
# =============================================================================

# ---------------------------------------------------------------------------
# filter_by_abundance_cases()
# ---------------------------------------------------------------------------

test_that("filter_by_abundance_cases(): exact filtered tibble snapshot", {
  count_long <- tibble::tribble(
    ~sample_id, ~taxon_id, ~count, ~library_size,
    "s1", "A", 500, 1000,
    "s2", "A", 600, 1000,
    "s1", "B",   1, 1000,
    "s2", "B",   2, 1000,
    "s1", "C", 100, 1000,
    "s2", "C", 200, 1000
  )
  out <- filter_by_abundance_cases(count_long, threshold = 0.05)
  expected <- tibble::tribble(
    ~sample_id, ~taxon_id, ~count, ~library_size,
    "s1", "A", 500, 1000,
    "s2", "A", 600, 1000,
    "s1", "C", 100, 1000,
    "s2", "C", 200, 1000
  )
  expect_equal(out, expected)
})

test_that("filter_by_abundance_cases(): threshold = 0 keeps every row (column-perfect)", {
  count_long <- tibble::tibble(
    sample_id = c("s1", "s1"), taxon_id = c("A", "B"),
    count = c(0, 0), library_size = c(100, 100)
  )
  out <- filter_by_abundance_cases(count_long, threshold = 0)
  expect_equal(out, count_long)
})

test_that("filter_by_abundance_cases(): drops every taxon with mean_prop < threshold", {
  count_long <- tibble::tribble(
    ~sample_id, ~taxon_id, ~count, ~library_size,
    "s1", "X",  1,  1000,    # 0.001
    "s2", "X",  2,  1000,    # 0.002 -> mean 0.0015 < 0.01
    "s1", "Y", 50,  1000,    # 0.050
    "s2", "Y", 80,  1000     # 0.080 -> mean 0.065 >= 0.01
  )
  out <- filter_by_abundance_cases(count_long, threshold = 0.01)
  expect_equal(unique(out$taxon_id), "Y")
})

test_that("filter_by_abundance_cases(): NULL threshold falls back to global abundance_threshold", {
  count_long <- tibble::tibble(
    sample_id = "s1", taxon_id = c("A", "B"),
    count = c(900, 1), library_size = c(1000, 1000)
  )
  assign("abundance_threshold", 0.1, envir = globalenv())
  on.exit(rm("abundance_threshold", envir = globalenv()), add = TRUE)
  out <- filter_by_abundance_cases(count_long)
  expect_equal(unique(out$taxon_id), "A")
})

# ---------------------------------------------------------------------------
# extract_job_params()
# ---------------------------------------------------------------------------

test_that("extract_job_params(): data.frame row -> exact list", {
  expect_identical(
    extract_job_params(data.frame(seed = 42, n_da_taxa = 5L)),
    list(seed = 42, n_da_taxa = 5L)
  )
})

test_that("extract_job_params(): list input -> exact list", {
  expect_identical(
    extract_job_params(list(seed = 11, n_da_taxa = 2L)),
    list(seed = 11, n_da_taxa = 2L)
  )
})

test_that("extract_job_params(): tibble row -> exact list", {
  expect_identical(
    extract_job_params(tibble::tibble(seed = 4, n_da_taxa = 1L)),
    list(seed = 4, n_da_taxa = 1L)
  )
})

test_that("extract_job_params(): rejects NULL / empty / missing columns", {
  expect_error(extract_job_params(NULL), regexp = "NULL")
  expect_error(extract_job_params(data.frame()), regexp = "empty data.frame")
  expect_error(extract_job_params(data.frame(other = 1)),
               regexp = "missing 'seed' column|missing 'n_da_taxa'")
})

test_that("extract_job_params(): rejects NA seed or n_da_taxa", {
  expect_error(extract_job_params(data.frame(seed = NA_real_, n_da_taxa = 1L)),
               regexp = "seed")
  expect_error(extract_job_params(data.frame(seed = 1, n_da_taxa = NA_integer_)),
               regexp = "n_da_taxa")
})

# ---------------------------------------------------------------------------
# simulate_case_realistic_taxa()
# ---------------------------------------------------------------------------

test_that("simulate_case_realistic_taxa(): vector dispersion -> expected sim list", {
  slope <- c(0.5, -0.5, 0, 0)
  mu <- c(-0.5, 0.5, -0.3, 0.3); mu <- mu - mean(mu)
  out <- .quietly({
    simulate_case_realistic_taxa(
      slope_vector = slope,
      mu_inv_softmax = mu,
      log_dispersion_assoc = 0.3,
      intercept_dispersion = rep(1.5, 4L),
      sd_log_overdispersion = 0.1,
      seed = 1L,
      n_taxa = 4L,
      n_samples_per_group = 5L,
      library_size_mean = 2000,
      library_size_sd = 200
    )
  })
  expect_equal(nrow(out$count_long), 10L * 4L)
  expect_equal(unname(rowSums(out$cohort_mu)), c(1, 1), tolerance = 1e-12)
  expect_equal(dim(out$cohort_mu), c(2L, 4L))
  # Two-group design built internally.
  expect_equal(out$n_cohorts, 2L)
  expect_equal(unname(out$unique_design[, "Group"]), c(0, 1))
})

test_that("simulate_case_realistic_taxa(): mutually exclusive dispersion inputs", {
  slope <- c(0.5, -0.5); mu <- c(0.5, -0.5)
  expect_error(
    .quietly({
      simulate_case_realistic_taxa(
        slope_vector = slope, mu_inv_softmax = mu,
        log_dispersion_assoc = 0.3,
        intercept_dispersion = NULL, intercept_dispersion_matrix = NULL,
        sd_log_overdispersion = 0.1,
        seed = 1L, n_taxa = 2L, n_samples_per_group = 4L,
        library_size_mean = 1000, library_size_sd = 100
      )
    }),
    regexp = "exactly one"
  )
  expect_error(
    .quietly({
      simulate_case_realistic_taxa(
        slope_vector = slope, mu_inv_softmax = mu,
        log_dispersion_assoc = 0.3,
        intercept_dispersion = c(1.5, 1.5),
        intercept_dispersion_matrix = matrix(1.5, nrow = 2, ncol = 2),
        sd_log_overdispersion = 0.1,
        seed = 1L, n_taxa = 2L, n_samples_per_group = 4L,
        library_size_mean = 1000, library_size_sd = 100
      )
    }),
    regexp = "exactly one"
  )
})

# ---------------------------------------------------------------------------
# run_simulation_job_brito()
# ---------------------------------------------------------------------------

test_that("run_simulation_job_brito(): seeded sim has exact counts and library sizes", {
  n_taxa <- 5L
  mu <- seq(-1, 1, length.out = n_taxa); mu <- mu - mean(mu)
  slopes_pool <- c(-0.6, -0.4, 0.4, 0.6, -0.2, 0.2)
  job_row <- data.frame(seed = 7L, n_da_taxa = 2L)
  sim_r <- .quietly({
    run_simulation_job_brito(
      job_row = job_row,
      composition_intercept_by_taxon = mu,
      dispersion_intercept_by_taxon = 1.5,
      mean_dispersion_assoc_slope = 0.3,
      slope_effects_distribution = slopes_pool,
      sd_log_overdispersion = 0.1,
      n_taxa = n_taxa,
      n_samples = 16L,
      n_samples_per_group = 8L,
      library_size_mean = 2000,
      library_size_sd = 200,
      group_levels = c("ctrl", "case")
    )
  })

  expect_equal(sim_r$seed, 7L)
  expect_equal(sim_r$n_da_taxa, 2L)
  expect_s3_class(sim_r$sample_metadata$group, "factor")
  expect_equal(levels(sim_r$sample_metadata$group), c("ctrl", "case"))
  expect_equal(nrow(sim_r$count_long), 16L * n_taxa)

  # Pinned snapshot from seed = 7 (count_long pivots by sample_id lexically).
  expect_equal(head(sim_r$count_long$count, 8L),
               c(0, 15, 0, 2029, 0, 0, 0, 2429))
  expect_equal(
    sim_r$sample_metadata$library_size,
    c(1811, 2150, 1977, 2031, 2438, 2071, 2543, 2456,
      2065, 2379, 2094, 1821, 1939, 1999, 2198, 2168)
  )
  # cohort_mu pinned.
  expect_equal(
    sim_r$cohort_mu[1, ],
    c(Taxon_1 = 0.05801222, Taxon_2 = 0.09564598,
      Taxon_3 = 0.15769358, Taxon_4 = 0.25999273,
      Taxon_5 = 0.42865549),
    tolerance = 1e-7
  )
  expect_equal(
    sim_r$cohort_mu[2, ],
    c(Taxon_1 = 0.08130150, Taxon_2 = 0.10974554,
      Taxon_3 = 0.18093977, Taxon_4 = 0.29831934,
      Taxon_5 = 0.32969384),
    tolerance = 1e-7
  )
})

test_that("run_simulation_job_brito(): no group_levels -> no factor populated", {
  n_taxa <- 4L
  mu <- rep(0, n_taxa)
  sim_r <- .quietly({
    run_simulation_job_brito(
      job_row = data.frame(seed = 3L, n_da_taxa = 0L),
      composition_intercept_by_taxon = mu,
      dispersion_intercept_by_taxon = 1.0,
      mean_dispersion_assoc_slope = 0.0,
      slope_effects_distribution = c(-0.1, 0.1),
      sd_log_overdispersion = 0.05,
      n_taxa = n_taxa, n_samples = 8L, n_samples_per_group = 4L,
      library_size_mean = 1000, library_size_sd = 100,
      group_levels = NULL
    )
  })
  expect_null(sim_r$sample_metadata$group)
  expect_equal(sim_r$seed, 3L)
  expect_equal(sim_r$n_da_taxa, 0L)
})

test_that("run_simulation_job_brito(): wrong-length composition intercept errors", {
  expect_error(
    .quietly({
      run_simulation_job_brito(
        job_row = data.frame(seed = 1L, n_da_taxa = 0L),
        composition_intercept_by_taxon = c(0, 0, 0),
        dispersion_intercept_by_taxon = 1.5,
        mean_dispersion_assoc_slope = 0.3,
        slope_effects_distribution = c(-0.3, 0.3),
        sd_log_overdispersion = 0.1,
        n_taxa = 4L, n_samples = 8L, n_samples_per_group = 4L,
        library_size_mean = 1000, library_size_sd = 100
      )
    }),
    regexp = "composition_intercept_by_taxon"
  )
})

test_that("run_simulation_job_brito(): empty slope_effects_distribution errors", {
  expect_error(
    .quietly({
      run_simulation_job_brito(
        job_row = data.frame(seed = 1L, n_da_taxa = 1L),
        composition_intercept_by_taxon = rep(0, 4L),
        dispersion_intercept_by_taxon = 1.5,
        mean_dispersion_assoc_slope = 0.3,
        slope_effects_distribution = numeric(0),
        sd_log_overdispersion = 0.1,
        n_taxa = 4L, n_samples = 8L, n_samples_per_group = 4L,
        library_size_mean = 1000, library_size_sd = 100
      )
    }),
    regexp = "empty"
  )
})

# ---------------------------------------------------------------------------
# run_analysis_job()
# ---------------------------------------------------------------------------

test_that("run_analysis_job(): default tibble snapshot with pinned p-values", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 8L, seed = 21L)
  set.seed(2026)
  res <- run_analysis_job(
    sim_r,
    max_perm = 999L, batch = 999L, target_se = 1,
    standard_distance_methods = c("bray", "jaccard")
  )
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("seed", "n_da_taxa", "method", "p_value"),
               ignore.order = FALSE)
  expect_equal(res$seed,      c(21L, 21L))
  expect_equal(res$n_da_taxa, c(2L, 2L))
  # Pretty method labels.
  expect_equal(res$method, c("Bray\u2013Curtis", "Jaccard"))
  # Pinned p-values (set.seed(2026) above + permutest internal RNG).
  expect_equal(res$p_value, c(0.865, 0.663), tolerance = 1e-3)
})

test_that("run_analysis_job(): return_full = TRUE -> full per-method list", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 8L, seed = 22L)
  set.seed(2027)
  res <- run_analysis_job(
    sim_r,
    max_perm = 999L, batch = 999L, target_se = 1,
    return_full = TRUE,
    standard_distance_methods = c("bray", "jaccard")
  )
  expect_named(res, c("method_results", "seed", "n_da_taxa"),
               ignore.order = FALSE)
  expect_equal(res$seed, 22L)
  expect_equal(res$n_da_taxa, 2L)
  expect_equal(names(res$method_results),
               c("permdisp_bray", "permdisp_jaccard"))
  one <- res$method_results[[1]]
  expect_named(one,
               c("method_key", "bd", "group1", "group2", "delta",
                 "p_value", "p_se", "n_perm", "distances", "group"),
               ignore.order = FALSE)
  expect_s3_class(one$bd, "betadisper")
  expect_equal(one$delta, one$group2 - one$group1, tolerance = 1e-12)
  expect_length(one$distances, 16L)
  expect_s3_class(one$group, "factor")
  expect_equal(levels(one$group), c("ctrl", "case"))
})

test_that("run_analysis_job(): missing group factor errors out", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 4L, seed = 1L)
  sim_r$sample_metadata$group <- NULL
  expect_error(run_analysis_job(sim_r), regexp = "factor column .group.")
})

test_that("run_analysis_job(): non-factor group errors out", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 4L, seed = 1L)
  sim_r$sample_metadata$group <- as.character(sim_r$sample_metadata$group)
  expect_error(run_analysis_job(sim_r), regexp = "must be a factor")
})

# ---------------------------------------------------------------------------
# lighten_sweep_analysis_result_for_fpr()
# ---------------------------------------------------------------------------

test_that("lighten_sweep_analysis_result_for_fpr(): exact trimmed structure", {
  fat <- list(
    method_results = list(
      m1 = list(method_key = "m1", bd = list("huge"), p_value = 0.04, p_se = 0.01,
                 n_perm = 999, distances = 1:50,
                 group = factor(rep(c("a", "b"), 25))),
      m2 = list(method_key = "m2", bd = list("huge"), p_value = 0.30, p_se = 0.02,
                 n_perm = 999, distances = 1:50,
                 group = factor(rep(c("a", "b"), 25)))
    ),
    seed = 99L,
    n_da_taxa = 3L,
    threshold = 0
  )
  lite <- lighten_sweep_analysis_result_for_fpr(fat)
  expect_named(lite, c("method_results", "seed", "n_da_taxa", "threshold"),
               ignore.order = FALSE)
  expect_equal(names(lite$method_results), c("m1", "m2"))
  for (mr in lite$method_results) {
    expect_named(mr, "p_value")
  }
  expect_equal(lite$method_results$m1, list(p_value = 0.04))
  expect_equal(lite$method_results$m2, list(p_value = 0.30))
  expect_equal(lite$seed, 99L)
  expect_equal(lite$n_da_taxa, 3L)
  expect_equal(lite$threshold, 0)
})

test_that("lighten_sweep_analysis_result_for_fpr(): omits threshold when absent", {
  fat <- list(
    method_results = list(m1 = list(p_value = 0.1)),
    seed = 1L,
    n_da_taxa = 0L
  )
  lite <- lighten_sweep_analysis_result_for_fpr(fat)
  expect_named(lite, c("method_results", "seed", "n_da_taxa"),
               ignore.order = FALSE)
  expect_false("threshold" %in% names(lite))
})

test_that("lighten_sweep_analysis_result_for_fpr(): pass-through for unrelated lists", {
  x <- list(a = 1, b = "hello")
  expect_identical(lighten_sweep_analysis_result_for_fpr(x), x)
})

# ---------------------------------------------------------------------------
# summarise_fpr()
# ---------------------------------------------------------------------------

test_that("summarise_fpr(): tibble input -> exact aggregated tibble", {
  results <- tibble::tribble(
    ~seed, ~n_da_taxa, ~threshold, ~method,            ~p_value,
       1L,         0L,          0, "Bray\u2013Curtis", 0.01,
       2L,         0L,          0, "Bray\u2013Curtis", 0.20,
       3L,         0L,          0, "Bray\u2013Curtis", 0.04,
       4L,         0L,          0, "Bray\u2013Curtis", 0.50,
       1L,         0L,          0, "Jaccard",          0.30,
       2L,         0L,          0, "Jaccard",          0.40
  )
  out <- summarise_fpr(results)
  expect_named(out, c("n_da_taxa", "method", "threshold", "n_sims", "n_fp", "FPR"),
               ignore.order = FALSE)
  bray <- out[out$method == "Bray\u2013Curtis", ]
  expect_equal(bray$n_sims, 4L)
  expect_equal(bray$n_fp,   2L)
  expect_equal(bray$FPR,    0.5)
  jac <- out[out$method == "Jaccard", ]
  expect_equal(jac$n_sims, 2L)
  expect_equal(jac$n_fp,   0L)
  expect_equal(jac$FPR,    0)
})

test_that("summarise_fpr(): list-of-jobs input -> exact aggregated tibble with NA threshold", {
  jobs <- list(
    list(method_results = list(m1 = list(p_value = 0.01),
                                m2 = list(p_value = 0.20)),
         seed = 1L, n_da_taxa = 0L),
    list(method_results = list(m1 = list(p_value = 0.04),
                                m2 = list(p_value = 0.40)),
         seed = 2L, n_da_taxa = 0L),
    list(method_results = list(m1 = list(p_value = 0.10),
                                m2 = list(p_value = 0.60)),
         seed = 3L, n_da_taxa = 0L)
  )
  out <- summarise_fpr(jobs)
  expect_named(out, c("n_da_taxa", "method", "threshold", "n_sims", "n_fp", "FPR"),
               ignore.order = FALSE)
  # threshold is NA when the per-job list lacks it.
  expect_true(all(is.na(out$threshold)))

  m1 <- out[out$method == "m1", ]
  expect_equal(m1$n_sims, 3L)
  expect_equal(m1$n_fp,   2L)
  expect_equal(m1$FPR,    2 / 3)

  m2 <- out[out$method == "m2", ]
  expect_equal(m2$n_sims, 3L)
  expect_equal(m2$n_fp,   0L)
  expect_equal(m2$FPR,    0)
})

test_that("summarise_fpr(): preserves threshold when present per-job", {
  jobs <- list(
    list(method_results = list(m1 = list(p_value = 0.01)),
         seed = 1L, n_da_taxa = 0L, threshold = 0.1),
    list(method_results = list(m1 = list(p_value = 0.20)),
         seed = 2L, n_da_taxa = 0L, threshold = 0.1)
  )
  out <- summarise_fpr(jobs)
  expect_equal(out$threshold, 0.1)
  expect_equal(out$n_sims, 2L)
  expect_equal(out$n_fp, 1L)
  expect_equal(out$FPR, 0.5)
})

test_that("summarise_fpr(): empty input errors out", {
  expect_error(summarise_fpr(list()), regexp = "empty")
})
