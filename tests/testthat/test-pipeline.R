# =============================================================================
# tests/testthat/test-pipeline.R - Tests for R/pipeline.R
# =============================================================================

# ---------------------------------------------------------------------------
# filter_by_abundance_cases
# ---------------------------------------------------------------------------

test_that("filter_by_abundance_cases(): threshold drops low-mean taxa", {
  count_long <- tibble::tribble(
    ~sample_id, ~taxon_id, ~count, ~library_size,
    "s1", "A", 500, 1000,   # mean prop = 0.5
    "s2", "A", 600, 1000,
    "s1", "B",   1, 1000,   # mean prop = 0.001
    "s2", "B",   2, 1000,
    "s1", "C", 100, 1000,   # mean prop = 0.1
    "s2", "C", 200, 1000
  )
  out <- filter_by_abundance_cases(count_long, threshold = 0.05)
  expect_setequal(unique(out$taxon_id), c("A", "C"))
})

test_that("filter_by_abundance_cases(): threshold = 0 keeps everything", {
  count_long <- tibble::tibble(
    sample_id = c("s1", "s1"), taxon_id = c("A", "B"),
    count = c(0, 0), library_size = c(100, 100)
  )
  out <- filter_by_abundance_cases(count_long, threshold = 0)
  # 0 / 100 = 0, and threshold = 0 keeps mean_prop >= 0 (both taxa).
  expect_setequal(unique(out$taxon_id), c("A", "B"))
})

test_that("filter_by_abundance_cases(): NULL threshold falls back to global `abundance_threshold`", {
  count_long <- tibble::tibble(
    sample_id = "s1", taxon_id = c("A", "B"),
    count = c(900, 1), library_size = c(1000, 1000)
  )
  # `filter_by_abundance_cases()` was defined in the global environment, so its
  # `abundance_threshold` lookup resolves via the search path (lexical scoping).
  # We set + clean up a global to match the convention used by the per-dataset
  # Quarto reports.
  assign("abundance_threshold", 0.1, envir = globalenv())
  on.exit(rm("abundance_threshold", envir = globalenv()), add = TRUE)

  out <- filter_by_abundance_cases(count_long)
  expect_setequal(unique(out$taxon_id), "A")
})

# ---------------------------------------------------------------------------
# extract_job_params
# ---------------------------------------------------------------------------

test_that("extract_job_params(): pulls (seed, n_da_taxa) from a data.frame row", {
  row <- data.frame(seed = 42, n_da_taxa = 5L)
  out <- extract_job_params(row)
  expect_equal(out, list(seed = 42, n_da_taxa = 5L))
})

test_that("extract_job_params(): handles nested column names like X1.seed", {
  row <- data.frame(X1.seed = 7, X2.n_da_taxa = 3L)
  out <- extract_job_params(row)
  expect_equal(out, list(seed = 7, n_da_taxa = 3L))
})

test_that("extract_job_params(): handles list-shaped input", {
  out <- extract_job_params(list(seed = 11, n_da_taxa = 2L))
  expect_equal(out, list(seed = 11, n_da_taxa = 2L))
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
# simulate_case_realistic_taxa
# ---------------------------------------------------------------------------

test_that("simulate_case_realistic_taxa(): builds expected sim_r with vector dispersion", {
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
  expect_true("count_long" %in% names(out))
  expect_equal(nrow(out$count_long), 10L * 4L)  # 2 groups * 5 samples * 4 taxa
  expect_equal(unname(rowSums(out$cohort_mu)), c(1, 1), tolerance = 1e-12)
})

test_that("simulate_case_realistic_taxa(): exactly one of vector / matrix dispersion required", {
  slope <- c(0.5, -0.5); mu <- c(0.5, -0.5)
  expect_error(
    .quietly({
      simulate_case_realistic_taxa(
        slope_vector = slope, mu_inv_softmax = mu,
        log_dispersion_assoc = 0.3,
        intercept_dispersion = NULL,
        intercept_dispersion_matrix = NULL,
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
# run_simulation_job_brito
# ---------------------------------------------------------------------------

test_that("run_simulation_job_brito(): runs end-to-end and tags seed/n_da_taxa", {
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
  expect_true(is.factor(sim_r$sample_metadata$group))
  expect_equal(levels(sim_r$sample_metadata$group), c("ctrl", "case"))
  expect_equal(nrow(sim_r$count_long), 16L * n_taxa)
})

test_that("run_simulation_job_brito(): wrong-length composition intercept errors", {
  expect_error(
    .quietly({
      run_simulation_job_brito(
        job_row = data.frame(seed = 1L, n_da_taxa = 0L),
        composition_intercept_by_taxon = c(0, 0, 0),  # length 3 != n_taxa = 4
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
# run_analysis_job
# ---------------------------------------------------------------------------

test_that("run_analysis_job(): default return -> long FPR-style tibble", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 8L, seed = 21L)
  methods <- c("bray", "jaccard")
  res <- run_analysis_job(
    sim_r,
    max_perm = 999L, batch = 999L, target_se = 1,
    standard_distance_methods = methods
  )
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("seed", "n_da_taxa", "method", "p_value"))
  expect_equal(nrow(res), length(methods))
  # Method labels are "pretty" (Bray-Curtis, Jaccard).
  expect_true(all(c("Bray\u2013Curtis", "Jaccard") %in% res$method))
})

test_that("run_analysis_job(): return_full = TRUE -> rich list with bd objects", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 8L, seed = 22L)
  methods <- c("bray", "jaccard")
  res <- run_analysis_job(
    sim_r,
    max_perm = 999L, batch = 999L, target_se = 1,
    return_full = TRUE,
    standard_distance_methods = methods
  )
  expect_named(res, c("method_results", "seed", "n_da_taxa"))
  expect_equal(length(res$method_results), length(methods))
  # Each method carries a betadisper fit and a p_value scalar.
  one <- res$method_results[[1]]
  expect_true(all(c("method_key", "bd", "p_value", "p_se", "distances", "group") %in% names(one)))
  expect_s3_class(one$bd, "betadisper")
})

test_that("run_analysis_job(): missing group factor errors out", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 4L, seed = 1L)
  sim_r$sample_metadata$group <- NULL
  expect_error(run_analysis_job(sim_r), regexp = "factor column .group.")
})

# ---------------------------------------------------------------------------
# lighten_sweep_analysis_result_for_fpr
# ---------------------------------------------------------------------------

test_that("lighten_sweep_analysis_result_for_fpr(): trims method_results to p_value-only", {
  fat <- list(
    method_results = list(
      m1 = list(method_key = "m1", bd = list("huge"), p_value = 0.04, p_se = 0.01,
                 n_perm = 999, distances = runif(50), group = factor(rep(c("a","b"), 25))),
      m2 = list(method_key = "m2", bd = list("huge"), p_value = 0.30, p_se = 0.02,
                 n_perm = 999, distances = runif(50), group = factor(rep(c("a","b"), 25)))
    ),
    seed = 99L,
    n_da_taxa = 3L,
    threshold = 0
  )
  lite <- lighten_sweep_analysis_result_for_fpr(fat)
  expect_named(lite, c("method_results", "seed", "n_da_taxa", "threshold"))
  expect_equal(names(lite$method_results), c("m1", "m2"))
  for (mr in lite$method_results) {
    expect_named(mr, "p_value")
  }
  expect_equal(lite$method_results$m1$p_value, 0.04)
  expect_equal(lite$method_results$m2$p_value, 0.30)
})

test_that("lighten_sweep_analysis_result_for_fpr(): pass-through for non-job lists", {
  x <- list(a = 1, b = "hello")
  expect_identical(lighten_sweep_analysis_result_for_fpr(x), x)
})

# ---------------------------------------------------------------------------
# summarise_fpr
# ---------------------------------------------------------------------------

test_that("summarise_fpr(): auto-bound tibble path", {
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
  expect_named(out, c("n_da_taxa", "method", "threshold", "n_sims", "n_fp", "FPR"))
  bray <- out[out$method == "Bray\u2013Curtis", ]
  expect_equal(bray$n_sims, 4L)
  expect_equal(bray$n_fp,   2L)  # 0.01 and 0.04
  expect_equal(bray$FPR,    0.5)
})

test_that("summarise_fpr(): list-of-jobs path", {
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
  m1 <- out[out$method == "m1", ]
  expect_equal(m1$n_sims, 3L)
  expect_equal(m1$n_fp,   2L)
  expect_equal(m1$FPR,    2 / 3)

  m2 <- out[out$method == "m2", ]
  expect_equal(m2$n_sims, 3L)
  expect_equal(m2$n_fp,   0L)
})

test_that("summarise_fpr(): empty input errors out", {
  expect_error(summarise_fpr(list()), regexp = "empty")
})
