# =============================================================================
# tests/testthat/test-permanova.R - Tests for R/permanova.R
# =============================================================================

# ---------------------------------------------------------------------------
# run_permanova_job
# ---------------------------------------------------------------------------

test_that("run_permanova_job(): returns one row per requested distance method", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 8L, seed = 51L)
  methods <- c("bray", "jaccard", "hellinger")
  res <- run_permanova_job(
    sim_r = sim_r,
    threshold = 0,
    seed = 51L,
    n_da_taxa = 2L,
    permutations = 49,
    standard_distance_methods = methods
  )
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("seed", "n_da_taxa", "threshold", "method", "p_value"))
  expect_equal(nrow(res), length(methods))
  expect_equal(res$method, paste0("permanova_", methods))
  expect_equal(unique(res$seed), 51)
  expect_equal(unique(res$n_da_taxa), 2L)
  expect_equal(unique(res$threshold), 0)
  expect_true(all(is.na(res$p_value) | (res$p_value >= 0 & res$p_value <= 1)))
})

test_that("run_permanova_job(): aggressive threshold leaves <3 taxa -> empty tibble", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 8L, seed = 51L)
  # Threshold = 1 forces all taxa to drop (mean proportion < 1 always).
  res <- run_permanova_job(
    sim_r = sim_r,
    threshold = 1,
    seed = 51L,
    n_da_taxa = 2L,
    permutations = 49
  )
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("seed", "n_da_taxa", "threshold", "method", "p_value"))
  expect_equal(nrow(res), 0L)
})

test_that("run_permanova_job(): non-factor group errors out", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 6L, seed = 1L)
  sim_r$sample_metadata$group <- as.character(sim_r$sample_metadata$group)
  expect_error(
    run_permanova_job(sim_r, threshold = 0, seed = 1L, n_da_taxa = 1L, permutations = 9),
    regexp = "must be a factor"
  )
})

# ---------------------------------------------------------------------------
# summarise_permanova_tpr
# ---------------------------------------------------------------------------

test_that("summarise_permanova_tpr(): computes TPR by (n_da_taxa, method, threshold)", {
  permanova_results <- tibble::tribble(
    ~seed, ~n_da_taxa, ~threshold, ~method,            ~p_value,
       1L,         2L,          0, "permanova_bray",   0.01,
       2L,         2L,          0, "permanova_bray",   0.20,
       3L,         2L,          0, "permanova_bray",   0.04,
       4L,         2L,          0, "permanova_bray",   NA_real_,
       1L,         2L,          0, "permanova_jaccard",0.30,
       2L,         2L,          0, "permanova_jaccard",0.40,
       1L,         5L,          0, "permanova_bray",   0.001,
       2L,         5L,          0, "permanova_bray",   0.002,
       # n_da_taxa == 0 rows must be EXCLUDED from TPR.
       1L,         0L,          0, "permanova_bray",   0.5,
       2L,         0L,          0, "permanova_bray",   0.5
  )
  out <- summarise_permanova_tpr(permanova_results)
  expect_named(out, c("n_da_taxa", "method", "threshold", "n_sims", "n_tp", "TPR"))
  # Null setting (n_da_taxa == 0) dropped.
  expect_false(0 %in% out$n_da_taxa)

  bray_2 <- out[out$n_da_taxa == 2 & out$method == "permanova_bray", ]
  expect_equal(bray_2$n_sims, 3L)  # NA excluded
  expect_equal(bray_2$n_tp,   2L)  # p < 0.05 -> 0.01 and 0.04
  expect_equal(bray_2$TPR,    2 / 3)

  jac_2 <- out[out$n_da_taxa == 2 & out$method == "permanova_jaccard", ]
  expect_equal(jac_2$n_sims, 2L)
  expect_equal(jac_2$n_tp,   0L)
  expect_equal(jac_2$TPR,    0)

  bray_5 <- out[out$n_da_taxa == 5 & out$method == "permanova_bray", ]
  expect_equal(bray_5$n_tp, 2L)
  expect_equal(bray_5$TPR,  1)
})

test_that("summarise_permanova_tpr(): n_sims == 0 yields NA TPR", {
  # Only NA p-values for one stratum.
  permanova_results <- tibble::tibble(
    seed = 1:2, n_da_taxa = c(3L, 3L), threshold = c(0, 0),
    method = c("permanova_bray", "permanova_bray"),
    p_value = c(NA_real_, NA_real_)
  )
  out <- summarise_permanova_tpr(permanova_results)
  expect_equal(out$n_sims, 0L)
  expect_true(is.na(out$TPR))
})
