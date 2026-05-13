# =============================================================================
# tests/testthat/test-permanova.R - Tests for R/permanova.R
# =============================================================================
#
# Strategy:
#   * Pin column names AND order of returned tibbles.
#   * Pin p-values via explicit set.seed() right before adonis2 is called.
#   * Pin summarise_permanova_tpr arithmetic via hand-calculated values.
# =============================================================================

# ---------------------------------------------------------------------------
# run_permanova_job()
# ---------------------------------------------------------------------------

test_that("run_permanova_job(): tibble shape, exact methods, deterministic p-values", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 8L, seed = 51L)
  methods <- c("bray", "jaccard", "hellinger")
  set.seed(2025)
  res <- run_permanova_job(
    sim_r = sim_r,
    threshold = 0,
    seed = 51L,
    n_da_taxa = 2L,
    permutations = 199,
    standard_distance_methods = methods
  )
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("seed", "n_da_taxa", "threshold", "method", "p_value"),
               ignore.order = FALSE)
  expect_equal(nrow(res), 3L)
  expect_equal(res$method, paste0("permanova_", methods))
  expect_equal(res$seed,     rep(51, 3))
  expect_equal(res$n_da_taxa, rep(2L, 3))
  expect_equal(res$threshold, rep(0, 3))
  # Pinned p-values (set.seed(2025) above + adonis2 internal RNG).
  expect_equal(res$p_value, c(0.445, 0.530, 0.570), tolerance = 1e-3)
})

test_that("run_permanova_job(): empty tibble (correct columns) when <3 taxa survive filtering", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 8L, seed = 51L)
  res <- run_permanova_job(
    sim_r = sim_r,
    threshold = 1,     # impossible -> drop all taxa
    seed = 51L,
    n_da_taxa = 2L,
    permutations = 49
  )
  expect_s3_class(res, "tbl_df")
  expect_named(res, c("seed", "n_da_taxa", "threshold", "method", "p_value"),
               ignore.order = FALSE)
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
# summarise_permanova_tpr()
# ---------------------------------------------------------------------------

test_that("summarise_permanova_tpr(): exact tibble, drops n_da_taxa == 0 rows", {
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
       1L,         0L,          0, "permanova_bray",   0.5,
       2L,         0L,          0, "permanova_bray",   0.5
  )
  out <- summarise_permanova_tpr(permanova_results)

  # Structure.
  expect_s3_class(out, "tbl_df")
  expect_named(out,
               c("n_da_taxa", "method", "threshold", "n_sims", "n_tp", "TPR"),
               ignore.order = FALSE)

  # Null rows excluded.
  expect_false(0 %in% out$n_da_taxa)

  # Exact rows (sorted by group_by keys: n_da_taxa, method, threshold).
  # Expected:
  #   n_da_taxa=2, bray:    n_sims=3, n_tp=2, TPR=2/3
  #   n_da_taxa=2, jaccard: n_sims=2, n_tp=0, TPR=0
  #   n_da_taxa=5, bray:    n_sims=2, n_tp=2, TPR=1
  expect_equal(nrow(out), 3L)

  bray_2 <- out[out$n_da_taxa == 2L & out$method == "permanova_bray", ]
  expect_equal(bray_2$n_sims, 3L)
  expect_equal(bray_2$n_tp,   2L)
  expect_equal(bray_2$TPR,    2 / 3)

  jac_2 <- out[out$n_da_taxa == 2L & out$method == "permanova_jaccard", ]
  expect_equal(jac_2$n_sims, 2L)
  expect_equal(jac_2$n_tp,   0L)
  expect_equal(jac_2$TPR,    0)

  bray_5 <- out[out$n_da_taxa == 5L & out$method == "permanova_bray", ]
  expect_equal(bray_5$n_sims, 2L)
  expect_equal(bray_5$n_tp,   2L)
  expect_equal(bray_5$TPR,    1)
})

test_that("summarise_permanova_tpr(): n_sims == 0 yields NA TPR", {
  permanova_results <- tibble::tibble(
    seed = 1:2, n_da_taxa = c(3L, 3L), threshold = c(0, 0),
    method = c("permanova_bray", "permanova_bray"),
    p_value = c(NA_real_, NA_real_)
  )
  out <- summarise_permanova_tpr(permanova_results)
  expect_equal(out$n_sims, 0L)
  expect_equal(out$n_tp, 0L)
  expect_true(is.na(out$TPR))
})

test_that("summarise_permanova_tpr(): preserves multiple thresholds and methods", {
  permanova_results <- tibble::tribble(
    ~seed, ~n_da_taxa, ~threshold, ~method,         ~p_value,
       1L,         3L,        0.0, "permanova_bray", 0.01,
       2L,         3L,        0.0, "permanova_bray", 0.02,
       1L,         3L,        0.1, "permanova_bray", 0.40,
       2L,         3L,        0.1, "permanova_bray", 0.50
  )
  out <- summarise_permanova_tpr(permanova_results)
  expect_equal(nrow(out), 2L)
  thr0 <- out[out$threshold == 0.0, ]
  thr1 <- out[out$threshold == 0.1, ]
  expect_equal(thr0$TPR, 1)
  expect_equal(thr1$TPR, 0)
})
