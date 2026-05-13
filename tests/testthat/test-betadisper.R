# =============================================================================
# tests/testthat/test-betadisper.R - Tests for R/betadisper.R
# =============================================================================

# ---------------------------------------------------------------------------
# build_standard_betadisper_inputs_paper
# ---------------------------------------------------------------------------

test_that("build_standard_betadisper_inputs_paper(): bray returns aligned (dist, group)", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 8L, seed = 3L)
  out <- build_standard_betadisper_inputs_paper(
    count_long = sim_r$count_long,
    sample_metadata = sim_r$sample_metadata,
    distance_method = "bray",
    group_levels = c("ctrl", "case")
  )
  expect_named(out, c("dist", "group"))
  expect_s3_class(out$dist, "dist")
  expect_s3_class(out$group, "factor")
  expect_equal(levels(out$group), c("ctrl", "case"))
  # Group length matches dist's underlying matrix dimension.
  expect_equal(length(out$group), attr(out$dist, "Size"))
  expect_true(all(is.finite(out$dist)))
})

test_that("build_standard_betadisper_inputs_paper(): all main distance methods work", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 6L, seed = 5L)
  for (dm in c("bray", "aitchison", "robust.aitchison", "jaccard", "hellinger")) {
    out <- build_standard_betadisper_inputs_paper(
      count_long = sim_r$count_long,
      sample_metadata = sim_r$sample_metadata,
      distance_method = dm,
      group_levels = c("ctrl", "case")
    )
    expect_s3_class(out$dist, "dist")
    expect_true(all(is.finite(out$dist)))
    expect_s3_class(out$group, "factor")
  }
})

test_that("build_standard_betadisper_inputs_paper(): missing required columns errors", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 4L, seed = 1L)
  bad_count_long <- sim_r$count_long[, setdiff(colnames(sim_r$count_long), "count"), drop = FALSE]
  expect_error(
    build_standard_betadisper_inputs_paper(
      count_long = bad_count_long,
      sample_metadata = sim_r$sample_metadata,
      distance_method = "bray",
      group_levels = c("ctrl", "case")
    ),
    regexp = "missing columns.*count"
  )

  bad_meta <- sim_r$sample_metadata[, "sample_id", drop = FALSE]
  expect_error(
    build_standard_betadisper_inputs_paper(
      count_long = sim_r$count_long,
      sample_metadata = bad_meta,
      distance_method = "bray",
      group_levels = c("ctrl", "case")
    ),
    regexp = "sample_metadata must contain"
  )
})

test_that("build_standard_betadisper_inputs_paper(): aitchison rejects bad pseudocount", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 4L, seed = 1L)
  expect_error(
    build_standard_betadisper_inputs_paper(
      count_long = sim_r$count_long,
      sample_metadata = sim_r$sample_metadata,
      distance_method = "aitchison",
      group_levels = c("ctrl", "case"),
      clr_pseudocount = -1
    ),
    regexp = "positive finite numeric scalar"
  )
})

# ---------------------------------------------------------------------------
# permutest_betadisper_fixed
# ---------------------------------------------------------------------------

test_that("permutest_betadisper_fixed(): returns finite p_value and echoes n_perm", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 10L, seed = 4L)
  inp <- build_standard_betadisper_inputs_paper(
    count_long = sim_r$count_long,
    sample_metadata = sim_r$sample_metadata,
    distance_method = "bray",
    group_levels = c("ctrl", "case")
  )
  bd <- vegan::betadisper(inp$dist, inp$group)
  set.seed(101)
  res <- permutest_betadisper_fixed(bd, permutations = 99)
  expect_named(res, c("p_value", "n_perm"))
  expect_equal(res$n_perm, 99)
  expect_true(is.finite(res$p_value))
  expect_true(res$p_value >= 0 && res$p_value <= 1)
})

# ---------------------------------------------------------------------------
# permutest_betadisper_stable
# ---------------------------------------------------------------------------

test_that("permutest_betadisper_stable(): returns p_value, p_se, n_perm", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 10L, seed = 4L)
  inp <- build_standard_betadisper_inputs_paper(
    count_long = sim_r$count_long,
    sample_metadata = sim_r$sample_metadata,
    distance_method = "bray",
    group_levels = c("ctrl", "case")
  )
  bd <- vegan::betadisper(inp$dist, inp$group)
  set.seed(202)
  # Use a small max_perm so the test is fast. batch must be >= 100.
  res <- permutest_betadisper_stable(bd, max_perm = 999, batch = 999, target_se = 1)
  expect_named(res, c("p_value", "p_se", "n_perm"))
  expect_true(is.finite(res$p_value))
  expect_true(res$p_value >= 0 && res$p_value <= 1)
  expect_true(is.finite(res$p_se) && res$p_se > 0)
  expect_equal(res$n_perm, 999L)
})

test_that("permutest_betadisper_stable(): rejects bad arguments", {
  fake_bd <- list()  # any object, validation happens before use
  expect_error(permutest_betadisper_stable(fake_bd, max_perm = 100), regexp = "max_perm")
  expect_error(permutest_betadisper_stable(fake_bd, batch = 50), regexp = "batch")
  expect_error(permutest_betadisper_stable(fake_bd, target_se = -1), regexp = "target_se")
})
