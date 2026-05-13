# =============================================================================
# tests/testthat/test-betadisper.R - Tests for R/betadisper.R
# =============================================================================
#
# Strategy: distance computation is deterministic given the (fixed-seed)
# simulation input, so we snapshot the exact dist values and group factor.
# permutest p-values use an explicit set.seed() right before the call so the
# pinned p-value is reproducible.
# =============================================================================

# ---------------------------------------------------------------------------
# build_standard_betadisper_inputs_paper()
# ---------------------------------------------------------------------------

test_that("build_standard_betadisper_inputs_paper(): bray output structure and exact values", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 10L, seed = 4L)
  out <- build_standard_betadisper_inputs_paper(
    count_long = sim_r$count_long,
    sample_metadata = sim_r$sample_metadata,
    distance_method = "bray",
    group_levels = c("ctrl", "case")
  )
  expect_named(out, c("dist", "group"), ignore.order = FALSE)
  expect_s3_class(out$dist, "dist")
  expect_equal(attr(out$dist, "Size"), 20L)
  expect_equal(attr(out$dist, "Diag"), FALSE)
  expect_equal(attr(out$dist, "Upper"), FALSE)
  expect_true(all(is.finite(out$dist)))

  # Snapshot of the first six pairwise distances (deterministic for seed = 4).
  expect_equal(
    head(as.numeric(out$dist), 6),
    c(0.8295621, 0.4358563, 0.4234590, 0.4553347, 0.4000611, 0.3760439),
    tolerance = 1e-6
  )

  expect_s3_class(out$group, "factor")
  expect_equal(levels(out$group), c("ctrl", "case"))
  # 10 ctrl then 10 case.
  expect_equal(as.character(out$group),
               c(rep("ctrl", 10), rep("case", 10)))
})

test_that("build_standard_betadisper_inputs_paper(): all main distance methods produce a (Size, group) consistent (dist, factor) pair", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 6L, seed = 5L)
  # With this seed Sample_6 has total count = 0 -> the function drops it, so
  # the expected Size is 11 (1 of 12 samples is removed). This is the
  # documented "silently drop zero-total samples" behaviour; pin it.
  for (dm in c("bray", "aitchison", "robust.aitchison", "jaccard", "hellinger")) {
    out <- build_standard_betadisper_inputs_paper(
      count_long = sim_r$count_long,
      sample_metadata = sim_r$sample_metadata,
      distance_method = dm,
      group_levels = c("ctrl", "case")
    )
    expect_s3_class(out$dist, "dist")
    expect_equal(attr(out$dist, "Size"), 11L)
    expect_equal(length(out$group), 11L)
    expect_true(all(is.finite(out$dist)))
    expect_s3_class(out$group, "factor")
    expect_equal(levels(out$group), c("ctrl", "case"))
  }
})

test_that("build_standard_betadisper_inputs_paper(): silently drops zero-total samples", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 6L, seed = 5L)
  # Sanity: confirm Sample_6 truly has 0 total counts in this fixture.
  totals <- tapply(sim_r$count_long$count, sim_r$count_long$sample_id, sum)
  zero_samples <- names(totals[totals == 0])
  expect_equal(zero_samples, "Sample_6")

  out <- build_standard_betadisper_inputs_paper(
    count_long = sim_r$count_long,
    sample_metadata = sim_r$sample_metadata,
    distance_method = "bray",
    group_levels = c("ctrl", "case")
  )
  # Dropped sample contributes neither to dist nor to group factor.
  expect_equal(attr(out$dist, "Size"), 11L)
  # The retained group factor should reflect the surviving samples in order.
  survivors <- sim_r$sample_metadata$sample_id[!sim_r$sample_metadata$sample_id %in% zero_samples]
  expected_group_chr <- as.character(
    sim_r$sample_metadata$group[match(survivors, sim_r$sample_metadata$sample_id)]
  )
  expect_equal(as.character(out$group), expected_group_chr)
})

test_that("build_standard_betadisper_inputs_paper(): missing count_long columns -> error", {
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
})

test_that("build_standard_betadisper_inputs_paper(): missing metadata columns -> error", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 4L, seed = 1L)
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
  expect_error(
    build_standard_betadisper_inputs_paper(
      count_long = sim_r$count_long,
      sample_metadata = sim_r$sample_metadata,
      distance_method = "aitchison",
      group_levels = c("ctrl", "case"),
      clr_pseudocount = c(1, 2)
    ),
    regexp = "positive finite numeric scalar"
  )
})

test_that("build_standard_betadisper_inputs_paper(): hellinger uses euclidean on hellinger-transformed counts", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 4L, seed = 17L)
  out <- build_standard_betadisper_inputs_paper(
    count_long = sim_r$count_long,
    sample_metadata = sim_r$sample_metadata,
    distance_method = "hellinger",
    group_levels = c("ctrl", "case")
  )

  # Reconstruct the expected hellinger-Euclidean distance directly.
  wide <- sim_r$count_long %>%
    dplyr::select(dplyr::all_of(c("sample_id", "taxon_id", "count"))) %>%
    tidyr::pivot_wider(names_from = "taxon_id", values_from = "count") %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()
  wide <- wide[rowSums(wide) > 0, , drop = FALSE]
  hel <- vegan::decostand(wide, method = "hellinger", MARGIN = 1)
  expected <- stats::dist(hel, method = "euclidean")
  expect_equal(as.numeric(out$dist), as.numeric(expected), tolerance = 1e-12)
})

# ---------------------------------------------------------------------------
# permutest_betadisper_fixed()
# ---------------------------------------------------------------------------

test_that("permutest_betadisper_fixed(): exact list shape and pinned p_value", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 10L, seed = 4L)
  inp <- build_standard_betadisper_inputs_paper(
    count_long = sim_r$count_long,
    sample_metadata = sim_r$sample_metadata,
    distance_method = "bray",
    group_levels = c("ctrl", "case")
  )
  set.seed(101)
  bd <- vegan::betadisper(inp$dist, inp$group)
  set.seed(101)
  res <- permutest_betadisper_fixed(bd, permutations = 99)

  expect_named(res, c("p_value", "n_perm"), ignore.order = FALSE)
  expect_type(res$p_value, "double")
  expect_equal(res$n_perm, 99)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
  # Pinned snapshot under the explicit seed above.
  expect_equal(res$p_value, 0.43, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# permutest_betadisper_stable()
# ---------------------------------------------------------------------------

test_that("permutest_betadisper_stable(): exact list shape and one-batch stop", {
  sim_r <- .make_tiny_sim(n_taxa = 5L, n_samples_per_group = 10L, seed = 4L)
  inp <- build_standard_betadisper_inputs_paper(
    count_long = sim_r$count_long,
    sample_metadata = sim_r$sample_metadata,
    distance_method = "bray",
    group_levels = c("ctrl", "case")
  )
  set.seed(202)
  bd <- vegan::betadisper(inp$dist, inp$group)
  set.seed(202)
  # target_se = 1 forces termination after one batch (SE is always < 1).
  res <- permutest_betadisper_stable(bd, max_perm = 999, batch = 999, target_se = 1)
  expect_named(res, c("p_value", "p_se", "n_perm"), ignore.order = FALSE)
  expect_equal(res$n_perm, 999L)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
  # SE is sqrt(p(1-p)/n_perm).
  expect_equal(res$p_se, sqrt(res$p_value * (1 - res$p_value) / 999),
               tolerance = 1e-12)
})

test_that("permutest_betadisper_stable(): rejects bad arguments", {
  fake_bd <- list()  # validation runs first
  expect_error(permutest_betadisper_stable(fake_bd, max_perm = 100), regexp = "max_perm")
  expect_error(permutest_betadisper_stable(fake_bd, batch = 50),     regexp = "batch")
  expect_error(permutest_betadisper_stable(fake_bd, target_se = -1), regexp = "target_se")
  expect_error(permutest_betadisper_stable(fake_bd, target_se = 0),  regexp = "target_se")
})
