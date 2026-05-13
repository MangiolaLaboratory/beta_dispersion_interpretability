# =============================================================================
# tests/testthat/test-alpha_diversity.R - Tests for R/alpha_diversity.R
# =============================================================================

# ---------------------------------------------------------------------------
# run_alpha_diversity_analysis (no rarefaction)
# ---------------------------------------------------------------------------

test_that("run_alpha_diversity_analysis(): basic happy path", {
  sim_r <- .make_tiny_sim(n_taxa = 6L, n_samples_per_group = 10L, seed = 12L)
  res <- run_alpha_diversity_analysis(
    sim_result = sim_r,
    case_label = "Test",
    group_levels = c("ctrl", "case"),
    group_colors = c(ctrl = "#377eb8", case = "#e41a1c")
  )

  expect_named(res, c("data", "summary", "tests", "plot"))
  expect_s3_class(res$data, "data.frame")
  expect_true(all(c("sample_id", "group", "richness", "shannon", "pielou_evenness") %in% colnames(res$data)))

  # One row per sample.
  expect_equal(nrow(res$data), 20L)
  # Richness in [0, n_taxa]; Shannon >= 0 and Pielou in [0, 1] (or NA when richness==1).
  expect_true(all(res$data$richness >= 0 & res$data$richness <= 6L))
  expect_true(all(res$data$shannon >= 0))
  ok <- is.finite(res$data$pielou_evenness)
  expect_true(all(res$data$pielou_evenness[ok] >= 0 - 1e-9))
  expect_true(all(res$data$pielou_evenness[ok] <= 1 + 1e-9))

  # Summary has one row per level.
  expect_equal(sort(as.character(res$summary$group)), sort(c("ctrl", "case")))

  # Tests carry t.test-shaped output.
  expect_true(all(c("richness", "shannon", "evenness") %in% names(res$tests)))

  # Plot is a ggplot.
  expect_s3_class(res$plot, "ggplot")
})

test_that("run_alpha_diversity_analysis(): missing group_levels/colors errors", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 4L, seed = 1L)
  expect_error(
    run_alpha_diversity_analysis(sim_r),
    regexp = "group_levels.*group_colors"
  )
})

test_that("run_alpha_diversity_analysis(): unnamed colors error", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 4L, seed = 1L)
  expect_error(
    run_alpha_diversity_analysis(
      sim_r,
      group_levels = c("ctrl", "case"),
      group_colors = c("#aaaaaa", "#bbbbbb")
    ),
    regexp = "must be named"
  )
})

test_that("run_alpha_diversity_analysis(): colors missing a level name -> error", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 4L, seed = 1L)
  expect_error(
    run_alpha_diversity_analysis(
      sim_r,
      group_levels = c("ctrl", "case"),
      group_colors = c(ctrl = "#aaaaaa")
    ),
    regexp = "missing names for"
  )
})

test_that("run_alpha_diversity_analysis(): missing count_long columns -> error", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 4L, seed = 1L)
  sim_r$count_long$count <- NULL
  expect_error(
    run_alpha_diversity_analysis(
      sim_r,
      group_levels = c("ctrl", "case"),
      group_colors = c(ctrl = "#a", case = "#b")
    ),
    regexp = "is missing columns"
  )
})

# ---------------------------------------------------------------------------
# run_alpha_diversity_analysis_normalised (rarefaction; needs mia)
# ---------------------------------------------------------------------------

test_that("run_alpha_diversity_analysis_normalised(): runs end-to-end", {
  sim_r <- .make_tiny_sim(n_taxa = 6L, n_samples_per_group = 8L, seed = 33L,
                          library_size_mean = 3000, library_size_sd = 200)
  res <- run_alpha_diversity_analysis_normalised(
    sim_result = sim_r,
    case_label = "Test",
    group_levels = c("ctrl", "case"),
    group_colors = c(ctrl = "#377eb8", case = "#e41a1c"),
    rarefy_niter = 5L
  )
  expect_named(res, c("data", "summary", "tests", "plot"))
  expect_equal(nrow(res$data), 16L)  # 8 per group * 2
  expect_s3_class(res$plot, "ggplot")
})

test_that("run_alpha_diversity_analysis_normalised(): input validation errors", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 4L, seed = 1L)
  expect_error(
    run_alpha_diversity_analysis_normalised(sim_r),
    regexp = "group_levels.*group_colors"
  )
  expect_error(
    run_alpha_diversity_analysis_normalised(
      sim_r,
      group_levels = c("ctrl", "case"),
      group_colors = c(ctrl = "#a", case = "#b"),
      rarefy_niter = 0L
    ),
    regexp = "rarefy_niter.*positive"
  )
})
