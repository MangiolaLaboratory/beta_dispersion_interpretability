# =============================================================================
# tests/testthat/test-alpha_diversity.R - Tests for R/alpha_diversity.R
# =============================================================================
#
# Strategy: simulation is deterministic via seed, so the per-sample richness
# / shannon / pielou_evenness values can be snapshotted exactly. We also pin
# the structure (column names, order, row counts) of the returned bundle.
# =============================================================================

# ---------------------------------------------------------------------------
# run_alpha_diversity_analysis()
# ---------------------------------------------------------------------------

test_that("run_alpha_diversity_analysis(): return-list names and shape are pinned", {
  sim_r <- .make_tiny_sim(n_taxa = 6L, n_samples_per_group = 10L, seed = 12L)
  res <- run_alpha_diversity_analysis(
    sim_result = sim_r,
    case_label = "Test",
    group_levels = c("ctrl", "case"),
    group_colors = c(ctrl = "#377eb8", case = "#e41a1c")
  )
  expect_named(res, c("data", "summary", "tests", "plot"), ignore.order = FALSE)
  expect_s3_class(res$data, "data.frame")
  expect_s3_class(res$plot, "ggplot")

  # data: one row per sample (20 here), exact columns in exact order.
  expect_equal(nrow(res$data), 20L)
  expect_named(res$data,
               c("sample_id", "group", "richness", "shannon", "pielou_evenness"),
               ignore.order = FALSE)
  expect_type(res$data$sample_id, "character")
  expect_s3_class(res$data$group, "factor")
  expect_equal(levels(res$data$group), c("ctrl", "case"))
  expect_type(res$data$richness, "integer")
  expect_type(res$data$shannon, "double")
  expect_type(res$data$pielou_evenness, "double")

  # tests: three named t-tests in the documented order.
  expect_named(res$tests, c("richness", "shannon", "evenness"), ignore.order = FALSE)

  # summary: one row per group level, exact columns in exact order.
  expect_named(res$summary,
               c("group", "mean_richness", "sd_richness",
                 "mean_shannon", "sd_shannon",
                 "mean_evenness", "sd_evenness"),
               ignore.order = FALSE)
  expect_equal(as.character(res$summary$group), c("ctrl", "case"))
})

test_that("run_alpha_diversity_analysis(): exact richness/shannon/pielou snapshot per sample", {
  sim_r <- .make_tiny_sim(n_taxa = 6L, n_samples_per_group = 10L, seed = 12L)
  res <- run_alpha_diversity_analysis(
    sim_result = sim_r,
    case_label = "Test",
    group_levels = c("ctrl", "case"),
    group_colors = c(ctrl = "#377eb8", case = "#e41a1c")
  )
  # Pinned snapshot for the first three samples. Rows are ordered by
  # sample_id from `dplyr::group_by()`; lexical order puts "Sample_10" /
  # "Sample_11" right after "Sample_1".
  expect_equal(
    head(res$data, 3),
    tibble::tibble(
      sample_id = c("Sample_1", "Sample_10", "Sample_11"),
      group = factor(c("ctrl", "ctrl", "case"), levels = c("ctrl", "case")),
      richness = c(3L, 1L, 2L),
      shannon = c(0.5359610497, 0.0000000000, 0.0511783527),
      pielou_evenness = c(0.4878527714, NaN, 0.0738347557)
    ),
    tolerance = 1e-7
  )
})

test_that("run_alpha_diversity_analysis(): per-group mean richness/shannon snapshot", {
  sim_r <- .make_tiny_sim(n_taxa = 6L, n_samples_per_group = 10L, seed = 12L)
  res <- run_alpha_diversity_analysis(
    sim_result = sim_r,
    case_label = "Test",
    group_levels = c("ctrl", "case"),
    group_colors = c(ctrl = "#377eb8", case = "#e41a1c")
  )
  # mean values come straight from .alpha_diversity_table_to_plot's group_by.
  expect_equal(res$summary$mean_richness, c(2.9, 2.6), tolerance = 1e-12)
  expect_equal(res$summary$sd_richness,
               c(stats::sd(res$data$richness[res$data$group == "ctrl"]),
                 stats::sd(res$data$richness[res$data$group == "case"])),
               tolerance = 1e-12)
})

test_that("run_alpha_diversity_analysis(): missing required args error", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 4L, seed = 1L)
  expect_error(run_alpha_diversity_analysis(sim_r),
               regexp = "group_levels.*group_colors")
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

test_that("run_alpha_diversity_analysis(): colors missing a level name errors", {
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

test_that("run_alpha_diversity_analysis(): missing required count_long columns errors", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 4L, seed = 1L)
  sim_r$count_long$count <- NULL
  expect_error(
    run_alpha_diversity_analysis(
      sim_r,
      group_levels = c("ctrl", "case"),
      group_colors = c(ctrl = "#a", case = "#b")
    ),
    regexp = "missing columns"
  )
})

# ---------------------------------------------------------------------------
# run_alpha_diversity_analysis_normalised()
# ---------------------------------------------------------------------------

test_that("run_alpha_diversity_analysis_normalised(): same bundle structure as raw variant", {
  sim_r <- .make_tiny_sim(n_taxa = 6L, n_samples_per_group = 8L, seed = 33L,
                          library_size_mean = 3000, library_size_sd = 200)
  res <- run_alpha_diversity_analysis_normalised(
    sim_result = sim_r,
    case_label = "Test",
    group_levels = c("ctrl", "case"),
    group_colors = c(ctrl = "#377eb8", case = "#e41a1c"),
    rarefy_niter = 5L
  )
  expect_named(res, c("data", "summary", "tests", "plot"), ignore.order = FALSE)
  expect_named(res$data,
               c("sample_id", "group", "richness", "shannon", "pielou_evenness"),
               ignore.order = FALSE)
  expect_equal(nrow(res$data), 16L)  # 8 per group * 2
  expect_named(res$tests, c("richness", "shannon", "evenness"), ignore.order = FALSE)
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
      group_colors = c("#a", "#b")
    ),
    regexp = "must be named"
  )
})
