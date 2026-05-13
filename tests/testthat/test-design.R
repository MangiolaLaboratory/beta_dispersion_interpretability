# =============================================================================
# tests/testthat/test-design.R - Tests for R/design.R
# =============================================================================
#
# Strategy: snapshot full design matrices and group factor values, not just
# shape, so a refactor that silently changes the row order or factor coding
# would be caught.
# =============================================================================

# ---------------------------------------------------------------------------
# design_matrix_from_groups()
# ---------------------------------------------------------------------------

test_that("design_matrix_from_groups(): balanced design matches exact expected matrix", {
  dm <- design_matrix_from_groups(5)
  expected <- cbind(
    Intercept = rep(1, 10),
    Group     = c(rep(0L, 5), rep(1L, 5))
  )
  expect_equal(dm, expected)
  expect_equal(colnames(dm), c("Intercept", "Group"))
  expect_equal(dim(dm), c(10L, 2L))
})

test_that("design_matrix_from_groups(): unbalanced (length-2) matches exact expected matrix", {
  dm <- design_matrix_from_groups(c(3, 7))
  expected <- cbind(
    Intercept = rep(1, 10),
    Group     = c(rep(0L, 3), rep(1L, 7))
  )
  expect_equal(dm, expected)
})

test_that("design_matrix_from_groups(): single-sample-per-group degenerate case", {
  dm <- design_matrix_from_groups(1L)
  expect_equal(dm,
               cbind(Intercept = c(1, 1), Group = c(0L, 1L)))
})

test_that("design_matrix_from_groups(): rejects bad lengths", {
  expect_error(design_matrix_from_groups(c(1, 2, 3)),
               regexp = "length 1.*length 2")
  expect_error(design_matrix_from_groups(integer(0)),
               regexp = "length 1.*length 2")
})

# ---------------------------------------------------------------------------
# synthetic_sim_add_two_group_factor()
# ---------------------------------------------------------------------------

test_that("synthetic_sim_add_two_group_factor(): adds factor with correct levels and full mapping", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 6L, seed = 42L,
                          group_levels = c("ctrl", "case"))
  # Strip the factor and re-add it to test the function in isolation.
  sim_r$sample_metadata$group <- NULL
  sim_r$count_long$group <- NULL

  out <- synthetic_sim_add_two_group_factor(sim_r, c("ctrl", "case"))

  # sample_metadata column added with right class and levels.
  expect_true("group" %in% colnames(out$sample_metadata))
  expect_s3_class(out$sample_metadata$group, "factor")
  expect_equal(levels(out$sample_metadata$group), c("ctrl", "case"))
  # 12 samples, 6 ctrl + 6 case.
  expect_equal(as.integer(table(out$sample_metadata$group)),
               c(6L, 6L))
  # Mapping: Group == 0L -> "ctrl"; Group == 1L -> "case" for EVERY row.
  expect_equal(
    as.character(out$sample_metadata$group),
    ifelse(out$sample_metadata$Group == 1L, "case", "ctrl")
  )

  # count_long column added with same levels.
  expect_s3_class(out$count_long$group, "factor")
  expect_equal(levels(out$count_long$group), c("ctrl", "case"))
  expect_equal(
    as.character(out$count_long$group),
    ifelse(out$count_long$Group == 1L, "case", "ctrl")
  )
})

test_that("synthetic_sim_add_two_group_factor(): preserves all other list elements unchanged", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 6L, seed = 42L,
                          group_levels = c("ctrl", "case"))
  before <- sim_r
  out <- synthetic_sim_add_two_group_factor(sim_r, c("ctrl", "case"))
  for (key in setdiff(names(before), c("sample_metadata", "count_long"))) {
    expect_identical(out[[key]], before[[key]], info = key)
  }
})

test_that("synthetic_sim_add_two_group_factor(): swapped labels map correctly", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 5L, seed = 1L)
  out <- synthetic_sim_add_two_group_factor(sim_r, c("case", "ctrl"))
  # First label now applies to Group == 0L.
  expect_equal(levels(out$sample_metadata$group), c("case", "ctrl"))
  expect_equal(
    as.character(out$sample_metadata$group),
    ifelse(out$sample_metadata$Group == 1L, "ctrl", "case")
  )
})

test_that("synthetic_sim_add_two_group_factor(): rejects non-length-2 labels", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 4L, seed = 1L)
  expect_error(synthetic_sim_add_two_group_factor(sim_r, "only_one"),
               regexp = "length 2")
  expect_error(synthetic_sim_add_two_group_factor(sim_r, c("a", "b", "c")),
               regexp = "length 2")
})
