# =============================================================================
# tests/testthat/test-design.R - Tests for R/design.R
# =============================================================================

test_that("design_matrix_from_groups(): balanced design has correct shape", {
  dm <- design_matrix_from_groups(5)
  expect_equal(dim(dm), c(10L, 2L))
  expect_equal(colnames(dm), c("Intercept", "Group"))
  expect_true(all(dm[, "Intercept"] == 1))
  # First half Group 0, second half Group 1.
  expect_equal(dm[, "Group"], c(rep(0L, 5), rep(1L, 5)))
})

test_that("design_matrix_from_groups(): unbalanced design (length-2 input)", {
  dm <- design_matrix_from_groups(c(3, 7))
  expect_equal(dim(dm), c(10L, 2L))
  expect_equal(dm[, "Group"], c(rep(0L, 3), rep(1L, 7)))
})

test_that("design_matrix_from_groups(): rejects bad length input", {
  expect_error(design_matrix_from_groups(c(1, 2, 3)), regexp = "length 1.*length 2")
  expect_error(design_matrix_from_groups(integer(0)), regexp = "length 1.*length 2")
})

test_that("synthetic_sim_add_two_group_factor(): adds factor with correct levels", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 6L, seed = 42L,
                          group_levels = c("ctrl", "case"))
  # Strip the group factor and re-add it to test the function in isolation.
  sim_r$sample_metadata$group <- NULL
  sim_r$count_long$group <- NULL

  out <- synthetic_sim_add_two_group_factor(sim_r, c("ctrl", "case"))

  expect_true("group" %in% names(out$sample_metadata))
  expect_true("group" %in% names(out$count_long))
  expect_s3_class(out$sample_metadata$group, "factor")
  expect_equal(levels(out$sample_metadata$group), c("ctrl", "case"))

  # Group == 1 -> "case", Group == 0 -> "ctrl".
  meta <- out$sample_metadata
  expect_equal(as.character(meta$group[meta$Group == 1L]),
               rep("case", sum(meta$Group == 1L)))
  expect_equal(as.character(meta$group[meta$Group == 0L]),
               rep("ctrl", sum(meta$Group == 0L)))
})

test_that("synthetic_sim_add_two_group_factor(): rejects non-length-2 labels", {
  sim_r <- .make_tiny_sim(n_taxa = 4L, n_samples_per_group = 4L, seed = 1L)
  expect_error(synthetic_sim_add_two_group_factor(sim_r, "only_one"),
               regexp = "length 2")
  expect_error(synthetic_sim_add_two_group_factor(sim_r, c("a", "b", "c")),
               regexp = "length 2")
})
