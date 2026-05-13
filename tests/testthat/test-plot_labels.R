# =============================================================================
# tests/testthat/test-plot_labels.R - Tests for R/plot_labels.R
# =============================================================================

# The en-dash literal used by `pretty_method_label("bray")`.
.endash <- "\u2013"

test_that("pretty_method_label(): known keys map to display strings", {
  expect_equal(pretty_method_label("bray"), paste0("Bray", .endash, "Curtis"))
  expect_equal(pretty_method_label("permdisp_bray"), paste0("Bray", .endash, "Curtis"))
  expect_equal(pretty_method_label("jaccard"), "Jaccard")
  expect_equal(pretty_method_label("permdisp_jaccard"), "Jaccard")
  expect_equal(pretty_method_label("hellinger"), "Hellinger")
  expect_equal(pretty_method_label("permdisp_hellinger"), "Hellinger")
  expect_equal(pretty_method_label("aitchison"), "Aitchison")
  expect_equal(pretty_method_label("permdisp_aitchison"), "Aitchison")
  expect_equal(pretty_method_label("robust.aitchison"), "robust Aitchison")
  expect_equal(pretty_method_label("permdisp_robust.aitchison"), "robust Aitchison")
  expect_equal(pretty_method_label("permdisp_euclidean"), "Euclidean")
})

test_that("pretty_method_label(): unknown keys are returned unchanged", {
  expect_equal(pretty_method_label("anything_else"), "anything_else")
  expect_equal(pretty_method_label(""), "")
})

test_that("pretty_method_label_roc() and _box() delegate to pretty_method_label()", {
  for (key in c("bray", "permdisp_jaccard", "hellinger", "unknown")) {
    expect_identical(pretty_method_label_roc(key), pretty_method_label(key))
    expect_identical(pretty_method_label_box(key), pretty_method_label(key))
  }
})
