# =============================================================================
# tests/testthat/test-plot_labels.R - Tests for R/plot_labels.R
# =============================================================================

# The en-dash literal that the function uses.
.endash <- "\u2013"

test_that("pretty_method_label(): full lookup table pinned exactly", {
  # Pin every documented (key -> label) mapping. Any new mapping is intentional
  # behaviour change and must update this table.
  expected <- list(
    "bray"                       = paste0("Bray", .endash, "Curtis"),
    "permdisp_bray"              = paste0("Bray", .endash, "Curtis"),
    "jaccard"                    = "Jaccard",
    "permdisp_jaccard"           = "Jaccard",
    "hellinger"                  = "Hellinger",
    "permdisp_hellinger"         = "Hellinger",
    "aitchison"                  = "Aitchison",
    "permdisp_aitchison"         = "Aitchison",
    "robust.aitchison"           = "robust Aitchison",
    "permdisp_robust.aitchison"  = "robust Aitchison",
    "permdisp_euclidean"         = "Euclidean"
  )
  for (key in names(expected)) {
    out <- pretty_method_label(key)
    expect_type(out, "character")
    expect_length(out, 1L)
    expect_identical(out, expected[[key]], info = key)
  }
})

test_that("pretty_method_label(): unknown keys are returned verbatim", {
  for (key in c("anything_else", "", "permdisp_unknown", "PERMDISP_BRAY")) {
    out <- pretty_method_label(key)
    expect_identical(out, key, info = key)
  }
})

test_that("pretty_method_label_roc() / _box() delegate to pretty_method_label()", {
  keys <- c("bray", "permdisp_jaccard", "hellinger", "aitchison",
            "robust.aitchison", "permdisp_euclidean", "unknown_key", "")
  for (key in keys) {
    expect_identical(pretty_method_label_roc(key), pretty_method_label(key),
                     info = paste("roc:", key))
    expect_identical(pretty_method_label_box(key), pretty_method_label(key),
                     info = paste("box:", key))
  }
})
