# =============================================================================
# tests/testthat/test-transforms.R - Tests for R/transforms.R
# =============================================================================
#
# Strategy: tests pin down (a) the mathematical contract (closed-form numerical
# equality, dimensions, dimnames) and (b) the exact behaviour on edge cases.
# Future refactors that preserve behaviour should leave all these passing.
# =============================================================================

# ---------------------------------------------------------------------------
# softmax()
# ---------------------------------------------------------------------------

test_that("softmax(): matches closed-form exp(x)/rowSums(exp(x)) exactly", {
  m <- matrix(c(0, 1, 2,
                3, 0, -3),
              nrow = 2, byrow = TRUE)
  expected <- exp(m) / rowSums(exp(m))
  out <- softmax(m)
  expect_equal(out, expected, tolerance = 1e-15)
  # Spot-check exact values (from the closed form).
  expect_equal(out[1, ], c(exp(0), exp(1), exp(2)) / (exp(0) + exp(1) + exp(2)),
               tolerance = 1e-15)
})

test_that("softmax(): output shape and rowSums == 1", {
  m <- matrix(c(0, 0, 0,
                1, 2, 3,
                -1, 0, 1), ncol = 3, byrow = TRUE)
  out <- softmax(m)
  expect_equal(dim(out), c(3L, 3L))
  expect_equal(unname(rowSums(out)), c(1, 1, 1), tolerance = 1e-15)
})

test_that("softmax(): equal inputs map to a uniform row", {
  out <- softmax(matrix(0, nrow = 2, ncol = 4))
  expect_equal(out, matrix(0.25, nrow = 2, ncol = 4), tolerance = 1e-15)
})

# ---------------------------------------------------------------------------
# .clamp01()
# ---------------------------------------------------------------------------

test_that(".clamp01(): eps = 0 clamps to [0, 1], preserves interior", {
  expect_equal(.clamp01(c(-1, 0, 0.5, 1, 2), eps = 0), c(0, 0, 0.5, 1, 1))
})

test_that(".clamp01(): eps > 0 clamps to [eps, 1 - eps]", {
  expect_equal(.clamp01(c(0, 0.05, 0.5, 0.95, 1), eps = 0.1),
               c(0.1, 0.1, 0.5, 0.9, 0.9))
})

test_that(".clamp01(): preserves matrix shape and applies element-wise", {
  m <- matrix(c(-1, 0.5, 2, 0.3), nrow = 2)
  out <- .clamp01(m, eps = 0)
  expect_equal(dim(out), dim(m))
  expect_equal(out, matrix(c(0, 0.5, 1, 0.3), nrow = 2))
})

# ---------------------------------------------------------------------------
# arcsin_sqrt()
# ---------------------------------------------------------------------------

test_that("arcsin_sqrt(): exact anchor values", {
  expect_equal(arcsin_sqrt(0),    0,     tolerance = 1e-15)
  expect_equal(arcsin_sqrt(1),    pi / 2, tolerance = 1e-15)
  expect_equal(arcsin_sqrt(0.5),  pi / 4, tolerance = 1e-15)
  expect_equal(arcsin_sqrt(0.25), pi / 6, tolerance = 1e-15)
  expect_equal(arcsin_sqrt(0.75), pi / 3, tolerance = 1e-15)
})

test_that("arcsin_sqrt(): vector input -> exact closed-form output", {
  p <- c(0, 0.25, 0.5, 0.75, 1)
  expect_equal(arcsin_sqrt(p), asin(sqrt(p)), tolerance = 1e-15)
})

test_that("arcsin_sqrt(): eps clamps boundary values before transform", {
  out <- arcsin_sqrt(c(0, 1), eps = 0.01)
  # After clamping: c(0.01, 0.99) -> asin(sqrt(...))
  expect_equal(out, asin(sqrt(c(0.01, 0.99))), tolerance = 1e-15)
  # And the boundary values are strictly inside the range.
  expect_true(all(out > 0 & out < pi / 2))
})

# ---------------------------------------------------------------------------
# inv_softmax()
# ---------------------------------------------------------------------------

test_that("inv_softmax(): exact snapshot for a vector input", {
  out <- inv_softmax(c(-1, 0, 1, 2))
  expect_equal(
    out,
    c(0.03205860, 0.08714432, 0.23688282, 0.64391426),
    tolerance = 1e-7
  )
  expect_equal(sum(out), 1, tolerance = 1e-15)
})

test_that("inv_softmax(): matches max-shifted closed form", {
  v <- c(-1, 0, 1, 2)
  shifted <- v - max(v)
  expected <- exp(shifted) / sum(exp(shifted))
  expect_equal(inv_softmax(v), expected, tolerance = 1e-15)
})

test_that("inv_softmax(): matrix input -> row softmax, dimnames preserved, exact values", {
  m <- matrix(c(0, 1, 2,
                3, 0, -3),
              nrow = 2, byrow = TRUE,
              dimnames = list(c("a", "b"), c("x", "y", "z")))
  out <- inv_softmax(m)
  expect_equal(dim(out), c(2L, 3L))
  expect_equal(dimnames(out), dimnames(m))
  # Row 1: numerically stable softmax of c(0, 1, 2).
  r1_shift <- c(0, 1, 2) - 2
  r1_expected <- exp(r1_shift) / sum(exp(r1_shift))
  expect_equal(unname(out["a", ]), r1_expected, tolerance = 1e-15)
  # Row 2: numerically stable softmax of c(3, 0, -3).
  r2_shift <- c(3, 0, -3) - 3
  r2_expected <- exp(r2_shift) / sum(exp(r2_shift))
  expect_equal(unname(out["b", ]), r2_expected, tolerance = 1e-15)
  expect_equal(unname(rowSums(out)), c(1, 1), tolerance = 1e-15)
})

test_that("inv_softmax(): numerically stable for huge inputs", {
  # softmax(c(1000, 1001)) overflows; the max-shifted form must stay finite.
  out <- inv_softmax(c(1000, 1001))
  expect_true(all(is.finite(out)))
  expect_equal(sum(out), 1, tolerance = 1e-15)
  # First entry uses shifted exp(-1), second exp(0) -> exact values.
  denom <- exp(-1) + 1
  expect_equal(out, c(exp(-1) / denom, 1 / denom), tolerance = 1e-15)
})

test_that("inv_softmax(): NA input errors with informative message", {
  expect_error(inv_softmax(c(1, NA)), regexp = "NA")
})

test_that("inv_softmax(): non-numeric input errors", {
  expect_error(inv_softmax("a"), regexp = "must be numeric")
  expect_error(inv_softmax(list(1, 2)), regexp = "must be numeric")
})

test_that("inv_softmax(): eps clamps outputs to [eps, 1 - eps]", {
  out <- inv_softmax(c(-50, 50), eps = 0.1)
  expect_equal(out, c(0.1, 0.9), tolerance = 1e-15)
})

test_that("inv_softmax(): scalar input -> probability 1", {
  expect_equal(inv_softmax(0), 1, tolerance = 1e-15)
  expect_equal(inv_softmax(42), 1, tolerance = 1e-15)
})
