# =============================================================================
# tests/testthat/test-transforms.R - Tests for R/transforms.R
# =============================================================================

test_that("softmax(): each row sums to 1 and preserves shape", {
  m <- matrix(c(0, 0, 0,
                1, 2, 3,
                -1, 0, 1), ncol = 3, byrow = TRUE)
  out <- softmax(m)
  expect_equal(dim(out), dim(m))
  expect_equal(unname(rowSums(out)), rep(1, nrow(m)), tolerance = 1e-12)
  # All probabilities strictly positive.
  expect_true(all(out > 0))
})

test_that("softmax(): equal inputs map to equal probabilities", {
  m <- matrix(0, nrow = 2, ncol = 4)
  out <- softmax(m)
  expect_equal(out, matrix(0.25, nrow = 2, ncol = 4))
})

test_that(".clamp01(): eps = 0 keeps interior values and clamps to [0, 1]", {
  x <- c(-1, 0, 0.5, 1, 2)
  expect_equal(.clamp01(x, eps = 0), c(0, 0, 0.5, 1, 1))
})

test_that(".clamp01(): eps > 0 clamps into [eps, 1 - eps]", {
  x <- c(0, 0.05, 0.5, 0.95, 1)
  expect_equal(.clamp01(x, eps = 0.1), c(0.1, 0.1, 0.5, 0.9, 0.9))
})

test_that(".clamp01(): preserves matrix shape", {
  m <- matrix(c(-1, 0.5, 2, 0.3), nrow = 2)
  out <- .clamp01(m, eps = 0)
  expect_equal(dim(out), dim(m))
  expect_true(all(out >= 0 & out <= 1))
})

test_that("arcsin_sqrt(): known anchor values", {
  expect_equal(arcsin_sqrt(0),   0)
  expect_equal(arcsin_sqrt(1),   pi / 2)
  expect_equal(arcsin_sqrt(0.5), pi / 4)
})

test_that("arcsin_sqrt(): eps clamps boundary values before transform", {
  # Without eps, arcsin_sqrt(c(0, 1)) is well-defined but at boundaries;
  # with eps > 0 the inputs are pulled inside the open interval.
  out <- arcsin_sqrt(c(0, 1), eps = 0.01)
  expect_true(all(out > 0 & out < pi / 2))
})

test_that("inv_softmax(): vector input -> probabilities summing to 1", {
  v <- c(-1, 0, 1, 2)
  out <- inv_softmax(v)
  expect_length(out, length(v))
  expect_equal(sum(out), 1, tolerance = 1e-12)
  expect_true(all(out > 0))
})

test_that("inv_softmax(): matrix input -> each row sums to 1 and dimnames preserved", {
  m <- matrix(c(0, 1, 2,
                3, 0, -3),
              nrow = 2, byrow = TRUE,
              dimnames = list(c("a", "b"), c("x", "y", "z")))
  out <- inv_softmax(m)
  expect_equal(dim(out), dim(m))
  expect_equal(dimnames(out), dimnames(m))
  expect_equal(unname(rowSums(out)), c(1, 1), tolerance = 1e-12)
})

test_that("inv_softmax(): numerically stable for huge inputs", {
  # softmax(c(1000, 1001)) overflows; inv_softmax must remain finite.
  out <- inv_softmax(c(1000, 1001))
  expect_true(all(is.finite(out)))
  expect_equal(sum(out), 1, tolerance = 1e-12)
  # Second component must dominate.
  expect_gt(out[2], out[1])
})

test_that("inv_softmax(): NA input errors clearly", {
  expect_error(inv_softmax(c(1, NA)), regexp = "NA")
})

test_that("inv_softmax(): non-numeric input errors", {
  expect_error(inv_softmax("a"), regexp = "numeric")
})

test_that("inv_softmax(): eps clamps outputs to [eps, 1 - eps]", {
  # A very lopsided input would normally produce probabilities near 0 / 1.
  v <- c(-50, 50)
  out <- inv_softmax(v, eps = 0.1)
  expect_true(all(out >= 0.1 - 1e-12))
  expect_true(all(out <= 0.9 + 1e-12))
})
