# =============================================================================
# tests/testthat.R - Top-level test runner
# =============================================================================
#
# Mirrors the layout used by R packages even though this project ships as a
# loose script + R/ module collection rather than as an installed package:
#
#   tests/
#     testthat.R              <- this file (entry point)
#     testthat/
#       helper-setup.R        <- shared fixtures, sourced before any test runs
#       test-<module>.R       <- one file per R/<module>.R
#
# Run from the project root with either:
#
#   Rscript tests/testthat.R
#
# or, interactively:
#
#   testthat::test_dir(here::here("tests", "testthat"))
#
# =============================================================================

library(testthat)

testthat::test_dir(
  here::here("tests", "testthat"),
  reporter = testthat::default_reporter(),
  stop_on_failure = TRUE
)
