# Standard testthat entry point for the betaDispersionInterpretability package.
#
# Run via:
#   R CMD check .                # full package check (installs, then runs)
#   Rscript -e 'devtools::test()' # in-place (no install, uses pkgload)
#   Rscript tests/testthat.R     # after installing the package

library(testthat)
library(betaDispersionInterpretability)

test_check("betaDispersionInterpretability")
