# =============================================================================
# functions.R - Simulation and Beta-Dispersion Analysis for Compositional Data
# =============================================================================
#
# Top-level loader. Sources every module under R/ in dependency order so that
# downstream scripts only need to call:
#
#     source(here::here("functions.R"))
#
# Supports the paper:
# "Beta dispersion lacks interpretability for differential stochastic dispersion
# analyses" (beta-dispersion interpretability project).
#
# =============================================================================

library(ggplot2)
library(dplyr)
library(purrr)
library(tibble)

source(here::here("R", "transforms.R"))
source(here::here("R", "simulation.R"))
source(here::here("R", "design.R"))
source(here::here("R", "sccomp.R"))
source(here::here("R", "betadisper.R"))
source(here::here("R", "alpha_diversity.R"))
source(here::here("R", "permanova.R"))
source(here::here("R", "plot_labels.R"))
source(here::here("R", "pipeline.R"))

