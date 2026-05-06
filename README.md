# Beta dispersion interpretability

This repository contains the code and analyses for the paper *"Beta dispersion lacks interpretability for differential stochastic dispersion analyses"*. Using a generative compositional Bayesian model, we benchmark the sensitivity and specificity of PERMDISP (`vegan::betadisper`) across common distance constructions (Bray–Curtis, Aitchison, robust Aitchison, Jaccard, Hellinger) under controlled differential abundance and differential dispersion scenarios.

## License

This project is released under the **GNU General Public License v2.0 (GPL-2.0)**, an [OSI-approved](https://opensource.org/licenses/) open-source license. The GPL is a strong copyleft license: you may use, modify, and distribute the code, provided derivatives are released under the same license. See [LICENSE](LICENSE) for the full text.

---

## Repository layout

| Location | Contents |
|---------|----------|
| **Repository root** | Shared R code [`functions.R`](functions.R); canonical **`{targets}`** pipeline script [`_targets.R`](_targets.R) for the main HMP2 benchmark; [`LICENSE`](LICENSE), this README |
| **`HMP_2019/`** | HMP longitudinal MIA-style benchmark (distinct filtering and cohort relative to Zenodo HMP2): [`benchmark_paper_analyses_HMP_2019.qmd`](HMP_2019/benchmark_paper_analyses_HMP_2019.qmd); supplementary [**alternative prevalence/detection**](HMP_2019/benchmark_paper_analyses_HMP_2019_SUPPLEMENTARY_ALTERNATIVE_FILTERING.qmd). Each report writes **`_targets_HMP_2019.R`** via `tar_script()` and builds store **`_targets_HMP_2019/`** |
| **`BritoIL_2016/`** | Same benchmark scaffold fit to the Brito *et al.* (IL 2016) cohort: [`benchmark_paper_analyses_BritoIL_2016.qmd`](BritoIL_2016/benchmark_paper_analyses_BritoIL_2016.qmd). Pipeline script **`_targets_BritoIL_2016.R`**, store **`_targets_BritoIL_2016/`** |
| **`ShaoY_2019/`** | Same benchmark scaffold for Shao *et al.* (2019): [`benchmark_paper_analyses_ShaoY_2019.qmd`](ShaoY_2019/benchmark_paper_analyses_ShaoY_2019.qmd). Pipeline script **`_targets_ShaoY_2019.R`**, store **`_targets_ShaoY_2019/`** (default location when running from this folder) |

Rendered HTML and figure chunks are written next to each `.qmd` (e.g. `*_paper_figures/` as set in the document YAML). There are no committed “final” `figure_*.png` exports at the repository root; open the rendered reports to view figures.

---

## Requirements

- **R** (≥ 4.0)
- **Quarto** (for `.qmd` reports)
- **Core R packages** (used across reports): `targets`, `dplyr`, `tidyr`, `purrr`, `tibble`, `stringr`, `ggplot2`, `vegan`, `parallelly`, `sccomp`, `patchwork`, `VGAM`, `knitr`
- **Additional packages** pulled in by specific reports: e.g. `phyloseq`, `readr` (HMP2 BIOM rebuild in `base_analyses/benchmark_paper_analyses.qmd`); `kableExtra`, `forcats`, `RColorBrewer`, `tidybulk`, `ggside`, `tidyHeatmap`, `ggrepel`, `lme4` (cohort-specific Quarto documents). Install missing packages when the render errors on `library(...)`.
- **Compute / memory**: the HMP 2019 MIA-style report states that a full render may need on the order of **~40 GB RAM**; sweeps are long-running. Only one `tar_make` should use a given store at a time (avoid concurrent renders that lock the same `_targets*` directory).

---

## Input data

### Main HMP2 benchmark (`base_analyses/`)

The pipeline expects **`HMP2_sccomp_results.rds`** in **`base_analyses/`** (same working directory as the Quarto execution for that report).

- **Rebuild from local inputs (default path in the report):** if `HMP2_sccomp_results.rds` is absent, `benchmark_paper_analyses.qmd` can rebuild a genus-level table and refit sccomp from **`taxonomic_profiles.biom`** and **`hmp2_metadata.csv`** placed in `base_analyses/`. That path is slower than reusing a precomputed fit.
- **Optional precomputed fit:** you can place `HMP2_sccomp_results.rds` in `base_analyses/` yourself. A copy is also available from [Zenodo (10.5281/zenodo.19144244)](https://zenodo.org/records/19144244) (large download). The automatic Zenodo `download.file()` fallback in the Quarto source is **commented out**; use `curl` or R if you want to fetch it manually.

### Cohort-specific benchmarks

Each of `HMP_2019/`, `BritoIL_2016/`, and `ShaoY_2019/` expects the corresponding **sccomp results `.rds`** and metadata paths referenced inside that `.qmd` (file names are defined in the setup chunks). Obtain or regenerate those inputs according to the methods section of the paper; they are not all shipped in this repository.

---

## How to reproduce the analyses

### 1. Install R and Quarto dependencies

```r
install.packages(c(
  "targets", "dplyr", "tidyr", "purrr", "tibble", "stringr",
  "ggplot2", "vegan", "parallelly", "sccomp", "patchwork", "VGAM",
  "knitr", "rmarkdown"
))
# As needed when rendering:
# install.packages(c("phyloseq", "readr", "quarto", "kableExtra", "forcats",
#                    "RColorBrewer", "tidybulk", "ggside", "tidyHeatmap",
#                    "ggrepel", "lme4"))
```

### 2. Link shared R sources from each analysis directory

Quarto sets the **code working directory to the folder that contains the `.qmd`**. Shared [`functions.R`](functions.R) and the main [`_targets.R`](_targets.R) live at the **repository root**, so from **`base_analyses/`** create symlinks (or copies) before the first render:

```bash
cd base_analyses
ln -sf ../functions.R functions.R
ln -sf ../_targets.R _targets.R
```

For **`HMP_2019/`**, **`BritoIL_2016/`**, and **`ShaoY_2019/`**, link only **`functions.R`** (each report generates its own `_targets_*.R` via `targets::tar_script()`):

```bash
(cd HMP_2019 && ln -sf ../functions.R functions.R)
(cd BritoIL_2016 && ln -sf ../functions.R functions.R)
(cd ShaoY_2019 && ln -sf ../functions.R functions.R)
```

The file [`HMP_2019/benchmark_paper_analyses_HMP_2019_SUPPLEMENTARY_ALTERNATIVE_FILTERING.qmd`](HMP_2019/benchmark_paper_analyses_HMP_2019_SUPPLEMENTARY_ALTERNATIVE_FILTERING.qmd) does **not** run `tar_make`; it reads from **`_targets_HMP_2019/`**. Build or render the main [`HMP_2019/benchmark_paper_analyses_HMP_2019.qmd`](HMP_2019/benchmark_paper_analyses_HMP_2019.qmd) first so that store exists.

### 3. Render Quarto reports

Run renders **from the repository root** with an explicit execution directory so chunk paths resolve correctly:

```bash
cd /path/to/beta_dispersion_interpretability

# Main HMP2 benchmark + supplementary figures (after symlinks in base_analyses/)
quarto render base_analyses/benchmark_paper_analyses.qmd --execute-dir=base_analyses
# Supplementary report only tar_read()s the default store; run after the main benchmark has built `_targets/`.
quarto render base_analyses/Supplementary_figure_report.qmd --execute-dir=base_analyses

# Cohort-specific benchmarks (after functions.R symlink in each folder)
quarto render HMP_2019/benchmark_paper_analyses_HMP_2019.qmd --execute-dir=HMP_2019
quarto render HMP_2019/benchmark_paper_analyses_HMP_2019_SUPPLEMENTARY_ALTERNATIVE_FILTERING.qmd --execute-dir=HMP_2019

quarto render BritoIL_2016/benchmark_paper_analyses_BritoIL_2016.qmd --execute-dir=BritoIL_2016
quarto render ShaoY_2019/benchmark_paper_analyses_ShaoY_2019.qmd --execute-dir=ShaoY_2019
```

The benchmark `.qmd` files call **`targets::tar_make()`** (or equivalent with a named `store =` for cohort pipelines). The first full build can take a long time and will populate the corresponding **`_targets/`** or **`_targets_*`** directory next to the executed report.

### 4. Optional: run `{targets}` only

If you prefer to build caches without Quarto:

```r
# From an R session whose working directory is base_analyses/ (with _targets.R present)
setwd("base_analyses")
targets::tar_make(script = "_targets.R")
```

Cohort scripts: use the path and `store` arguments matching the relevant `.qmd` (e.g. `store = "_targets_HMP_2019"` for the HMP 2019 report).

On Linux clusters, the generated pipeline scripts may be configured for heavier parallel backends (e.g. `crew` / `crew.cluster`); adjust `_targets.R` or the cohort-specific script if your environment differs.

---

## Citation

If you use this code, please cite the associated paper when published.

---

## Contact

Mangiola Laboratory — [https://github.com/MangiolaLaboratory/beta_dispersion_interpretability](https://github.com/MangiolaLaboratory/beta_dispersion_interpretability)
