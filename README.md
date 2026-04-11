# Beta dispersion interpretability

This repository contains the code and analyses for the paper *"Beta dispersion lacks interpretability for differential stochastic dispersion analyses"*. Using a generative compositional Bayesian model, we benchmark the sensitivity and specificity of PERMDISP (`vegan::betadisper`) across common distance constructions (Bray–Curtis, Aitchison, robust Aitchison, Jaccard, Hellinger) under controlled differential abundance and differential dispersion scenarios.

## License

This project is released under the **GNU General Public License v2.0 (GPL-2.0)**, an [OSI-approved](https://opensource.org/licenses/) open-source license. The GPL is a strong copyleft license: you may use, modify, and distribute the code, provided derivatives are released under the same license. See [LICENSE](LICENSE) for the full text.

## Requirements

- **R** (≥ 4.0)
- **R packages**: `targets`, `dplyr`, `tidyr`, `purrr`, `tibble`, `stringr`, `ggplot2`, `vegan`, `parallelly`, `sccomp`, `patchwork`, `VGAM`, `Quarto` (for reports)
- **Input data**: The HMP2 sccomp fit (`HMP2_sccomp_results.rds`) is required for the main pipeline. This file is available from [Zenodo (10.5281/zenodo.19144244)](https://zenodo.org/records/19144244) and is downloaded automatically when rendering the reports.

---

## File descriptions

Below is a brief summary of the files present in this repository.

| File | Description |
|------|-------------|
| `benchmark_paper_analyses.qmd` | Main Quarto report producing paper figures: compares PERMDISP distance inputs (Bray–Curtis, Aitchison, robust Aitchison, Jaccard, Hellinger), plots sensitivity under pure differential overdispersion, specificity under pure differential abundance, and sensitivity–specificity trade-offs |
| `Supplementary_figure_report.qmd` | Supplementary figures and diagnostics from targets outputs |
| `figure_1.png` | Main paper Figure 1 |
| `figure_2.png` | Main paper Figure 2 |
| `LICENSE` | GNU General Public License v2.0 |

---

## How to reproduce the results

### 1. Obtain the sccomp fit

The pipeline requires `HMP2_sccomp_results.rds` in the project root, an sccomp fit to the HMP2 gut microbiome dataset (Healthy vs IBD).

**Recommended (automatic)**: When you render `benchmark_paper_analyses.qmd`, the report automatically downloads this file from [Zenodo (10.5281/zenodo.19144244)](https://zenodo.org/records/19144244) if it is not already present.

**Manual download**: You can also download it directly to the project root:

```bash
curl -L -o HMP2_sccomp_results.rds "https://zenodo.org/records/19144244/files/HMP2_sccomp_results.rds?download=1"
```

Or in R:

```r
download.file("https://zenodo.org/records/19144244/files/HMP2_sccomp_results.rds?download=1",
              "HMP2_sccomp_results.rds", mode = "wb")
```

### 2. Install R dependencies

```r
install.packages(c("targets", "dplyr", "tidyr", "purrr", "tibble", "stringr", 
                   "ggplot2", "vegan", "parallelly", "sccomp", "patchwork", "VGAM", 
                   "quarto"))
# Optional, for SLURM/crew parallelization:
# install.packages("crew")
# install.packages("crew.cluster")
```

### 3. Run the targets pipeline

The full reproduction pipeline uses `_targets.R` (not in this repository). If you have access to it, from the project root:

```r
targets::tar_make()
```

This will:
- Load `HMP2_sccomp_results.rds`
- Extract simulation parameters
- Run DA and DO sweeps across distance methods
- Compute false positive and false negative rates
- Populate `_targets/` with cached objects

On macOS, the pipeline uses local parallelism (8 workers). On Linux with SLURM, it uses `crew.cluster` controllers (when `_targets.R` is available).

### 4. Render the main and supplementary reports

```bash
quarto render benchmark_paper_analyses.qmd
quarto render Supplementary_figure_report.qmd
```

These reports read from the targets store and produce the paper figures and supplementary diagnostics.

---

## Citation

If you use this code, please cite the associated paper when published.

---

## Contact

Mangiola Laboratory — [https://github.com/MangiolaLaboratory/beta_dispersion_interpretability](https://github.com/MangiolaLaboratory/beta_dispersion_interpretability)
