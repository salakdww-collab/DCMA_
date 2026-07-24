# DCMA code

This repository contains the R implementation used for DCMA and the analysis scripts for the reported synthetic, semi-synthetic, and NHANES experiments.

## Structure

- `model/`: DCMA fitting and reconstruction code.
- `utils/`: neural network and effect-summary utilities.
- `experiment/`: experiment-specific code.

Data files are not included. Run scripts from the repository root. Scripts source the required files from `model/`, `utils/`, and their experiment folders.

## Main scripts

- Synthetic data generation: `experiment/synthetic/export_shared_data.R`
- Synthetic DCMA-ES: `experiment/synthetic/run_dcma_es.R`
- Synthetic DCMA-WGR: `experiment/synthetic/run_dcma_wgr.R`
- Synthetic linear benchmark: `experiment/synthetic/run_linear.R`
- Synthetic result collection: `experiment/synthetic/collect_results.R`
- IHDP preparation: `experiment/ihdp_semisynthetic/prepare_ihdp_semisynth_3m.R`
- IHDP DCMA: `experiment/ihdp_semisynthetic/run_ihdp_dcma.R`
- NHANES preparation: `experiment/nhanes_liver/prepare_nhanes_liver_2017_2018.R`
- NHANES bootstrap: `experiment/nhanes_liver/run_nhanes_liver.R`

The required R packages are listed in the `library()` calls at the top of each script.
