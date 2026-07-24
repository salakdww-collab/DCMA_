# DCMA code

This repository contains the R implementation used for DCMA and the analysis scripts for the reported synthetic, semi-synthetic, and NHANES experiments.

## Structure

- `model/`: DCMA fitting and reconstruction code.
- `utils/`: neural network and effect-summary utilities.
- `experiment/`: experiment-specific code.

Run scripts from the repository root. Data files are not included.

## Run examples

Synthetic benchmark:

```sh
BENCH_SCENARIOS=S1,S2 BENCH_R_REP=100 Rscript experiment/synthetic/export_shared_data.R
BENCH_SCENARIOS=S1,S2 Rscript experiment/synthetic/run_dcma_es.R
BENCH_SCENARIOS=S1,S2 Rscript experiment/synthetic/run_dcma_wgr.R
BENCH_SCENARIOS=S1,S2 Rscript experiment/synthetic/run_linear.R
Rscript experiment/synthetic/collect_results.R
```

S1 outcome-noise ablation:

```sh
S1_NOISE_R_REP=100 Rscript experiment/ablation/run_s1_noise_ablation.R
```

S2 joint-mediator ablation:

```sh
S2_ABLAT_R_REP=100 Rscript experiment/ablation/run_s2_joint_ablation.R
```

IHDP semi-synthetic experiment:

```sh
IHDP_BASE_CSV=/path/to/ihdp_npci_1.csv \
IHDP_SEMISYNTH_OUT_DIR=data/ihdp_covtail_skew_tnnls_v7d_n5000_rep100 \
IHDP_R_REP=100 \
IHDP_SCENARIO_LIST=IHDPMechanismThreeChannelSmoothTail \
Rscript experiment/ihdp_semisynthetic/prepare_ihdp_semisynth_3m.R

TNNLS_COVTAIL_REP100_MANIFEST=data/ihdp_covtail_skew_tnnls_v7d_n5000_rep100/manifest.csv \
Rscript experiment/ihdp_semisynthetic/run_ihdp_dcma.R
```

NHANES liver experiment:

```sh
Rscript experiment/nhanes_liver/prepare_nhanes_liver_2017_2018.R

NHANES_LIVER_IN_CSV=data/nhanes_2017_2018_liver/nhanes_2017_2018_liver_analysis_fasting.csv \
NHANES_LIVER_FULLBOOT_N=100 \
Rscript experiment/nhanes_liver/run_nhanes_liver.R
```
