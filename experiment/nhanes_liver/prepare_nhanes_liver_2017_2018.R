#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(readr)
  library(tibble)
})

base_dir <- Sys.getenv("DCMA_DATA_ROOT", ".")
out_dir <- file.path(base_dir, "data", "nhanes_2017_2018_liver")
raw_dir <- file.path(out_dir, "raw_xpt")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

base_url <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles"
file_map <- c(
  DEMO = "DEMO_J.XPT",
  BMX = "BMX_J.XPT",
  LUX = "LUX_J.XPT",
  TRIGLY = "TRIGLY_J.XPT",
  HDL = "HDL_J.XPT",
  GLU = "GLU_J.XPT",
  INS = "INS_J.XPT",
  SMQ = "SMQ_J.XPT",
  ALQ = "ALQ_J.XPT",
  PAQ = "PAQ_J.XPT"
)

download_one <- function(fname) {
  dest <- file.path(raw_dir, fname)
  if (!file.exists(dest)) {
    url <- sprintf("%s/%s", base_url, fname)
    download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
  }
  dest
}

read_one <- function(tag, cols) {
  path <- download_one(file_map[[tag]])
  read_xpt(path) |>
    as_tibble() |>
    select(any_of(cols))
}

options(timeout = 180)

demo <- read_one(
  "DEMO",
  c(
    "SEQN", "SDMVSTRA", "SDMVPSU", "WTMEC2YR", "RIAGENDR", "RIDAGEYR", "RIDRETH3",
    "DMDEDUC2", "INDFMPIR"
  )
)
bmx <- read_one("BMX", c("SEQN", "BMXBMI", "BMXWAIST", "BMXHT"))
lux <- read_one("LUX", c("SEQN", "LUAXSTAT", "LUXSMED", "LUXSIQR", "LUXSIQRM", "LUXCAPM", "LUXCPIQR"))
trigly <- read_one("TRIGLY", c("SEQN", "WTSAF2YR", "LBXTR"))
hdl <- read_one("HDL", c("SEQN", "LBDHDD"))
glu <- read_one("GLU", c("SEQN", "LBXGLU"))
ins <- read_one("INS", c("SEQN", "LBXIN"))
smq <- read_one("SMQ", c("SEQN", "SMQ020", "SMQ040"))
alq <- read_one("ALQ", c("SEQN", "ALQ111", "ALQ121", "ALQ151"))
paq <- read_one("PAQ", c("SEQN", "PAQ605", "PAQ620", "PAQ650", "PAQ665"))

valid_code <- function(x) {
  if_else(x %in% c(7, 9, 77, 99, 777, 999, 7777, 9999), NA_real_, as.numeric(x))
}

merged <- demo |>
  left_join(bmx, by = "SEQN") |>
  left_join(lux, by = "SEQN") |>
  left_join(trigly, by = "SEQN") |>
  left_join(hdl, by = "SEQN") |>
  left_join(glu, by = "SEQN") |>
  left_join(ins, by = "SEQN") |>
  left_join(smq, by = "SEQN") |>
  left_join(alq, by = "SEQN") |>
  left_join(paq, by = "SEQN") |>
  mutate(
    sex = case_when(
      RIAGENDR == 1 ~ "Male",
      RIAGENDR == 2 ~ "Female",
      TRUE ~ NA_character_
    ),
    race_eth = case_when(
      RIDRETH3 == 1 ~ "Mexican American",
      RIDRETH3 == 2 ~ "Other Hispanic",
      RIDRETH3 == 3 ~ "Non-Hispanic White",
      RIDRETH3 == 4 ~ "Non-Hispanic Black",
      RIDRETH3 == 6 ~ "Non-Hispanic Asian",
      RIDRETH3 == 7 ~ "Other/Multiracial",
      TRUE ~ NA_character_
    ),
    education = case_when(
      valid_code(DMDEDUC2) == 1 ~ "Less than 9th grade",
      valid_code(DMDEDUC2) == 2 ~ "9-11th grade",
      valid_code(DMDEDUC2) == 3 ~ "High school/GED",
      valid_code(DMDEDUC2) == 4 ~ "Some college/AA",
      valid_code(DMDEDUC2) == 5 ~ "College graduate+",
      TRUE ~ NA_character_
    ),
    smoking_status = case_when(
      valid_code(SMQ020) == 2 ~ "Never",
      valid_code(SMQ020) == 1 & valid_code(SMQ040) %in% c(1, 2) ~ "Current",
      valid_code(SMQ020) == 1 & valid_code(SMQ040) == 3 ~ "Former",
      TRUE ~ NA_character_
    ),
    alcohol_status = case_when(
      valid_code(ALQ111) == 2 ~ "No lifetime alcohol",
      valid_code(ALQ111) == 1 & valid_code(ALQ121) == 0 ~ "No past-year alcohol",
      valid_code(ALQ111) == 1 & !is.na(valid_code(ALQ121)) & valid_code(ALQ121) > 0 ~ "Past-year alcohol",
      TRUE ~ NA_character_
    ),
    binge_drinking = case_when(
      alcohol_status %in% c("No lifetime alcohol", "No past-year alcohol") ~ "No",
      valid_code(ALQ151) == 1 ~ "Yes",
      valid_code(ALQ151) == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    any_physical_activity = case_when(
      valid_code(PAQ605) == 1 | valid_code(PAQ620) == 1 |
        valid_code(PAQ650) == 1 | valid_code(PAQ665) == 1 ~ "Yes",
      valid_code(PAQ605) == 2 & valid_code(PAQ620) == 2 &
        valid_code(PAQ650) == 2 & valid_code(PAQ665) == 2 ~ "No",
      TRUE ~ NA_character_
    ),
    obesity_bmi30 = if_else(!is.na(BMXBMI), as.integer(BMXBMI >= 30), NA_integer_),
    whtr = if_else(!is.na(BMXWAIST) & !is.na(BMXHT) & BMXHT > 0, BMXWAIST / BMXHT, NA_real_),
    central_obesity = case_when(
      RIAGENDR == 1 & !is.na(BMXWAIST) ~ as.integer(BMXWAIST >= 102),
      RIAGENDR == 2 & !is.na(BMXWAIST) ~ as.integer(BMXWAIST >= 88),
      TRUE ~ NA_integer_
    ),
    homa_ir = if_else(!is.na(LBXGLU) & !is.na(LBXIN), LBXGLU * LBXIN / 405, NA_real_),
    log_homa_ir = if_else(!is.na(homa_ir) & homa_ir > 0, log(homa_ir), NA_real_),
    tg_hdl_ratio = if_else(!is.na(LBXTR) & !is.na(LBDHDD) & LBDHDD > 0, LBXTR / LBDHDD, NA_real_),
    log_tg_hdl_ratio = if_else(!is.na(tg_hdl_ratio) & tg_hdl_ratio > 0, log(tg_hdl_ratio), NA_real_),
    te_complete = as.integer(LUAXSTAT == 1),
    te_valid_lsm = as.integer(!is.na(LUXSMED) & !is.na(LUXSIQRM) & LUXSIQRM <= 30),
    te_valid_cap = as.integer(!is.na(LUXCAPM)),
    te_valid_both = as.integer(te_complete == 1 & te_valid_lsm == 1 & te_valid_cap == 1)
  )

analysis_fasting <- merged |>
  filter(te_valid_both == 1) |>
  filter(!is.na(obesity_bmi30), !is.na(RIDAGEYR), !is.na(RIAGENDR), !is.na(RIDRETH3)) |>
  filter(!is.na(LUXCAPM), !is.na(LUXSMED), !is.na(LBXTR), !is.na(LBDHDD), !is.na(LBXGLU), !is.na(LBXIN), !is.na(homa_ir)) |>
  filter(
    !is.na(education), !is.na(INDFMPIR), !is.na(smoking_status),
    !is.na(alcohol_status), !is.na(binge_drinking), !is.na(any_physical_activity)
  ) |>
  mutate(
    cov_age = as.numeric(RIDAGEYR),
    cov_sex = sex,
    cov_race = race_eth,
    cov_education = education,
    cov_pir = as.numeric(INDFMPIR),
    cov_smoking = smoking_status,
    cov_alcohol = alcohol_status,
    cov_binge = binge_drinking,
    cov_pa_any = any_physical_activity
  ) |>
  select(
    SEQN, WTMEC2YR, WTSAF2YR, SDMVSTRA, SDMVPSU,
    RIDAGEYR, RIAGENDR, RIDRETH3, sex, race_eth,
    DMDEDUC2, INDFMPIR, education, smoking_status, alcohol_status, binge_drinking, any_physical_activity,
    cov_age, cov_sex, cov_race, cov_education, cov_pir,
    cov_smoking, cov_alcohol, cov_binge, cov_pa_any,
    BMXBMI, BMXWAIST, BMXHT, obesity_bmi30, central_obesity, whtr,
    LBXTR, LBDHDD, LBXGLU, LBXIN, homa_ir, log_homa_ir, tg_hdl_ratio, log_tg_hdl_ratio,
    LUXCAPM, LUXSMED, LUXSIQRM, LUXCPIQR
  )

summary_df <- tibble(
  n_demo = nrow(demo),
  n_with_te = sum(!is.na(merged$LUXCAPM) | !is.na(merged$LUXSMED)),
  n_te_valid_both = sum(merged$te_valid_both == 1, na.rm = TRUE),
  n_analysis_fasting = nrow(analysis_fasting),
  obesity_rate = mean(analysis_fasting$obesity_bmi30, na.rm = TRUE),
  cap_mean = mean(analysis_fasting$LUXCAPM, na.rm = TRUE),
  cap_sd = sd(analysis_fasting$LUXCAPM, na.rm = TRUE),
  lsm_mean = mean(analysis_fasting$LUXSMED, na.rm = TRUE),
  lsm_sd = sd(analysis_fasting$LUXSMED, na.rm = TRUE),
  homa_ir_median = median(analysis_fasting$homa_ir, na.rm = TRUE),
  tg_hdl_ratio_median = median(analysis_fasting$tg_hdl_ratio, na.rm = TRUE)
)

cov_summary_df <- tibble(
  variable = c(
    "cov_education", "cov_pir_raw", "cov_smoking", "cov_alcohol", "cov_binge", "cov_pa_any"
  ),
  n_missing_after_complete_case_filter = c(
    sum(is.na(analysis_fasting$cov_education)),
    sum(is.na(analysis_fasting$cov_pir)),
    sum(is.na(analysis_fasting$cov_smoking)),
    sum(is.na(analysis_fasting$cov_alcohol)),
    sum(is.na(analysis_fasting$cov_binge)),
    sum(is.na(analysis_fasting$cov_pa_any))
  )
)

dict_df <- tribble(
  ~variable, ~role, ~description,
  "obesity_bmi30", "A (exposure)", "Binary obesity indicator (BMI >= 30)",
  "homa_ir", "M (mediator)", "Insulin resistance proxy, glucose*insulin/405",
  "LBXTR", "M (mediator)", "Triglycerides (mg/dL)",
  "LBDHDD", "M (mediator)", "HDL cholesterol (mg/dL)",
  "LUXCAPM", "Y (outcome)", "Median CAP from FibroScan (dB/m)",
  "LUXSMED", "Y (outcome)", "Median liver stiffness from FibroScan (kPa)",
  "RIDAGEYR", "C (covariate)", "Age in years",
  "RIAGENDR", "C (covariate)", "Sex",
  "RIDRETH3", "C (covariate)", "Race/ethnicity category",
  "cov_education", "C (covariate)", "Educational attainment",
  "cov_pir", "C (covariate)", "Family income-to-poverty ratio",
  "cov_smoking", "C (covariate)", "Smoking status: never, former, or current",
  "cov_alcohol", "C (covariate)", "Alcohol-use status: no lifetime, no past-year, or past-year",
  "cov_binge", "C (covariate)", "Binge drinking indicator",
  "cov_pa_any", "C (covariate)", "Any moderate/vigorous work or recreational physical activity",
  "WTMEC2YR", "weight", "MEC exam weight (2-year)",
  "WTSAF2YR", "weight", "Fasting subsample weight (2-year)"
)

out_merged <- file.path(out_dir, "nhanes_2017_2018_liver_merged.csv")
out_analysis <- file.path(out_dir, "nhanes_2017_2018_liver_analysis_fasting.csv")
out_summary <- file.path(out_dir, "nhanes_2017_2018_liver_summary.csv")
out_cov_summary <- file.path(out_dir, "nhanes_2017_2018_liver_covariate_missingness.csv")
out_dict <- file.path(out_dir, "nhanes_2017_2018_liver_dictionary.csv")

write_csv(merged, out_merged, na = "")
write_csv(analysis_fasting, out_analysis, na = "")
write_csv(summary_df, out_summary, na = "")
write_csv(cov_summary_df, out_cov_summary, na = "")
write_csv(dict_df, out_dict, na = "")

cat("Wrote files:\n")
cat(" - ", out_merged, "\n", sep = "")
cat(" - ", out_analysis, "\n", sep = "")
cat(" - ", out_summary, "\n", sep = "")
cat(" - ", out_cov_summary, "\n", sep = "")
cat(" - ", out_dict, "\n", sep = "")
print(summary_df)
print(cov_summary_df)
