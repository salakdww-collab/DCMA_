nhanes_liver_covariate_vars <- function() {
  c(
    "C_age", "C_pir",
    "C_sex", "C_race", "C_education",
    "C_smoking", "C_alcohol", "C_binge", "C_pa_any"
  )
}

nhanes_liver_factor_covariate_vars <- function() {
  c("C_sex", "C_race", "C_education", "C_smoking", "C_alcohol", "C_binge", "C_pa_any")
}

nhanes_liver_make_analysis_dat <- function(dat0, y_var) {
  out <- dat0 |>
    dplyr::transmute(
      A = as.integer(obesity_bmi30),
      Y = as.numeric(.data[[y_var]]),
      M1 = as.numeric(log_homa_ir),
      M2 = as.numeric(log_tg_hdl_ratio),
      C_age = as.numeric(cov_age),
      C_pir = as.numeric(cov_pir),
      C_sex = factor(cov_sex, levels = c("Male", "Female")),
      C_race = factor(
        cov_race,
        levels = c(
          "Mexican American", "Other Hispanic", "Non-Hispanic White",
          "Non-Hispanic Black", "Non-Hispanic Asian", "Other/Multiracial"
        )
      ),
      C_education = factor(
        cov_education,
        levels = c(
          "Less than 9th grade", "9-11th grade", "High school/GED",
          "Some college/AA", "College graduate+"
        )
      ),
      C_smoking = factor(cov_smoking, levels = c("Never", "Former", "Current")),
      C_alcohol = factor(
        cov_alcohol,
        levels = c("No lifetime alcohol", "No past-year alcohol", "Past-year alcohol")
      ),
      C_binge = factor(cov_binge, levels = c("No", "Yes")),
      C_pa_any = factor(cov_pa_any, levels = c("No", "Yes"))
    ) |>
    dplyr::filter(dplyr::if_all(dplyr::everything(), ~ !is.na(.x)))

  out
}

nhanes_liver_make_dcma_dat <- function(dat_base) {
  fac_vars <- nhanes_liver_factor_covariate_vars()
  xcat <- stats::model.matrix(stats::reformulate(fac_vars), data = dat_base)[, -1, drop = FALSE]
  colnames(xcat) <- make.names(colnames(xcat))

  cbind(
    X = as.numeric(dat_base$A),
    Y = as.numeric(dat_base$Y),
    M1 = as.numeric(dat_base$M1),
    M2 = as.numeric(dat_base$M2),
    C_age = as.numeric(dat_base$C_age),
    C_pir = as.numeric(dat_base$C_pir),
    as.data.frame(xcat)
  ) |>
    as.data.frame()
}
