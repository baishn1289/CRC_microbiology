library(tidyverse)
library(MMUPHin)
library(metafor)

code_dir <- file.path(
  base_symbol,
  "downstream_analysis/src/github_code_availability/submission_version_261102/FigS2_sensitivity_analysis"
)
sensitivity_root <- file.path(
  base_symbol,
  "downstream_analysis/output/submission_version_261102/FigS2_sensitivity_analysis"
)
source(file.path(code_dir, "00.sensitivity_functions.R"))

results_dir <- file.path(sensitivity_root, "REML_sensitivity_analysis/results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Reuse cohort fits; only the random-effects estimator changes.
dir.create(
  file.path(sensitivity_root, "models/REML_EO"),
  recursive = TRUE,
  showWarnings = FALSE
)
dir.create(
  file.path(sensitivity_root, "models/REML_LO"),
  recursive = TRUE,
  showWarnings = FALSE
)
fit_meta.crc_EO_REML <- fit_meta.crc_EO
fit_meta.crc_LO_REML <- fit_meta.crc_LO
fit_meta.crc_EO_REML$control$rma_method <- "REML"
fit_meta.crc_LO_REML$control$rma_method <- "REML"
fit_meta.crc_EO_REML$meta_fits <- MMUPHin:::rma_wrapper(
  fit_meta.crc_EO$maaslin_fits,
  method = "REML",
  output = file.path(sensitivity_root, "models/REML_EO"),
  forest_plot = NULL
)
fit_meta.crc_LO_REML$meta_fits <- MMUPHin:::rma_wrapper(
  fit_meta.crc_LO$maaslin_fits,
  method = "REML",
  output = file.path(sensitivity_root, "models/REML_LO"),
  forest_plot = NULL
)

# Apply shared meta-QC and add estimator-specific prediction intervals.
reml_qc <- function(fit, estimator) {
  sensitivity_qc(fit) %>%
    mutate(
      estimator = estimator,
      Species = str_remove(str_extract(feature, "s__[A-Za-z0-9_]+"), "^s__"),
      pi_se = if_else(k >= 2 & !is.na(tau2), sqrt(stderr^2 + pmax(tau2, 0)), NA_real_),
      pi_lower_95 = coef - qnorm(0.975) * pi_se,
      pi_upper_95 = coef + qnorm(0.975) * pi_se,
      prediction_interval_excludes_zero = !is.na(pi_lower_95) & (pi_lower_95 > 0 | pi_upper_95 < 0)
    )
}

# Pair EO and LO estimates among taxa passing QC in both strata.
make_shared_qc_table <- function(results_eo, results_lo, estimator) {
  eo <- results_eo %>%
    filter(pass_qc) %>%
    select(
      feature,
      Species,
      coef_EO = coef,
      stderr_EO = stderr,
      pval_EO = pval,
      qval_EO = qval.fdr,
      tau2_EO = tau2,
      I2_EO = I2,
      max_weight_EO = max_weight,
      consistency_EO = consistency_score,
      significant_EO = significant,
      pi_lower_EO = pi_lower_95,
      pi_upper_EO = pi_upper_95,
      pi_excludes_zero_EO = prediction_interval_excludes_zero
    )
  lo <- results_lo %>%
    filter(pass_qc) %>%
    select(
      feature,
      coef_LO = coef,
      stderr_LO = stderr,
      pval_LO = pval,
      qval_LO = qval.fdr,
      tau2_LO = tau2,
      I2_LO = I2,
      max_weight_LO = max_weight,
      consistency_LO = consistency_score,
      significant_LO = significant,
      pi_lower_LO = pi_lower_95,
      pi_upper_LO = pi_upper_95,
      pi_excludes_zero_LO = prediction_interval_excludes_zero
    )
  inner_join(eo, lo, by = "feature") %>%
    mutate(
      estimator = estimator,
      same_direction = sign(coef_EO) == sign(coef_LO),
      both_significant = significant_EO & significant_LO,
      abs_coef_EO = abs(coef_EO),
      abs_coef_LO = abs(coef_LO),
      delta_abs = abs_coef_EO - abs_coef_LO
    )
}

# Calculate the original conditional taxon-bootstrap median interval.
reml_median_ci <- function(delta, seed, n_boot = 10000L) {
  set.seed(seed)
  boot_medians <- replicate(n_boot, median(sample(
    delta,
    size = length(delta),
    replace = TRUE
  )))
  unname(quantile(boot_medians, probs = c(0.025, 0.975), names = FALSE))
}

# Summarise each direction group without changing its taxon set.
summarise_effect_magnitude <- function(shared_df, set_name, estimator, seed_base) {
  groups <- list(
    all = rep(TRUE, nrow(shared_df)),
    positive = shared_df$coef_EO > 0,
    negative = shared_df$coef_EO <
      0
  )
  bind_rows(lapply(
    seq_along(groups),
    function(i) {
      group_name <- names(groups)[i]
      dat <- shared_df[groups[[i]], , drop = FALSE]
      ci <- reml_median_ci(dat$delta_abs, seed = seed_base + i)
      p_wilcox <- if (nrow(dat) >= 2 && any(dat$delta_abs != 0)) {
        as.numeric(wilcox.test(
          abs(dat$coef_EO),
          abs(dat$coef_LO),
          paired = TRUE,
          alternative = "greater",
          exact = TRUE,
          correct = FALSE
        )$p.value)
      } else {
        NA_real_
      }

      tibble(
        set_name = set_name,
        estimator = estimator,
        direction_group = group_name,
        n_species = nrow(dat),
        n_delta_positive = sum(dat$delta_abs > 0),
        median_delta_abs = median(dat$delta_abs),
        bootstrap_ci_lower_95 = ci[1],
        bootstrap_ci_upper_95 = ci[2],
        paired_exact_wilcoxon_one_sided_p = p_wilcox,
        bootstrap_definition = "10,000 taxon-level resamples; conditional on the fixed species set"
      )
    }
  ))
}

hs_eo <- reml_qc(fit_meta.crc_EO, "HS")
hs_lo <- reml_qc(fit_meta.crc_LO, "HS")
reml_eo <- reml_qc(fit_meta.crc_EO_REML, "REML")
reml_lo <- reml_qc(fit_meta.crc_LO_REML, "REML")
hs_shared_qc <- make_shared_qc_table(hs_eo, hs_lo, "HS")
reml_shared_qc <- make_shared_qc_table(reml_eo, reml_lo, "REML")

hs_original27 <- hs_shared_qc %>%
  filter(both_significant, same_direction)

reml_defined_core <- reml_shared_qc %>%
  filter(both_significant, same_direction)

original27_comparison <- hs_original27 %>%
  select(
    feature,
    Species,
    hs_coef_EO = coef_EO,
    hs_stderr_EO = stderr_EO,
    hs_qval_EO = qval_EO,
    hs_tau2_EO = tau2_EO,
    hs_I2_EO = I2_EO,
    hs_pi_lower_EO = pi_lower_EO,
    hs_pi_upper_EO = pi_upper_EO,
    hs_coef_LO = coef_LO,
    hs_stderr_LO = stderr_LO,
    hs_qval_LO = qval_LO,
    hs_tau2_LO = tau2_LO,
    hs_I2_LO = I2_LO,
    hs_pi_lower_LO = pi_lower_LO,
    hs_pi_upper_LO = pi_upper_LO,
    hs_delta_abs = delta_abs
  ) %>%
  left_join(
    reml_eo %>%
      select(
        feature,
        reml_coef_EO = coef,
        reml_stderr_EO = stderr,
        reml_pval_EO = pval,
        reml_qval_EO = qval.fdr,
        reml_tau2_EO = tau2,
        reml_I2_EO = I2,
        reml_max_weight_EO = max_weight,
        reml_consistency_EO = consistency_score,
        reml_pass_qc_EO = pass_qc,
        reml_significant_EO = significant,
        reml_pi_lower_EO = pi_lower_95,
        reml_pi_upper_EO = pi_upper_95,
        reml_pi_excludes_zero_EO = prediction_interval_excludes_zero
      ),
    by = "feature"
  ) %>%
  left_join(
    reml_lo %>%
      select(
        feature,
        reml_coef_LO = coef,
        reml_stderr_LO = stderr,
        reml_pval_LO = pval,
        reml_qval_LO = qval.fdr,
        reml_tau2_LO = tau2,
        reml_I2_LO = I2,
        reml_max_weight_LO = max_weight,
        reml_consistency_LO = consistency_score,
        reml_pass_qc_LO = pass_qc,
        reml_significant_LO = significant,
        reml_pi_lower_LO = pi_lower_95,
        reml_pi_upper_LO = pi_upper_95,
        reml_pi_excludes_zero_LO = prediction_interval_excludes_zero
      ),
    by = "feature"
  ) %>%
  mutate(
    reml_delta_abs = abs(reml_coef_EO) -
      abs(reml_coef_LO),
    delta_abs_change_REML_minus_HS = reml_delta_abs - hs_delta_abs,
    same_direction_EO_HS_vs_REML = sign(hs_coef_EO) ==
      sign(reml_coef_EO),
    same_direction_LO_HS_vs_REML = sign(hs_coef_LO) == sign(reml_coef_LO),
    same_direction_EO_vs_LO_under_REML = sign(reml_coef_EO) ==
      sign(reml_coef_LO),
    retained_both_significant_under_REML = reml_significant_EO & reml_significant_LO,
    retained_shared_core_under_REML = retained_both_significant_under_REML &
      same_direction_EO_vs_LO_under_REML
  ) %>%
  arrange(desc(reml_delta_abs))

original27_hs_for_summary <- original27_comparison %>%
  transmute(
    Species,
    coef_EO = hs_coef_EO,
    coef_LO = hs_coef_LO,
    delta_abs = hs_delta_abs
  )

original27_reml_for_summary <- original27_comparison %>%
  transmute(
    Species,
    coef_EO = reml_coef_EO,
    coef_LO = reml_coef_LO,
    delta_abs = reml_delta_abs
  )

effect_magnitude_summary <- bind_rows(
  summarise_effect_magnitude(
    original27_hs_for_summary,
    "original_HS_27_fixed",
    "HS",
    202608200L
  ),
  summarise_effect_magnitude(
    original27_reml_for_summary,
    "original_HS_27_fixed",
    "REML",
    202608210L
  ),
  summarise_effect_magnitude(
    reml_defined_core,
    "REML_defined_shared_significant",
    "REML",
    202608220L
  )
)

write_csv(
  original27_comparison,
  file.path(results_dir, "01_original27_HS_REML_comparison.csv"),
  na = ""
)

write_csv(
  effect_magnitude_summary,
  file.path(results_dir, "02_effect_magnitude_summary.csv"),
  na = ""
)

write_csv(
  reml_defined_core,
  file.path(results_dir, "03_REML_defined_shared_significant_species.csv"),
  na = ""
)
