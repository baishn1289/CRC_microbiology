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

age_dir <- file.path(sensitivity_root, "Age_exclusion_analysis/results")
loco_dir <- file.path(sensitivity_root, "leave_out_cohort_analysis/results")
dir.create(age_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(loco_dir, recursive = TRUE, showWarnings = FALSE)

# Summarise a fixed taxon set using the original age-exclusion tests.
summarise_age_effect_set <- function(values, analysis_set, repetitions) {
  values <- values[is.finite(values)]
  nonzero_values <- values[values != 0]
  exact_possible <- length(values) > 0 &&
    length(values) <= 50 &&
    length(nonzero_values) == length(values) &&
    !anyDuplicated(abs(nonzero_values))
  wilcoxon_one_sided <- if (length(nonzero_values) > 0) {
    wilcox.test(
      values,
      mu = 0,
      alternative = "greater",
      exact = exact_possible,
      correct = FALSE
    )
  } else {
    NULL
  }

  wilcoxon_two_sided <- if (length(nonzero_values) > 0) {
    wilcox.test(
      values,
      mu = 0,
      alternative = "two.sided",
      exact = exact_possible,
      correct = FALSE
    )
  } else {
    NULL
  }

  sign_test <- if (length(nonzero_values) > 0) {
    binom.test(
      sum(nonzero_values > 0),
      length(nonzero_values),
      p = 0.5,
      alternative = "greater"
    )
  } else {
    NULL
  }

  bootstrap_ci <- if (length(values) > 1) {
    bootstrap_medians <- replicate(
      repetitions,
      median(sample(values, length(values), replace = TRUE))
    )
    as.numeric(quantile(
      bootstrap_medians,
      probs = c(0.025, 0.975),
      names = FALSE,
      type = 6
    ))
  } else if (length(values) == 1) {
    rep(values, 2)
  } else {
    c(NA_real_, NA_real_)
  }

  tibble(
    analysis_set = analysis_set,
    n_species = length(values),
    n_delta_positive = sum(values > 0),
    proportion_delta_positive = if_else(
      length(values) >
        0,
      mean(values > 0),
      NA_real_
    ),
    median_delta_abs = if_else(length(values) > 0, median(values), NA_real_),
    mean_delta_abs = if_else(length(values) > 0, mean(values), NA_real_),
    conditional_taxon_bootstrap_ci_low = bootstrap_ci[1],
    conditional_taxon_bootstrap_ci_high = bootstrap_ci[2],
    bootstrap_repetitions = repetitions,
    wilcoxon_p_one_sided = if (is.null(wilcoxon_one_sided)) {
      NA_real_
    } else {
      wilcoxon_one_sided$p.value
    },
    wilcoxon_p_two_sided = if (is.null(wilcoxon_two_sided)) {
      NA_real_
    } else {
      wilcoxon_two_sided$p.value
    },
    wilcoxon_method = if_else(exact_possible, "exact", "asymptotic_no_cc"),
    sign_test_p_one_sided = if (is.null(sign_test))
      NA_real_
    else sign_test$p.value
  )
}

# Pair significant EO and LO taxa with concordant directions.
select_strict_core <- function(processed_eo, processed_lo, beta_prefix) {
  significant_eo <- processed_eo %>%
    filter(pass_qc, qval_fdr < 0.1) %>%
    transmute(
      feature,
      Species = str_remove(feature, ".*s__"),
      beta_EO = coef,
      se_EO = stderr,
      p_EO = pval,
      q_EO = qval_fdr,
      k_EO = k,
      I2_EO = I2,
      max_weight_EO = max_weight,
      consistency_EO = consistency_score,
      qc_status_EO = qc_status
    )
  significant_lo <- processed_lo %>%
    filter(pass_qc, qval_fdr < 0.1) %>%
    transmute(
      feature,
      beta_LO = coef,
      se_LO = stderr,
      p_LO = pval,
      q_LO = qval_fdr,
      k_LO = k,
      I2_LO = I2,
      max_weight_LO = max_weight,
      consistency_LO = consistency_score,
      qc_status_LO = qc_status
    )
  inner_join(significant_eo, significant_lo, by = "feature") %>%
    filter(sign(beta_EO) == sign(beta_LO)) %>%
    mutate(
      abs_beta_EO = abs(beta_EO),
      abs_beta_LO = abs(beta_LO),
      delta_abs = abs_beta_EO - abs_beta_LO,
      effect_direction = if_else(
        beta_EO >
          0,
        "CRC-enriched",
        "CRC-depleted"
      ),
      analysis_source = beta_prefix
    ) %>%
    arrange(Species)
}

main_eo <- sensitivity_qc(fit_meta.crc_EO) %>%
  mutate(qval_fdr = qval.fdr)
main_lo <- sensitivity_qc(fit_meta.crc_LO) %>%
  mutate(qval_fdr = qval.fdr)
main_core <- select_strict_core(main_eo, main_lo, "main")

# Exclude ages 48-51 in both disease groups, then refit and reselect.
age_metadata <- metadata_allcohorts %>%
  filter(Group %in% c("Control", "CRC"), Age <= 47 | Age >= 52) %>%
  mutate(Age_class = if_else(Age <= 47, "EO", "LO")) %>%
  sensitivity_metadata()
eo <- age_metadata %>%
  filter(Age_class == "EO")
lo <- age_metadata %>%
  filter(Age_class == "LO")
abundance <- run_consensus_EO_LO_combine(
  meta_list_metadata_intersect_samples,
  eo,
  lo,
  min_cohorts = 3,
  min_prevalence = 0.05,
  MIN_STEP1_COUNT = 3,
  abund_threshold = 0.001
)
age_fit_eo <- sensitivity_fit(eo, abundance$EO, file.path(sensitivity_root, "models/age_EO"))
age_fit_lo <- sensitivity_fit(lo, abundance$LO, file.path(sensitivity_root, "models/age_LO"))
age_eo <- sensitivity_qc(age_fit_eo) %>%
  mutate(qval_fdr = qval.fdr)
age_lo <- sensitivity_qc(age_fit_lo) %>%
  mutate(qval_fdr = qval.fdr)
age_core <- select_strict_core(age_eo, age_lo, "age_exclusion")
fixed_core <- main_core %>%
  select(feature, Species) %>%
  left_join(age_eo %>%
    select(feature, beta_EO = coef), by = "feature") %>%
  left_join(age_lo %>%
    select(feature, beta_LO = coef), by = "feature") %>%
  mutate(delta_abs = abs(beta_EO) - abs(beta_LO))
set.seed(20260828L)
age_summary <- bind_rows(
  summarise_age_effect_set(main_core$delta_abs, "main_reselected_core_age_lt50_vs_ge50", 10000L),
  summarise_age_effect_set(
    age_core$delta_abs,
    "sensitivity_reselected_core_age_le47_vs_ge52",
    10000L
  ),
  summarise_age_effect_set(fixed_core$delta_abs, "sensitivity_fixed_original_27_species", 10000L)
)
write_csv(age_summary, file.path(age_dir, "00_key_summary.csv"))
write_csv(age_core, file.path(age_dir, "reselected_species.csv"))
write_csv(fixed_core, file.path(age_dir, "fixed_species.csv"))

# Each shared cohort is removed from both age strata; the 27-species target stays fixed.
md <- sensitivity_metadata(metadata_allcohorts)
shared_cohorts <- intersect(
  names(fit_meta.crc_EO$maaslin_fits),
  names(fit_meta.crc_LO$maaslin_fits)
)
cohort_rows <- list()

for (cohort in shared_cohorts) {
  eo <- md %>%
    filter(Age_class == "EO", Cohort != cohort)
  lo <- md %>%
    filter(Age_class == "LO", Cohort != cohort)
  abundance <- run_consensus_EO_LO_combine(
    meta_list_metadata_intersect_samples,
    eo,
    lo,
    min_cohorts = 3,
    min_prevalence = 0.05,
    MIN_STEP1_COUNT = 3,
    abund_threshold = 0.001
  )
  eo_fit <- sensitivity_fit(
    eo,
    abundance$EO,
    file.path(sensitivity_root, "models/LOCO", cohort, "EO")
  )
  lo_fit <- sensitivity_fit(
    lo,
    abundance$LO,
    file.path(sensitivity_root, "models/LOCO", cohort, "LO")
  )
  fixed <- main_core %>%
    select(feature) %>%
    left_join(eo_fit$meta_fits %>%
      select(feature, coef_EO = coef), by = "feature") %>%
    left_join(lo_fit$meta_fits %>%
      select(feature, coef_LO = coef), by = "feature")
  stopifnot(nrow(fixed) == 27L, all(is.finite(fixed$coef_EO)), all(is.finite(fixed$coef_LO)))
  delta <- median(abs(fixed$coef_EO) - abs(fixed$coef_LO))
  cohort_rows[[cohort]] <- tibble(
    omitted_cohort = cohort,
    n_taxa = nrow(fixed),
    median_delta_abs = delta,
    change_from_baseline = delta - median(main_core$delta_abs),
    absolute_change_from_baseline = abs(delta - median(main_core$delta_abs))
  )
}

write_csv(
  bind_rows(cohort_rows) %>%
    arrange(desc(absolute_change_from_baseline)),
  file.path(loco_dir, "fixed27_cohort_influence_ranking.csv")
)
