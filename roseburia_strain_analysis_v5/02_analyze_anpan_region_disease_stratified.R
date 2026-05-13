#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(anpan)
  library(ape)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(vegan)
  library(permute)
})

env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

env_or_default_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.integer(default))
  }
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed)) {
    stop(sprintf("Environment variable %s must be an integer, got '%s'.", name, value), call. = FALSE)
  }
  parsed
}

env_or_default_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(as.numeric(default))
  }
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed)) {
    stop(sprintf("Environment variable %s must be numeric, got '%s'.", name, value), call. = FALSE)
  }
  parsed
}

env_or_default_bool <- function(name, default) {
  value <- tolower(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) {
    return(isTRUE(default))
  }
  if (value %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (value %in% c("0", "false", "f", "no", "n")) return(FALSE)
  stop(sprintf("Environment variable %s must be boolean, got '%s'.", name, value), call. = FALSE)
}

split_env_values <- function(name, default_values = character()) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) {
    return(default_values)
  }
  values <- trimws(unlist(strsplit(raw, ",", fixed = TRUE)))
  values[nzchar(values)]
}

path_basename <- function(path) {
  basename(gsub("[/\\\\]+$", "", path))
}

find_required_file <- function(directory, pattern, label) {
  hits <- list.files(directory, pattern = pattern, full.names = TRUE)
  if (length(hits) != 1) {
    stop(
      sprintf("Expected exactly one %s file in %s, found %d.", label, directory, length(hits)),
      call. = FALSE
    )
  }
  hits[[1]]
}

compute_within_group_center <- function(values, groups) {
  out <- rep(NA_real_, length(values))
  keep <- !is.na(values) & !is.na(groups)
  if (!any(keep)) {
    return(out)
  }
  group_means <- tapply(values[keep], groups[keep], mean, na.rm = TRUE)
  out[keep] <- values[keep] - unname(group_means[as.character(groups[keep])])
  out
}

drop_unused_factor_levels <- function(df) {
  df |>
    mutate(across(where(is.factor), droplevels))
}

screen_covariates <- function(model_metadata, covariates) {
  keep <- character()
  dropped <- character()

  for (covariate in covariates) {
    if (!covariate %in% names(model_metadata)) {
      dropped <- c(dropped, paste0(covariate, ":missing"))
      next
    }

    values <- model_metadata[[covariate]]
    non_missing <- values[!is.na(values)]

    if (length(non_missing) == 0) {
      dropped <- c(dropped, paste0(covariate, ":all_na"))
      next
    }

    n_unique <- dplyr::n_distinct(non_missing)
    if (n_unique < 2) {
      dropped <- c(dropped, paste0(covariate, ":single_level"))
      next
    }

    keep <- c(keep, covariate)
  }

  list(keep = keep, dropped = dropped)
}

build_formula <- function(response, terms) {
  rhs <- if (length(terms) == 0) "1" else paste(terms, collapse = " + ")
  as.formula(sprintf("%s ~ %s", response, rhs))
}

combine_notes <- function(...) {
  parts <- unlist(list(...), use.names = FALSE)
  parts <- parts[!is.na(parts) & nzchar(parts)]
  if (length(parts) == 0) {
    return(NA_character_)
  }
  paste(parts, collapse = "; ")
}

sanitize_for_tsv <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x[[1]])) {
    return(NA_character_)
  }
  str_squish(gsub("[\r\n\t]+", " ", as.character(x[[1]])))
}

to_tsv_matrix <- function(mat, row_id = "sample_id") {
  as_tibble(
    data.frame(
      setNames(list(rownames(mat)), row_id),
      mat,
      check.names = FALSE
    )
  )
}

extract_param_summary <- function(fit_obj) {
  if (is.null(fit_obj)) {
    return(tibble())
  }
  fit_obj$summary() |>
    as_tibble()
}

null_default <- function(value, default) {
  if (is.null(value)) default else value
}

extract_max_rhat <- function(param_summary) {
  if (nrow(param_summary) == 0 || !"rhat" %in% names(param_summary)) {
    return(NA_real_)
  }
  rhat_vals <- suppressWarnings(as.numeric(param_summary$rhat))
  rhat_vals <- rhat_vals[is.finite(rhat_vals)]
  if (length(rhat_vals) == 0) {
    return(NA_real_)
  }
  max(rhat_vals)
}

build_anpan_diagnostic_note <- function(pglmm_max_rhat, base_max_rhat, threshold = 1.05) {
  notes <- character()
  if (!is.na(pglmm_max_rhat) && pglmm_max_rhat > threshold) {
    notes <- c(notes, sprintf("pglmm_max_rhat_gt_%.2f", threshold))
  }
  if (!is.na(base_max_rhat) && base_max_rhat > threshold) {
    notes <- c(notes, sprintf("base_max_rhat_gt_%.2f", threshold))
  }
  if (length(notes) == 0) {
    NA_character_
  } else {
    paste(notes, collapse = ";")
  }
}

extract_first_matching_value <- function(tbl, patterns, row = 1) {
  if (is.null(tbl) || NROW(tbl) < row) {
    return(NA_real_)
  }
  tbl_df <- as.data.frame(tbl)
  if (nrow(tbl_df) < row || is.null(colnames(tbl_df))) {
    return(NA_real_)
  }
  match_idx <- which(vapply(
    colnames(tbl_df),
    function(nm) any(vapply(patterns, function(pattern) grepl(pattern, nm, ignore.case = TRUE, perl = TRUE), logical(1))),
    logical(1)
  ))
  if (length(match_idx) == 0) {
    return(NA_real_)
  }
  parsed <- suppressWarnings(as.numeric(tbl_df[row, match_idx[[1]], drop = TRUE]))
  if (length(parsed) == 0 || is.na(parsed)) NA_real_ else parsed
}

run_dispersion_test <- function(spec, model_metadata, dist_obj, perm_scheme, test_dir) {
  group_vars <- null_default(spec$dispersion_group_vars, character())

  if (length(group_vars) == 0) {
    summary_row <- tibble(
      test_id = spec$id,
      question = spec$question,
      focus_term = spec$focus_term,
      dispersion_grouping = NA_character_,
      n_used_samples = nrow(model_metadata),
      n_groups = NA_integer_,
      min_group_n = NA_integer_,
      anova_f = NA_real_,
      permutest_f = NA_real_,
      permutest_p = NA_real_,
      note = "dispersion_not_requested"
    )
    write_tsv(summary_row, file.path(test_dir, "dispersion_test_summary.tsv"))
    return(summary_row)
  }

  group_counts <- model_metadata |>
    group_by(across(all_of(group_vars))) |>
    summarise(n = n(), .groups = "drop")
  write_tsv(group_counts, file.path(test_dir, "dispersion_group_counts.tsv"))

  min_group_n <- if (nrow(group_counts) == 0) 0L else min(group_counts$n)
  n_groups <- nrow(group_counts)

  if (n_groups < 2 || min_group_n < 2) {
    summary_row <- tibble(
      test_id = spec$id,
      question = spec$question,
      focus_term = spec$focus_term,
      dispersion_grouping = paste(group_vars, collapse = " x "),
      n_used_samples = nrow(model_metadata),
      n_groups = n_groups,
      min_group_n = min_group_n,
      anova_f = NA_real_,
      permutest_f = NA_real_,
      permutest_p = NA_real_,
      note = "skipped_insufficient_dispersion_groups"
    )
    write_tsv(summary_row, file.path(test_dir, "dispersion_test_summary.tsv"))
    return(summary_row)
  }

  group_factor <- do.call(
    interaction,
    c(as.list(model_metadata[, group_vars, drop = FALSE]), list(drop = TRUE, sep = " | "))
  )

  bd_fit <- tryCatch(
    vegan::betadisper(dist_obj, group_factor, type = "median"),
    error = function(e) e
  )

  if (inherits(bd_fit, "error")) {
    summary_row <- tibble(
      test_id = spec$id,
      question = spec$question,
      focus_term = spec$focus_term,
      dispersion_grouping = paste(group_vars, collapse = " x "),
      n_used_samples = nrow(model_metadata),
      n_groups = n_groups,
      min_group_n = min_group_n,
      anova_f = NA_real_,
      permutest_f = NA_real_,
      permutest_p = NA_real_,
      note = paste("dispersion_fit_failed:", conditionMessage(bd_fit))
    )
    writeLines(conditionMessage(bd_fit), file.path(test_dir, "dispersion_test_error.txt"))
    write_tsv(summary_row, file.path(test_dir, "dispersion_test_summary.tsv"))
    return(summary_row)
  }

  anova_fit <- tryCatch(anova(bd_fit), error = function(e) e)
  perm_fit <- tryCatch(vegan::permutest(bd_fit, permutations = perm_scheme), error = function(e) e)

  if (!inherits(anova_fit, "error")) {
    write_tsv(as_tibble(as.data.frame(anova_fit), rownames = "term"), file.path(test_dir, "dispersion_anova.tsv"))
  }
  if (!inherits(perm_fit, "error") && !is.null(perm_fit$tab)) {
    write_tsv(as_tibble(as.data.frame(perm_fit$tab), rownames = "term"), file.path(test_dir, "dispersion_permutest.tsv"))
  }

  if (inherits(anova_fit, "error") || inherits(perm_fit, "error")) {
    note_parts <- character()
    if (inherits(anova_fit, "error")) {
      note_parts <- c(note_parts, paste("dispersion_anova_failed:", conditionMessage(anova_fit)))
    }
    if (inherits(perm_fit, "error")) {
      note_parts <- c(note_parts, paste("dispersion_permutest_failed:", conditionMessage(perm_fit)))
    }
    summary_row <- tibble(
      test_id = spec$id,
      question = spec$question,
      focus_term = spec$focus_term,
      dispersion_grouping = paste(group_vars, collapse = " x "),
      n_used_samples = nrow(model_metadata),
      n_groups = n_groups,
      min_group_n = min_group_n,
      anova_f = if (!inherits(anova_fit, "error")) extract_first_matching_value(anova_fit, c("^F$", "^F value$")) else NA_real_,
      permutest_f = NA_real_,
      permutest_p = NA_real_,
      note = paste(note_parts, collapse = "; ")
    )
    write_tsv(summary_row, file.path(test_dir, "dispersion_test_summary.tsv"))
    return(summary_row)
  }

  perm_tab <- as.data.frame(perm_fit$tab)
  summary_row <- tibble(
    test_id = spec$id,
    question = spec$question,
    focus_term = spec$focus_term,
    dispersion_grouping = paste(group_vars, collapse = " x "),
    n_used_samples = nrow(model_metadata),
    n_groups = n_groups,
    min_group_n = min_group_n,
    anova_f = extract_first_matching_value(anova_fit, c("^F$", "^F value$")),
    permutest_f = extract_first_matching_value(perm_tab, c("^F$", "^F value$")),
    permutest_p = extract_first_matching_value(perm_tab, c("Pr\\(>F\\)", "^P$")),
    note = NA_character_
  )

  write_tsv(summary_row, file.path(test_dir, "dispersion_test_summary.tsv"))
  summary_row
}

extract_loo_table <- function(result_obj) {
  comparison <- result_obj$loo$comparison
  if (is.null(comparison)) {
    return(tibble())
  }
  as.data.frame(comparison) |>
    rownames_to_column("candidate_model") |>
    as_tibble()
}

derive_phylogeny_gap <- function(loo_tbl) {
  if (nrow(loo_tbl) < 2 || !"candidate_model" %in% names(loo_tbl) || !"elpd_diff" %in% names(loo_tbl)) {
    return(tibble(
      best_model = NA_character_,
      phylogeny_best = NA,
      pglmm_minus_base_elpd = NA_real_,
      pglmm_minus_base_se = NA_real_,
      interpretation = NA_character_
    ))
  }

  pglmm_row <- loo_tbl |> filter(candidate_model == "pglmm_fit")
  base_row <- loo_tbl |> filter(candidate_model == "base_fit")

  if (nrow(pglmm_row) != 1 || nrow(base_row) != 1) {
    return(tibble(
      best_model = loo_tbl$candidate_model[[1]],
      phylogeny_best = loo_tbl$candidate_model[[1]] == "pglmm_fit",
      pglmm_minus_base_elpd = NA_real_,
      pglmm_minus_base_se = NA_real_,
      interpretation = "loo_comparison_unexpected_shape"
    ))
  }

  best_model <- loo_tbl$candidate_model[[1]]
  phylogeny_best <- identical(best_model, "pglmm_fit")
  nonbest_row <- if (phylogeny_best) base_row else pglmm_row
  gap <- if (phylogeny_best) abs(nonbest_row$elpd_diff[[1]]) else -abs(nonbest_row$elpd_diff[[1]])
  se_gap <- abs(nonbest_row$se_diff[[1]])

  tibble(
    best_model = best_model,
    phylogeny_best = phylogeny_best,
    pglmm_minus_base_elpd = gap,
    pglmm_minus_base_se = se_gap,
    interpretation = "pending_effect_size_threshold"
  )
}

classify_phylogeny_gap <- function(gap, se_gap, n_used_samples,
                                   z_threshold = 2,
                                   min_abs_gap = 2,
                                   min_per_sample_gap = 0.01) {
  if (is.na(gap)) {
    return(NA_character_)
  }
  if (is.na(se_gap) || se_gap == 0) {
    return("loo_comparison_no_se")
  }
  if (is.na(n_used_samples) || n_used_samples <= 0) {
    return("loo_comparison_no_sample_size")
  }

  z_like <- abs(gap / se_gap)
  per_sample_gap <- abs(gap) / n_used_samples
  statistically_supported <- z_like > z_threshold
  effect_size_supported <- abs(gap) >= min_abs_gap && per_sample_gap >= min_per_sample_gap

  if (gap > 0) {
    if (statistically_supported && effect_size_supported) {
      "phylogeny_better_clear"
    } else if (statistically_supported) {
      "phylogeny_better_statistically_detectable_but_tiny"
    } else {
      "phylogeny_better_uncertain"
    }
  } else if (gap < 0) {
    if (statistically_supported && effect_size_supported) {
      "base_better_clear"
    } else if (statistically_supported) {
      "base_better_statistically_detectable_but_tiny"
    } else {
      "base_better_uncertain"
    }
  } else {
    "models_tied"
  }
}

build_anpan_model_specs <- function(config) {
  within_region_covariates <- c("Age_within_region_ageclass", "Sex", "BMI", "Cohort")

  list(
    list(
      id = "a01_eu_eo_cc",
      question = "phylogeny ~ CRC_vs_Control within Europe EO",
      outcome = "Diagnosis_crc",
      family = "binomial",
      covariates = within_region_covariates,
      subset_fn = function(meta) filter(meta, region == config$europe_label, Age_class == "EO"),
      outcome_note = "1 = Europe EOCRC, 0 = Europe EOControl"
    ),
    list(
      id = "a02_ea_eo_cc",
      question = "phylogeny ~ CRC_vs_Control within East Asia EO",
      outcome = "Diagnosis_crc",
      family = "binomial",
      covariates = within_region_covariates,
      subset_fn = function(meta) filter(meta, region == config$east_asia_label, Age_class == "EO"),
      outcome_note = "1 = East Asia EOCRC, 0 = East Asia EOControl"
    ),
    list(
      id = "a03_eu_lo_cc",
      question = "phylogeny ~ CRC_vs_Control within Europe LO",
      outcome = "Diagnosis_crc",
      family = "binomial",
      covariates = within_region_covariates,
      subset_fn = function(meta) filter(meta, region == config$europe_label, Age_class == "LO"),
      outcome_note = "1 = Europe LOCRC, 0 = Europe LOControl"
    ),
    list(
      id = "a04_ea_lo_cc",
      question = "phylogeny ~ CRC_vs_Control within East Asia LO",
      outcome = "Diagnosis_crc",
      family = "binomial",
      covariates = within_region_covariates,
      subset_fn = function(meta) filter(meta, region == config$east_asia_label, Age_class == "LO"),
      outcome_note = "1 = East Asia LOCRC, 0 = East Asia LOControl"
    ),
    list(
      id = "a05_eu_crc_eolo",
      question = "phylogeny ~ EO_vs_LO within Europe CRC",
      outcome = "EO_bin",
      family = "binomial",
      covariates = within_region_covariates,
      subset_fn = function(meta) filter(meta, region == config$europe_label, Diagnosis == "CRC"),
      outcome_note = "1 = Europe EOCRC, 0 = Europe LOCRC"
    ),
    list(
      id = "a06_ea_crc_eolo",
      question = "phylogeny ~ EO_vs_LO within East Asia CRC",
      outcome = "EO_bin",
      family = "binomial",
      covariates = within_region_covariates,
      subset_fn = function(meta) filter(meta, region == config$east_asia_label, Diagnosis == "CRC"),
      outcome_note = "1 = East Asia EOCRC, 0 = East Asia LOCRC"
    ),
    list(
      id = "a07_eu_ctrl_eolo",
      question = "phylogeny ~ EO_vs_LO within Europe Control",
      outcome = "EO_bin",
      family = "binomial",
      covariates = within_region_covariates,
      subset_fn = function(meta) filter(meta, region == config$europe_label, Diagnosis == "Control"),
      outcome_note = "1 = Europe EOControl, 0 = Europe LOControl"
    ),
    list(
      id = "a08_ea_ctrl_eolo",
      question = "phylogeny ~ EO_vs_LO within East Asia Control",
      outcome = "EO_bin",
      family = "binomial",
      covariates = within_region_covariates,
      subset_fn = function(meta) filter(meta, region == config$east_asia_label, Diagnosis == "Control"),
      outcome_note = "1 = East Asia EOControl, 0 = East Asia LOControl"
    )
  )
}

build_distance_model_specs <- function(config) {
  list(
    list(
      id = "01_interaction_eo",
      question = "EO patristic distance ~ region * Diagnosis",
      subset_fn = function(meta) filter(meta, Age_class == "EO"),
      formula_terms = c("region * Diagnosis", "Sex", "BMI", "Age_within_region_ageclass"),
      test_type = "interaction",
      focus_term = "region:Diagnosis",
      count_by_vars = c("region", "Diagnosis"),
      expected_n_groups = 4L,
      min_group_n = 3L,
      min_total_n = 12L,
      dispersion_group_vars = c("region", "Diagnosis")
    ),
    list(
      id = "02_interaction_lo",
      question = "LO patristic distance ~ region * Diagnosis",
      subset_fn = function(meta) filter(meta, Age_class == "LO"),
      formula_terms = c("region * Diagnosis", "Sex", "BMI", "Age_within_region_ageclass"),
      test_type = "interaction",
      focus_term = "region:Diagnosis",
      count_by_vars = c("region", "Diagnosis"),
      expected_n_groups = 4L,
      min_group_n = 3L,
      min_total_n = 12L,
      dispersion_group_vars = c("region", "Diagnosis")
    ),
    list(
      id = "03_case_control_europe_eo",
      question = "Europe EO patristic distance ~ Diagnosis",
      subset_fn = function(meta) filter(meta, region == config$europe_label, Age_class == "EO"),
      formula_terms = c("Diagnosis", "Sex", "BMI", "Age_within_region_ageclass", "Cohort"),
      test_type = "case_control",
      focus_term = "Diagnosis",
      count_by_vars = c("Diagnosis"),
      expected_n_groups = 2L,
      min_group_n = 3L,
      min_total_n = 10L,
      dispersion_group_vars = c("Diagnosis")
    ),
    list(
      id = "04_case_control_eastasia_eo",
      question = "East Asia EO patristic distance ~ Diagnosis",
      subset_fn = function(meta) filter(meta, region == config$east_asia_label, Age_class == "EO"),
      formula_terms = c("Diagnosis", "Sex", "BMI", "Age_within_region_ageclass", "Cohort"),
      test_type = "case_control",
      focus_term = "Diagnosis",
      count_by_vars = c("Diagnosis"),
      expected_n_groups = 2L,
      min_group_n = 3L,
      min_total_n = 10L,
      dispersion_group_vars = c("Diagnosis")
    ),
    list(
      id = "05_case_control_europe_lo",
      question = "Europe LO patristic distance ~ Diagnosis",
      subset_fn = function(meta) filter(meta, region == config$europe_label, Age_class == "LO"),
      formula_terms = c("Diagnosis", "Sex", "BMI", "Age_within_region_ageclass", "Cohort"),
      test_type = "case_control",
      focus_term = "Diagnosis",
      count_by_vars = c("Diagnosis"),
      expected_n_groups = 2L,
      min_group_n = 3L,
      min_total_n = 10L,
      dispersion_group_vars = c("Diagnosis")
    ),
    list(
      id = "06_case_control_eastasia_lo",
      question = "East Asia LO patristic distance ~ Diagnosis",
      subset_fn = function(meta) filter(meta, region == config$east_asia_label, Age_class == "LO"),
      formula_terms = c("Diagnosis", "Sex", "BMI", "Age_within_region_ageclass", "Cohort"),
      test_type = "case_control",
      focus_term = "Diagnosis",
      count_by_vars = c("Diagnosis"),
      expected_n_groups = 2L,
      min_group_n = 3L,
      min_total_n = 10L,
      dispersion_group_vars = c("Diagnosis")
    ),
    list(
      id = "07_interaction_europe_age",
      question = "Europe patristic distance ~ Age_class * Diagnosis",
      subset_fn = function(meta) filter(meta, region == config$europe_label),
      formula_terms = c("Age_class * Diagnosis", "Sex", "BMI", "Age_within_region_ageclass", "Cohort"),
      test_type = "formal_age_contrast",
      focus_term = "Age_class:Diagnosis",
      count_by_vars = c("Age_class", "Diagnosis"),
      expected_n_groups = 4L,
      min_group_n = 3L,
      min_total_n = 12L,
      dispersion_group_vars = c("Age_class", "Diagnosis")
    ),
    list(
      id = "08_interaction_eastasia_age",
      question = "East Asia patristic distance ~ Age_class * Diagnosis",
      subset_fn = function(meta) filter(meta, region == config$east_asia_label),
      formula_terms = c("Age_class * Diagnosis", "Sex", "BMI", "Age_within_region_ageclass", "Cohort"),
      test_type = "formal_age_contrast",
      focus_term = "Age_class:Diagnosis",
      count_by_vars = c("Age_class", "Diagnosis"),
      expected_n_groups = 4L,
      min_group_n = 3L,
      min_total_n = 12L,
      dispersion_group_vars = c("Age_class", "Diagnosis")
    ),
    list(
      id = "09_interaction_3way_global",
      question = "Global patristic distance ~ region * Age_class * Diagnosis",
      subset_fn = function(meta) meta,
      formula_terms = c("region * Age_class * Diagnosis", "Sex", "BMI", "Age_within_region_ageclass"),
      test_type = "formal_age_contrast",
      focus_term = "region:Age_class:Diagnosis",
      count_by_vars = c("region", "Age_class", "Diagnosis"),
      expected_n_groups = 8L,
      min_group_n = 3L,
      min_total_n = 24L,
      dispersion_group_vars = c("region", "Age_class", "Diagnosis")
    ),
    list(
      id = "10_eolo_europe_crc",
      question = "Europe CRC patristic distance ~ Age_class",
      subset_fn = function(meta) filter(meta, region == config$europe_label, Diagnosis == "CRC"),
      formula_terms = c("Age_class", "Sex", "BMI", "Age_within_region_ageclass", "Cohort"),
      test_type = "eo_lo_contrast",
      focus_term = "Age_class",
      count_by_vars = c("Age_class"),
      expected_n_groups = 2L,
      min_group_n = 3L,
      min_total_n = 10L,
      dispersion_group_vars = c("Age_class")
    ),
    list(
      id = "11_eolo_eastasia_crc",
      question = "East Asia CRC patristic distance ~ Age_class",
      subset_fn = function(meta) filter(meta, region == config$east_asia_label, Diagnosis == "CRC"),
      formula_terms = c("Age_class", "Sex", "BMI", "Age_within_region_ageclass", "Cohort"),
      test_type = "eo_lo_contrast",
      focus_term = "Age_class",
      count_by_vars = c("Age_class"),
      expected_n_groups = 2L,
      min_group_n = 3L,
      min_total_n = 10L,
      dispersion_group_vars = c("Age_class")
    ),
    list(
      id = "12_eolo_europe_control",
      question = "Europe Control patristic distance ~ Age_class",
      subset_fn = function(meta) filter(meta, region == config$europe_label, Diagnosis == "Control"),
      formula_terms = c("Age_class", "Sex", "BMI", "Age_within_region_ageclass", "Cohort"),
      test_type = "eo_lo_contrast",
      focus_term = "Age_class",
      count_by_vars = c("Age_class"),
      expected_n_groups = 2L,
      min_group_n = 3L,
      min_total_n = 10L,
      dispersion_group_vars = c("Age_class")
    ),
    list(
      id = "13_eolo_eastasia_control",
      question = "East Asia Control patristic distance ~ Age_class",
      subset_fn = function(meta) filter(meta, region == config$east_asia_label, Diagnosis == "Control"),
      formula_terms = c("Age_class", "Sex", "BMI", "Age_within_region_ageclass", "Cohort"),
      test_type = "eo_lo_contrast",
      focus_term = "Age_class",
      count_by_vars = c("Age_class"),
      expected_n_groups = 2L,
      min_group_n = 3L,
      min_total_n = 10L,
      dispersion_group_vars = c("Age_class")
    )
  )
}

prepare_region_metadata <- function(config, tree_file) {
  meta_env <- new.env(parent = emptyenv())
  load(config$metadata_rdata, envir = meta_env)

  if (!exists(config$metadata_object, envir = meta_env, inherits = FALSE)) {
    stop(sprintf("metadata.Rdata does not contain %s.", config$metadata_object), call. = FALSE)
  }

  meta_raw <- get(config$metadata_object, envir = meta_env)
  if (!is.data.frame(meta_raw)) {
    stop(sprintf("%s is not a data.frame/tibble.", config$metadata_object), call. = FALSE)
  }

  required_cols <- c(
    config$sample_id_col,
    config$age_col,
    config$age_class_col,
    config$sex_col,
    config$bmi_col,
    config$country_col,
    config$continent_col,
    config$diagnosis_col,
    config$group_col,
    config$cohort_col
  )
  missing_cols <- setdiff(required_cols, names(meta_raw))
  if (length(missing_cols) > 0) {
    stop(
      sprintf("Metadata object %s is missing required columns: %s", config$metadata_object, paste(missing_cols, collapse = ", ")),
      call. = FALSE
    )
  }

  tree <- ape::read.tree(tree_file)
  tip_tbl <- tibble(tree_tip_label = tree$tip.label) |>
    mutate(
      is_reference_tip = str_detect(tree_tip_label, config$reference_tip_pattern),
      sample_id = str_remove(tree_tip_label, config$trim_pattern)
    )

  duplicate_tip_ids <- tip_tbl |>
    filter(!is_reference_tip, !is.na(sample_id), nzchar(sample_id)) |>
    count(sample_id, name = "n") |>
    filter(n > 1)
  if (nrow(duplicate_tip_ids) > 0) {
    stop(
      sprintf(
        "Tree tips are not uniquely identifiable after trim_pattern. Example duplicated Sample_IDs: %s",
        paste(head(duplicate_tip_ids$sample_id, 10), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  meta_tbl_raw <- as_tibble(meta_raw) |>
    transmute(
      sample_id = as.character(.data[[config$sample_id_col]]),
      Age = as.numeric(.data[[config$age_col]]),
      Age_class = as.character(.data[[config$age_class_col]]),
      Sex = as.character(.data[[config$sex_col]]),
      BMI = as.numeric(.data[[config$bmi_col]]),
      Country = as.character(.data[[config$country_col]]),
      Continent = as.character(.data[[config$continent_col]]),
      Diagnosis = as.character(.data[[config$diagnosis_col]]),
      Group = as.character(.data[[config$group_col]]),
      Cohort = as.character(.data[[config$cohort_col]])
    )

  duplicate_meta_ids <- meta_tbl_raw |>
    filter(!is.na(sample_id), nzchar(sample_id)) |>
    count(sample_id, name = "n") |>
    filter(n > 1)
  if (nrow(duplicate_meta_ids) > 0) {
    stop(
      sprintf(
        "Metadata Sample_ID is not globally unique. Use cohort-qualified sample IDs before joining tree tips. Example duplicated IDs: %s",
        paste(head(duplicate_meta_ids$sample_id, 10), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  meta_tbl <- meta_tbl_raw |>
    filter(!is.na(sample_id), nzchar(sample_id))

  audit_tbl <- tip_tbl |>
    filter(!is_reference_tip) |>
    left_join(meta_tbl, by = "sample_id") |>
    mutate(
      has_metadata = !is.na(Country),
      region = case_when(
        Country %in% config$east_asia_countries ~ config$east_asia_label,
        length(config$europe_countries) > 0 & Country %in% config$europe_countries ~ config$europe_label,
        Continent %in% config$europe_continents ~ config$europe_label,
        TRUE ~ NA_character_
      ),
      is_target_region = !is.na(region),
      is_target_diagnosis = Diagnosis %in% config$included_diagnoses,
      complete_covariates = !is.na(Age) & !is.na(BMI) & nzchar(Sex) & Age_class %in% c("EO", "LO"),
      exclusion_reason = case_when(
        !has_metadata ~ "no_metadata",
        !is_target_region ~ "non_target_region",
        !is_target_diagnosis ~ "excluded_diagnosis",
        !complete_covariates ~ "missing_covariates",
        TRUE ~ "retained"
      )
    )

  retained_tbl <- audit_tbl |>
    filter(exclusion_reason == "retained") |>
    mutate(
      region = factor(region, levels = c(config$europe_label, config$east_asia_label)),
      Diagnosis = factor(Diagnosis, levels = c("Control", "CRC")),
      Age_class = factor(Age_class, levels = c("LO", "EO")),
      Sex = factor(Sex),
      Cohort = factor(Cohort),
      Age_within_region = compute_within_group_center(Age, region),
      Age_within_region_ageclass = compute_within_group_center(Age, interaction(region, Age_class, drop = TRUE)),
      EastAsia_bin = as.integer(region == config$east_asia_label),
      Diagnosis_crc = as.integer(Diagnosis == "CRC"),
      EO_bin = as.integer(Age_class == "EO"),
      Sample_ID = sample_id
    ) |>
    arrange(sample_id)

  list(
    tree = tree,
    tip_tbl = tip_tbl,
    audit_tbl = audit_tbl,
    retained_tbl = retained_tbl
  )
}

build_retained_tree <- function(tree, retained_labels, reference_tip_pattern) {
  sample_tree <- drop.tip(tree, tree$tip.label[str_detect(tree$tip.label, reference_tip_pattern)])
  drop.tip(sample_tree, setdiff(sample_tree$tip.label, retained_labels)) |>
    ladderize(right = FALSE)
}

collapse_lines <- function(tbl, empty_line = "- none") {
  if (nrow(tbl) == 0 || !"line" %in% names(tbl)) {
    return(empty_line)
  }
  paste(tbl$line, collapse = "\n")
}

run_single_anpan_model <- function(spec, tree_file, base_metadata, config) {
  model_dir <- file.path(config$anpan_dir, spec$id)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

  model_metadata <- spec$subset_fn(base_metadata) |>
    arrange(sample_id) |>
    drop_unused_factor_levels()

  covariate_screen <- screen_covariates(model_metadata, spec$covariates)
  used_covariates <- covariate_screen$keep
  dropped_covariates <- covariate_screen$dropped

  write_tsv(model_metadata, file.path(model_dir, "input_metadata_pre_overlap.tsv"))

  result <- tryCatch(
    anpan::anpan_pglmm(
      meta_file = model_metadata,
      tree_file = tree_file,
      outcome = spec$outcome,
      covariates = used_covariates,
      out_dir = model_dir,
      trim_pattern = config$trim_pattern,
      bug_name = paste(config$bug_name, spec$id, sep = "__"),
      omit_na = config$omit_na,
      ladderize = config$ladderize,
      family = spec$family,
      show_plot_cor_mat = config$show_plot_cor_mat,
      show_plot_tree = config$show_plot_tree,
      show_post = config$show_post,
      show_yrep = config$show_yrep,
      save_object = config$save_object,
      verbose = config$verbose,
      loo_comparison = config$loo_comparison,
      run_diagnostics = config$run_diagnostics,
      reg_noise = FALSE,
      reg_gamma_params = config$reg_gamma_params,
      plot_ext = config$plot_ext,
      int_prior_scale = config$int_prior_scale,
      sigma_phylo_scale = config$sigma_phylo_scale,
      parallel_chains = config$parallel_chains,
      iter_warmup = config$iter_warmup,
      iter_sampling = config$iter_sampling,
      seed = config$seed,
      refresh = config$refresh,
      adapt_delta = config$adapt_delta,
      max_treedepth = config$max_treedepth
    ),
    error = function(e) e
  )

  if (inherits(result, "error")) {
    error_tbl <- tibble(
      model_id = spec$id,
      question = spec$question,
      family = spec$family,
      outcome = spec$outcome,
      outcome_note = spec$outcome_note,
      candidate_covariates = paste(spec$covariates, collapse = " + "),
      used_covariates = paste(used_covariates, collapse = " + "),
      dropped_covariates = paste(dropped_covariates, collapse = "; "),
      n_input_samples = nrow(model_metadata),
      n_used_samples = NA_integer_,
      n_tree_tips_used = NA_integer_,
      n_cohorts = dplyr::n_distinct(model_metadata$Cohort),
      error_message = sanitize_for_tsv(conditionMessage(result)),
      best_model = NA_character_,
      phylogeny_best = NA,
      pglmm_minus_base_elpd = NA_real_,
      pglmm_minus_base_elpd_per_sample = NA_real_,
      pglmm_minus_base_se = NA_real_,
      interpretation = "model_failed",
      pglmm_max_rhat = NA_real_,
      base_max_rhat = NA_real_,
      diagnostic_note = NA_character_
    )
    writeLines(conditionMessage(result), file.path(model_dir, "model_error.txt"))
    write_tsv(error_tbl, file.path(model_dir, "model_summary.tsv"))
    return(error_tbl)
  }

  loo_tbl <- extract_loo_table(result)
  loo_gap <- derive_phylogeny_gap(loo_tbl)
  pglmm_summary <- extract_param_summary(result$pglmm_fit)
  base_summary <- extract_param_summary(result$base_fit)
  pglmm_max_rhat <- extract_max_rhat(pglmm_summary)
  base_max_rhat <- extract_max_rhat(base_summary)
  diagnostic_note <- build_anpan_diagnostic_note(pglmm_max_rhat, base_max_rhat)

  write_tsv(as_tibble(result$model_input), file.path(model_dir, "input_metadata_used.tsv"))
  write_tsv(to_tsv_matrix(result$cor_mat), file.path(model_dir, "phylo_correlation_matrix.tsv"))
  if (nrow(loo_tbl) > 0) {
    write_tsv(loo_tbl, file.path(model_dir, "loo_comparison.tsv"))
  }
  write_tsv(pglmm_summary, file.path(model_dir, "pglmm_parameter_summary.tsv"))
  if (nrow(base_summary) > 0) {
    write_tsv(base_summary, file.path(model_dir, "base_parameter_summary.tsv"))
  }
  saveRDS(result, file.path(model_dir, "anpan_result_object.rds"))

  summary_row <- tibble(
    model_id = spec$id,
    question = spec$question,
    family = spec$family,
    outcome = spec$outcome,
    outcome_note = spec$outcome_note,
    candidate_covariates = paste(spec$covariates, collapse = " + "),
    used_covariates = paste(used_covariates, collapse = " + "),
    dropped_covariates = paste(dropped_covariates, collapse = "; "),
    n_input_samples = nrow(model_metadata),
    n_used_samples = nrow(result$model_input),
    n_tree_tips_used = nrow(result$cor_mat),
    n_cohorts = dplyr::n_distinct(model_metadata$Cohort),
    error_message = NA_character_
  ) |>
    bind_cols(loo_gap) |>
    mutate(
      pglmm_minus_base_elpd_per_sample = if_else(
        !is.na(pglmm_minus_base_elpd) & n_used_samples > 0,
        pglmm_minus_base_elpd / n_used_samples,
        NA_real_
      ),
      pglmm_minus_base_z = if_else(
        !is.na(pglmm_minus_base_elpd) & !is.na(pglmm_minus_base_se) & pglmm_minus_base_se > 0,
        abs(pglmm_minus_base_elpd / pglmm_minus_base_se),
        NA_real_
      ),
      interpretation = mapply(
        classify_phylogeny_gap,
        pglmm_minus_base_elpd,
        pglmm_minus_base_se,
        n_used_samples,
        USE.NAMES = FALSE
      ),
      loo_effect_threshold_note = "clear requires |ELPD gap| >= 2, |ELPD gap|/sample >= 0.01, and |gap/SE| > 2",
      pglmm_max_rhat = pglmm_max_rhat,
      base_max_rhat = base_max_rhat,
      diagnostic_note = diagnostic_note
    )

  write_tsv(summary_row, file.path(model_dir, "model_summary.tsv"))
  summary_row
}

run_single_distance_test <- function(spec, base_metadata, cophenetic_mat, config) {
  test_dir <- file.path(config$distance_dir, spec$id)
  dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)

  model_metadata <- spec$subset_fn(base_metadata) |>
    arrange(sample_id) |>
    drop_unused_factor_levels()

  formula_terms <- spec$formula_terms
  dropped_formula_terms <- character()
  if ("Cohort" %in% formula_terms && dplyr::n_distinct(model_metadata$Cohort) < 2) {
    formula_terms <- setdiff(formula_terms, "Cohort")
    dropped_formula_terms <- c(dropped_formula_terms, "Cohort:single_level")
  }
  formula_obj <- build_formula("dist_obj", formula_terms)
  formula_text <- paste(deparse(formula_obj), collapse = " ")
  note_prefix <- if (length(dropped_formula_terms) == 0) NA_character_ else paste(dropped_formula_terms, collapse = "; ")

  write_tsv(model_metadata, file.path(test_dir, "input_metadata_used.tsv"))

  count_by_vars <- null_default(spec$count_by_vars, c("region", "Diagnosis"))
  grouping_label <- paste(count_by_vars, collapse = " x ")
  group_counts <- model_metadata |>
    group_by(across(all_of(count_by_vars))) |>
    summarise(n = n(), .groups = "drop")
  write_tsv(group_counts, file.path(test_dir, "group_counts.tsv"))

  min_group_n <- if (nrow(group_counts) == 0) 0L else min(group_counts$n)
  n_groups_observed <- nrow(group_counts)
  expected_n_groups <- null_default(spec$expected_n_groups, n_groups_observed)
  min_group_required <- null_default(spec$min_group_n, 2L)
  min_total_n <- null_default(spec$min_total_n, 8L)
  cohort_n <- dplyr::n_distinct(model_metadata$Cohort)
  failed_checks <- character()

  if (nrow(model_metadata) < min_total_n) {
    failed_checks <- c(failed_checks, sprintf("n_lt_%d", min_total_n))
  }
  if (n_groups_observed < expected_n_groups) {
    failed_checks <- c(failed_checks, sprintf("observed_groups_lt_%d", expected_n_groups))
  }
  if (min_group_n < min_group_required) {
    failed_checks <- c(failed_checks, sprintf("min_group_n_lt_%d", min_group_required))
  }
  if (cohort_n < 2) {
    failed_checks <- c(failed_checks, "cohorts_lt_2")
  }

  if (length(failed_checks) > 0) {
    summary_row <- tibble(
      test_id = spec$id,
      question = spec$question,
      test_type = spec$test_type,
      focus_term = spec$focus_term,
      n_used_samples = nrow(model_metadata),
      n_cohorts = cohort_n,
      grouping_vars = grouping_label,
      expected_n_groups = expected_n_groups,
      n_groups_observed = n_groups_observed,
      min_group_n = min_group_n,
      formula = formula_text,
      term = NA_character_,
      df = NA_real_,
      sum_of_sqs = NA_real_,
      r2 = NA_real_,
      f_stat = NA_real_,
      p_value = NA_real_,
      note = combine_notes(note_prefix, "skipped_insufficient_samples", paste(failed_checks, collapse = "; "))
    )
    dispersion_row <- tibble(
      test_id = spec$id,
      question = spec$question,
      focus_term = spec$focus_term,
      dispersion_grouping = paste(null_default(spec$dispersion_group_vars, character()), collapse = " x "),
      n_used_samples = nrow(model_metadata),
      n_groups = n_groups_observed,
      min_group_n = min_group_n,
      anova_f = NA_real_,
      permutest_f = NA_real_,
      permutest_p = NA_real_,
      note = "dispersion_not_run_distance_skipped"
    )
    write_tsv(summary_row, file.path(test_dir, "distance_test_summary.tsv"))
    write_tsv(dispersion_row, file.path(test_dir, "dispersion_test_summary.tsv"))
    return(list(distance_summary = summary_row, dispersion_summary = dispersion_row))
  }

  labels <- model_metadata$tree_tip_label
  dist_obj <- as.dist(cophenetic_mat[labels, labels, drop = FALSE])
  perm_scheme <- permute::how(nperm = config$distance_n_permutations, blocks = model_metadata$Cohort)

  adonis_fit <- tryCatch(
    vegan::adonis2(
      formula = formula_obj,
      data = model_metadata,
      permutations = perm_scheme,
      by = "margin"
    ),
    error = function(e) e
  )

  if (inherits(adonis_fit, "error")) {
    summary_row <- tibble(
      test_id = spec$id,
      question = spec$question,
      test_type = spec$test_type,
      focus_term = spec$focus_term,
      n_used_samples = nrow(model_metadata),
      n_cohorts = cohort_n,
      grouping_vars = grouping_label,
      expected_n_groups = expected_n_groups,
      n_groups_observed = n_groups_observed,
      min_group_n = min_group_n,
      formula = formula_text,
      term = NA_character_,
      df = NA_real_,
      sum_of_sqs = NA_real_,
      r2 = NA_real_,
      f_stat = NA_real_,
      p_value = NA_real_,
      note = combine_notes(note_prefix, paste("distance_test_failed:", conditionMessage(adonis_fit)))
    )
    dispersion_row <- tibble(
      test_id = spec$id,
      question = spec$question,
      focus_term = spec$focus_term,
      dispersion_grouping = paste(null_default(spec$dispersion_group_vars, character()), collapse = " x "),
      n_used_samples = nrow(model_metadata),
      n_groups = n_groups_observed,
      min_group_n = min_group_n,
      anova_f = NA_real_,
      permutest_f = NA_real_,
      permutest_p = NA_real_,
      note = "dispersion_not_run_distance_failed"
    )
    writeLines(conditionMessage(adonis_fit), file.path(test_dir, "distance_test_error.txt"))
    write_tsv(summary_row, file.path(test_dir, "distance_test_summary.tsv"))
    write_tsv(dispersion_row, file.path(test_dir, "dispersion_test_summary.tsv"))
    return(list(distance_summary = summary_row, dispersion_summary = dispersion_row))
  }

  adonis_tbl <- as.data.frame(adonis_fit) |>
    rownames_to_column("term") |>
    as_tibble() |>
    rename(
      df = Df,
      sum_of_sqs = SumOfSqs,
      r2 = R2,
      f_stat = F,
      p_value = `Pr(>F)`
    )

  write_tsv(adonis_tbl, file.path(test_dir, "distance_test_table.tsv"))

  summary_tbl <- adonis_tbl |>
    filter(term != "Total") |>
    mutate(
      test_id = spec$id,
      question = spec$question,
      test_type = spec$test_type,
      focus_term = spec$focus_term,
      n_used_samples = nrow(model_metadata),
      n_cohorts = cohort_n,
      grouping_vars = grouping_label,
      expected_n_groups = expected_n_groups,
      n_groups_observed = n_groups_observed,
      min_group_n = min_group_n,
      formula = formula_text,
      note = note_prefix
    ) |>
    select(test_id, question, test_type, focus_term, n_used_samples, n_cohorts, grouping_vars, expected_n_groups, n_groups_observed, min_group_n, formula, term, df, sum_of_sqs, r2, f_stat, p_value, note)

  write_tsv(summary_tbl, file.path(test_dir, "distance_test_summary.tsv"))
  dispersion_tbl <- run_dispersion_test(spec, model_metadata, dist_obj, perm_scheme, test_dir)

  list(distance_summary = summary_tbl, dispersion_summary = dispersion_tbl)
}

main <- function() {
  project_root <- normalizePath(env_or_default("PROJECT_ROOT", getwd()), winslash = "/", mustWork = FALSE)
  default_strain_dir <- file.path(project_root, "data/strainphlan/Europe_Asia_t__SGB4925")
  default_output_root <- file.path(project_root, "results/strain_analysis_v5/rf4925_fig3")

  config <- list(
    strain_dir = env_or_default("STRAINPHLAN_STRAIN_DIR", default_strain_dir),
    tree_file = env_or_default("STRAINPHLAN_TREE_FILE", ""),
    metadata_rdata = env_or_default("STRAINPHLAN_METADATA_RDATA", file.path(project_root, "metadata.Rdata")),
    metadata_object = env_or_default("STRAINPHLAN_REGION_METADATA_OBJECT", "metadata_allcohorts"),
    sample_id_col = env_or_default("STRAINPHLAN_REGION_SAMPLE_ID_COL", "Sample_ID"),
    age_col = env_or_default("STRAINPHLAN_REGION_AGE_COL", "Age"),
    age_class_col = env_or_default("STRAINPHLAN_REGION_AGE_CLASS_COL", "Age_class"),
    sex_col = env_or_default("STRAINPHLAN_REGION_SEX_COL", "Sex"),
    bmi_col = env_or_default("STRAINPHLAN_REGION_BMI_COL", "BMI"),
    country_col = env_or_default("STRAINPHLAN_REGION_COUNTRY_COL", "Country"),
    continent_col = env_or_default("STRAINPHLAN_REGION_CONTINENT_COL", "Continent"),
    diagnosis_col = env_or_default("STRAINPHLAN_REGION_DIAGNOSIS_COL", "Diagnosis"),
    group_col = env_or_default("STRAINPHLAN_REGION_GROUP_COL", "Group"),
    cohort_col = env_or_default("STRAINPHLAN_REGION_COHORT_COL", "Cohort"),
    east_asia_countries = split_env_values("STRAINPHLAN_EAST_ASIA_COUNTRIES", c("China", "Japan")),
    europe_countries = split_env_values("STRAINPHLAN_EUROPE_COUNTRIES", character()),
    europe_continents = split_env_values("STRAINPHLAN_EUROPE_CONTINENTS", c("Europe")),
    included_diagnoses = split_env_values("STRAINPHLAN_INCLUDED_DIAGNOSES", c("CRC", "Control")),
    east_asia_label = env_or_default("STRAINPHLAN_EAST_ASIA_LABEL", "East_Asia"),
    europe_label = env_or_default("STRAINPHLAN_EUROPE_LABEL", "Europe"),
    output_root = env_or_default("STRAINPHLAN_REGION_DISEASE_OUTPUT_ROOT", default_output_root),
    trim_pattern = env_or_default("STRAINPHLAN_TRIM_PATTERN", "\\.temp$|\\.json\\.bz2$|\\.json$|\\.bz2$"),
    reference_tip_pattern = env_or_default("STRAINPHLAN_REFERENCE_TIP_PATTERN", "^(GCF_|GCA_)"),
    distance_n_permutations = env_or_default_int("STRAINPHLAN_DISTANCE_N_PERM", 999),
    omit_na = env_or_default_bool("STRAINPHLAN_ANPAN_OMIT_NA", TRUE),
    ladderize = env_or_default_bool("STRAINPHLAN_ANPAN_LADDERIZE", TRUE),
    show_plot_cor_mat = env_or_default_bool("STRAINPHLAN_ANPAN_SHOW_COR_MAT", FALSE),
    show_plot_tree = env_or_default_bool("STRAINPHLAN_ANPAN_SHOW_TREE", FALSE),
    show_post = env_or_default_bool("STRAINPHLAN_ANPAN_SHOW_POST", FALSE),
    show_yrep = env_or_default_bool("STRAINPHLAN_ANPAN_SHOW_YREP", FALSE),
    save_object = env_or_default_bool("STRAINPHLAN_ANPAN_SAVE_OBJECT", TRUE),
    verbose = env_or_default_bool("STRAINPHLAN_ANPAN_VERBOSE", TRUE),
    loo_comparison = env_or_default_bool("STRAINPHLAN_ANPAN_LOO", TRUE),
    run_diagnostics = env_or_default_bool("STRAINPHLAN_ANPAN_DIAGNOSTICS", TRUE),
    reg_gamma_params = c(
      env_or_default_num("STRAINPHLAN_ANPAN_REG_GAMMA_SHAPE", 1),
      env_or_default_num("STRAINPHLAN_ANPAN_REG_GAMMA_RATE", 2)
    ),
    plot_ext = env_or_default("STRAINPHLAN_ANPAN_PLOT_EXT", "pdf"),
    int_prior_scale = env_or_default_num("STRAINPHLAN_ANPAN_INT_PRIOR_SCALE", 1),
    sigma_phylo_scale = env_or_default_num("STRAINPHLAN_ANPAN_SIGMA_PHYLO_SCALE", 0.333),
    parallel_chains = env_or_default_int("STRAINPHLAN_ANPAN_PARALLEL_CHAINS", 4),
    iter_warmup = env_or_default_int("STRAINPHLAN_ANPAN_ITER_WARMUP", 1000),
    iter_sampling = env_or_default_int("STRAINPHLAN_ANPAN_ITER_SAMPLING", 1000),
    seed = env_or_default_int("STRAINPHLAN_ANPAN_SEED", 20260421),
    refresh = env_or_default_int("STRAINPHLAN_ANPAN_REFRESH", 100),
    adapt_delta = env_or_default_num("STRAINPHLAN_ANPAN_ADAPT_DELTA", 0.95),
    max_treedepth = env_or_default_int("STRAINPHLAN_ANPAN_MAX_TREEDEPTH", 12)
  )

  config$species_id <- path_basename(config$strain_dir)
  extracted_short <- str_extract(config$species_id, "t__SGB[0-9]+")
  config$species_short <- env_or_default(
    "STRAINPHLAN_SPECIES_SHORT_ID",
    ifelse(is.na(extracted_short) || !nzchar(extracted_short), config$species_id, extracted_short)
  )
  config$bug_name <- env_or_default("STRAINPHLAN_ANPAN_BUG_NAME", "rf4925")
  config$output_dir <- file.path(config$output_root, config$species_short)
  config$anpan_dir <- file.path(config$output_dir, "a")
  config$distance_dir <- file.path(config$output_dir, "distance_tests")

  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(config$anpan_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(config$distance_dir, recursive = TRUE, showWarnings = FALSE)

  tree_file <- if (nzchar(config$tree_file)) {
    config$tree_file
  } else {
    find_required_file(config$strain_dir, "^RAxML_bestTree\\..*\\.tre$", "RAxML best tree")
  }

  prep <- prepare_region_metadata(config, tree_file)
  retained_tree <- build_retained_tree(prep$tree, prep$retained_tbl$tree_tip_label, config$reference_tip_pattern)
  cophenetic_mat <- cophenetic.phylo(retained_tree)

  anpan_specs <- build_anpan_model_specs(config)
  distance_specs <- build_distance_model_specs(config)

  exclusion_summary <- prep$audit_tbl |>
    count(exclusion_reason, sort = TRUE)
  retained_region_diagnosis <- prep$retained_tbl |>
    count(region, Diagnosis, sort = FALSE)
  retained_region_diagnosis_ageclass <- prep$retained_tbl |>
    count(region, Diagnosis, Age_class, sort = FALSE)
  retained_country_counts <- prep$retained_tbl |>
    count(region, Country, sort = TRUE)

  config_tbl <- tibble(
    strain_dir = config$strain_dir,
    tree_file = tree_file,
    metadata_rdata = config$metadata_rdata,
    metadata_object = config$metadata_object,
    east_asia_countries = paste(config$east_asia_countries, collapse = ","),
    europe_countries = paste(config$europe_countries, collapse = ","),
    europe_continents = paste(config$europe_continents, collapse = ","),
    included_diagnoses = paste(config$included_diagnoses, collapse = ","),
    output_dir = config$output_dir,
    trim_pattern = config$trim_pattern,
    reference_tip_pattern = config$reference_tip_pattern,
    distance_n_permutations = config$distance_n_permutations,
    cmdstan_path = tryCatch(cmdstanr::cmdstan_path(), error = function(e) paste("ERROR:", e$message)),
    parallel_chains = config$parallel_chains,
    iter_warmup = config$iter_warmup,
    iter_sampling = config$iter_sampling,
    seed = config$seed,
    adapt_delta = config$adapt_delta,
    max_treedepth = config$max_treedepth
  )

  write_tsv(config_tbl, file.path(config$output_dir, "run_config.tsv"))
  write_tsv(prep$audit_tbl, file.path(config$output_dir, "sample_filter_audit.tsv"))
  write_tsv(exclusion_summary, file.path(config$output_dir, "sample_exclusion_summary.tsv"))
  write_tsv(prep$retained_tbl, file.path(config$output_dir, "retained_metadata.tsv"))
  write_tsv(retained_region_diagnosis, file.path(config$output_dir, "retained_region_diagnosis_counts.tsv"))
  write_tsv(retained_region_diagnosis_ageclass, file.path(config$output_dir, "retained_region_diagnosis_ageclass_counts.tsv"))
  write_tsv(retained_country_counts, file.path(config$output_dir, "retained_country_counts.tsv"))
  write_tsv(to_tsv_matrix(cophenetic_mat, row_id = "tree_tip_label"), file.path(config$output_dir, "retained_cophenetic_matrix.tsv"))

  distance_results <- lapply(
    distance_specs,
    run_single_distance_test,
    base_metadata = prep$retained_tbl,
    cophenetic_mat = cophenetic_mat,
    config = config
  )
  distance_summaries <- bind_rows(lapply(distance_results, `[[`, "distance_summary"))
  dispersion_summaries <- bind_rows(lapply(distance_results, `[[`, "dispersion_summary"))
  write_tsv(distance_summaries, file.path(config$output_dir, "distance_analysis_summary.tsv"))
  write_tsv(dispersion_summaries, file.path(config$output_dir, "dispersion_analysis_summary.tsv"))

  anpan_summaries <- bind_rows(lapply(anpan_specs, run_single_anpan_model, tree_file = tree_file, base_metadata = prep$retained_tbl, config = config))
  write_tsv(anpan_summaries, file.path(config$output_dir, "anpan_analysis_summary.tsv"))

  interaction_summary <- distance_summaries |>
    filter(test_id %in% c("01_interaction_eo", "02_interaction_lo"), term == focus_term)

  formal_age_contrast_summary <- distance_summaries |>
    filter(test_id %in% c("07_interaction_europe_age", "08_interaction_eastasia_age", "09_interaction_3way_global"), term == focus_term)

  case_control_distance_summary <- distance_summaries |>
    filter(test_type == "case_control", term == focus_term)

  eo_lo_distance_summary <- distance_summaries |>
    filter(test_type == "eo_lo_contrast", term == focus_term)

  eo_lo_anpan_summary <- anpan_summaries |>
    filter(model_id %in% c("a05_eu_crc_eolo", "a06_ea_crc_eolo", "a07_eu_ctrl_eolo", "a08_ea_ctrl_eolo"))

  write_tsv(interaction_summary, file.path(config$output_dir, "eo_lo_region_interaction_summary.tsv"))
  write_tsv(formal_age_contrast_summary, file.path(config$output_dir, "formal_age_contrast_summary.tsv"))
  write_tsv(case_control_distance_summary, file.path(config$output_dir, "case_control_distance_focus_summary.tsv"))
  write_tsv(eo_lo_distance_summary, file.path(config$output_dir, "eo_lo_distance_focus_summary.tsv"))
  write_tsv(eo_lo_anpan_summary, file.path(config$output_dir, "eo_lo_anpan_focus_summary.tsv"))

  summary_lines <- c(
    sprintf("Species directory: %s", config$strain_dir),
    sprintf("Tree file: %s", tree_file),
    sprintf("Tree sample tips (non-reference): %d", sum(!prep$tip_tbl$is_reference_tip)),
    sprintf("Tree reference tips: %d", sum(prep$tip_tbl$is_reference_tip)),
    sprintf("Matched to metadata: %d", sum(prep$audit_tbl$has_metadata)),
    sprintf("Retained after cleaning: %d", nrow(prep$retained_tbl)),
    "",
    "Retained counts by region, diagnosis and age class:",
    collapse_lines(
      retained_region_diagnosis_ageclass |>
        mutate(line = sprintf("- %s | %s | %s | n=%d", region, Diagnosis, Age_class, n)) |>
        select(line)
    ),
    "",
    "EO/LO interaction tests on phylogenetic distance (focus on region:Diagnosis term):",
    collapse_lines(
      interaction_summary |>
        mutate(
          line = sprintf(
            "- %s | focus=%s | n=%d | cohorts=%d | R2=%s | F=%s | p=%s | note=%s",
            test_id,
            focus_term,
            n_used_samples,
            n_cohorts,
            ifelse(is.na(r2), "NA", sprintf("%.4f", r2)),
            ifelse(is.na(f_stat), "NA", sprintf("%.4f", f_stat)),
            ifelse(is.na(p_value), "NA", sprintf("%.4g", p_value)),
            ifelse(is.na(note), "NA", note)
          )
        ) |>
        select(line)
    ),
    "",
    "Formal EO-vs-LO contrast tests (focus on Age_class:Diagnosis or region:Age_class:Diagnosis):",
    collapse_lines(
      formal_age_contrast_summary |>
        mutate(
          line = sprintf(
            "- %s | focus=%s | n=%d | cohorts=%d | R2=%s | F=%s | p=%s | note=%s",
            test_id,
            focus_term,
            n_used_samples,
            n_cohorts,
            ifelse(is.na(r2), "NA", sprintf("%.4f", r2)),
            ifelse(is.na(f_stat), "NA", sprintf("%.4f", f_stat)),
            ifelse(is.na(p_value), "NA", sprintf("%.4g", p_value)),
            ifelse(is.na(note), "NA", note)
          )
        ) |>
        select(line)
    ),
    "",
    "Region-specific case-control distance tests (Diagnosis term):",
    collapse_lines(
      case_control_distance_summary |>
        mutate(
          line = sprintf(
            "- %s | focus=%s | n=%d | cohorts=%d | R2=%s | F=%s | p=%s | note=%s",
            test_id,
            focus_term,
            n_used_samples,
            n_cohorts,
            ifelse(is.na(r2), "NA", sprintf("%.4f", r2)),
            ifelse(is.na(f_stat), "NA", sprintf("%.4f", f_stat)),
            ifelse(is.na(p_value), "NA", sprintf("%.4g", p_value)),
            ifelse(is.na(note), "NA", note)
          )
        ) |>
        select(line)
    ),
    "",
    "Region-specific EO/LO distance tests (Age_class term within diagnosis strata):",
    collapse_lines(
      eo_lo_distance_summary |>
        mutate(
          line = sprintf(
            "- %s | focus=%s | n=%d | cohorts=%d | R2=%s | F=%s | p=%s | note=%s",
            test_id,
            focus_term,
            n_used_samples,
            n_cohorts,
            ifelse(is.na(r2), "NA", sprintf("%.4f", r2)),
            ifelse(is.na(f_stat), "NA", sprintf("%.4f", f_stat)),
            ifelse(is.na(p_value), "NA", sprintf("%.4g", p_value)),
            ifelse(is.na(note), "NA", note)
          )
        ) |>
        select(line)
    ),
    "",
    "Dispersion checks paired with distance models (betadisper/permutest):",
    collapse_lines(
      dispersion_summaries |>
        mutate(
          line = sprintf(
            "- %s | grouping=%s | n=%d | groups=%s | min_group_n=%s | perm_p=%s | note=%s",
            test_id,
            ifelse(is.na(dispersion_grouping) | !nzchar(dispersion_grouping), "NA", dispersion_grouping),
            n_used_samples,
            ifelse(is.na(n_groups), "NA", as.character(n_groups)),
            ifelse(is.na(min_group_n), "NA", as.character(min_group_n)),
            ifelse(is.na(permutest_p), "NA", sprintf("%.4g", permutest_p)),
            ifelse(is.na(note), "NA", note)
          )
        ) |>
        select(line)
    ),
    "",
    "Region-specific case-control Anpan summary:",
    collapse_lines(
      anpan_summaries |>
        mutate(
          line = sprintf(
            "- %s | used_n=%s | cohorts=%s | covariates=%s | best=%s | phylogeny_best=%s | elpd_gap=%s | elpd_gap_per_sample=%s | pglmm_max_rhat=%s | base_max_rhat=%s | interpretation=%s | note=%s",
            model_id,
            ifelse(is.na(n_used_samples), "NA", as.character(n_used_samples)),
            ifelse(is.na(n_cohorts), "NA", as.character(n_cohorts)),
            ifelse(is.na(used_covariates) | !nzchar(used_covariates), "1", used_covariates),
            ifelse(is.na(best_model), "NA", best_model),
            ifelse(is.na(phylogeny_best), "NA", as.character(phylogeny_best)),
            ifelse(is.na(pglmm_minus_base_elpd), "NA", sprintf("%.3f", pglmm_minus_base_elpd)),
            ifelse(is.na(pglmm_minus_base_elpd_per_sample), "NA", sprintf("%.4f", pglmm_minus_base_elpd_per_sample)),
            ifelse(is.na(pglmm_max_rhat), "NA", sprintf("%.3f", pglmm_max_rhat)),
            ifelse(is.na(base_max_rhat), "NA", sprintf("%.3f", base_max_rhat)),
            ifelse(is.na(interpretation), "NA", interpretation),
            ifelse(is.na(diagnostic_note), "NA", diagnostic_note)
          )
        ) |>
        select(line)
    ),
    "",
    "Region-specific EO/LO Anpan summary:",
    collapse_lines(
      eo_lo_anpan_summary |>
        mutate(
          line = sprintf(
            "- %s | used_n=%s | cohorts=%s | covariates=%s | best=%s | phylogeny_best=%s | elpd_gap=%s | elpd_gap_per_sample=%s | pglmm_max_rhat=%s | base_max_rhat=%s | interpretation=%s | note=%s",
            model_id,
            ifelse(is.na(n_used_samples), "NA", as.character(n_used_samples)),
            ifelse(is.na(n_cohorts), "NA", as.character(n_cohorts)),
            ifelse(is.na(used_covariates) | !nzchar(used_covariates), "1", used_covariates),
            ifelse(is.na(best_model), "NA", best_model),
            ifelse(is.na(phylogeny_best), "NA", as.character(phylogeny_best)),
            ifelse(is.na(pglmm_minus_base_elpd), "NA", sprintf("%.3f", pglmm_minus_base_elpd)),
            ifelse(is.na(pglmm_minus_base_elpd_per_sample), "NA", sprintf("%.4f", pglmm_minus_base_elpd_per_sample)),
            ifelse(is.na(pglmm_max_rhat), "NA", sprintf("%.3f", pglmm_max_rhat)),
            ifelse(is.na(base_max_rhat), "NA", sprintf("%.3f", base_max_rhat)),
            ifelse(is.na(interpretation), "NA", interpretation),
            ifelse(is.na(diagnostic_note), "NA", diagnostic_note)
          )
        ) |>
        select(line)
    ),
    "",
    "Interpretation note:",
    "- Within-region models now include Cohort as an explicit fixed effect where estimable; global cross-region interaction models still avoid direct cohort adjustment because region and cohort are partially nested."
  )

  writeLines(summary_lines, file.path(config$output_dir, "analysis_summary.txt"))
  message("Saved Fig3-oriented region-disease strain results to: ", config$output_dir)
}

main()
