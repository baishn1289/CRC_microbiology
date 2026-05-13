#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(anpan)
  library(ape)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
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

sanitize_for_tsv <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x[[1]])) {
    return(NA_character_)
  }
  str_squish(gsub("[\r\n\t]+", " ", as.character(x[[1]])))
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

build_model_specs <- function() {
  list(
    list(
      id = "01_phylogeny_region_overall",
      question = "phylogeny ~ Europe_vs_EastAsia overall",
      outcome = "EastAsia_bin",
      family = "binomial",
      covariates = c("Diagnosis_crc", "EO_bin", "Age_within_region_ageclass", "Sex", "BMI"),
      subset_fn = function(meta) meta,
      outcome_note = "1 = East_Asia, 0 = Europe"
    ),
    list(
      id = "02_phylogeny_region_within_crc",
      question = "phylogeny ~ Europe_vs_EastAsia within CRC",
      outcome = "EastAsia_bin",
      family = "binomial",
      covariates = c("EO_bin", "Age_within_region_ageclass", "Sex", "BMI"),
      subset_fn = function(meta) filter(meta, Diagnosis == "CRC"),
      outcome_note = "1 = East_Asia CRC, 0 = Europe CRC"
    ),
    list(
      id = "03_phylogeny_region_within_control",
      question = "phylogeny ~ Europe_vs_EastAsia within Control",
      outcome = "EastAsia_bin",
      family = "binomial",
      covariates = c("EO_bin", "Age_within_region_ageclass", "Sex", "BMI"),
      subset_fn = function(meta) filter(meta, Diagnosis == "Control"),
      outcome_note = "1 = East_Asia Control, 0 = Europe Control"
    ),
    list(
      id = "04_phylogeny_region_within_eocrc",
      question = "phylogeny ~ Europe_vs_EastAsia within EOCRC",
      outcome = "EastAsia_bin",
      family = "binomial",
      covariates = c("Age_within_region_ageclass", "Sex", "BMI"),
      subset_fn = function(meta) filter(meta, Diagnosis == "CRC", Age_class == "EO"),
      outcome_note = "1 = East_Asia EOCRC, 0 = Europe EOCRC"
    ),
    list(
      id = "05_phylogeny_region_within_locrc",
      question = "phylogeny ~ Europe_vs_EastAsia within LOCRC",
      outcome = "EastAsia_bin",
      family = "binomial",
      covariates = c("Age_within_region_ageclass", "Sex", "BMI"),
      subset_fn = function(meta) filter(meta, Diagnosis == "CRC", Age_class == "LO"),
      outcome_note = "1 = East_Asia LOCRC, 0 = Europe LOCRC"
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

run_single_model <- function(spec, tree_file, base_metadata, config) {
  model_dir <- file.path(config$output_dir, spec$id)
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

  model_metadata <- spec$subset_fn(base_metadata) |>
    arrange(sample_id)

  write_tsv(model_metadata, file.path(model_dir, "input_metadata_pre_overlap.tsv"))

  result <- tryCatch(
    anpan::anpan_pglmm(
      meta_file = model_metadata,
      tree_file = tree_file,
      outcome = spec$outcome,
      covariates = spec$covariates,
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
      covariates = paste(spec$covariates, collapse = " + "),
      n_input_samples = nrow(model_metadata),
      n_used_samples = NA_integer_,
      n_tree_tips_used = NA_integer_,
      error_message = sanitize_for_tsv(conditionMessage(result)),
      best_model = NA_character_,
      phylogeny_best = NA,
      pglmm_minus_base_elpd = NA_real_,
      pglmm_minus_base_se = NA_real_,
      interpretation = "model_failed"
    )
    writeLines(conditionMessage(result), file.path(model_dir, "model_error.txt"))
    write_tsv(error_tbl, file.path(model_dir, "model_summary.tsv"))
    return(error_tbl)
  }

  loo_tbl <- extract_loo_table(result)
  loo_gap <- derive_phylogeny_gap(loo_tbl)
  pglmm_summary <- extract_param_summary(result$pglmm_fit)
  base_summary <- extract_param_summary(result$base_fit)

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
    covariates = paste(spec$covariates, collapse = " + "),
    n_input_samples = nrow(model_metadata),
    n_used_samples = nrow(result$model_input),
    n_tree_tips_used = nrow(result$cor_mat),
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
      loo_effect_threshold_note = "clear requires |ELPD gap| >= 2, |ELPD gap|/sample >= 0.01, and |gap/SE| > 2"
    )

  write_tsv(summary_row, file.path(model_dir, "model_summary.tsv"))
  summary_row
}

main <- function() {
  project_root <- normalizePath(env_or_default("PROJECT_ROOT", getwd()), winslash = "/", mustWork = FALSE)
  default_strain_dir <- file.path(project_root, "data/strainphlan/Europe_Asia_t__SGB4925")
  default_output_root <- file.path(project_root, "results/strain_analysis_v5/anpan_region_comparison")

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
    output_root = env_or_default("STRAINPHLAN_REGION_OUTPUT_ROOT", default_output_root),
    trim_pattern = env_or_default("STRAINPHLAN_TRIM_PATTERN", "\\.temp$|\\.json\\.bz2$|\\.json$|\\.bz2$"),
    reference_tip_pattern = env_or_default("STRAINPHLAN_REFERENCE_TIP_PATTERN", "^(GCF_|GCA_)"),
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
  config$bug_name <- env_or_default("STRAINPHLAN_ANPAN_BUG_NAME", config$species_id)
  config$output_dir <- file.path(config$output_root, config$species_id)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  tree_file <- if (nzchar(config$tree_file)) {
    config$tree_file
  } else {
    find_required_file(config$strain_dir, "^RAxML_bestTree\\..*\\.tre$", "RAxML best tree")
  }

  prep <- prepare_region_metadata(config, tree_file)
  model_specs <- build_model_specs()
  selected_model_ids <- split_env_values("STRAINPHLAN_REGION_MODEL_IDS", character())
  if (length(selected_model_ids) > 0) {
    model_specs <- Filter(function(spec) spec$id %in% selected_model_ids, model_specs)
    if (length(model_specs) == 0) {
      stop("STRAINPHLAN_REGION_MODEL_IDS was set, but none of the requested ids matched.", call. = FALSE)
    }
  }

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

  model_summaries <- bind_rows(lapply(model_specs, run_single_model, tree_file = tree_file, base_metadata = prep$retained_tbl, config = config))
  write_tsv(model_summaries, file.path(config$output_dir, "analysis_summary.tsv"))

  summary_lines <- c(
    sprintf("Species directory: %s", config$strain_dir),
    sprintf("Tree file: %s", tree_file),
    sprintf("Tree sample tips (non-reference): %d", sum(!prep$tip_tbl$is_reference_tip)),
    sprintf("Tree reference tips: %d", sum(prep$tip_tbl$is_reference_tip)),
    sprintf("Matched to metadata: %d", sum(prep$audit_tbl$has_metadata)),
    sprintf("Retained after cleaning: %d", nrow(prep$retained_tbl)),
    sprintf("East Asia countries: %s", paste(config$east_asia_countries, collapse = ", ")),
    sprintf("Included diagnoses: %s", paste(config$included_diagnoses, collapse = ", ")),
    "",
    "Retained counts by region, diagnosis and age class:",
    paste(
      retained_region_diagnosis_ageclass |>
        mutate(line = sprintf("- %s | %s | %s | n=%d", region, Diagnosis, Age_class, n)) |>
        pull(line),
      collapse = "\n"
    ),
    "",
    "LOO-based summary:",
    paste(
      model_summaries |>
        mutate(
          line = sprintf(
            "- %s | used_n=%s | best=%s | phylogeny_best=%s | pglmm_minus_base_elpd=%s | interpretation=%s",
            model_id,
            ifelse(is.na(n_used_samples), "NA", as.character(n_used_samples)),
            ifelse(is.na(best_model), "NA", best_model),
            ifelse(is.na(phylogeny_best), "NA", as.character(phylogeny_best)),
            ifelse(is.na(pglmm_minus_base_elpd), "NA", sprintf("%.3f", pglmm_minus_base_elpd)),
            ifelse(is.na(interpretation), "NA", interpretation)
          )
        ) |>
        pull(line),
      collapse = "\n"
    )
  )

  writeLines(summary_lines, file.path(config$output_dir, "analysis_summary.txt"))
  message("Saved region-comparison Anpan results to: ", config$output_dir)
}

main()
