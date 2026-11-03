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

analysis_dir <- file.path(sensitivity_root, "cross_fitting")

source(file.path(code_dir, "05.cross_fitting_functions.R"))

results_dir <- file.path(analysis_dir, "results")
mmuphin_output_dir <- file.path(analysis_dir, "models")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

bootstrap_repetitions <- 10000L
random_seed <- 20260820L

abundance_sample_ids <- unique(unlist(lapply(
  meta_list_metadata_intersect_samples,
  function(cohort_data) setdiff(colnames(cohort_data), "clade_name")
)))

base_metadata <- metadata_allcohorts %>%
  filter(
    Group %in% c("CRC", "Control"),
    Age_class %in% c("EO", "LO"),
    Sample_ID %in% abundance_sample_ids
  ) %>%
  mutate(
    Age = suppressWarnings(as.numeric(as.character(Age))),
    BMI = suppressWarnings(as.numeric(as.character(BMI))),
    Sex = as.factor(Sex),
    Group = factor(as.character(Group), levels = c("Control", "CRC")),
    covariate_complete = is.finite(Age) &
      !is.na(Sex) & is.finite(BMI)
  )

if (anyDuplicated(base_metadata$Sample_ID)) {
  stop("metadata 中存在重复 Sample_ID，不能继续 cross-fitting。")
}

fold_assignment_base <- tribble(
  ~Cohort,                 ~fold, ~source_group,
  "FUSCC-SHSD",            1L,    "FUSCC-SHSD",
  "Feng_2015",             1L,    "Feng_2015",
  "ThomasAM_2018b",        1L,    "ThomasAM_2018b",
  "Zeller_2014",           1L,    "Zeller_2014",
  "Gupta_2019",            2L,    "Gupta_2019",
  "Hannigan_2018",         2L,    "Hannigan_2018",
  "Liu_2022",              2L,    "Liu_2022",
  "YachidaS_2019",         2L,    "YachidaS_2019",
  "Piccinno_2025_IIGM_CZ", 3L,    "Piccinno_2025_IIGM",
  "Piccinno_2025_IIGM_IT", 3L,    "Piccinno_2025_IIGM",
  "Piccinno_2025_IIGM_TU", 3L,    "Piccinno_2025_IIGM",
  "Vogtmann_2016",         3L,    "Vogtmann_2016",
  "Wirbel_2019",           3L,    "Wirbel_2019",
  "Yang_2020",             3L,    "Yang_2020",
  "Yu_2017",               3L,    "Yu_2017"
)

if (!setequal(fold_assignment_base$Cohort, unique(base_metadata$Cohort))) {
  stop("固定 fold 与输入数据中的 cohort 名称不完全一致。")
}

if (anyDuplicated(fold_assignment_base$Cohort)) {
  stop("固定 fold 中存在重复 cohort。")
}

fold_counts <- base_metadata %>%
  count(Cohort, Age_class, Group, name = "n_raw") %>%
  pivot_wider(
    names_from = c(Age_class, Group),
    values_from = n_raw,
    values_fill = 0,
    names_glue = "{Age_class}{Group}_n"
  )

complete_counts <- base_metadata %>%
  filter(covariate_complete) %>%
  count(Cohort, Age_class, Group, name = "n_complete") %>%
  pivot_wider(
    names_from = c(Age_class, Group),
    values_from = n_complete,
    values_fill = 0,
    names_glue = "{Age_class}{Group}_complete_n"
  )

geography <- base_metadata %>%
  group_by(Cohort) %>%
  summarise(
    country = paste(sort(unique(Country)), collapse = ";"),
    continent = paste(sort(unique(Continent)), collapse = ";"),
    .groups = "drop"
  )

fold_assignment <- fold_assignment_base %>%
  left_join(fold_counts, by = "Cohort") %>%
  left_join(complete_counts, by = "Cohort") %>%
  left_join(geography, by = "Cohort") %>%
  mutate(
    sequencing_platform = NA_character_,
    platform_note = "not available in the supplied Rdata objects"
  ) %>%
  arrange(fold, Cohort)

write_tsv(fold_assignment, file.path(results_dir, "fold_assignment.tsv"), na = "NA")

analysis_designs <- c("main15", "common12")
all_effect_results <- list()
all_tier_effect_results <- list()
all_round_metadata <- list()
all_fit_objects <- list()

for (design_name in analysis_designs) {
  for (validation_fold in 1:3) {
    round_label <- paste0(design_name, "_fold", validation_fold)
    validation_cohorts <- fold_assignment_base %>%
      filter(fold == validation_fold) %>%
      pull(Cohort)
    discovery_cohorts <- setdiff(fold_assignment_base$Cohort, validation_cohorts)
    discovery_metadata_all <- base_metadata %>%
      filter(Cohort %in% discovery_cohorts)
    validation_metadata_all <- base_metadata %>%
      filter(Cohort %in% validation_cohorts)
    discovery_eligibility <- get_eligible_cohorts(discovery_metadata_all)
    validation_eligibility <- get_eligible_cohorts(validation_metadata_all)
    discovery_eo_cohorts <- discovery_eligibility$eo
    validation_eo_cohorts <- validation_eligibility$eo
    if (design_name == "main15") {
      discovery_lo_cohorts <- discovery_eligibility$lo
      validation_lo_cohorts <- validation_eligibility$lo
    } else {
      discovery_lo_cohorts <- intersect(discovery_eligibility$lo, discovery_eo_cohorts)
      validation_lo_cohorts <- intersect(validation_eligibility$lo, validation_eo_cohorts)
    }

    discovery_eo <- discovery_metadata_all %>%
      filter(Age_class == "EO", Cohort %in% discovery_eo_cohorts)
    discovery_lo <- discovery_metadata_all %>%
      filter(Age_class == "LO", Cohort %in% discovery_lo_cohorts)
    validation_eo <- validation_metadata_all %>%
      filter(Age_class == "EO", Cohort %in% validation_eo_cohorts)
    validation_lo <- validation_metadata_all %>%
      filter(Age_class == "LO", Cohort %in% validation_lo_cohorts)
    consensus <- run_discovery_consensus(
      meta_list = meta_list_metadata_intersect_samples,
      discovery_eo = discovery_eo,
      discovery_lo = discovery_lo,
      min_prevalence = 0.05,
      min_cohorts = 3,
      min_detected_samples = 3,
      abundance_threshold = 0.001
    )
    if (nrow(consensus$merged) == 0) {
      stop("Round ", round_label, " discovery consensus 为空。")
    }

    consensus_matrix <- consensus$merged %>%
      column_to_rownames("clade_name") %>%
      as.matrix()
    discovery_fit_eo <- fit_mmuphin_meta(
      metadata = discovery_eo,
      abundance_matrix = consensus_matrix,
      output_label = paste0(round_label, "_discovery_EO")
    )
    discovery_fit_lo <- fit_mmuphin_meta(
      metadata = discovery_lo,
      abundance_matrix = consensus_matrix,
      output_label = paste0(round_label, "_discovery_LO")
    )
    discovery_processed_eo <- process_discovery_qc(discovery_fit_eo)
    discovery_processed_lo <- process_discovery_qc(discovery_fit_lo)
    tier_species <- build_discovery_tier_species(discovery_processed_eo, discovery_processed_lo) %>%
      mutate(
        design = design_name,
        validation_fold = validation_fold,
        .before = 1
      )
    selected_species <- tier_species %>%
      filter(tier_strict_core) %>%
      select(-starts_with("tier_"))
    write_tsv(
      discovery_processed_eo,
      file.path(results_dir, paste0(round_label, "_discovery_EO_meta.tsv")),
      na = "NA"
    )
    write_tsv(
      discovery_processed_lo,
      file.path(results_dir, paste0(round_label, "_discovery_LO_meta.tsv")),
      na = "NA"
    )
    write_tsv(
      selected_species,
      file.path(results_dir, paste0(round_label, "_selected_species.tsv")),
      na = "NA"
    )
    write_tsv(
      tier_species,
      file.path(results_dir, paste0(round_label, "_three_tier_selected_species.tsv")),
      na = "NA"
    )
    if (nrow(tier_species) == 0) {
      effect_results <- tibble(
        design = character(),
        validation_fold = integer(),
        feature = character(),
        Species = character(),
        discovery_direction = double(),
        discovery_direction_label = character(),
        beta_EO = double(),
        beta_LO = double(),
        paired_estimate_available = logical(),
        delta_abs = double(),
        delta_aligned = double()
      )
      tier_effect_results <- tibble()
      validation_meta_eo <- tibble()
      validation_meta_lo <- tibble()
      validation_meta_eo_all <- tibble()
      validation_meta_lo_all <- tibble()
      local_validation <- tibble()
      local_validation_all <- tibble()
      validation_fit_eo <- NULL
      validation_fit_lo <- NULL
    } else {
      selected_matrix <- consensus_matrix[tier_species$feature, , drop = FALSE]
      validation_fit_eo <- fit_mmuphin_meta(
        metadata = validation_eo,
        abundance_matrix = selected_matrix,
        output_label = paste0(round_label, "_validation_EO")
      )
      validation_fit_lo <- fit_mmuphin_meta(
        metadata = validation_lo,
        abundance_matrix = selected_matrix,
        output_label = paste0(round_label, "_validation_LO")
      )
      validation_meta_eo_all <- extract_validation_meta(validation_fit_eo, "EO")
      validation_meta_lo_all <- extract_validation_meta(validation_fit_lo, "LO")
      local_validation_all <- bind_rows(
        extract_local_validation(validation_fit_eo, "EO"),
        extract_local_validation(validation_fit_lo, "LO")
      ) %>%
        semi_join(tier_species, by = "feature") %>%
        mutate(
          design = design_name,
          validation_fold = validation_fold,
          .before = 1
        )
      effect_results_base <- tier_species %>%
        left_join(validation_meta_eo_all, by = "feature") %>%
        left_join(validation_meta_lo_all, by = "feature") %>%
        rename(beta_EO = coef_EO, beta_LO = coef_LO) %>%
        mutate(
          paired_estimate_available = is.finite(beta_EO) &
            is.finite(beta_LO),
          low_contributing_cohort_flag = paired_estimate_available & (k_EO < 3 | k_LO <
            3),
          delta_abs = if_else(paired_estimate_available, abs(beta_EO) - abs(beta_LO), NA_real_),
          delta_aligned = if_else(
            paired_estimate_available,
            discovery_direction * (beta_EO - beta_LO),
            NA_real_
          )
        )
      tier_effect_results <- bind_rows(
        effect_results_base %>%
          filter(tier_strict_core) %>%
          mutate(selection_tier = "strict_core"),
        effect_results_base %>%
          filter(tier_one_side_significant) %>%
          mutate(selection_tier = "one_side_significant"),
        effect_results_base %>%
          filter(tier_all_directional_qc_pass) %>%
          mutate(selection_tier = "all_directional_qc_pass")
      ) %>%
        relocate(selection_tier, .after = validation_fold)
      effect_results <- tier_effect_results %>%
        filter(selection_tier == "strict_core") %>%
        select(-selection_tier, -starts_with("tier_"))
      validation_meta_eo <- validation_meta_eo_all %>%
        semi_join(selected_species, by = "feature")
      validation_meta_lo <- validation_meta_lo_all %>%
        semi_join(selected_species, by = "feature")
      local_validation <- local_validation_all %>%
        semi_join(selected_species, by = "feature")
    }

    write_tsv(
      validation_meta_eo,
      file.path(results_dir, paste0(round_label, "_validation_EO_meta.tsv")),
      na = "NA"
    )
    write_tsv(
      validation_meta_lo,
      file.path(results_dir, paste0(round_label, "_validation_LO_meta.tsv")),
      na = "NA"
    )
    write_tsv(
      local_validation,
      file.path(results_dir, paste0(round_label, "_validation_local_effects.tsv")),
      na = "NA"
    )
    write_tsv(
      effect_results,
      file.path(results_dir, paste0(round_label, "_validation_effect_difference.tsv")),
      na = "NA"
    )
    write_tsv(
      validation_meta_eo_all,
      file.path(results_dir, paste0(round_label, "_three_tier_validation_EO_meta.tsv")),
      na = "NA"
    )
    write_tsv(
      validation_meta_lo_all,
      file.path(results_dir, paste0(round_label, "_three_tier_validation_LO_meta.tsv")),
      na = "NA"
    )
    write_tsv(
      local_validation_all,
      file.path(results_dir, paste0(round_label, "_three_tier_validation_local_effects.tsv")),
      na = "NA"
    )
    write_tsv(
      tier_effect_results,
      file.path(results_dir, paste0(round_label, "_three_tier_validation_effects.tsv")),
      na = "NA"
    )
    all_effect_results[[round_label]] <- effect_results
    all_tier_effect_results[[round_label]] <- tier_effect_results
    all_round_metadata[[round_label]] <- tibble(
      design = design_name,
      validation_fold = validation_fold,
      discovery_cohorts = paste(sort(discovery_cohorts), collapse = ";"),
      validation_cohorts = paste(sort(validation_cohorts), collapse = ";"),
      discovery_eo_cohorts = paste(sort(discovery_eo_cohorts), collapse = ";"),
      discovery_lo_cohorts = paste(sort(discovery_lo_cohorts), collapse = ";"),
      validation_eo_cohorts = paste(sort(validation_eo_cohorts), collapse = ";"),
      validation_lo_cohorts = paste(sort(validation_lo_cohorts), collapse = ";"),
      discovery_eo_n = nrow(discovery_eo),
      discovery_lo_n = nrow(discovery_lo),
      validation_eo_n = nrow(validation_eo),
      validation_lo_n = nrow(validation_lo),
      discovery_consensus_n = nrow(consensus_matrix),
      discovery_selected_n = nrow(selected_species),
      discovery_one_side_selected_n = sum(tier_species$tier_one_side_significant),
      discovery_all_directional_selected_n = sum(tier_species$tier_all_directional_qc_pass),
      validation_estimable_n = sum(effect_results$paired_estimate_available),
      validation_one_side_estimable_n = sum(tier_effect_results$selection_tier == "one_side_significant" &
        tier_effect_results$paired_estimate_available),
      validation_all_directional_estimable_n = sum(tier_effect_results$selection_tier ==
        "all_directional_qc_pass" & tier_effect_results$paired_estimate_available)
    )
    all_fit_objects[[round_label]] <- list(
      discovery_EO = discovery_fit_eo,
      discovery_LO = discovery_fit_lo,
      validation_EO = validation_fit_eo,
      validation_LO = validation_fit_lo,
      selected_species = selected_species,
      tier_species = tier_species
    )
  }
}

effect_results_all <- bind_rows(all_effect_results)
tier_effect_results_all <- bind_rows(all_tier_effect_results)
round_metadata <- bind_rows(all_round_metadata)

write_tsv(
  effect_results_all,
  file.path(results_dir, "crossfit_all_validation_effects.tsv"),
  na = "NA"
)

write_tsv(round_metadata, file.path(results_dir, "crossfit_round_metadata.tsv"), na = "NA")

write_tsv(
  tier_effect_results_all,
  file.path(results_dir, "crossfit_three_tier_validation_effects.tsv"),
  na = "NA"
)

set.seed(random_seed)

fold_summaries <- list()

for (design_name in analysis_designs) {
  for (validation_fold in 1:3) {
    fold_data <- effect_results_all %>%
      filter(
        design == design_name,
        validation_fold == !!validation_fold,
        paired_estimate_available
      )

    for (metric_name in c("delta_abs", "delta_aligned")) {
      metric_values <- fold_data[[metric_name]]
      metric_test <- summarise_test_statistics(metric_values)
      metric_ci <- bootstrap_median_ci(metric_values, bootstrap_repetitions)
      fold_summaries[[paste(design_name, validation_fold, metric_name)]] <- metric_test %>%
        mutate(
          design = design_name,
          scope = paste0("validation_fold", validation_fold),
          metric = metric_name,
          ci_lower = metric_ci[1],
          ci_upper = metric_ci[2],
          ci_method = paste0(
            "taxon bootstrap percentile, B=",
            bootstrap_repetitions,
            "; conditional on discovery selection"
          ),
          .before = 1
        )
    }
  }
}

fold_summary <- bind_rows(fold_summaries)

species_balanced <- effect_results_all %>%
  group_by(design, feature, Species) %>%
  summarise(
    n_selected_folds = n(),
    n_estimable_folds = sum(paired_estimate_available),
    selected_folds = paste(sort(unique(validation_fold)), collapse = ";"),
    estimable_folds = paste(
      sort(unique(validation_fold[paired_estimate_available])),
      collapse = ";"
    ),
    crossfit_delta_abs = if (any(paired_estimate_available)) {
      mean(delta_abs[paired_estimate_available])
    } else {
      NA_real_
    },
    crossfit_delta_aligned = if (any(paired_estimate_available)) {
      mean(delta_aligned[paired_estimate_available])
    } else {
      NA_real_
    },
    .groups = "drop"
  )

write_tsv(
  species_balanced,
  file.path(results_dir, "crossfit_species_balanced_effects.tsv"),
  na = "NA"
)

overall_summaries <- list()

for (design_name in analysis_designs) {
  design_species <- species_balanced %>%
    filter(design == design_name)

  for (metric_name in c("crossfit_delta_abs", "crossfit_delta_aligned")) {
    metric_values <- design_species[[metric_name]]
    metric_test <- summarise_test_statistics(metric_values)
    metric_ci <- bootstrap_median_ci(metric_values, bootstrap_repetitions)
    overall_summaries[[paste(design_name, metric_name)]] <- metric_test %>%
      mutate(
        design = design_name,
        scope = "overall_species_balanced",
        metric = str_remove(metric_name, "^crossfit_"),
        ci_lower = metric_ci[1],
        ci_upper = metric_ci[2],
        ci_method = paste0(
          "unique-taxon bootstrap percentile, B=",
          bootstrap_repetitions,
          "; repeated folds averaged before resampling"
        ),
        .before = 1
      )
  }
}

overall_summary <- bind_rows(overall_summaries)
crossfit_summary <- bind_rows(fold_summary, overall_summary)

fold_balanced_summary <- fold_summary %>%
  group_by(design, metric) %>%
  summarise(
    n_folds = n(),
    n_positive_fold_medians = sum(
      median >
        0,
      na.rm = TRUE
    ),
    mean_of_fold_medians = mean(median, na.rm = TRUE),
    median_of_fold_medians = median(median, na.rm = TRUE),
    minimum_fold_median = min(median, na.rm = TRUE),
    maximum_fold_median = max(median, na.rm = TRUE),
    .groups = "drop"
  )

k3_summary <- effect_results_all %>%
  filter(paired_estimate_available, k_EO >= 3, k_LO >= 3) %>%
  group_by(design, validation_fold) %>%
  group_modify(function(data, keys) {
    bind_rows(
      summarise_test_statistics(data$delta_abs) %>%
        mutate(metric = "delta_abs"),
      summarise_test_statistics(data$delta_aligned) %>%
        mutate(metric = "delta_aligned")
    )
  }) %>%
  ungroup() %>%
  mutate(scope = paste0("validation_fold", validation_fold, "_k_EO_LO_ge_3")) %>%
  select(-validation_fold)

write_tsv(crossfit_summary, file.path(results_dir, "crossfit_summary.tsv"), na = "NA")

write_tsv(
  fold_balanced_summary,
  file.path(results_dir, "crossfit_fold_balanced_summary.tsv"),
  na = "NA"
)

write_tsv(k3_summary, file.path(results_dir, "crossfit_k_ge_3_sensitivity.tsv"), na = "NA")

selection_tiers <- c("strict_core", "one_side_significant", "all_directional_qc_pass")

set.seed(random_seed + 1L)

three_tier_fold_summaries <- list()

for (design_name in analysis_designs) {
  for (tier_name in selection_tiers) {
    for (validation_fold in 1:3) {
      fold_data <- tier_effect_results_all %>%
        filter(
          design == design_name,
          selection_tier == tier_name,
          validation_fold == !!validation_fold,
          paired_estimate_available
        )

      for (metric_name in c("delta_abs", "delta_aligned")) {
        metric_values <- fold_data[[metric_name]]
        metric_test <- summarise_test_statistics(metric_values)
        metric_ci <- bootstrap_median_ci(metric_values, bootstrap_repetitions)
        result_name <- paste(design_name, tier_name, validation_fold, metric_name)
        three_tier_fold_summaries[[result_name]] <- metric_test %>%
          mutate(
            design = design_name,
            selection_tier = tier_name,
            scope = paste0("validation_fold", validation_fold),
            metric = metric_name,
            ci_lower = metric_ci[1],
            ci_upper = metric_ci[2],
            ci_method = paste0(
              "taxon bootstrap percentile, B=",
              bootstrap_repetitions,
              "; conditional on discovery tier"
            ),
            .before = 1
          )
      }
    }
  }
}

three_tier_fold_summary <- bind_rows(three_tier_fold_summaries)

three_tier_species_balanced <- tier_effect_results_all %>%
  group_by(design, selection_tier, feature, Species) %>%
  summarise(
    n_selected_folds = n(),
    n_estimable_folds = sum(paired_estimate_available),
    selected_folds = paste(sort(unique(validation_fold)), collapse = ";"),
    estimable_folds = paste(
      sort(unique(validation_fold[paired_estimate_available])),
      collapse = ";"
    ),
    crossfit_delta_abs = if (any(paired_estimate_available)) {
      mean(delta_abs[paired_estimate_available])
    } else {
      NA_real_
    },
    crossfit_delta_aligned = if (any(paired_estimate_available)) {
      mean(delta_aligned[paired_estimate_available])
    } else {
      NA_real_
    },
    .groups = "drop"
  )

three_tier_overall_summaries <- list()

for (design_name in analysis_designs) {
  for (tier_name in selection_tiers) {
    tier_species <- three_tier_species_balanced %>%
      filter(design == design_name, selection_tier == tier_name)

    for (metric_name in c("crossfit_delta_abs", "crossfit_delta_aligned")) {
      metric_values <- tier_species[[metric_name]]
      metric_test <- summarise_test_statistics(metric_values)
      metric_ci <- bootstrap_median_ci(metric_values, bootstrap_repetitions)
      result_name <- paste(design_name, tier_name, metric_name)
      three_tier_overall_summaries[[result_name]] <- metric_test %>%
        mutate(
          design = design_name,
          selection_tier = tier_name,
          scope = "overall_species_balanced",
          metric = str_remove(metric_name, "^crossfit_"),
          ci_lower = metric_ci[1],
          ci_upper = metric_ci[2],
          ci_method = paste0(
            "unique-taxon bootstrap percentile, B=",
            bootstrap_repetitions,
            "; repeated folds averaged before resampling"
          ),
          .before = 1
        )
    }
  }
}

three_tier_overall_summary <- bind_rows(three_tier_overall_summaries)
three_tier_summary <- bind_rows(three_tier_fold_summary, three_tier_overall_summary)

strict_summary_for_three_tier <- crossfit_summary %>%
  mutate(selection_tier = "strict_core") %>%
  select(design, selection_tier, everything())

three_tier_summary <- three_tier_summary %>%
  filter(selection_tier != "strict_core") %>%
  bind_rows(strict_summary_for_three_tier) %>%
  arrange(design, selection_tier, scope, metric)

three_tier_fold_balanced_summary <- three_tier_fold_summary %>%
  group_by(design, selection_tier, metric) %>%
  summarise(
    n_folds = n(),
    n_positive_fold_medians = sum(median > 0, na.rm = TRUE),
    mean_of_fold_medians = mean(median, na.rm = TRUE),
    median_of_fold_medians = median(median, na.rm = TRUE),
    minimum_fold_median = min(median, na.rm = TRUE),
    maximum_fold_median = max(median, na.rm = TRUE),
    .groups = "drop"
  )

three_tier_k3_summary <- tier_effect_results_all %>%
  filter(paired_estimate_available, k_EO >= 3, k_LO >= 3) %>%
  group_by(design, selection_tier, validation_fold) %>%
  group_modify(function(data, keys) {
    bind_rows(
      summarise_test_statistics(data$delta_abs) %>%
        mutate(metric = "delta_abs"),
      summarise_test_statistics(data$delta_aligned) %>%
        mutate(metric = "delta_aligned")
    )
  }) %>%
  ungroup() %>%
  mutate(scope = paste0("validation_fold", validation_fold, "_k_EO_LO_ge_3")) %>%
  select(-validation_fold)

write_tsv(
  three_tier_summary,
  file.path(results_dir, "crossfit_three_tier_summary.tsv"),
  na = "NA"
)

write_tsv(
  three_tier_species_balanced,
  file.path(results_dir, "crossfit_three_tier_species_balanced_effects.tsv"),
  na = "NA"
)

write_tsv(
  three_tier_fold_balanced_summary,
  file.path(results_dir, "crossfit_three_tier_fold_balanced_summary.tsv"),
  na = "NA"
)

write_tsv(
  three_tier_k3_summary,
  file.path(results_dir, "crossfit_three_tier_k_ge_3_sensitivity.tsv"),
  na = "NA"
)

crossfit_models <- list(
  fold_assignment = fold_assignment,
  round_metadata = round_metadata,
  fits = all_fit_objects,
  validation_effects = effect_results_all,
  three_tier_validation_effects = tier_effect_results_all,
  three_tier_summary = three_tier_summary,
  species_balanced = species_balanced,
  summary = crossfit_summary
)

saveRDS(crossfit_models, file.path(results_dir, "cross_fitting_analysis_objects.rds"))
