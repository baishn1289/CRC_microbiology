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

table_dir <- file.path(sensitivity_root, "taxon")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
out_table_dir <- table_dir

# Genus nominations use the same thresholds with a genus-specific row filter.
run_consensus_EO_LO_combine_rank <- function(
  meta_list,
  meta_EO,
  meta_LO,
  rank_pattern,
  min_prevalence = 0.05,
  min_cohorts = 3,
  min_step1_count = 3,
  abund_threshold = 0.001,
  renormalize = TRUE
) {
  rank_meta_list <- lapply(meta_list, extract_rank_table, pattern = rank_pattern)
  names(rank_meta_list) <- names(meta_list)
  get_passed_taxa <- function(cohort_name, full_df, target_samples) {
    sample_cols <- full_df %>%
      column_to_rownames("clade_name") %>%
      select(any_of(target_samples))
    n_samples <- ncol(sample_cols)
    if (n_samples == 0) {
      message(sprintf("Cohort %s: no target samples.", cohort_name))
      return(character(0))
    }

    n_detected <- rowSums(sample_cols > abund_threshold)
    passed_taxa <- rownames(sample_cols)[n_detected >= min_step1_count &
      n_detected / n_samples >= min_prevalence]
    message(sprintf(
      "Cohort %s (N=%d): %d genus-level clades passed step 1.",
      cohort_name,
      n_samples,
      length(passed_taxa)
    ))
    passed_taxa
  }

  cohort_votes_EO <- lapply(
    names(rank_meta_list),
    function(cn) {
      get_passed_taxa(cn, rank_meta_list[[cn]], meta_EO$Sample_ID)
    }
  )
  cohort_votes_LO <- lapply(
    names(rank_meta_list),
    function(cn) {
      get_passed_taxa(cn, rank_meta_list[[cn]], meta_LO$Sample_ID)
    }
  )
  keep_from_votes <- function(votes) {
    votes <- unlist(votes, use.names = FALSE)
    if (length(votes) == 0) {
      return(character(0))
    }

    tab <- table(votes)
    names(tab)[tab >= min_cohorts]
  }

  keep_taxa_EO <- keep_from_votes(cohort_votes_EO)
  keep_taxa_LO <- keep_from_votes(cohort_votes_LO)
  keep_taxa_final <- unique(c(keep_taxa_EO, keep_taxa_LO))
  if (length(keep_taxa_final) == 0) {
    stop("No genus-level clades passed the consensus filter.")
  }

  message(sprintf(
    "Consensus result: %d genus-level clades retained (EO=%d, LO=%d).",
    length(keep_taxa_final),
    length(keep_taxa_EO),
    length(keep_taxa_LO)
  ))
  merged_df <- reduce(rank_meta_list, full_join, by = "clade_name") %>%
    filter(clade_name %in% keep_taxa_final)
  taxa_names <- merged_df$clade_name
  mat <- merged_df %>%
    select(-clade_name) %>%
    as.matrix()
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  if (renormalize) {
    col_sums <- colSums(mat)
    col_sums[col_sums == 0] <- 1
    mat <- sweep(mat, 2, col_sums, "/") * 100
  }

  merged_df <- as.data.frame(mat, check.names = FALSE) %>%
    mutate(clade_name = taxa_names, .before = 1)

  list(
    merged = merged_df,
    EO = merged_df %>%
      select(clade_name, any_of(meta_EO$Sample_ID)),
    LO = merged_df %>%
      select(clade_name, any_of(meta_LO$Sample_ID)),
    keep_taxa_EO = keep_taxa_EO,
    keep_taxa_LO = keep_taxa_LO
  )
}

clean_taxon_label <- function(feature, rank) {
  if (rank == "species") {
    str_replace(feature, ".*s__", "") %>%
      str_replace_all("_", " ")
  } else if (rank == "genus") {
    str_replace(feature, ".*g__", "") %>%
      str_replace_all("_", " ")
  } else {
    feature
  }
}

safe_wilcox <- function(x, y, alternative = "greater") {
  if (length(x) < 2 || length(y) < 2 || all((x - y) == 0)) {
    return(NA_real_)
  }

  tryCatch(
    as.numeric(wilcox.test(
      x,
      y,
      paired = TRUE,
      alternative = alternative,
      exact = TRUE,
      correct = FALSE
    )$p.value),
    error = function(e) NA_real_
  )
}

safe_sign_test <- function(delta) {
  n_positive <- sum(delta > 0, na.rm = TRUE)
  n_nonzero <- sum(delta != 0 & !is.na(delta))
  if (n_nonzero == 0)
    return(NA_real_)
  binom.test(
    n_positive,
    n_nonzero,
    p = 0.5,
    alternative = "greater"
  )$p.value
}

make_shared_table <- function(results_eo, results_lo, rank) {
  rank_name <- rank
  pass_eo <- results_eo %>%
    filter(pass_qc) %>%
    transmute(
      feature,
      coef_EO = coef,
      pval_EO = pval,
      qval_EO = qval.fdr,
      I2_EO = I2,
      max_weight_EO = max_weight,
      consistency_EO = consistency_score,
      significant_EO = significant,
      qc_status_EO = qc_status
    )
  pass_lo <- results_lo %>%
    filter(pass_qc) %>%
    transmute(
      feature,
      coef_LO = coef,
      pval_LO = pval,
      qval_LO = qval.fdr,
      I2_LO = I2,
      max_weight_LO = max_weight,
      consistency_LO = consistency_score,
      significant_LO = significant,
      qc_status_LO = qc_status
    )
  inner_join(pass_eo, pass_lo, by = "feature") %>%
    mutate(
      rank = rank_name,
      taxon = clean_taxon_label(feature, rank_name),
      abs_coef_EO = abs(coef_EO),
      abs_coef_LO = abs(coef_LO),
      delta_abs = abs_coef_EO - abs_coef_LO,
      same_direction = sign(coef_EO) == sign(coef_LO),
      both_significant = significant_EO & significant_LO
    ) %>%
    select(
      rank,
      feature,
      taxon,
      coef_EO,
      coef_LO,
      abs_coef_EO,
      abs_coef_LO,
      delta_abs,
      same_direction,
      both_significant,
      everything()
    )
}

build_main_species_core_table <- function(results_eo, results_lo) {
  main_sig_eo <- results_eo %>%
    filter(significant) %>%
    mutate(Species = str_remove(str_extract(feature, "s__[A-Za-z0-9_]+"), "^s__"))
  main_sig_lo <- results_lo %>%
    filter(significant) %>%
    mutate(Species = str_remove(str_extract(feature, "s__[A-Za-z0-9_]+"), "^s__"))
  species_core <- full_join(
    main_sig_eo %>%
      transmute(
        feature_EO = feature,
        Species,
        coef_EO = coef,
        pval_EO = pval,
        qval_EO = qval.fdr,
        I2_EO = I2,
        max_weight_EO = max_weight,
        consistency_EO = consistency_score,
        qc_status_EO = qc_status
      ),
    main_sig_lo %>%
      transmute(
        feature_LO = feature,
        Species,
        coef_LO = coef,
        pval_LO = pval,
        qval_LO = qval.fdr,
        I2_LO = I2,
        max_weight_LO = max_weight,
        consistency_LO = consistency_score,
        qc_status_LO = qc_status
      ),
    by = "Species"
  ) %>%
    filter(!is.na(coef_EO), !is.na(coef_LO), sign(coef_EO) == sign(coef_LO)) %>%
    mutate(
      rank = "species",
      feature = feature_EO,
      taxon = str_replace_all(Species, "_", " "),
      abs_coef_EO = abs(coef_EO),
      abs_coef_LO = abs(coef_LO),
      delta_abs = abs_coef_EO - abs_coef_LO,
      same_direction = TRUE,
      both_significant = TRUE
    ) %>%
    arrange(match(feature, fit_meta.crc_EO$meta_fits$feature)) %>%
    select(
      rank,
      feature,
      taxon,
      coef_EO,
      coef_LO,
      abs_coef_EO,
      abs_coef_LO,
      delta_abs,
      same_direction,
      both_significant,
      everything()
    )
  if (nrow(species_core) != 27) {
    stop("Expected 27 main-analysis shared core species, got ", nrow(species_core))
  }

  if (!identical(species_core$feature_EO, species_core$feature_LO)) {
    stop("EO and LO species feature identifiers differ in the main shared core set.")
  }

  species_core
}

add_effect_direction <- function(shared_df) {
  shared_df %>%
    mutate(effect_direction = case_when(
      coef_EO > 0 & coef_LO > 0 ~ "positive",
      coef_EO < 0 & coef_LO <
        0 ~ "negative",
      TRUE ~ "mixed_or_zero"
    ))
}

summarise_shared <- function(shared_df, set_name, sign_group = "all_signs") {
  if (nrow(shared_df) == 0) {
    return(tibble())
  }

  if (!("effect_direction" %in% colnames(shared_df))) {
    shared_df <- add_effect_direction(shared_df)
  }

  shared_df %>%
    summarise(
      rank = first(rank),
      set = set_name,
      sign_group = sign_group,
      n_taxa = n(),
      n_positive = sum(
        effect_direction ==
          "positive",
        na.rm = TRUE
      ),
      n_negative = sum(effect_direction == "negative", na.rm = TRUE),
      median_abs_EO = median(abs_coef_EO, na.rm = TRUE),
      median_abs_LO = median(abs_coef_LO, na.rm = TRUE),
      median_delta_abs = median(delta_abs, na.rm = TRUE),
      mean_delta_abs = mean(delta_abs, na.rm = TRUE),
      n_EO_gt_LO = sum(delta_abs > 0, na.rm = TRUE),
      prop_EO_gt_LO = mean(delta_abs > 0, na.rm = TRUE),
      wilcox_p_greater = safe_wilcox(abs_coef_EO, abs_coef_LO, "greater"),
      wilcox_p_two_sided = safe_wilcox(abs_coef_EO, abs_coef_LO, "two.sided"),
      sign_test_p_greater = safe_sign_test(delta_abs),
      n_same_direction = sum(same_direction, na.rm = TRUE),
      n_both_significant = sum(both_significant, na.rm = TRUE),
      .groups = "drop"
    )
}

summarise_with_sign_splits <- function(shared_df, set_name) {
  shared_df <- add_effect_direction(shared_df)
  bind_rows(
    summarise_shared(shared_df, set_name, "all_signs"),
    summarise_shared(
      shared_df %>%
        filter(effect_direction ==
          "positive"),
      set_name,
      "positive"
    ),
    summarise_shared(
      shared_df %>%
        filter(effect_direction == "negative"),
      set_name,
      "negative"
    )
  )
}

make_one_sided_significant_table <- function(shared_df) {
  shared_df %>%
    mutate(
      significant_EO = coalesce(significant_EO, FALSE),
      significant_LO = coalesce(significant_LO, FALSE),
      significant_side = case_when(
        significant_EO & !significant_LO ~ "EOCRC_only",
        !significant_EO &
          significant_LO ~ "LOCRC_only",
        TRUE ~ "not_one_sided"
      )
    ) %>%
    filter(
      same_direction,
      significant_side %in% c("EOCRC_only", "LOCRC_only")
    )
}

make_at_least_one_side_significant_table <- function(shared_df) {
  shared_df %>%
    mutate(
      significant_EO = coalesce(significant_EO, FALSE),
      significant_LO = coalesce(significant_LO, FALSE),
      significant_side = case_when(
        significant_EO & significant_LO ~ "both_EOCRC_LOCRC",
        significant_EO &
          !significant_LO ~ "EOCRC_only",
        !significant_EO & significant_LO ~ "LOCRC_only",
        TRUE ~ "not_significant"
      )
    ) %>%
    filter(same_direction, significant_EO | significant_LO)
}

md <- sensitivity_metadata(metadata_allcohorts)
world_meta_EO <- md %>%
  filter(Age_class == "EO")
world_meta_LO <- md %>%
  filter(Age_class == "LO")
genus_abundance <- run_consensus_EO_LO_combine_rank(
  meta_list_metadata_intersect_samples,
  world_meta_EO,
  world_meta_LO,
  rank_pattern = "\\|g__[^|]+$",
  min_cohorts = 3,
  min_prevalence = 0.05,
  min_step1_count = 3,
  abund_threshold = 0.001,
  renormalize = TRUE
)
fit_meta.genus_EO <- sensitivity_fit(
  world_meta_EO,
  genus_abundance$EO,
  file.path(sensitivity_root, "models/genus_EO")
)
fit_meta.genus_LO <- sensitivity_fit(
  world_meta_LO,
  genus_abundance$LO,
  file.path(sensitivity_root, "models/genus_LO")
)
species_results_EO <- sensitivity_qc(fit_meta.crc_EO) %>%
  mutate(rank = "species")

species_results_LO <- sensitivity_qc(fit_meta.crc_LO) %>%
  mutate(rank = "species")

genus_results_EO_main_qc <- sensitivity_qc(fit_meta.genus_EO) %>%
  mutate(rank = "genus")

genus_results_LO_main_qc <- sensitivity_qc(fit_meta.genus_LO) %>%
  mutate(rank = "genus")

species_shared_all <- make_shared_table(species_results_EO, species_results_LO, "species")

species_shared_dir <- species_shared_all %>%
  filter(same_direction)

species_shared_core_dir <- build_main_species_core_table(species_results_EO, species_results_LO)
genus_shared_all <- make_shared_table(genus_results_EO_main_qc, genus_results_LO_main_qc, "genus")

genus_shared_dir <- genus_shared_all %>%
  filter(same_direction)

genus_shared_core_dir <- genus_shared_all %>%
  filter(both_significant, same_direction)

species_one_sided_sig_dir <- make_one_sided_significant_table(species_shared_all)
genus_one_sided_sig_dir <- make_one_sided_significant_table(genus_shared_all)
species_at_least_one_side_sig_dir <- make_at_least_one_side_significant_table(species_shared_all)
genus_at_least_one_side_sig_dir <- make_at_least_one_side_significant_table(genus_shared_all)
shared_all <- bind_rows(species_shared_all, genus_shared_all)

summary_table <- bind_rows(
  summarise_with_sign_splits(species_shared_dir, "direction_consistent_QC_pass"),
  summarise_with_sign_splits(
    species_shared_core_dir,
    "shared_core_significant_direction_consistent"
  ),
  summarise_with_sign_splits(
    species_at_least_one_side_sig_dir,
    "at_least_one_side_significant_direction_consistent"
  ),
  summarise_with_sign_splits(
    species_one_sided_sig_dir,
    "one_sided_significant_direction_consistent"
  ),
  summarise_with_sign_splits(
    species_one_sided_sig_dir %>%
      filter(significant_side ==
        "EOCRC_only"),
    "EOCRC_only_significant_direction_consistent"
  ),
  summarise_with_sign_splits(
    species_one_sided_sig_dir %>%
      filter(significant_side == "LOCRC_only"),
    "LOCRC_only_significant_direction_consistent"
  ),
  summarise_with_sign_splits(genus_shared_dir, "direction_consistent_QC_pass"),
  summarise_with_sign_splits(genus_shared_core_dir, "shared_core_significant_direction_consistent"),
  summarise_with_sign_splits(
    genus_at_least_one_side_sig_dir,
    "at_least_one_side_significant_direction_consistent"
  ),
  summarise_with_sign_splits(genus_one_sided_sig_dir, "one_sided_significant_direction_consistent"),
  summarise_with_sign_splits(
    genus_one_sided_sig_dir %>%
      filter(significant_side == "EOCRC_only"),
    "EOCRC_only_significant_direction_consistent"
  ),
  summarise_with_sign_splits(
    genus_one_sided_sig_dir %>%
      filter(significant_side == "LOCRC_only"),
    "LOCRC_only_significant_direction_consistent"
  )
) %>%
  mutate(across(where(is.numeric), ~ signif(.x, 6)))

write.csv(
  species_shared_all,
  file.path(out_table_dir, "5.3.all_qc_species_effect_size_sensitivity.csv"),
  row.names = FALSE
)

write.csv(
  genus_shared_all,
  file.path(out_table_dir, "5.3.all_qc_genus_effect_size_sensitivity.csv"),
  row.names = FALSE
)

write.csv(
  species_shared_core_dir,
  file.path(out_table_dir, "5.3.main_shared_core_species_effect_size_sensitivity.csv"),
  row.names = FALSE
)

write.csv(
  genus_shared_core_dir,
  file.path(out_table_dir, "5.3.main_shared_core_genus_effect_size_sensitivity.csv"),
  row.names = FALSE
)

write.csv(
  species_one_sided_sig_dir,
  file.path(out_table_dir, "5.3.one_sided_significant_species_effect_size_sensitivity.csv"),
  row.names = FALSE
)

write.csv(
  genus_one_sided_sig_dir,
  file.path(out_table_dir, "5.3.one_sided_significant_genus_effect_size_sensitivity.csv"),
  row.names = FALSE
)

write.csv(
  species_at_least_one_side_sig_dir,
  file.path(
    out_table_dir,
    "5.3.at_least_one_side_significant_species_effect_size_sensitivity.csv"
  ),
  row.names = FALSE
)

write.csv(
  genus_at_least_one_side_sig_dir,
  file.path(out_table_dir, "5.3.at_least_one_side_significant_genus_effect_size_sensitivity.csv"),
  row.names = FALSE
)

write.csv(
  summary_table,
  file.path(out_table_dir, "5.3.all_qc_taxa_effect_size_sensitivity_summary.csv"),
  row.names = FALSE
)

# Exact paired tests retain the original one-sided direction.
direction_summary <- species_shared_core_dir %>%
  mutate(direction = if_else(coef_EO > 0, "positive", "negative")) %>%
  group_by(direction) %>%
  summarise(
    rank = "species",
    n_taxa = n(),
    wilcox_p_greater_abs_EO_vs_LO = wilcox.test(
      abs_coef_EO,
      abs_coef_LO,
      paired = TRUE,
      alternative = "greater",
      exact = TRUE,
      correct = FALSE
    )$p.value,
    .groups = "drop"
  )
write_csv(direction_summary, file.path(table_dir, "5.3.species_effect_size_by_direction.csv"))
