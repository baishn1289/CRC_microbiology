
suppressPackageStartupMessages({
  library(MMUPHin)
  library(metafor)
  library(tidyverse)
  library(data.table)
  library(furrr)
  library(future)
  library(glue)
  library(readr)
})


disk_symbol <- "J:/"
base_symbol <- "/home/global_EOCRC_microbiology/"
project_dir <- base_symbol
downstream_dir <- file.path(project_dir, "downstream_analysis")
generated_dir <- file.path(downstream_dir, "generated_data")
results_dir <- file.path(downstream_dir, "output", "table", "permutation_test")
picture_dir <- file.path(downstream_dir, "output", "picture", "permutation_test")

dir.create(generated_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(picture_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  hit <- match(flag, args)
  if (is.na(hit) || hit == length(args)) {
    return(default)
  }
  args[[hit + 1]]
}

n_permutations <- as.integer(get_arg("--n-permutations", "10000"))
n_workers <- as.integer(get_arg("--workers", "1"))
seed_base <- as.integer(get_arg("--seed", "20260620"))
keep_mmuphin_output <- isTRUE(as.logical(get_arg("--keep-mmuphin-output", "FALSE")))
coef_tolerance <- as.numeric(get_arg("--coef-tolerance", "1e-8"))
sample_mode <- get_arg("--sample-mode", "main-compatible")
allow_failures <- isTRUE(as.logical(get_arg("--allow-failures", "FALSE")))

if (is.na(n_permutations) || n_permutations < 0) {
  stop("--n-permutations must be a non-negative integer")
}
if (is.na(n_workers) || n_workers < 1) {
  stop("--workers must be a positive integer")
}
if (!sample_mode %in% c("main-compatible", "complete-case")) {
  stop("--sample-mode must be either 'main-compatible' or 'complete-case'")
}

analysis_tag <- sprintf("fixed27_B%d_seed%d", n_permutations, seed_base)
output_prefix <- file.path(results_dir, analysis_tag)
plot_prefix <- file.path(picture_dir, analysis_tag)
mmuphin_tmp_dir <- file.path(results_dir, "mmuphin_tmp", analysis_tag)
dir.create(mmuphin_tmp_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading saved analysis objects...")

process_meta_with_rescue <- function(meta_res_obj,
                                     group_name = "CRC",
                                     i2_tier1 = 50,
                                     i2_tier2 = 75,
                                     weight_limit = 40,
                                     base_consistency = 0.6,
                                     rescue_consistency = 0.75) {
  df_meta <- meta_res_obj$meta_fits %>%
    filter(!is.na(coef)) %>%
    rownames_to_column("Pathway")

  df_individual <- bind_rows(meta_res_obj$maaslin_fits) %>%
    select(feature, coef, batch) %>%
    rename(Pathway = feature, cohort_coef = coef, cohort_id = batch) %>%
    filter(!is.na(cohort_coef))

  consistency_stats <- df_meta %>%
    select(Pathway, global_coef = coef) %>%
    left_join(df_individual, by = "Pathway") %>%
    group_by(Pathway) %>%
    filter(cohort_coef != 0) %>%
    summarise(
      n_total_cohorts = n_distinct(cohort_id),
      n_same_dir = sum(sign(cohort_coef) == sign(global_coef), na.rm = TRUE),
      consistency_score = n_same_dir / n_total_cohorts,
      n_opp_dir = n_total_cohorts - n_same_dir,
      pass_min_cohort = n_total_cohorts >= 3,
      .groups = "drop"
    )

  df_meta %>%
    left_join(consistency_stats, by = "Pathway") %>%
    rowwise() %>%
    mutate(
      max_weight = ifelse(
        all(is.na(c_across(starts_with("weight_")))),
        NA,
        max(c_across(starts_with("weight_")), na.rm = TRUE)
      ),
      pass_weight = !is.na(max_weight) & max_weight < weight_limit,
      is_tier1 = I2 < i2_tier1,
      is_in_rescue_zone = (I2 >= i2_tier1) & (I2 < i2_tier2),
      pass_base_consistency = !is.na(consistency_score) &
        (consistency_score >= base_consistency),
      pass_rescue_consistency = !is.na(consistency_score) &
        (consistency_score >= rescue_consistency),
      pass_qc = case_when(
        !pass_weight ~ FALSE,
        !pass_min_cohort ~ FALSE,
        is_tier1 & pass_base_consistency ~ TRUE,
        is_tier1 & !pass_base_consistency ~ FALSE,
        is_in_rescue_zone & pass_rescue_consistency ~ TRUE,
        TRUE ~ FALSE
      ),
      qc_status = case_when(
        !pass_weight ~ "Fail_High_Weight",
        !pass_min_cohort ~ "Fail_Too_Few_Cohorts",
        is_tier1 & pass_base_consistency ~ "Pass_Tier1_Robust",
        is_tier1 & !pass_base_consistency ~ "Fail_Tier1_Inconsistent_Direction",
        is_in_rescue_zone & pass_rescue_consistency ~ "Pass_Tier2_Rescued",
        is_in_rescue_zone & !pass_rescue_consistency ~ "Fail_Tier2_Inconsistent",
        I2 >= i2_tier2 ~ "Fail_High_Heterogeneity",
        TRUE ~ "Fail_Other"
      )
    ) %>%
    ungroup()
}

apply_step2_masking <- function(feature_abd, meta_data, batch_col = "Cohort") {
  min_model_count <- 5
  feature_abd_masked <- feature_abd
  cohorts <- unique(meta_data[[batch_col]])

  for (cohort in cohorts) {
    sample_ids <- rownames(meta_data)[meta_data[[batch_col]] == cohort]
    sample_ids <- intersect(sample_ids, colnames(feature_abd_masked))
    if (length(sample_ids) == 0) {
      next
    }

    sub_dat <- feature_abd_masked[, sample_ids, drop = FALSE]
    n_nonzero <- rowSums(sub_dat > 0)
    sparse_taxa <- names(n_nonzero)[n_nonzero > 0 & n_nonzero < min_model_count]

    if (length(sparse_taxa) > 0) {
      feature_abd_masked[sparse_taxa, sample_ids] <- 0
    }
  }

  feature_abd_masked
}

format_meta_for_model <- function(meta_df) {
  meta_model <- meta_df %>%
    filter(Group %in% c("CRC", "Control")) %>%
    mutate(
      Group = factor(as.character(Group), levels = c("Control", "CRC")),
      Cohort = as.factor(Cohort),
      Age = suppressWarnings(as.numeric(as.character(Age))),
      Sex = as.factor(Sex),
      BMI = suppressWarnings(as.numeric(as.character(BMI)))
    )

  if (anyDuplicated(meta_model$Sample_ID)) {
    stop("Duplicated Sample_ID values in model metadata")
  }

  if (sample_mode == "complete-case") {
    meta_model <- meta_model %>%
      filter(is.finite(Age), !is.na(Sex), is.finite(BMI))
  }

  meta_model %>% column_to_rownames("Sample_ID")
}

prepare_abundance_for_meta <- function(meta_df, abundance_matrix) {
  meta_model <- format_meta_for_model(meta_df)
  sample_ids <- rownames(meta_model)
  missing_samples <- setdiff(sample_ids, colnames(abundance_matrix))
  if (length(missing_samples) > 0) {
    stop("Samples missing from abundance matrix: ", paste(missing_samples, collapse = ", "))
  }

  feature_abd <- abundance_matrix[, sample_ids, drop = FALSE]
  feature_abd <- apply_step2_masking(
    feature_abd = feature_abd,
    meta_data = meta_model,
    batch_col = "Cohort"
  )

  list(meta = meta_model, abundance = as.data.frame(feature_abd))
}

compute_meta_coefs <- function(meta_df,
                               abundance_matrix,
                               age_group,
                               perm_id,
                               seed_value,
                               output_root) {
  prepared <- prepare_abundance_for_meta(meta_df, abundance_matrix)
  output_dir <- file.path(
    output_root,
    sprintf("%s_seed_%d_pid_%s", age_group, seed_value, Sys.getpid())
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  fit <- suppressMessages(
    MMUPHin::lm_meta(
      feature_abd = prepared$abundance,
      exposure = "Group",
      batch = "Cohort",
      covariates = c("Age", "Sex", "BMI"),
      data = prepared$meta,
      control = list(
        rma_method = "HS",
        transform = "LOG",
        normalization = "NONE",
        min_prevalence = 0,
        min_abundance = 0.001,
        forest_plot = NULL,
        verbose = FALSE,
        output = output_dir
      )
    )
  )

  if (!"meta_fits" %in% names(fit)) {
    stop("MMUPHin::lm_meta() output does not contain meta_fits")
  }

  fit$meta_fits %>%
    as_tibble() %>%
    select(feature, coef) %>%
    mutate(coef = as.numeric(coef))
}

calculate_statistic <- function(res_eo, res_lo, target_features) {
  joined <- inner_join(res_eo, res_lo, by = "feature", suffix = c("_EO", "_LO")) %>%
    arrange(match(feature, target_features)) %>%
    mutate(delta_abs = abs(coef_EO) - abs(coef_LO))

  if (nrow(joined) != length(target_features)) {
    stop("Expected ", length(target_features), " joined features, got ", nrow(joined))
  }
  if (!identical(joined$feature, target_features)) {
    stop("Joined features do not match target feature order")
  }
  if (!all(is.finite(joined$coef_EO)) || !all(is.finite(joined$coef_LO))) {
    stop("Non-finite EO or LO coefficients")
  }

  list(
    statistic = median(joined$delta_abs),
    joined = joined
  )
}

build_main_results <- function(fit_obj) {
  process_meta_with_rescue(
    meta_res_obj = fit_obj,
    i2_tier1 = 50,
    i2_tier2 = 80,
    weight_limit = 40,
    base_consistency = 0.6,
    rescue_consistency = 0.75
  ) %>%
    filter(str_starts(qc_status, "Pass")) %>%
    mutate(qval.fdr = p.adjust(pval, method = "fdr")) %>%
    filter(qval.fdr < 0.1) %>%
    mutate(Species = str_remove(str_extract(feature, "s__[A-Za-z0-9_]+"), "^s__"))
}

message("Reconstructing the 27 target species from main-analysis results...")
meta_results_EO <- build_main_results(fit_meta.crc_EO)
meta_results_LO <- build_main_results(fit_meta.crc_LO)

target_27_table <- full_join(
  meta_results_EO %>% select(feature_EO = feature, Species, main_coef_EO = coef),
  meta_results_LO %>% select(feature_LO = feature, Species, main_coef_LO = coef),
  by = "Species"
) %>%
  mutate(
    Status = case_when(
      is.na(main_coef_LO) ~ "EO Unique",
      is.na(main_coef_EO) ~ "LO Unique",
      sign(main_coef_EO) == sign(main_coef_LO) ~ "Shared (Consistent)",
      TRUE ~ "Shared (Inconsistent)"
    )
  ) %>%
  filter(Status == "Shared (Consistent)") %>%
  mutate(feature = feature_EO) %>%
  arrange(match(feature, fit_meta.crc_EO$meta_fits$feature))

if (nrow(target_27_table) != 27) {
  stop("Expected exactly 27 shared direction-consistent species, got ", nrow(target_27_table))
}
if (!identical(target_27_table$feature_EO, target_27_table$feature_LO)) {
  stop("EO and LO feature names differ for at least one target species")
}
if (anyDuplicated(target_27_table$feature) || anyDuplicated(target_27_table$Species)) {
  stop("Duplicated target feature or species names")
}

target_features <- target_27_table$feature
write_csv(target_27_table, paste0(output_prefix, "_target_27_species.csv"))

message("Reconstructing main-analysis cohort eligibility...")
cohort_counts <- metadata_allcohorts %>%
  filter(Group %in% c("CRC", "Control")) %>%
  count(Cohort, Age_class, Group) %>%
  pivot_wider(
    names_from = c(Age_class, Group),
    values_from = n,
    values_fill = 0,
    names_glue = "{Age_class}{Group}"
  )

eo_cohorts <- cohort_counts %>%
  filter(EOControl >= 3, EOCRC >= 3) %>%
  pull(Cohort) %>%
  as.character()

lo_cohorts <- cohort_counts %>%
  filter(LOControl >= 3, LOCRC >= 3) %>%
  pull(Cohort) %>%
  as.character()

world_meta_EO <- metadata_allcohorts %>%
  filter(Age_class == "EO", Group %in% c("CRC", "Control"), Cohort %in% eo_cohorts)

world_meta_LO <- metadata_allcohorts %>%
  filter(Age_class == "LO", Group %in% c("CRC", "Control"), Cohort %in% lo_cohorts)

message("Reconstructing the main-analysis consensus abundance matrix...")
abundance_world_EO_LO_combine <- run_consensus_EO_LO_combine(
  meta_list = meta_list_metadata_intersect_samples,
  meta_EO = world_meta_EO,
  meta_LO = world_meta_LO,
  min_cohorts = 3,
  min_prevalence = 0.05,
  abund_threshold = 0.001
)

merged_df <- abundance_world_EO_LO_combine$merged
if (!setequal(rownames(fit_meta.crc_EO$meta_fits), merged_df$clade_name)) {
  stop("Reconstructed abundance features do not match fit_meta.crc_EO$meta_fits")
}
if (!setequal(rownames(fit_meta.crc_LO$meta_fits), merged_df$clade_name)) {
  stop("Reconstructed abundance features do not match fit_meta.crc_LO$meta_fits")
}

target_abundance_matrix <- merged_df %>%
  filter(clade_name %in% target_features) %>%
  arrange(match(clade_name, target_features)) %>%
  column_to_rownames("clade_name") %>%
  as.matrix()

if (nrow(target_abundance_matrix) != 27 ||
    !identical(rownames(target_abundance_matrix), target_features)) {
  stop("Target abundance matrix is not an exact 27-feature match")
}

message("Fixing the analysis sample set before permutation...")
fixed_meta_all <- metadata_allcohorts %>%
  filter(
    Group %in% c("CRC", "Control"),
    Age_class %in% c("EO", "LO"),
    Cohort %in% union(eo_cohorts, lo_cohorts),
    Sample_ID %in% colnames(target_abundance_matrix)
  ) %>%
  mutate(
    Age = suppressWarnings(as.numeric(as.character(Age))),
    BMI = suppressWarnings(as.numeric(as.character(BMI))),
    Sex = as.factor(Sex),
    covariate_complete = is.finite(Age) & !is.na(Sex) & is.finite(BMI)
  )

if (sample_mode == "complete-case") {
  fixed_meta_all <- fixed_meta_all %>%
    filter(covariate_complete)
}

fixed_sample_qc <- fixed_meta_all %>%
  group_by(Cohort, Group, Age_class) %>%
  summarise(
    n = n(),
    complete_covariate_n = sum(covariate_complete),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = c(Age_class, Group),
    values_from = c(n, complete_covariate_n),
    values_fill = 0,
    names_glue = "{.value}_{Age_class}{Group}"
  )
write_csv(fixed_sample_qc, paste0(output_prefix, "_fixed_sample_qc.csv"))

message("Computing observed coefficients with the fixed 27-feature workflow...")
observed_tmp_root <- file.path(mmuphin_tmp_dir, "observed")
dir.create(observed_tmp_root, recursive = TRUE, showWarnings = FALSE)
res_EO_obs <- compute_meta_coefs(
  meta_df = fixed_meta_all %>% filter(Age_class == "EO", Cohort %in% eo_cohorts),
  abundance_matrix = target_abundance_matrix,
  age_group = "EO_observed",
  perm_id = 0,
  seed_value = seed_base,
  output_root = observed_tmp_root
)
res_LO_obs <- compute_meta_coefs(
  meta_df = fixed_meta_all %>% filter(Age_class == "LO", Cohort %in% lo_cohorts),
  abundance_matrix = target_abundance_matrix,
  age_group = "LO_observed",
  perm_id = 0,
  seed_value = seed_base,
  output_root = observed_tmp_root
)

obs_stat <- calculate_statistic(res_EO_obs, res_LO_obs, target_features)
observed_delta <- obs_stat$statistic

observed_coef_compare <- obs_stat$joined %>%
  left_join(
    target_27_table %>% select(feature, Species, main_coef_EO, main_coef_LO),
    by = "feature"
  ) %>%
  mutate(
    abs_diff_EO = abs(coef_EO - main_coef_EO),
    abs_diff_LO = abs(coef_LO - main_coef_LO)
  )
write_csv(observed_coef_compare, paste0(output_prefix, "_observed_coef_compare.csv"))

max_abs_diff_EO <- max(observed_coef_compare$abs_diff_EO)
max_abs_diff_LO <- max(observed_coef_compare$abs_diff_LO)
cor_EO <- cor(observed_coef_compare$coef_EO, observed_coef_compare$main_coef_EO)
cor_LO <- cor(observed_coef_compare$coef_LO, observed_coef_compare$main_coef_LO)

message(glue(
  "Observed statistic = {round(observed_delta, 6)}; ",
  "max abs coef diff EO = {signif(max_abs_diff_EO, 4)}, ",
  "LO = {signif(max_abs_diff_LO, 4)}"
))

if (max_abs_diff_EO > coef_tolerance || max_abs_diff_LO > coef_tolerance) {
  stop(
    "Observed coefficients do not match main analysis within tolerance. ",
    "Inspect ", paste0(output_prefix, "_observed_coef_compare.csv")
  )
}
run_single_permutation <- function(perm_id) {
  seed_value <- seed_base + perm_id
  perm_tmp_root <- file.path(
    mmuphin_tmp_dir,
    sprintf("perm_%05d_seed_%d_pid_%s", perm_id, seed_value, Sys.getpid())
  )
  dir.create(perm_tmp_root, recursive = TRUE, showWarnings = FALSE)

  tryCatch({
    set.seed(seed_value)

    meta_perm <- fixed_meta_all %>%
      group_by(Cohort, Group, covariate_complete) %>%
      mutate(Age_class_perm = sample(Age_class)) %>%
      ungroup()

    counts_before <- fixed_meta_all %>%
      count(Cohort, Group, covariate_complete, Age_class) %>%
      arrange(Cohort, Group, covariate_complete, Age_class)
    counts_after <- meta_perm %>%
      count(Cohort, Group, covariate_complete, Age_class = Age_class_perm) %>%
      arrange(Cohort, Group, covariate_complete, Age_class)
    if (!identical(counts_before, counts_after)) {
      stop("Age-stratum counts changed after permutation")
    }

    meta_EO_perm <- meta_perm %>%
      filter(Age_class_perm == "EO", Cohort %in% eo_cohorts) %>%
      select(-Age_class) %>%
      rename(Age_class = Age_class_perm)

    meta_LO_perm <- meta_perm %>%
      filter(Age_class_perm == "LO", Cohort %in% lo_cohorts) %>%
      select(-Age_class) %>%
      rename(Age_class = Age_class_perm)

    res_EO <- compute_meta_coefs(
      meta_df = meta_EO_perm,
      abundance_matrix = target_abundance_matrix,
      age_group = "EO",
      perm_id = perm_id,
      seed_value = seed_value,
      output_root = perm_tmp_root
    )

    res_LO <- compute_meta_coefs(
      meta_df = meta_LO_perm,
      abundance_matrix = target_abundance_matrix,
      age_group = "LO",
      perm_id = perm_id,
      seed_value = seed_value,
      output_root = perm_tmp_root
    )

    stat <- calculate_statistic(res_EO, res_LO, target_features)

    tibble(
      perm_id = perm_id,
      seed = seed_value,
      success = TRUE,
      statistic = stat$statistic,
      failure_reason = NA_character_,
      n_features = nrow(stat$joined)
    )
  }, error = function(e) {
    tibble(
      perm_id = perm_id,
      seed = seed_value,
      success = FALSE,
      statistic = NA_real_,
      failure_reason = conditionMessage(e),
      n_features = NA_integer_
    )
  })
}

if (n_permutations > 0) {
  message(glue("Running {n_permutations} permutations with {n_workers} worker(s)..."))
  if (n_workers > 1) {
    future::plan(multisession, workers = n_workers)
    on.exit(future::plan(sequential), add = TRUE)
  } else {
    future::plan(sequential)
  }

  permutation_results <- future_map_dfr(
    seq_len(n_permutations),
    run_single_permutation,
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  )
} else {
  permutation_results <- tibble(
    perm_id = integer(),
    seed = integer(),
    success = logical(),
    statistic = double(),
    failure_reason = character(),
    n_features = integer()
  )
}

write_csv(permutation_results, paste0(output_prefix, "_permutation_results.csv"))

success_stats <- permutation_results %>%
  filter(success, is.finite(statistic)) %>%
  pull(statistic)

B <- length(success_stats)
n_failed <- nrow(permutation_results) - B
b_extreme <- sum(success_stats >= observed_delta)
empirical_p <- if (B > 0) (b_extreme + 1) / (B + 1) else NA_real_
mc_se <- if (B > 0) sqrt(empirical_p * (1 - empirical_p) / B) else NA_real_
binom_ci <- if (B > 0) {
  stats::binom.test(b_extreme, B)$conf.int
} else {
  c(NA_real_, NA_real_)
}

summary_tbl <- tibble(
  n_requested = n_permutations,
  n_success = B,
  n_failed = n_failed,
  target_feature_n = length(target_features),
  observed_statistic = observed_delta,
  null_median = if (B > 0) median(success_stats) else NA_real_,
  null_q025 = if (B > 0) unname(quantile(success_stats, 0.025)) else NA_real_,
  null_q975 = if (B > 0) unname(quantile(success_stats, 0.975)) else NA_real_,
  b_extreme = b_extreme,
  empirical_p_plus1 = empirical_p,
  mc_se = mc_se,
  binom_ci_low_raw = binom_ci[[1]],
  binom_ci_high_raw = binom_ci[[2]],
  max_abs_diff_EO_vs_main = max_abs_diff_EO,
  max_abs_diff_LO_vs_main = max_abs_diff_LO,
  cor_EO_vs_main = cor_EO,
  cor_LO_vs_main = cor_LO,
  seed_base = seed_base,
  n_workers = n_workers,
  sample_mode = sample_mode,
  permutation_strata = "Cohort + Group + covariate_complete"
)
write_csv(summary_tbl, paste0(output_prefix, "_summary.csv"))

failure_summary <- permutation_results %>%
  filter(!success) %>%
  count(failure_reason, sort = TRUE)
write_csv(failure_summary, paste0(output_prefix, "_failure_summary.csv"))

if (n_failed > 0 && !allow_failures) {
  stop(
    n_failed,
    " permutations failed. Inspect ",
    paste0(output_prefix, "_failure_summary.csv"),
    " before inference."
  )
}

if (B > 0) {
  p_null <- tibble(statistic = success_stats) %>%
    ggplot(aes(x = statistic)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 60,
      fill = "grey85",
      color = "grey40",
      linewidth = 0.2
    ) +
    geom_density(color = "#2b8cbe", linewidth = 0.8) +
    geom_vline(xintercept = observed_delta, color = "#de2d26",
               linetype = "dashed", linewidth = 0.9) +
    geom_vline(xintercept = median(success_stats), color = "black",
               linewidth = 0.5) +
    labs(
      title = "Cohort-stratified permutation null distribution",
      subtitle = glue(
        "Fixed 27 species; B = {B}; empirical P = {signif(empirical_p, 4)}"
      ),
      x = expression("Median(" * abs(beta[EO]) - abs(beta[LO]) * ")"),
      y = "Density"
    ) +
    theme_classic(base_size = 12)

  ggsave(
    filename = paste0(plot_prefix, "_null_distribution.pdf"),
    plot = p_null,
    width = 8,
    height = 5,
    device = grDevices::pdf
  )
}

if (!keep_mmuphin_output) {
  unlink(mmuphin_tmp_dir, recursive = TRUE, force = TRUE)
}

message("Permutation analysis finished.")
message("Summary: ", paste0(output_prefix, "_summary.csv"))
message("Figure: ", paste0(plot_prefix, "_null_distribution.pdf"))
