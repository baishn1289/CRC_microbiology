
A30_metagenomic_id <- read_xlsx(meta_path, sheet = 2)
full_meta <- read_xlsx(meta_path, sheet = 1)

sample_id_col <- names(full_meta)[10]
batch_col <- names(full_meta)[11]
stage_col <- names(full_meta)[20]
tumor_site_col <- names(full_meta)[21]

expected_metadata_columns <- c(
  sample_id_col = "\u5b8f\u57fa\u56e0\u7ec4.\u4ee3\u8c22\u7f16\u53f7",
  batch_col = "Batch",
  stage_col = "TNM\u5206\u671f",
  tumor_site_col = "Location.0.Right.hemicolon.1..Left.hemicolon.2.Rectum."
)
actual_metadata_columns <- c(sample_id_col, batch_col, stage_col, tumor_site_col)
if (!identical(unname(actual_metadata_columns), unname(expected_metadata_columns))) {
  stop(
    "Unexpected metadata column layout. Expected: ",
    paste(expected_metadata_columns, collapse = " | "),
    "; observed: ",
    paste(actual_metadata_columns, collapse = " | ")
  )
}

FUSCC_cl2 <- full_meta %>%
  filter(ID %in% A30_metagenomic_id$`*ID`) %>%
  transmute(
    Patient_ID = ID,
    Sex = factor(gender),
    Age = as.numeric(age),
    BMI = as.numeric(BMI),
    Group = Group,
    Sample_ID = .data[[sample_id_col]],
    Batch = factor(.data[[batch_col]]),
    Stage_raw = .data[[stage_col]],
    Tumor_Site_raw = .data[[tumor_site_col]],
    Cohort = "A30_FUSCC",
    Country = "China",
    Diagnosis = Group
  ) %>%
  mutate(
    Group = ifelse(Group == 'CTRL',"Control",Group ),
    Disease = factor(Group, levels = c("Control", "CRC")),
    Age_group = factor(ifelse(Age < 50, "EO", "LO"), levels = c("LO", "EO")),
    Group4 = factor(paste0(Age_group, Disease), levels = group_levels),
    Age_within_group = Age - ave(Age, Age_group, FUN = function(x) mean(x, na.rm = TRUE)),
    Tumor_Site3 = factor(case_when(
      Tumor_Site_raw == 0 ~ "right_colon",
      Tumor_Site_raw == 1 ~ "left_colon",
      Tumor_Site_raw == 2 ~ "rectum",
      TRUE ~ "missing"
    ), levels = c("right_colon", "left_colon", "rectum", "missing")),
    Tumor_Site2 = factor(case_when(
      Tumor_Site3 %in% c("right_colon", "left_colon") ~ "colon",
      Tumor_Site3 == "rectum" ~ "rectum",
      TRUE ~ "missing"
    ), levels = c("colon", "rectum", "missing")),
    Stage_chr = str_trim(as.character(Stage_raw)),
    Stage_num = suppressWarnings(readr::parse_number(Stage_chr)),
    Stage_collapsed = factor(case_when(
      !is.na(Stage_num) & Stage_num <= 2 ~ "Stage0_2",
      !is.na(Stage_num) & Stage_num >= 3 ~ "Stage3_4",
      TRUE ~ "Missing_or_pCR"
    ), levels = c("Stage0_2", "Stage3_4", "Missing_or_pCR"))
  )

stopifnot(
  nrow(FUSCC_cl2) > 0,
  all(!is.na(FUSCC_cl2$Sample_ID)),
  !any(duplicated(FUSCC_cl2$Sample_ID))
)

write.csv(
  FUSCC_cl2 %>%
    group_by(Group4) %>%
    summarise(
      n = n(),
      age_mean = mean(Age, na.rm = TRUE),
      age_sd = sd(Age, na.rm = TRUE),
      bmi_mean = mean(BMI, na.rm = TRUE),
      bmi_sd = sd(BMI, na.rm = TRUE),
      female_n = sum(Sex == "female", na.rm = TRUE),
      male_n = sum(Sex == "male", na.rm = TRUE),
      .groups = "drop"
    ),
  file.path(table_dir, "clinical_group_summary_v29.csv"),
  row.names = FALSE
)

write.csv(
  FUSCC_cl2 %>%
    filter(Disease == "CRC") %>%
    count(Age_group, Stage_collapsed) %>%
    tidyr::pivot_wider(names_from = Stage_collapsed, values_from = n, values_fill = 0),
  file.path(table_dir, "crc_stage_by_age_group_v29.csv"),
  row.names = FALSE
)

write.csv(
  FUSCC_cl2 %>%
    filter(Disease == "CRC") %>%
    count(Age_group, Tumor_Site3) %>%
    tidyr::pivot_wider(names_from = Tumor_Site3, values_from = n, values_fill = 0),
  file.path(table_dir, "crc_tumor_site_by_age_group_v29.csv"),
  row.names = FALSE
)

LC_MS_negative <- read_mode_file(file.path(analysis_dir, "OMIX006518-09.csv"), "negative")
LC_MS_positive <- read_mode_file(file.path(analysis_dir, "OMIX006518-08.csv"), "positive")
mode_dedup_scored <- score_mode_candidates(bind_rows(LC_MS_negative, LC_MS_positive))
LC_MS_selected <- deduplicate_modes(mode_dedup_scored)

write.csv(
  mode_dedup_scored %>% arrange(metabolite_identification, selection_rank),
  file.path(table_dir, "targeted_mode_dedup_decisions.csv"),
  row.names = FALSE
)

sample_cols <- intersect(names(LC_MS_selected), FUSCC_cl2$Sample_ID)
missing_lcms_samples <- setdiff(FUSCC_cl2$Sample_ID, names(LC_MS_selected))
if (length(missing_lcms_samples) > 0) {
  stop("LC-MS table is missing metadata samples: ", paste(missing_lcms_samples, collapse = ", "))
}

LC_long <- LC_MS_selected %>%
  select(all_of(c("metabolite_identification", sample_cols))) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample_ID",
    values_to = "intensity"
  ) %>%
  left_join(FUSCC_cl2, by = "Sample_ID")

metab_qc <- LC_long %>%
  group_by(metabolite_identification) %>%
  summarise(
    overall_missing_rate = mean(is.na(intensity)),
    overall_detect_rate = mean(!is.na(intensity) & intensity > 0),
    .groups = "drop"
  ) %>%
  left_join(
    LC_long %>%
      group_by(metabolite_identification, Group4) %>%
      summarise(group_detect_rate = mean(!is.na(intensity) & intensity > 0), .groups = "drop") %>%
      group_by(metabolite_identification) %>%
      summarise(max_group_detect_rate = max(group_detect_rate), .groups = "drop"),
    by = "metabolite_identification"
  ) %>%
  mutate(
    keep_feature = overall_missing_rate <= max_missing_fraction &
      max_group_detect_rate >= min_present_fraction
  )

write.csv(metab_qc, file.path(table_dir, "metabolite_qc_summary.csv"), row.names = FALSE)

keep_metabs <- metab_qc %>%
  filter(keep_feature) %>%
  pull(metabolite_identification)

LC_long2 <- LC_long %>%
  filter(metabolite_identification %in% keep_metabs) %>%
  group_by(metabolite_identification) %>%
  mutate(
    min_pos = suppressWarnings(min(intensity[intensity > 0], na.rm = TRUE)),
    min_pos = ifelse(is.infinite(min_pos), NA_real_, min_pos),
    intensity_imp = case_when(
      !is.na(intensity) & intensity > 0 ~ intensity,
      is.na(min_pos) ~ NA_real_,
      TRUE ~ min_pos / 2
    ),
    log2_intensity = log2(intensity_imp)
  ) %>%
  ungroup()

feature_mat <- LC_long2 %>%
  select(metabolite_identification, Sample_ID, log2_intensity) %>%
  distinct() %>%
  pivot_wider(names_from = Sample_ID, values_from = log2_intensity) %>%
  as.data.frame()

rownames(feature_mat) <- feature_mat$metabolite_identification
feature_mat$metabolite_identification <- NULL
feature_mat <- as.matrix(feature_mat)
feature_mat <- feature_mat[, FUSCC_cl2$Sample_ID, drop = FALSE]

sample_mat_log2 <- t(feature_mat)
sample_mat_pareto <- pareto_scale_matrix(sample_mat_log2)

pca_obj <- prcomp(sample_mat_pareto, center = FALSE, scale. = FALSE)
pca_var <- round(100 * summary(pca_obj)$importance[2, 1:2], 1)

pca_df <- as.data.frame(pca_obj$x[, 1:5]) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(FUSCC_cl2, by = "Sample_ID")

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group4, shape = Disease)) +
  geom_point(size = 2.0, alpha = 0.75) +
  stat_ellipse(aes(group = Group4), linewidth = 0.6, alpha = 0.5) +
  scale_color_manual(values = group_palette) +
  labs(
    title = "Pareto-scaled PCA",
    x = paste0("PC1 (", pca_var[1], "%)"),
    y = paste0("PC2 (", pca_var[2], "%)")
  ) +
  theme(panel.grid.minor = element_blank())

save_gg(pca_plot, file.path(plot_dir, "pca_group4.pdf"), width = 8.5, height = 6.5)

adonis_group4_df <- as.data.frame(
  vegan::adonis2(
    sample_mat_pareto ~ Group4 + Age_within_group + Sex + BMI,
    data = FUSCC_cl2,
    method = "euclidean",
    by = "margin",
    permutations = 999
  )
) %>%
  rownames_to_column("term")

adonis_interaction_df <- as.data.frame(
  vegan::adonis2(
    # by = "terms" is sequential, so covariates must be entered before the
    # disease-by-age interaction to avoid attributing covariate variance to it.
    sample_mat_pareto ~ Age_within_group + Sex + BMI + Disease * Age_group,
    data = FUSCC_cl2,
    method = "euclidean",
    by = "terms",
    permutations = 999
  )
) %>%
  rownames_to_column("term")

write.csv(adonis_group4_df, file.path(table_dir, "adonis2_group4_overall.csv"), row.names = FALSE)
write.csv(
  adonis_interaction_df,
  file.path(table_dir, "adonis2_disease_age_interaction_terms.csv"),
  row.names = FALSE
)

dispersion_obj <- vegan::betadisper(dist(sample_mat_pareto), group = FUSCC_cl2$Group4)
dispersion_anova <- as.data.frame(anova(dispersion_obj)) %>%
  rownames_to_column("term")
dispersion_perm <- as.data.frame(run_betadisper_permutest(dispersion_obj, permutations = 999)$tab) %>%
  rownames_to_column("term")

write.csv(dispersion_anova, file.path(table_dir, "betadisper_anova.csv"), row.names = FALSE)
write.csv(dispersion_perm, file.path(table_dir, "betadisper_permutest.csv"), row.names = FALSE)

dispersion_df <- FUSCC_cl2 %>%
  mutate(distance_to_centroid = dispersion_obj$distances)

dispersion_plot <- ggplot(dispersion_df, aes(x = Group4, y = distance_to_centroid, fill = Group4)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.85) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 0.8, color = "black") +
  scale_fill_manual(values = group_palette) +
  labs(
    title = "Multivariate dispersion by four-group status",
    x = NULL,
    y = "Distance to centroid"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank()
  )

save_gg(dispersion_plot, file.path(plot_dir, "group_dispersion_boxplot.pdf"), width = 8, height = 6)

panel_membership <- build_targeted_panel_membership(rownames(feature_mat))
write.csv(panel_membership, file.path(table_dir, "targeted_panel_membership.csv"), row.names = FALSE)
write.csv(
  panel_membership %>%
    distinct(panel, panel_label, panel_tier, panel_tier_label, canonical_metabolite) %>%
    arrange(panel, canonical_metabolite),
  file.path(table_dir, "targeted_panel_catalog_v29.csv"),
  row.names = FALSE
)
write.csv(
  build_focus_metabolite_review(rownames(feature_mat)),
  file.path(table_dir, "panel_design_review_v29.csv"),
  row.names = FALSE
)
write.csv(
  mode_dedup_scored %>%
    filter(selected_mode) %>%
    semi_join(
      panel_membership %>%
        filter(panel_tier == "strict_manuscript_aligned") %>%
        distinct(metabolite_identification),
      by = "metabolite_identification"
    ) %>%
    arrange(metabolite_identification, mode) %>%
    select(
      metabolite_identification,
      mode,
      prefer_mode,
      mode_rank,
      selection_rank,
      selected_mode,
      selection_reason,
      n_detected,
      missing_rate,
      median_signal
    ),
  file.path(table_dir, "targeted_mode_dedup_strict_metabolites.csv"),
  row.names = FALSE
)

panel_score_coverage_by_sample <- build_panel_score_coverage_table(
  raw_long_df = LC_long %>% filter(metabolite_identification %in% rownames(feature_mat)),
  panel_membership = panel_membership
)
write.csv(
  panel_score_coverage_by_sample,
  file.path(table_dir, "panel_score_coverage_by_sample.csv"),
  row.names = FALSE
)
write.csv(
  summarise_panel_score_coverage(panel_score_coverage_by_sample),
  file.path(table_dir, "panel_score_coverage_summary_by_group.csv"),
  row.names = FALSE
)

eo_meta <- FUSCC_cl2 %>%
  filter(Age_group == "EO") %>%
  mutate(Age_centered = Age - mean(Age, na.rm = TRUE))

lo_meta <- FUSCC_cl2 %>%
  filter(Age_group == "LO") %>%
  mutate(Age_centered = Age - mean(Age, na.rm = TRUE))

crc_meta <- FUSCC_cl2 %>%
  filter(Disease == "CRC") %>%
  group_by(Age_group) %>%
  mutate(Age_within_crc_group = Age - mean(Age, na.rm = TRUE)) %>%
  ungroup()

write.csv(
  tibble(
    model_label = c(
      "EO_stratified_main",
      "LO_stratified_main",
      "pooled_interaction_ageadj",
      "CRC_only_sensitivity"
    ),
    n_samples = c(nrow(eo_meta), nrow(lo_meta), nrow(FUSCC_cl2), nrow(crc_meta))
  ),
  file.path(table_dir, "model_sample_sizes_v29.csv"),
  row.names = FALSE
)

eo_feature_fit <- fit_limma_contrasts(
  expr_mat = feature_mat[, eo_meta$Sample_ID, drop = FALSE],
  sample_meta = eo_meta,
  design_formula = ~ Disease + Age_centered + Sex + BMI,
  contrast_definitions = c(EOCRC_vs_EOControl_ageadj = "DiseaseCRC"),
  model_label = "EO_stratified_main",
  test_type = "treat",
  lfc_threshold = effect_cutoff
)

lo_feature_fit <- fit_limma_contrasts(
  expr_mat = feature_mat[, lo_meta$Sample_ID, drop = FALSE],
  sample_meta = lo_meta,
  design_formula = ~ Disease + Age_centered + Sex + BMI,
  contrast_definitions = c(LOCRC_vs_LOControl_ageadj = "DiseaseCRC"),
  model_label = "LO_stratified_main",
  test_type = "treat",
  lfc_threshold = effect_cutoff
)

pooled_feature_fit <- fit_limma_contrasts(
  expr_mat = feature_mat,
  sample_meta = FUSCC_cl2,
  design_formula = ~ Disease * Age_group + Age_within_group + Sex + BMI,
  contrast_definitions = c(
    LOCRC_vs_LOControl_pooled_ageadj = "DiseaseCRC",
    EOCRC_vs_EOControl_pooled_ageadj = "DiseaseCRC + DiseaseCRC.Age_groupEO",
    Disease_age_interaction_pooled_ageadj = "DiseaseCRC.Age_groupEO"
  ),
  model_label = "pooled_interaction_ageadj",
  test_type = "treat",
  lfc_threshold = effect_cutoff
)

crc_feature_fit <- fit_limma_contrasts(
  expr_mat = feature_mat[, crc_meta$Sample_ID, drop = FALSE],
  sample_meta = crc_meta,
  design_formula = ~ Age_group + Age_within_crc_group + Sex + BMI + Tumor_Site3 + Stage_collapsed,
  contrast_definitions = c(EOCRC_vs_LOCRC_crc_age_stage_site_adj = "Age_groupEO"),
  model_label = "CRC_only_sensitivity",
  test_type = "treat",
  lfc_threshold = effect_cutoff
)

feature_results_display <- bind_rows(
  eo_feature_fit$results,
  lo_feature_fit$results,
  crc_feature_fit$results
)

feature_results_pooled <- pooled_feature_fit$results
feature_results_all <- bind_rows(feature_results_display, feature_results_pooled)

write.csv(
  feature_results_display,
  file.path(table_dir, "feature_results_display_models_treat.csv"),
  row.names = FALSE
)
write.csv(
  feature_results_pooled,
  file.path(table_dir, "feature_results_pooled_interaction_treat.csv"),
  row.names = FALSE
)
write.csv(
  feature_results_all,
  file.path(table_dir, "feature_results_all_models_treat.csv"),
  row.names = FALSE
)

write.csv(
  feature_results_display %>% filter(significant_fdr),
  file.path(table_dir, "feature_results_display_models_treat_fdr_lt_0_05.csv"),
  row.names = FALSE
)

targeted_results_display <- feature_results_display %>%
  inner_join(
    panel_membership,
    by = "metabolite_identification",
    relationship = "many-to-many"
  ) %>%
  arrange(contrast, adj.P.Val, canonical_metabolite, metabolite_identification)

targeted_results_pooled <- feature_results_pooled %>%
  inner_join(
    panel_membership,
    by = "metabolite_identification",
    relationship = "many-to-many"
  ) %>%
  arrange(contrast, adj.P.Val, canonical_metabolite, metabolite_identification)

targeted_results_all <- bind_rows(targeted_results_display, targeted_results_pooled)

write.csv(
  targeted_results_display,
  file.path(table_dir, "targeted_panel_feature_results_display_models.csv"),
  row.names = FALSE
)
write.csv(
  targeted_results_pooled,
  file.path(table_dir, "targeted_panel_feature_results_pooled_interaction.csv"),
  row.names = FALSE
)
write.csv(
  targeted_results_all,
  file.path(table_dir, "targeted_panel_feature_results_all_models.csv"),
  row.names = FALSE
)

significant_targeted_display <- targeted_results_display %>%
  filter(significant_fdr)

canonical_targeted_display <- targeted_results_display %>%
  group_by(model_label, contrast, panel, panel_label, panel_tier, panel_tier_label, canonical_metabolite) %>%
  summarise(
    best_adj.P.Val = min(adj.P.Val, na.rm = TRUE),
    best_P.Value = min(P.Value, na.rm = TRUE),
    any_significant_fdr = any(significant_fdr),
    any_significant_fdr_fc = any(significant_fdr_fc),
    mean_logFC = mean(logFC, na.rm = TRUE),
    n_raw_features = n_distinct(metabolite_identification),
    raw_features = paste(sort(unique(metabolite_identification)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(direction = ifelse(mean_logFC >= 0, "higher", "lower"))

canonical_targeted_pooled <- targeted_results_pooled %>%
  group_by(model_label, contrast, panel, panel_label, panel_tier, panel_tier_label, canonical_metabolite) %>%
  summarise(
    best_adj.P.Val = min(adj.P.Val, na.rm = TRUE),
    best_P.Value = min(P.Value, na.rm = TRUE),
    any_significant_fdr = any(significant_fdr),
    any_significant_fdr_fc = any(significant_fdr_fc),
    mean_logFC = mean(logFC, na.rm = TRUE),
    n_raw_features = n_distinct(metabolite_identification),
    raw_features = paste(sort(unique(metabolite_identification)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(direction = ifelse(mean_logFC >= 0, "higher", "lower"))

significant_targeted_display_canonical <- canonical_targeted_display %>%
  filter(any_significant_fdr)

write.csv(
  significant_targeted_display,
  file.path(table_dir, "targeted_panel_feature_results_display_models_fdr_lt_0_05.csv"),
  row.names = FALSE
)
write.csv(
  canonical_targeted_display,
  file.path(table_dir, "targeted_panel_feature_results_display_models_canonical_summary.csv"),
  row.names = FALSE
)
write.csv(
  significant_targeted_display_canonical,
  file.path(table_dir, "targeted_panel_feature_results_display_models_canonical_fdr_lt_0_05.csv"),
  row.names = FALSE
)
write.csv(
  canonical_targeted_display %>% filter(any_significant_fdr),
  file.path(table_dir, "targeted_panel_feature_results_display_models_canonical_fdr_lt_0_05_only.csv"),
  row.names = FALSE
)
write.csv(
  canonical_targeted_pooled,
  file.path(table_dir, "targeted_panel_feature_results_pooled_interaction_canonical_summary.csv"),
  row.names = FALSE
)

panel_score_mat <- build_panel_score_matrix(feature_mat, panel_membership)

eo_panel_fit <- fit_limma_contrasts(
  expr_mat = panel_score_mat[, eo_meta$Sample_ID, drop = FALSE],
  sample_meta = eo_meta,
  design_formula = ~ Disease + Age_centered + Sex + BMI,
  contrast_definitions = c(EOCRC_vs_EOControl_ageadj = "DiseaseCRC"),
  model_label = "EO_stratified_main",
  test_type = "ebayes",
  lfc_threshold = 0
)

lo_panel_fit <- fit_limma_contrasts(
  expr_mat = panel_score_mat[, lo_meta$Sample_ID, drop = FALSE],
  sample_meta = lo_meta,
  design_formula = ~ Disease + Age_centered + Sex + BMI,
  contrast_definitions = c(LOCRC_vs_LOControl_ageadj = "DiseaseCRC"),
  model_label = "LO_stratified_main",
  test_type = "ebayes",
  lfc_threshold = 0
)

pooled_panel_fit <- fit_limma_contrasts(
  expr_mat = panel_score_mat,
  sample_meta = FUSCC_cl2,
  design_formula = ~ Disease * Age_group + Age_within_group + Sex + BMI,
  contrast_definitions = c(
    LOCRC_vs_LOControl_pooled_ageadj = "DiseaseCRC",
    EOCRC_vs_EOControl_pooled_ageadj = "DiseaseCRC + DiseaseCRC.Age_groupEO",
    Disease_age_interaction_pooled_ageadj = "DiseaseCRC.Age_groupEO"
  ),
  model_label = "pooled_interaction_ageadj",
  test_type = "ebayes",
  lfc_threshold = 0
)

crc_panel_fit <- fit_limma_contrasts(
  expr_mat = panel_score_mat[, crc_meta$Sample_ID, drop = FALSE],
  sample_meta = crc_meta,
  design_formula = ~ Age_group + Age_within_crc_group + Sex + BMI + Tumor_Site3 + Stage_collapsed,
  contrast_definitions = c(EOCRC_vs_LOCRC_crc_age_stage_site_adj = "Age_groupEO"),
  model_label = "CRC_only_sensitivity",
  test_type = "ebayes",
  lfc_threshold = 0
)

panel_results_display <- bind_rows(
  eo_panel_fit$results,
  lo_panel_fit$results,
  crc_panel_fit$results
) %>%
  rename(panel = metabolite_identification) %>%
  mutate(
    panel_label = recode(panel, !!!panel_labels),
    panel_tier = recode(panel, !!!panel_tier_map),
    panel_tier_label = recode(panel_tier, !!!panel_tier_labels)
  )

panel_results_pooled <- pooled_panel_fit$results %>%
  rename(panel = metabolite_identification) %>%
  mutate(
    panel_label = recode(panel, !!!panel_labels),
    panel_tier = recode(panel, !!!panel_tier_map),
    panel_tier_label = recode(panel_tier, !!!panel_tier_labels)
  )

panel_results_all <- bind_rows(panel_results_display, panel_results_pooled)

write.csv(
  panel_results_display,
  file.path(table_dir, "targeted_panel_score_results_display_models.csv"),
  row.names = FALSE
)
write.csv(
  panel_results_pooled,
  file.path(table_dir, "targeted_panel_score_results_pooled_interaction.csv"),
  row.names = FALSE
)
write.csv(
  panel_results_all,
  file.path(table_dir, "targeted_panel_score_results_all_models.csv"),
  row.names = FALSE
)

panel_score_long <- as.data.frame(t(panel_score_mat)) %>%
  rownames_to_column("Sample_ID") %>%
  pivot_longer(
    cols = -Sample_ID,
    names_to = "panel",
    values_to = "panel_score"
  ) %>%
  mutate(
    panel = factor(panel, levels = panel_levels),
    panel_label = recode(panel, !!!panel_labels),
    panel_tier = recode(panel, !!!panel_tier_map),
    panel_tier_label = recode(panel_tier, !!!panel_tier_labels)
  ) %>%
  left_join(FUSCC_cl2, by = "Sample_ID")

write.csv(panel_score_long, file.path(table_dir, "targeted_panel_scores_by_sample.csv"), row.names = FALSE)

panel_score_annotations <- build_panel_score_annotations(panel_score_long, panel_results_display)
write.csv(
  panel_score_annotations,
  file.path(table_dir, "targeted_panel_score_annotations_v29.csv"),
  row.names = FALSE
)

panel_score_plot <- plot_panel_scores(panel_score_long, 
                                      panel_score_annotations, 
                                      label_style = "beta_sig",
                                      ncol = 3,eolocr_lift = 0.15,
                                      text_size = 12)
save_gg(panel_score_plot, file.path(plot_dir, "targeted_panel_scores.pdf"), 
        width = 15, height = 13)

strict_panel_score_plot <- plot_panel_scores(
  panel_score_long %>% filter(panel_tier == "strict_manuscript_aligned"),
  panel_score_annotations %>% filter(panel_tier == "strict_manuscript_aligned"),
  label_style = "no_sig_only",beta_offset = 0.03,sig_offset = 0.04,
  ncol = 1
)
save_gg(
  strict_panel_score_plot,
  file.path(plot_dir, "targeted_panel_scores_strict.pdf"),
  width = 4.5,
  height = 6.5
)

extended_panel_score_plot <- plot_panel_scores(
  panel_score_long %>% filter(panel_tier == "extended_class_based"),
  panel_score_annotations %>% filter(panel_tier == "extended_class_based"),
  label_style = "beta_sig",ncol = 5,
  beta_offset = 0.08,sig_offset = 0.08,eolocr_lift = 0.2,
  text_size = 10
)
save_gg(
  extended_panel_score_plot,
  file.path(plot_dir, "targeted_panel_scores_extended.pdf"),
  width = 12,
  height = 7
)

volcano_plots <- map(
  contrast_display_order,
  ~ plot_volcano(
    df = feature_results_display %>% filter(contrast == .x),
    contrast_name = .x,
    targeted_membership = panel_membership
  )
)

volcano_panel <- wrap_plots(volcano_plots, ncol = 2) +
  plot_annotation(title = "Differential metabolite volcano plots for display models")

save_gg(volcano_panel, file.path(plot_dir, "volcano_plots_display_models.pdf"), width = 14, height = 10)

targeted_canonical_long <- build_targeted_canonical_long(LC_long2, panel_membership)
write.csv(
  targeted_canonical_long,
  file.path(table_dir, "targeted_panel_canonical_scores_by_sample.csv"),
  row.names = FALSE
)

targeted_heatmap_plot <- plot_targeted_heatmap(
  long_df = targeted_canonical_long,
  significant_targeted_panel_df = significant_targeted_display_canonical
)

if (!is.null(targeted_heatmap_plot)) {
  save_gg(
    targeted_heatmap_plot,
    file.path(plot_dir, "significant_targeted_metabolite_heatmap.pdf"),
    width = 8.5,
    height = 10
  )
}

all_panel_member_heatmap <- plot_panel_member_heatmap(
  long_df = targeted_canonical_long,
  panel_lookup = panel_membership,
  feature_results = canonical_targeted_display
)

if (!is.null(all_panel_member_heatmap)) {
  save_gg(
    all_panel_member_heatmap,
    file.path(plot_dir, "all_panel_members_heatmap.pdf"),
    width = 11,
    height = 14
  )
}

strict_targeted_heatmap_plot <- plot_targeted_heatmap(
  long_df = targeted_canonical_long,
  significant_targeted_panel_df = significant_targeted_display_canonical %>%
    filter(panel_tier == "strict_manuscript_aligned")
)

if (!is.null(strict_targeted_heatmap_plot)) {
  save_gg(
    strict_targeted_heatmap_plot,
    file.path(plot_dir, "significant_targeted_metabolite_heatmap_strict.pdf"),
    width = 8.2,
    height = 5.8
  )
}

extended_targeted_heatmap_plot <- plot_targeted_heatmap(
  long_df = targeted_canonical_long,
  significant_targeted_panel_df = significant_targeted_display_canonical %>%
    filter(panel_tier == "extended_class_based")
)

if (!is.null(extended_targeted_heatmap_plot)) {
  save_gg(
    extended_targeted_heatmap_plot,
    file.path(plot_dir, "significant_targeted_metabolite_heatmap_extended.pdf"),
    width = 8.8,
    height = 10.5
  )
}

strict_all_panel_member_heatmap <- plot_panel_member_heatmap(
  long_df = targeted_canonical_long,
  panel_lookup = panel_membership %>% filter(panel_tier == "strict_manuscript_aligned"),
  feature_results = canonical_targeted_display %>% filter(panel_tier == "strict_manuscript_aligned")
)

if (!is.null(strict_all_panel_member_heatmap)) {
  save_gg(
    strict_all_panel_member_heatmap,
    file.path(plot_dir, "all_panel_members_heatmap_strict.pdf"),
    width = 8.4,
    height = 6.2
  )
}

extended_all_panel_member_heatmap <- plot_panel_member_heatmap(
  long_df = targeted_canonical_long,
  panel_lookup = panel_membership %>% filter(panel_tier == "extended_class_based"),
  feature_results = canonical_targeted_display %>% filter(panel_tier == "extended_class_based")
)

if (!is.null(extended_all_panel_member_heatmap)) {
  save_gg(
    extended_all_panel_member_heatmap,
    file.path(plot_dir, "all_panel_members_heatmap_extended.pdf"),
    width = 9.2,
    height = 11
  )
}

predefined_focus <- c(
  "5-Aminopentanoic acid",
  "N-Acetylputrescine",
  "Citrulline",
  "D-Ornithine",
  "N-Acetylornithine",
  "Indole-3-acetaldehyde",
  "Indole-3-carboxyaldehyde",
  "Indoxyl sulfate",
  "4-Hydroxyphenylacetylglutamine",
  "Phenylacetic acid",
  "N-Acetylputrescine",
  "N-Acetylcadaverine",
  "Histamine",
  "Tyramine",
  "p-Cresol",
  "p-Cresol sulfate",
  "Phenol",
  "Phenol sulphate",
  "Indole",
  "Indole-3-acetamide",
  "Arginine"
)

representative_feature_map <- panel_membership %>%
  group_by(canonical_metabolite) %>%
  summarise(
    representative_feature = sort(unique(metabolite_identification))[1],
    .groups = "drop"
  )

top_targeted_features <- significant_targeted_display_canonical %>%
  group_by(canonical_metabolite) %>%
  summarise(best_fdr = min(best_adj.P.Val), .groups = "drop") %>%
  arrange(best_fdr) %>%
  slice_head(n = 12) %>%
  left_join(representative_feature_map, by = "canonical_metabolite") %>%
  pull(representative_feature)

focus_features <- unique(c(predefined_focus, top_targeted_features))
focus_features <- focus_features[focus_features %in% rownames(feature_mat)]
focus_features <- focus_features[seq_len(min(length(focus_features), 15))]

focus_boxplot <- plot_focus_boxplots(LC_long2, focus_features)

if (!is.null(focus_boxplot)) {
  save_gg(
    focus_boxplot,
    file.path(plot_dir, "focus_metabolite_boxplots.pdf"),
    width = 13,
    height = 13
  )
}

main_figure_panel_effects <- build_main_figure_panel_effects(
  panel_score_long = panel_score_long,
  panel_results_display = panel_results_display
)

main_figure_canonical_effects <- build_main_figure_canonical_effects(
  canonical_targeted_pooled = canonical_targeted_pooled
)

main_figure_leave_one_out <- run_main_figure_leave_one_out(
  feature_mat = feature_mat,
  panel_membership = panel_membership,
  eo_meta = eo_meta
)

write.csv(
  main_figure_panel_effects,
  file.path(table_dir, "main_figure_panel_effects_v29.csv"),
  row.names = FALSE
)
write.csv(
  main_figure_canonical_effects %>%
    select(-text_color, -max_abs),
  file.path(table_dir, "main_figure_canonical_effects_v29.csv"),
  row.names = FALSE
)
write.csv(
  main_figure_leave_one_out,
  file.path(table_dir, "main_figure_leave_one_out_v29.csv"),
  row.names = FALSE
)

main_panel_effect_plot <- plot_main_figure_panel_effects(main_figure_panel_effects %>% 
                                                           mutate(panel_label =  factor(panel_label,
                                                                                        levels = c('Arginine-ornithine routing axis',
                                                                                                   'Amino-acid putrefaction core')  )))
main_canonical_effect_plot <- plot_main_figure_canonical_effects(main_figure_canonical_effects)
main_leave_one_out_plot <- plot_main_figure_leave_one_out(main_figure_leave_one_out)

save_gg(
  main_panel_effect_plot,
  file.path(plot_dir, "main_figure_panel_effects_v29.pdf"),
  width = 8,
  height = 8
)
save_gg(
  main_canonical_effect_plot,
  file.path(plot_dir, "main_figure_canonical_effects_v29.pdf"),
  width = 10,
  height = 8
)
save_gg(
  main_leave_one_out_plot,
  file.path(plot_dir, "main_figure_leave_one_out_v29.pdf"),
  width = 11,
  height = 8
)

main_figure_composite <- ((main_panel_effect_plot / main_leave_one_out_plot) | main_canonical_effect_plot) +
  plot_layout(widths = c(1.05, 1.2), guides = "collect") +
  plot_annotation(
    title = "EOCRC nitrogen diversion axis: effect size, metabolite order, and robustness",
    subtitle = paste(
      "LO effect = LOCRC - LOControl; EO effect = EOCRC - EOControl; EO > LO = disease-by-age interaction.",
      "Leave-one-out sensitivity uses the EO age-adjusted model."
    )
  ) &
  theme(legend.position = "bottom")

contrast_summary <- feature_results_display %>%
  group_by(model_label, contrast) %>%
  summarise(
    n_fdr_sig = sum(significant_fdr),
    n_fdr_fc_sig = sum(significant_fdr_fc),
    .groups = "drop"
  )

pooled_interaction_summary <- feature_results_pooled %>%
  group_by(model_label, contrast) %>%
  summarise(
    n_fdr_sig = sum(significant_fdr),
    n_fdr_fc_sig = sum(significant_fdr_fc),
    .groups = "drop"
  )

targeted_contrast_summary <- targeted_results_display %>%
  group_by(model_label, contrast, panel, panel_label, panel_tier, panel_tier_label) %>%
  summarise(
    n_fdr_sig = n_distinct(canonical_metabolite[significant_fdr]),
    n_fdr_fc_sig = n_distinct(canonical_metabolite[significant_fdr_fc]),
    n_panel_members = n_distinct(canonical_metabolite),
    .groups = "drop"
  )

write.csv(contrast_summary, file.path(table_dir, "contrast_level_hit_counts_display_models.csv"), row.names = FALSE)
write.csv(
  pooled_interaction_summary,
  file.path(table_dir, "contrast_level_hit_counts_pooled_interaction.csv"),
  row.names = FALSE
)
write.csv(
  targeted_contrast_summary,
  file.path(table_dir, "targeted_panel_hit_counts_by_contrast.csv"),
  row.names = FALSE
)

writeLines(
  c(
    paste0("Formal cross-omics framework (", analysis_version, ")"),
    "1. Primary discovery is restricted to the EO stratum (EOCRC vs EOControl), adjusted for continuous age, sex, and BMI.",
    "2. Candidate species must first associate with the target metabolite panel after covariate adjustment, rather than being selected only by crude abundance differences.",
    "3. EO disease association is evaluated separately so that panel-linked species can be distinguished from background panel correlates.",
    "4. EOCRC-only sensitivity models add tumor site and stage to test whether panel-species links persist within EO cancers.",
    "5. The former injurious_luminal_panel is split into biogenic amine, parent phenol, and host-processed aromatic toxin subpanels to avoid chemistry-driven signal cancellation.",
    "6. Pooled Disease x Age_group interaction is retained as secondary EO-specificity support, not as the primary discovery screen.",
    "7. HUMAnN MetaCyc/EC anchors are evaluated at the community level as mechanistic support, not as proof that a specific species encodes the function.",
    "8. Species prevalence filtering is evaluated separately in EO and LO using identical criteria; the union defines the common CLR feature space.",
    "9. EO models are interpreted for EO-eligible species, LO models for LO-eligible species, and pooled interaction models for species eligible in both strata.",
    "10. Candidate support labels are descriptive (EO-supported only, LO-supported only, supported in both strata) and do not by themselves prove age specificity.",
    "11. Formal significance is defined only by FDR < 0.05 in this version; the former FDR < 0.10 exploratory layer is disabled."
  ),
  file.path(output_dir, "crossomics_analysis_framework_v29.txt")
)
