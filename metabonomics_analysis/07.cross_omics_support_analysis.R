
loaded_rdata_objects <- load(file.path(analysis_dir, "FUSCC_metabonomics.Rdata"))
required_rdata_objects <- c("meta_list_species", "metacyc_list", "ec_list")
missing_rdata_objects <- setdiff(required_rdata_objects, loaded_rdata_objects)
if (length(missing_rdata_objects) > 0) {
  stop("FUSCC_metabonomics.Rdata is missing required objects: ", paste(missing_rdata_objects, collapse = ", "))
}
if (!"FUSCC-SHSD" %in% names(meta_list_species)) {
  stop("meta_list_species does not contain the required FUSCC-SHSD entry.")
}
if (!"A30" %in% names(metacyc_list)) {
  stop("metacyc_list does not contain the required A30 entry.")
}
if (!"A30" %in% names(ec_list)) {
  stop("ec_list does not contain the required A30 entry.")
}

omics_sample_map <- FUSCC_cl2 %>%
  select(Patient_ID, Sample_ID) %>%
  distinct()

eo_crc_only_meta <- FUSCC_cl2 %>%
  filter(Age_group == "EO", Disease == "CRC") %>%
  mutate(Age_centered_crc = Age - mean(Age, na.rm = TRUE))

lo_crc_only_meta <- FUSCC_cl2 %>%
  filter(Age_group == "LO", Disease == "CRC") %>%
  mutate(Age_centered_crc = Age - mean(Age, na.rm = TRUE))

species_raw_mat <- prepare_crossomic_matrix(
  raw_df = meta_list_species[["FUSCC-SHSD"]],
  feature_col = "clade_name",
  sample_map = omics_sample_map
)
species_raw_mat <- species_raw_mat[, FUSCC_cl2$Sample_ID, drop = FALSE]

species_qc_eo <- summarise_crossomic_qc(
  feature_mat = species_raw_mat,
  sample_meta = FUSCC_cl2,
  subset_mask = FUSCC_cl2$Age_group == "EO",
  group_var = "Disease",
  min_prevalence = species_prevalence_min,
  min_group_prevalence = species_group_prevalence_min
)

species_qc_lo <- summarise_crossomic_qc(
  feature_mat = species_raw_mat,
  sample_meta = FUSCC_cl2,
  subset_mask = FUSCC_cl2$Age_group == "LO",
  group_var = "Disease",
  min_prevalence = species_prevalence_min,
  min_group_prevalence = species_group_prevalence_min
)

keep_species_eo <- species_qc_eo %>%
  filter(keep_feature) %>%
  pull(feature_id)

keep_species_lo <- species_qc_lo %>%
  filter(keep_feature) %>%
  pull(feature_id)

keep_species_union <- union(keep_species_eo, keep_species_lo)
keep_species_both_strata <- intersect(keep_species_eo, keep_species_lo)

species_qc <- full_join(
  species_qc_eo %>%
    rename(
      eo_overall_prevalence = overall_prevalence,
      eo_max_group_prevalence = max_group_prevalence,
      eo_keep_feature = keep_feature
    ),
  species_qc_lo %>%
    rename(
      lo_overall_prevalence = overall_prevalence,
      lo_max_group_prevalence = max_group_prevalence,
      lo_keep_feature = keep_feature
    ),
  by = "feature_id"
) %>%
  mutate(
    eo_keep_feature = coalesce(eo_keep_feature, FALSE),
    lo_keep_feature = coalesce(lo_keep_feature, FALSE),
    keep_feature = eo_keep_feature | lo_keep_feature,
    keep_feature_both_strata = eo_keep_feature & lo_keep_feature,
    overall_prevalence = pmax(
      coalesce(eo_overall_prevalence, 0),
      coalesce(lo_overall_prevalence, 0)
    ),
    max_group_prevalence = pmax(
      coalesce(eo_max_group_prevalence, 0),
      coalesce(lo_max_group_prevalence, 0)
    )
  ) %>%
  arrange(desc(keep_feature_both_strata), desc(keep_feature), desc(overall_prevalence), desc(max_group_prevalence))

write.csv(species_qc_eo, file.path(table_dir, "species_qc_summary_eo_v29.csv"), row.names = FALSE)
write.csv(species_qc_lo, file.path(table_dir, "species_qc_summary_lo_v29.csv"), row.names = FALSE)
write.csv(species_qc, file.path(table_dir, "species_qc_summary_v29.csv"), row.names = FALSE)
write.csv(
  tibble(
    species_set = c("EO_eligible", "LO_eligible", "union_for_clr", "eligible_in_both_strata"),
    n_species = c(
      length(keep_species_eo),
      length(keep_species_lo),
      length(keep_species_union),
      length(keep_species_both_strata)
    )
  ),
  file.path(table_dir, "species_qc_eligible_set_counts_v29.csv"),
  row.names = FALSE
)

species_clr_mat <- clr_transform_matrix(species_raw_mat[keep_species_union, , drop = FALSE])
species_annotation <- build_species_taxonomy_annotation(rownames(species_clr_mat))

species_eo_clr_mat <- species_clr_mat[intersect(keep_species_eo, rownames(species_clr_mat)), , drop = FALSE]
species_lo_clr_mat <- species_clr_mat[intersect(keep_species_lo, rownames(species_clr_mat)), , drop = FALSE]
species_interaction_clr_mat <- species_clr_mat[intersect(keep_species_both_strata, rownames(species_clr_mat)), , drop = FALSE]

species_eo_disease_fit <- fit_single_coefficient_limma(
  expr_mat = species_eo_clr_mat[, eo_meta$Sample_ID, drop = FALSE],
  sample_meta = eo_meta,
  design_formula = ~ Disease + Age_centered + Sex + BMI,
  coefficient_expression = "DiseaseCRC",
  contrast_name = "EOCRC_vs_EOControl_ageadj",
  model_label = "species_EO_stratified_main",
  test_type = "ebayes",
  lfc_threshold = 0
)

species_lo_disease_fit <- fit_single_coefficient_limma(
  expr_mat = species_lo_clr_mat[, lo_meta$Sample_ID, drop = FALSE],
  sample_meta = lo_meta,
  design_formula = ~ Disease + Age_centered + Sex + BMI,
  coefficient_expression = "DiseaseCRC",
  contrast_name = "LOCRC_vs_LOControl_ageadj",
  model_label = "species_LO_stratified_main",
  test_type = "ebayes",
  lfc_threshold = 0
)

species_pooled_interaction_fit <- fit_single_coefficient_limma(
  expr_mat = species_interaction_clr_mat,
  sample_meta = FUSCC_cl2,
  design_formula = ~ Disease * Age_group + Age_within_group + Sex + BMI,
  coefficient_expression = "DiseaseCRC.Age_groupEO",
  contrast_name = "Disease_age_interaction_pooled_ageadj",
  model_label = "species_pooled_interaction_ageadj",
  test_type = "ebayes",
  lfc_threshold = 0
)

species_panel_assoc_results <- run_panel_association_models(
  expr_mat = species_eo_clr_mat,
  sample_meta_base = eo_meta,
  panel_score_long = panel_score_long,
  target_panels = eo_target_panels,
  design_formula = ~ panel_score + Disease + Age_centered + Sex + BMI,
  coefficient_expression = "panel_score",
  model_label = "species_EO_panel_association"
)

species_lo_panel_assoc_results <- run_panel_association_models(
  expr_mat = species_lo_clr_mat,
  sample_meta_base = lo_meta,
  panel_score_long = panel_score_long,
  target_panels = eo_target_panels,
  design_formula = ~ panel_score + Disease + Age_centered + Sex + BMI,
  coefficient_expression = "panel_score",
  model_label = "species_LO_panel_association"
)

species_panel_sensitivity_results <- run_panel_association_models(
  expr_mat = species_eo_clr_mat,
  sample_meta_base = eo_crc_only_meta,
  panel_score_long = panel_score_long,
  target_panels = eo_target_panels,
  design_formula = ~ panel_score + Age_centered_crc + Sex + BMI + Tumor_Site3 + Stage_collapsed,
  coefficient_expression = "panel_score",
  model_label = "species_EOCRC_panel_sensitivity"
)

species_lo_panel_sensitivity_results <- run_panel_association_models(
  expr_mat = species_lo_clr_mat,
  sample_meta_base = lo_crc_only_meta,
  panel_score_long = panel_score_long,
  target_panels = eo_target_panels,
  design_formula = ~ panel_score + Age_centered_crc + Sex + BMI + Tumor_Site3 + Stage_collapsed,
  coefficient_expression = "panel_score",
  model_label = "species_LOCRC_panel_sensitivity"
)

write.csv(
  species_eo_disease_fit$results %>%
    left_join(species_annotation, by = c("metabolite_identification" = "feature_id")),
  file.path(table_dir, "species_eo_disease_results_v29.csv"),
  row.names = FALSE
)
write.csv(
  species_lo_disease_fit$results %>%
    left_join(species_annotation, by = c("metabolite_identification" = "feature_id")),
  file.path(table_dir, "species_lo_disease_results_v29.csv"),
  row.names = FALSE
)
write.csv(
  species_pooled_interaction_fit$results %>%
    left_join(species_annotation, by = c("metabolite_identification" = "feature_id")),
  file.path(table_dir, "species_pooled_interaction_results_v29.csv"),
  row.names = FALSE
)
write.csv(
  species_panel_assoc_results %>%
    left_join(species_annotation, by = c("metabolite_identification" = "feature_id")),
  file.path(table_dir, "species_panel_association_results_v29.csv"),
  row.names = FALSE
)
write.csv(
  species_lo_panel_assoc_results %>%
    left_join(species_annotation, by = c("metabolite_identification" = "feature_id")),
  file.path(table_dir, "species_lo_panel_association_results_v29.csv"),
  row.names = FALSE
)
write.csv(
  species_panel_sensitivity_results %>%
    left_join(species_annotation, by = c("metabolite_identification" = "feature_id")),
  file.path(table_dir, "species_panel_sensitivity_results_v29.csv"),
  row.names = FALSE
)
write.csv(
  species_lo_panel_sensitivity_results %>%
    left_join(species_annotation, by = c("metabolite_identification" = "feature_id")),
  file.path(table_dir, "species_lo_panel_sensitivity_results_v29.csv"),
  row.names = FALSE
)

species_panel_assoc_annotated <- build_species_association_table(
  panel_assoc_results = species_panel_assoc_results,
  feature_annotation = species_annotation,
  qc_tbl = species_qc
)

write.csv(
  species_panel_assoc_annotated,
  file.path(table_dir, "species_panel_association_annotated_v29.csv"),
  row.names = FALSE
)

write.csv(
  species_panel_assoc_annotated %>%
    filter(adj.P.Val < fdr_cutoff) %>%
    arrange(panel, adj.P.Val, desc(abs(logFC))),
  file.path(table_dir, "species_panel_top_hits_primary_v29.csv"),
  row.names = FALSE
)

top_species_plot <- plot_top_associated_species(
  species_panel_assoc_annotated,
  fdr_threshold = fdr_cutoff,
  top_n = top_species_hits_per_direction
)
if (!is.null(top_species_plot)) {
  save_gg(
    top_species_plot,
    file.path(plot_dir, "species_top_associated_forest_v29.pdf"),
    width = 12.5,
    height = 9.5
  )
}

species_cladogram_primary_data <- build_taxonomic_cladogram_data(
  species_panel_assoc_annotated,
  fdr_threshold = fdr_cutoff
)
if (!is.null(species_cladogram_primary_data)) {
  primary_cladogram_panels <- unique(as.character(species_cladogram_primary_data$species_points$panel))
  primary_cladogram_panels <- primary_cladogram_panels[!is.na(primary_cladogram_panels)]

  if (length(primary_cladogram_panels) == 1) {
    primary_plot <- plot_taxonomic_cladogram(
      species_cladogram_primary_data,
      panel_name = primary_cladogram_panels[[1]]
    )
    if (!is.null(primary_plot)) {
      save_gg(primary_plot, file.path(plot_dir, "species_taxonomic_cladogram_primary_v29.pdf"), width = 12, height = 9)
    }
  } else {
    walk(primary_cladogram_panels, function(panel_name) {
      panel_plot <- plot_taxonomic_cladogram(species_cladogram_primary_data, panel_name = panel_name)
      if (!is.null(panel_plot)) {
        output_name <- paste0(
          "species_taxonomic_cladogram_primary_",
          safe_file_tag(panel_name),
          "_v29.pdf"
        )
        save_gg(panel_plot, file.path(plot_dir, output_name), width = 12, height = 9)
      }
    })
  }
}

species_taxonomic_enrichment_primary <- run_taxonomic_enrichment(
  association_tbl = species_panel_assoc_annotated,
  direction_var = "association_direction",
  significance_layer = "FDR<0.05"
)

write.csv(
  species_taxonomic_enrichment_primary,
  file.path(table_dir, "species_taxonomic_enrichment_primary_v29.csv"),
  row.names = FALSE
)

species_taxonomic_enrichment_primary_plot <- plot_taxonomic_enrichment(
  species_taxonomic_enrichment_primary,
  adj_threshold = fdr_cutoff
)
if (!is.null(species_taxonomic_enrichment_primary_plot)) {
  save_gg(
    species_taxonomic_enrichment_primary_plot,
    file.path(plot_dir, "species_taxonomic_enrichment_primary_v29.pdf"),
    width = 14,
    height = 8
  )
}

species_candidates <- extract_crossomic_candidates(
  panel_assoc_results = species_panel_assoc_results,
  disease_results = species_eo_disease_fit$results,
  interaction_results = species_pooled_interaction_fit$results,
  sensitivity_results = species_panel_sensitivity_results,
  lo_sensitivity_results = species_lo_panel_sensitivity_results,
  lo_disease_results = species_lo_disease_fit$results,
  lo_panel_assoc_results = species_lo_panel_assoc_results,
  feature_annotation = species_annotation
) %>%
  mutate(
    target_panel_label = recode(panel, !!!panel_labels),
    eo_specificity_supported = interaction_fdr < fdr_cutoff & sign(interaction_beta) == sign(disease_beta)
  )

write.csv(species_candidates, file.path(table_dir, "species_candidate_summary_v29.csv"), row.names = FALSE)
write.csv(
  species_candidates %>%
    count(panel, target_panel_label, candidate_tier, eo_specificity_supported),
  file.path(table_dir, "species_candidate_counts_v29.csv"),
  row.names = FALSE
)
write.csv(
  species_candidates %>%
    filter(candidate_tier != "not_shortlisted"),
  file.path(table_dir, "species_candidate_shortlist_v29.csv"),
  row.names = FALSE
)

candidate_species_heatmap <- plot_candidate_species_heatmap(species_candidates)
if (!is.null(candidate_species_heatmap)) {
  save_gg(
    candidate_species_heatmap,
    file.path(plot_dir, "species_candidate_heatmap_v29.pdf"),
    width = 11.5,
    height = 9.5
  )
}

metacyc_raw_mat <- prepare_crossomic_matrix(
  raw_df = metacyc_list[["A30"]],
  feature_col = names(metacyc_list[["A30"]])[1],
  sample_map = omics_sample_map
)
metacyc_raw_mat <- metacyc_raw_mat[, FUSCC_cl2$Sample_ID, drop = FALSE]
metacyc_anchor_raw <- extract_anchor_matrix(
  feature_mat = metacyc_raw_mat,
  anchor_ids = manuscript_metacyc_anchor_ids,
  id_type = "pathway"
)

ec_raw_mat <- prepare_crossomic_matrix(
  raw_df = ec_list[["A30"]],
  feature_col = names(ec_list[["A30"]])[1],
  sample_map = omics_sample_map
)
ec_raw_mat <- ec_raw_mat[, FUSCC_cl2$Sample_ID, drop = FALSE]
ec_anchor_raw <- extract_anchor_matrix(
  feature_mat = ec_raw_mat,
  anchor_ids = manuscript_ec_anchor_ids,
  id_type = "ec"
)

metacyc_anchor_qc <- summarise_crossomic_qc(
  feature_mat = metacyc_anchor_raw,
  sample_meta = FUSCC_cl2,
  subset_mask = FUSCC_cl2$Age_group == "EO",
  group_var = "Disease",
  min_prevalence = function_prevalence_min,
  min_group_prevalence = function_group_prevalence_min
)
ec_anchor_qc <- summarise_crossomic_qc(
  feature_mat = ec_anchor_raw,
  sample_meta = FUSCC_cl2,
  subset_mask = FUSCC_cl2$Age_group == "EO",
  group_var = "Disease",
  min_prevalence = function_prevalence_min,
  min_group_prevalence = function_group_prevalence_min
)

write.csv(metacyc_anchor_qc, file.path(table_dir, "metacyc_anchor_qc_v29.csv"), row.names = FALSE)
write.csv(ec_anchor_qc, file.path(table_dir, "ec_anchor_qc_v29.csv"), row.names = FALSE)

metacyc_anchor_clr <- clr_transform_matrix(
  metacyc_anchor_raw[metacyc_anchor_qc$feature_id[metacyc_anchor_qc$keep_feature], , drop = FALSE]
)
ec_anchor_clr <- clr_transform_matrix(
  ec_anchor_raw[ec_anchor_qc$feature_id[ec_anchor_qc$keep_feature], , drop = FALSE]
)

metacyc_anchor_disease_fit <- fit_single_coefficient_limma(
  expr_mat = metacyc_anchor_clr[, eo_meta$Sample_ID, drop = FALSE],
  sample_meta = eo_meta,
  design_formula = ~ Disease + Age_centered + Sex + BMI,
  coefficient_expression = "DiseaseCRC",
  contrast_name = "EOCRC_vs_EOControl_ageadj",
  model_label = "metacyc_anchor_EO_stratified_main",
  test_type = "ebayes",
  lfc_threshold = 0
)
ec_anchor_disease_fit <- fit_single_coefficient_limma(
  expr_mat = ec_anchor_clr[, eo_meta$Sample_ID, drop = FALSE],
  sample_meta = eo_meta,
  design_formula = ~ Disease + Age_centered + Sex + BMI,
  coefficient_expression = "DiseaseCRC",
  contrast_name = "EOCRC_vs_EOControl_ageadj",
  model_label = "ec_anchor_EO_stratified_main",
  test_type = "ebayes",
  lfc_threshold = 0
)

metacyc_anchor_panel_results <- run_panel_association_models(
  expr_mat = metacyc_anchor_clr,
  sample_meta_base = eo_meta,
  panel_score_long = panel_score_long,
  target_panels = eo_target_panels,
  design_formula = ~ panel_score + Disease + Age_centered + Sex + BMI,
  coefficient_expression = "panel_score",
  model_label = "metacyc_anchor_EO_panel_association"
)
ec_anchor_panel_results <- run_panel_association_models(
  expr_mat = ec_anchor_clr,
  sample_meta_base = eo_meta,
  panel_score_long = panel_score_long,
  target_panels = eo_target_panels,
  design_formula = ~ panel_score + Disease + Age_centered + Sex + BMI,
  coefficient_expression = "panel_score",
  model_label = "ec_anchor_EO_panel_association"
)

metacyc_anchor_results_combined <- bind_rows(
  metacyc_anchor_disease_fit$results %>% mutate(anchor_type = "MetaCyc_pathway"),
  metacyc_anchor_panel_results %>% mutate(anchor_type = "MetaCyc_pathway")
)

ec_anchor_results_combined <- bind_rows(
  ec_anchor_disease_fit$results %>% mutate(anchor_type = "EC"),
  ec_anchor_panel_results %>% mutate(anchor_type = "EC")
)

write.csv(
  metacyc_anchor_results_combined,
  file.path(table_dir, "metacyc_anchor_results_v29.csv"),
  row.names = FALSE
)
write.csv(
  ec_anchor_results_combined,
  file.path(table_dir, "ec_anchor_results_v29.csv"),
  row.names = FALSE
)

candidate_species_ids <- species_candidates %>%
  filter(candidate_tier != "not_shortlisted") %>%
  pull(feature_id) %>%
  unique()

if (length(candidate_species_ids) > 0) {
  candidate_anchor_support <- bind_rows(
    map_dfr(eo_target_panels, function(panel_name) {
      panel_species_ids <- species_candidates %>%
        filter(panel == panel_name, candidate_tier != "not_shortlisted") %>%
        pull(feature_id) %>%
        unique()

      if (length(panel_species_ids) == 0 || nrow(metacyc_anchor_clr) == 0) {
        return(tibble())
      }

      map_dfr(panel_species_ids, function(species_id) {
        map_dfr(rownames(metacyc_anchor_clr), function(anchor_id) {
          lm_df <- eo_meta %>%
            mutate(
              species_value = species_clr_mat[species_id, Sample_ID],
              anchor_value = metacyc_anchor_clr[anchor_id, Sample_ID]
            )

          tidy(lm(anchor_value ~ species_value + Disease + Age_centered + Sex + BMI, data = lm_df)) %>%
            filter(term == "species_value") %>%
            transmute(
              panel = panel_name,
              target_panel_label = unname(panel_labels[panel_name]),
              anchor_type = "MetaCyc_pathway",
              anchor_id = anchor_id,
              feature_id = species_id,
              estimate,
              statistic,
              p.value
            )
        })
      })
    }),
    map_dfr(eo_target_panels, function(panel_name) {
      panel_species_ids <- species_candidates %>%
        filter(panel == panel_name, candidate_tier != "not_shortlisted") %>%
        pull(feature_id) %>%
        unique()

      if (length(panel_species_ids) == 0 || nrow(ec_anchor_clr) == 0) {
        return(tibble())
      }

      map_dfr(panel_species_ids, function(species_id) {
        map_dfr(rownames(ec_anchor_clr), function(anchor_id) {
          lm_df <- eo_meta %>%
            mutate(
              species_value = species_clr_mat[species_id, Sample_ID],
              anchor_value = ec_anchor_clr[anchor_id, Sample_ID]
            )

          tidy(lm(anchor_value ~ species_value + Disease + Age_centered + Sex + BMI, data = lm_df)) %>%
            filter(term == "species_value") %>%
            transmute(
              panel = panel_name,
              target_panel_label = unname(panel_labels[panel_name]),
              anchor_type = "EC",
              anchor_id = anchor_id,
              feature_id = species_id,
              estimate,
              statistic,
              p.value
            )
        })
      })
    })
  ) %>%
    group_by(panel, anchor_type) %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "BH")) %>%
    ungroup() %>%
    left_join(species_annotation, by = "feature_id")

  write.csv(
    candidate_anchor_support,
    file.path(table_dir, "candidate_species_anchor_support_v29.csv"),
    row.names = FALSE
  )
} else {
  write.csv(
    tibble(),
    file.path(table_dir, "candidate_species_anchor_support_v29.csv"),
    row.names = FALSE
  )
}

main_figure_support_strip <- build_species_support_strip_data(species_candidates = species_candidates)

write.csv(
  main_figure_support_strip %>% select(-text_color, -max_abs),
  file.path(table_dir, "main_figure_support_strip_species_v29.csv"),
  row.names = FALSE
)

support_strip_plot <- plot_species_support_strip(main_figure_support_strip)
save_gg(
  support_strip_plot,
  file.path(plot_dir, "main_figure_support_strip_species_v29.pdf"),
  width =8,
  height = 5
)
