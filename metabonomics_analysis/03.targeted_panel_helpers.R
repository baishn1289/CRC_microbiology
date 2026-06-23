
build_targeted_panel_membership <- function(feature_names) {
  feature_tbl <- tibble(metabolite_identification = feature_names)

  canonical_map <- tribble(
    ~metabolite_identification, ~canonical_metabolite,
    "Arginine", "Arginine",
    "Citrulline", "Citrulline",
    "L-Citrulline", "Citrulline",
    "D-Ornithine", "Ornithine_related",
    "DL-Ornithine", "Ornithine_related",
    "N-.alpha.-Acetyl-L-ornithine", "N-Acetylornithine",
    "N-Acetylornithine", "N-Acetylornithine",
    "N-Acetylputrescine", "N-Acetylputrescine",
    "5-Aminopentanoic acid", "5-Aminopentanoate",
    "5-Aminovaleric acid", "5-Aminopentanoate",
    "5-Aminopentanal", "5-Aminopentanal",
    "Indole-3-acetaldehyde", "Indole-3-acetaldehyde",
    "Indoleacetaldehyde", "Indole-3-acetaldehyde",
    "Indole-3-carboxyaldehyde", "Indole-3-carboxyaldehyde",
    "3-Formylindole", "Indole-3-carboxyaldehyde",
    "N-Carbamoylputrescine", "N-Carbamoylputrescine",
    "N-Acetylcadaverine", "N-Acetylcadaverine",
    "Lysine", "Lysine",
    "Tryptophan", "Tryptophan",
    "Urocanic acid", "Urocanic acid",
    "Gamma-Aminobutyric acid", "Gamma-Aminobutyric acid",
    "Histamine", "Histamine",
    "N-Acetylhistamine", "N-Acetylhistamine",
    "Tyramine", "Tyramine",
    "Tyramine-O-sulfate", "Tyramine-O-sulfate",
    "N-Methyltyramine", "N-Methyltyramine",
    "Indole", "Indole",
    "Indoleacetic acid", "Indole-3-acetic acid",
    "3-Indoleacetic acid", "Indole-3-acetic acid",
    "Indolelactic acid", "Indole-3-lactate",
    "Indole-3-acetamide", "Indole-3-acetamide",
    "Indole-3-carboxylic acid", "Indole-3-carboxylic acid",
    "1H-Indole-3-carboxylic acid", "Indole-3-carboxylic acid",
    "p-Cresol", "p-Cresol",
    "p-Cresol sulfate", "p-Cresol sulfate",
    "Phenol", "Phenol",
    "Phenol sulphate", "Phenol sulphate",
    "Indoxyl sulfate", "Indoxyl sulfate",
    "Phenylacetic acid", "Phenylacetic acid",
    "Phenylacetyl-L-glutamine", "Phenylacetyl-L-glutamine",
    "N-Phenylacetylglutamic acid", "N-Phenylacetylglutamic acid",
    "4-Hydroxyphenylacetylglutamine", "4-Hydroxyphenylacetylglutamine"
  )

  panel_catalog <- tribble(
    ~panel, ~canonical_metabolite,
    "arginine_ornithine_routing_panel", "Arginine",
    "arginine_ornithine_routing_panel", "Citrulline",
    "arginine_ornithine_routing_panel", "Ornithine_related",
    "arginine_ornithine_routing_panel", "N-Acetylornithine",
    "amino_acid_putrefaction_core_panel", "5-Aminopentanoate",
    "amino_acid_putrefaction_core_panel", "N-Acetylputrescine",
    "amino_acid_putrefaction_core_panel", "N-Carbamoylputrescine",
    "amino_acid_putrefaction_auxiliary_panel", "N-Acetylcadaverine",
    "amino_acid_putrefaction_auxiliary_panel", "5-Aminopentanal",
    "amino_acid_putrefaction_auxiliary_panel", "Gamma-Aminobutyric acid",
    "amino_acid_putrefaction_auxiliary_panel", "Urocanic acid",
    "tryptophan_indole_branch_panel", "Indole",
    "tryptophan_indole_branch_panel", "Indole-3-acetaldehyde",
    "tryptophan_indole_branch_panel", "Indole-3-carboxyaldehyde",
    "tryptophan_indole_branch_panel", "Indole-3-acetic acid",
    "tryptophan_indole_branch_panel", "Indole-3-lactate",
    "polyamine_lysine_putrefaction_panel", "N-Acetylputrescine",
    "polyamine_lysine_putrefaction_panel", "N-Carbamoylputrescine",
    "polyamine_lysine_putrefaction_panel", "N-Acetylcadaverine",
    "polyamine_lysine_putrefaction_panel", "5-Aminopentanoate",
    "biogenic_amine_panel", "Histamine",
    "biogenic_amine_panel", "N-Acetylhistamine",
    "biogenic_amine_panel", "Tyramine",
    "biogenic_amine_panel", "Tyramine-O-sulfate",
    "biogenic_amine_panel", "N-Methyltyramine",
    "indole_catabolism_panel", "Indole",
    "indole_catabolism_panel", "Indole-3-acetaldehyde",
    "indole_catabolism_panel", "Indole-3-acetamide",
    "indole_catabolism_panel", "Indole-3-carboxyaldehyde",
    "indole_catabolism_panel", "Indole-3-acetic acid",
    "indole_catabolism_panel", "Indole-3-lactate",
    "indole_catabolism_panel", "Indole-3-carboxylic acid",
    "aromatic_host_cometabolite_panel", "Phenylacetic acid",
    "aromatic_host_cometabolite_panel", "Phenylacetyl-L-glutamine",
    "aromatic_host_cometabolite_panel", "N-Phenylacetylglutamic acid",
    "aromatic_host_cometabolite_panel", "4-Hydroxyphenylacetylglutamine",
    "aromatic_host_cometabolite_panel", "Indoxyl sulfate",
    "injurious_biogenic_amine_subpanel", "Histamine",
    "injurious_biogenic_amine_subpanel", "Tyramine",
    "injurious_biogenic_amine_subpanel", "N-Methyltyramine",
    "injurious_parent_phenol_subpanel", "p-Cresol",
    "injurious_parent_phenol_subpanel", "Phenol",
    "injurious_host_processed_aromatic_toxin_subpanel", "p-Cresol sulfate",
    "injurious_host_processed_aromatic_toxin_subpanel", "Phenol sulphate",
    "injurious_host_processed_aromatic_toxin_subpanel", "Indoxyl sulfate"
  )

  feature_tbl %>%
    inner_join(canonical_map, by = "metabolite_identification") %>%
    inner_join(panel_catalog, by = "canonical_metabolite", relationship = "many-to-many") %>%
    distinct(metabolite_identification, canonical_metabolite, panel, .keep_all = TRUE) %>%
    mutate(
      panel = factor(panel, levels = panel_levels),
      panel_label = recode(panel, !!!panel_labels),
      panel_tier = recode(panel, !!!panel_tier_map),
      panel_tier_label = recode(panel_tier, !!!panel_tier_labels)
    ) %>%
    arrange(panel, canonical_metabolite, metabolite_identification)
}

build_focus_metabolite_review <- function(feature_names) {
  tibble(
    metabolite_identification = c(
      "5-Aminopentanoic acid",
      "5-Aminovaleric acid",
      "N-Acetylputrescine",
      "N-Carbamoylputrescine",
      "N-Acetylcadaverine",
      "5-Aminopentanal",
      "Gamma-Aminobutyric acid",
      "Urocanic acid",
      "Indole",
      "Indole-3-acetaldehyde",
      "Indole-3-carboxyaldehyde",
      "3-Formylindole",
      "Indoleacetic acid",
      "3-Indoleacetic acid",
      "Indolelactic acid",
      "Indoxyl sulfate",
      "4-Hydroxyphenylacetylglutamine",
      "Indole-3-carbinol",
      "5-Hydroxyindole-3-acetic acid",
      "L-Histidine",
      "L-Phenylalanine"
    ),
    final_decision = c(
      "keep_as_canonical_5_aminopentanoate",
      "collapse_into_5_aminopentanoate",
      "keep_putrefaction_core_tier1",
      "keep_putrefaction_core_tier1",
      "keep_putrefaction_auxiliary_tier2",
      "keep_putrefaction_auxiliary_tier2",
      "keep_putrefaction_auxiliary_tier2",
      "keep_putrefaction_auxiliary_tier2",
      "keep_indole_branch_supplementary",
      "keep_indole_branch_supplementary",
      "keep_indole_branch_supplementary",
      "collapse_into_indole_3_carboxyaldehyde",
      "keep_indole_branch_supplementary",
      "collapse_into_indole_3_acetic_acid",
      "keep_indole_branch_supplementary",
      "keep_extended_only",
      "keep_extended_only",
      "exclude_from_curated_panels",
      "exclude_from_curated_panels",
      "exclude_from_curated_panels",
      "exclude_from_curated_panels"
    ),
    representative_panel = c(
      "amino_acid_putrefaction_core_panel",
      "amino_acid_putrefaction_core_panel",
      "amino_acid_putrefaction_core_panel",
      "amino_acid_putrefaction_core_panel",
      "amino_acid_putrefaction_auxiliary_panel",
      "amino_acid_putrefaction_auxiliary_panel",
      "amino_acid_putrefaction_auxiliary_panel",
      "amino_acid_putrefaction_auxiliary_panel",
      "tryptophan_indole_branch_panel",
      "tryptophan_indole_branch_panel",
      "tryptophan_indole_branch_panel",
      "tryptophan_indole_branch_panel",
      "tryptophan_indole_branch_panel",
      "tryptophan_indole_branch_panel",
      "tryptophan_indole_branch_panel",
      "injurious_host_processed_aromatic_toxin_subpanel",
      "aromatic_host_cometabolite_panel",
      NA,
      NA,
      NA,
      NA
    ),
    biochemical_comment = c(
      "Primary lysine/cadaverine branch marker aligned with the amino-acid putrefaction mainline; retained as a first-tier core feature.",
      "Synonym-level duplicate of 5-aminopentanoic acid; retained in raw results but collapsed for panel scoring and heatmaps.",
      "Detectable putrescine-related spillover marker aligned with the arginine-ornithine-polyamine axis; retained as a first-tier core feature.",
      "Putrescine/agmatine-related intermediate retained as a same-axis support marker in the strict putrefaction core.",
      "Cadaverine itself is not detected, so this acetylated cadaverine signal is retained only as a second-tier auxiliary lysine-branch marker.",
      "Near-downstream lysine/cadaverine branch aldehyde retained only as a second-tier auxiliary marker.",
      "Broad amino-acid decarboxylation/putrefaction-adjacent readout retained only as a second-tier auxiliary marker because of lower pathway specificity.",
      "Direct product of histidine ammonia-lyase (EC 4.3.1.3); retained as an auxiliary histidine-derived ammonia-release marker rather than a core feature.",
      "Direct product of tryptophan indole-lyase (EC 4.1.99.1); retained only in the supplementary tryptophan-indole branch.",
      "Tryptophan-derived indole aldehyde marker retained only in the supplementary indole branch rather than the main amino-acid putrefaction core.",
      "Tryptophan-derived indole aldehyde marker retained only in the supplementary indole branch rather than the main amino-acid putrefaction core.",
      "Synonym-level duplicate of indole-3-carboxyaldehyde rather than an independent metabolite signal.",
      "Tryptophan-derived indole acid marker retained only in the supplementary indole branch.",
      "Synonym-level duplicate of indoleacetic acid and collapsed into indole-3-acetic acid.",
      "Tryptophan-derived indole lactate marker retained only in the supplementary indole branch.",
      "Strong gut microbiota-host co-metabolite toxin signal, but it extends beyond the exact manuscript EC/pathway anchors, so it is kept only in the extended panels.",
      "Aromatic host-microbe conjugate consistent with broadened aromatic fermentation support, but not specific enough for strict manuscript-aligned validation.",
      "Likely dominated by dietary cruciferous origin rather than a robust fecal putrefaction readout; excluded to avoid overclaiming.",
      "Predominantly a host serotonin catabolite rather than a clean microbial putrefaction marker; excluded.",
      "Represents broad amino-acid pool status rather than the arginine-ornithine precursor axis described in the manuscript; excluded from curated panels.",
      "Represents a broad amino-acid reservoir rather than a direct nitrogen-flux-diversion anchor; excluded from curated panels."
    )
  ) %>%
    mutate(
      present_in_dataset = metabolite_identification %in% feature_names
    )
}

build_panel_score_matrix <- function(feature_mat, panel_membership) {
  z_mat <- z_score_features(feature_mat)

  canonical_feature_map <- panel_membership %>%
    distinct(metabolite_identification, canonical_metabolite)

  canonical_feature_list <- split(
    canonical_feature_map$metabolite_identification,
    canonical_feature_map$canonical_metabolite
  )

  canonical_mat <- vapply(
    canonical_feature_list,
    function(mets) safe_colmeans_available(z_mat[unique(mets), , drop = FALSE]),
    numeric(ncol(z_mat))
  )
  canonical_mat <- t(canonical_mat)
  rownames(canonical_mat) <- names(canonical_feature_list)
  colnames(canonical_mat) <- colnames(z_mat)

  panel_list <- split(panel_membership$canonical_metabolite, panel_membership$panel)
  panel_list <- panel_list[lengths(panel_list) > 0]

  panel_mat <- vapply(
    panel_list,
    function(mets) safe_colmeans_available(canonical_mat[unique(mets), , drop = FALSE]),
    numeric(ncol(canonical_mat))
  )
  panel_mat <- t(panel_mat)
  rownames(panel_mat) <- names(panel_list)
  colnames(panel_mat) <- colnames(canonical_mat)
  panel_mat
}

build_targeted_canonical_long <- function(long_df, panel_membership) {
  canonical_map <- panel_membership %>%
    distinct(metabolite_identification, canonical_metabolite)

  long_df %>%
    inner_join(canonical_map, by = "metabolite_identification") %>%
    group_by(canonical_metabolite, Sample_ID, Group4, Disease, Age_group, Sex, BMI) %>%
    summarise(log2_intensity = mean(log2_intensity, na.rm = TRUE), .groups = "drop")
}

build_panel_score_coverage_table <- function(raw_long_df, panel_membership) {
  canonical_map <- panel_membership %>%
    distinct(panel, panel_label, panel_tier, panel_tier_label, metabolite_identification, canonical_metabolite)

  raw_long_df %>%
    inner_join(canonical_map, by = "metabolite_identification", relationship = "many-to-many") %>%
    mutate(has_value = !is.na(intensity) & intensity > 0) %>%
    group_by(Sample_ID, Group4, Disease, Age_group, Sex, BMI, panel, panel_label, panel_tier, panel_tier_label, canonical_metabolite) %>%
    summarise(has_value = any(has_value), .groups = "drop") %>%
    group_by(Sample_ID, Group4, Disease, Age_group, Sex, BMI, panel, panel_label, panel_tier, panel_tier_label) %>%
    summarise(
      n_observed = sum(has_value),
      n_total = dplyr::n(),
      observed_fraction = n_observed / n_total,
      .groups = "drop"
    )
}

summarise_panel_score_coverage <- function(panel_coverage_tbl) {
  panel_coverage_tbl %>%
    group_by(panel, panel_label, panel_tier, panel_tier_label, Group4, Disease, Age_group) %>%
    summarise(
      n_samples = dplyr::n(),
      mean_observed_fraction = mean(observed_fraction, na.rm = TRUE),
      median_observed_fraction = median(observed_fraction, na.rm = TRUE),
      min_observed_fraction = min(observed_fraction, na.rm = TRUE),
      complete_coverage_rate = mean(observed_fraction >= 1, na.rm = TRUE),
      .groups = "drop"
    )
}

# final canonicalized heatmap implementation
plot_targeted_heatmap <- function(long_df, significant_targeted_panel_df, id_col = "canonical_metabolite") {
  if (nrow(significant_targeted_panel_df) == 0) {
    return(NULL)
  }

  order_tbl <- significant_targeted_panel_df %>%
    distinct(panel, panel_label, panel_tier, panel_tier_label, .data[[id_col]]) %>%
    arrange(panel, .data[[id_col]]) %>%
    mutate(panel_feature = paste(panel, .data[[id_col]], sep = "||"))

  heatmap_df <- long_df %>%
    filter(.data[[id_col]] %in% order_tbl[[id_col]]) %>%
    group_by(.data[[id_col]], Group4) %>%
    summarise(mean_log2 = mean(log2_intensity, na.rm = TRUE), .groups = "drop") %>%
    group_by(.data[[id_col]]) %>%
    mutate(group_mean_z = as.numeric(scale(mean_log2))) %>%
    ungroup() %>%
    inner_join(order_tbl, by = id_col, relationship = "many-to-many") %>%
    mutate(
      Group4 = factor(Group4, levels = group_levels),
      panel_feature = factor(panel_feature, levels = rev(order_tbl$panel_feature))
    )

  ggplot(heatmap_df, aes(x = Group4, y = panel_feature, fill = group_mean_z)) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(panel_label ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(
      low = "#1D4E89",
      mid = "white",
      high = "#B22222",
      midpoint = 0
    ) +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) +
    labs(
      title = "Group-mean z-scores for significant curated metabolites",
      x = NULL,
      y = NULL,
      fill = "z-score"
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.text.y = element_text(angle = 0)
    )
}

plot_panel_member_heatmap <- function(long_df, panel_lookup, feature_results, id_col = "canonical_metabolite") {
  order_tbl <- panel_lookup %>%
    distinct(panel, panel_label, panel_tier, panel_tier_label, .data[[id_col]]) %>%
    arrange(panel, .data[[id_col]]) %>%
    mutate(panel_feature = paste(panel, .data[[id_col]], sep = "||"))

  if (nrow(order_tbl) == 0) {
    return(NULL)
  }

  sig_col <- if ("significant_fdr" %in% names(feature_results)) {
    "significant_fdr"
  } else if ("any_significant_fdr" %in% names(feature_results)) {
    "any_significant_fdr"
  } else {
    stop("feature_results must contain either significant_fdr or any_significant_fdr.")
  }

  feature_sig <- feature_results %>%
    group_by(.data[[id_col]]) %>%
    summarise(any_fdr_sig = any(.data[[sig_col]], na.rm = TRUE), .groups = "drop")

  heatmap_df <- long_df %>%
    filter(.data[[id_col]] %in% order_tbl[[id_col]]) %>%
    group_by(.data[[id_col]], Group4) %>%
    summarise(mean_log2 = mean(log2_intensity, na.rm = TRUE), .groups = "drop") %>%
    group_by(.data[[id_col]]) %>%
    mutate(group_mean_z = as.numeric(scale(mean_log2))) %>%
    ungroup() %>%
    inner_join(order_tbl, by = id_col, relationship = "many-to-many") %>%
    left_join(feature_sig, by = id_col) %>%
    mutate(
      Group4 = factor(Group4, levels = group_levels),
      panel_feature = factor(panel_feature, levels = rev(order_tbl$panel_feature)),
      sig_flag = ifelse(any_fdr_sig, "FDR<0.05 in >=1 contrast", "Not FDR-significant")
    )

  ggplot(heatmap_df, aes(x = Group4, y = panel_feature, fill = group_mean_z#, alpha = sig_flag
                         )) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(
       panel_label~.,
      scales = "free_y",
      space = "free_y",
      labeller = as_labeller(function(x) {
        stringr::str_replace(x, "^Extended: ", "Extended panel\n")
      })
    )+
    # facet_grid(panel_label ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(
      low = "#1D4E89",
      mid = "white",
      high = "#B22222",
      midpoint = 0,
      guide = guide_colorbar(
        barwidth = unit(6, "cm"),
        barheight = unit(0.45, "cm"),
        title.position = "left"
        # ticks.colour = "black",
        # frame.colour = "black"
      )
    ) +
    # scale_alpha_manual(values = c("FDR<0.05 in >=1 contrast" = 1, "Not FDR-significant" = 0.55)) +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) +
    labs(
      title = "Group-mean z-scores for all curated metabolites",
      subtitle = "Canonicalized metabolite labels remove synonym-level double counting",
      x = NULL,
      y = NULL,
      fill = "z-score"#,
      # alpha = NULL
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size =10,angle = 30, hjust = 1),
      axis.text.y = element_text(size =10),
      strip.text.y = element_text(angle = 0,size =10),
      title = element_text(size =12),
      plot.subtitle = element_text(size =10),
      # strip.background = element_blank(),
      legend.position = "bottom"
    )
}
