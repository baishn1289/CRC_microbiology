
run_panel_association_models <- function(
    expr_mat,
    sample_meta_base,
    panel_score_long,
    target_panels,
    design_formula,
    coefficient_expression,
    model_label) {
  map_dfr(target_panels, function(panel_name) {
    panel_meta <- sample_meta_base %>%
      left_join(
        panel_score_long %>%
          filter(panel == panel_name) %>%
          select(Sample_ID, panel_score),
        by = "Sample_ID"
      )

    fit_single_coefficient_limma(
      expr_mat = expr_mat[, panel_meta$Sample_ID, drop = FALSE],
      sample_meta = panel_meta,
      design_formula = design_formula,
      coefficient_expression = coefficient_expression,
      contrast_name = paste0(panel_name, "_association"),
      model_label = model_label,
      test_type = "ebayes",
      lfc_threshold = 0
    )$results %>%
      mutate(
        target_panel = panel_name,
        target_panel_label = unname(panel_labels[panel_name])
      )
  })
}

build_panel_score_annotations <- function(panel_score_long, panel_results) {
  contrast_positions <- tibble(
    contrast = contrast_display_order,
    x_start = c(1, 3, 2),
    x_end = c(2, 4, 4)
  )

  panel_ranges <- panel_score_long %>%
    group_by(panel, panel_label) %>%
    summarise(
      y_max = max(panel_score, na.rm = TRUE),
      y_min = min(panel_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y_span = pmax(y_max - y_min, 0.35))

  panel_results %>%
    filter(contrast %in% contrast_display_order) %>%
    left_join(contrast_positions, by = "contrast") %>%
    left_join(panel_ranges, by = c("panel", "panel_label")) %>%
    mutate(
      contrast = factor(contrast, levels = contrast_display_order),
      rank_within_panel = as.numeric(contrast),
      y_position = y_max + y_span * (0.12 + 0.12 * rank_within_panel),
      x_mid = (x_start + x_end) / 2,
      sig_label = significance_label(adj.P.Val),
      beta_label = if_else(is.na(logFC), "NA", sprintf("%+.2f", logFC)),
      annot_label = paste0("\u03b2=", beta_label, "\n", sig_label)
    )
}

plot_panel_scores <- function(panel_score_long,
                              panel_annotations,
                              label_style = c("no_sig_only", "beta_sig"),
                              ncol = 4,
                              beta_offset = 0.08,
                              sig_offset  = 0.045,
                              bracket_depth = 0.03,
                              eolocr_lift = 0.12,
                              text_size = 16
                          
) {
  label_style <- match.arg(label_style)
  
  group_palette <- c(
    "LOControl" = "#B8C7DA",
    "LOCRC"     = "#3F5F8F",
    "EOControl" = "#E8B8B0",
    "EOCRC"     = "#C83E4D"
  )
  
  panel_score_long <- panel_score_long %>%
    mutate(
      Group4 = factor(
        Group4,
        levels = c("LOControl", "LOCRC", "EOControl", "EOCRC")
      )
    )
  
  ann <- panel_annotations
  
  if (!"beta_label" %in% names(ann)) {
    ann$beta_label <- if ("annot_label" %in% names(ann)) {
      sub("\\n.*$", "", ann$annot_label)
    } else {
      NA_character_
    }
  }
  
  if (!"sig_label_plot" %in% names(ann)) {
    if ("sig_label" %in% names(ann)) {
      ann$sig_label_plot <- ann$sig_label
    } else if ("annot_label" %in% names(ann)) {
      ann$sig_label_plot <- ifelse(
        grepl("\\n", ann$annot_label),
        sub("^.*\\n", "", ann$annot_label),
        NA_character_
      )
    } else {
      ann$sig_label_plot <- NA_character_
    }
  }
  
  ann <- ann %>%
    mutate(
      is_LOCRC_vs_EOCRC = 
        (pmin(x_start, x_end) == 2 & pmax(x_start, x_end) == 4)
    )
  # Lift the LOCRC vs EOCRC bracket when beta/significance labels would overlap.
  if (identical(label_style, "beta_sig")) {
    ann <- ann %>%
      mutate(
        y_position = ifelse(
          is_LOCRC_vs_EOCRC,
          y_position + eolocr_lift * y_span,
          y_position
        )
      )
  }
  # Place beta and significance labels around the bracket tick.
  ann$y_tick <- ann$y_position - bracket_depth * ann$y_span
  ann$y_beta <- ann$y_position + beta_offset * ann$y_span
  ann$y_sig  <- ann$y_position - sig_offset  * ann$y_span
  
  p <- ggplot(panel_score_long, aes(x = Group4, y = panel_score, fill = Group4)) +
    geom_jitter(
      aes(color = Group4),
      width = 0.15,
      alpha = 0.8,
      size = 0.55,
      show.legend = FALSE
    ) + geom_boxplot(
      outlier.shape = NA,
      width = 0.65,
      alpha = 0.88,
      linewidth = 0.35,
      color = "grey20"
    ) +

    scale_color_manual(
      values = group_palette,
      breaks = c("LOControl", "LOCRC", "EOControl", "EOCRC"),
      drop = FALSE
    ) +
    # facet_wrap(~ panel_label, scales = "free_y", ncol = ncol) +
    facet_wrap(
      ~ panel_label,
      scales = "free_y",
      ncol = ncol,
      labeller = as_labeller(function(x) {
        stringr::str_replace(x, "^Extended: ", "Extended panel\n")
      })
    )+
    scale_fill_manual(
      values = group_palette,
      breaks = c("LOControl", "LOCRC", "EOControl", "EOCRC"),
      drop = FALSE
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
    labs(
      x = NULL,
      y = "Mean feature z-score"
    ) +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(hjust = 1,angle = 25, size = text_size),
      axis.text.y = element_text(hjust = 0.5, size = text_size),
      panel.grid = element_blank(),
      plot.margin = margin(5.5, 12, 5.5, 5.5),
      strip.background = element_blank(),
      strip.text = element_text(size = text_size),
      axis.title.y = element_text(size = text_size)
    )
  
  if (identical(label_style, "beta_sig")) {
    p <- p +
      geom_segment(
        data = ann,
        aes(x = x_start, xend = x_end, y = y_position, yend = y_position),
        inherit.aes = FALSE,
        linewidth = 0.45
      ) +
      geom_segment(
        data = ann,
        aes(x = x_start, xend = x_start, y = y_tick, yend = y_position),
        inherit.aes = FALSE,
        linewidth = 0.45
      ) +
      geom_segment(
        data = ann,
        aes(x = x_end, xend = x_end, y = y_tick, yend = y_position),
        inherit.aes = FALSE,
        linewidth = 0.45
      ) +
      geom_text(
        data = ann,
        aes(x = x_mid, y = y_beta, label = beta_label),
        inherit.aes = FALSE,
        size = 4,
        lineheight = 0.92
      ) +
      geom_text(
        data = ann,
        aes(x = x_mid, y = y_sig, label = sig_label_plot),
        inherit.aes = FALSE,
        size = 4,
        lineheight = 0.92
      )
  }
  
  p
}

plot_focus_boxplots <- function(long_df, feature_order) {
  if (length(feature_order) == 0) {
    return(NULL)
  }

  plot_df <- long_df %>%
    filter(metabolite_identification %in% feature_order) %>%
    mutate(
      Group4 = factor(Group4, levels = group_levels),
      metabolite_identification = factor(metabolite_identification, levels = feature_order)
    )

  ggplot(plot_df, aes(x = Group4, y = log2_intensity, fill = Group4)) +
    geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.85) +
    geom_jitter(width = 0.15, alpha = 0.3, size = 0.6, color = "black") +
    facet_wrap(~ metabolite_identification, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = group_palette) +
    labs(
      title = "Boxplots for selected nitrogen / putrefaction-related metabolites",
      x = NULL,
      y = "log2 intensity"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

extract_lm_term <- function(fit, term_name, contrast_name) {
  if (!(term_name %in% names(coef(fit)))) {
    stop("Requested term was not found in the fitted model: ", term_name)
  }

  estimate <- unname(coef(fit)[term_name])
  std_error <- sqrt(diag(vcov(fit)))[term_name]
  df_resid <- df.residual(fit)
  statistic <- estimate / std_error
  p_value <- 2 * pt(-abs(statistic), df = df_resid)
  crit <- qt(0.975, df = df_resid)

  tibble(
    contrast = contrast_name,
    estimate = estimate,
    std.error = std_error,
    statistic = statistic,
    p.value = p_value,
    conf.low = estimate - crit * std_error,
    conf.high = estimate + crit * std_error
  )
}

compare_independent_lm_terms <- function(fit_a, fit_b, term_name, contrast_name) {
  if (!(term_name %in% names(coef(fit_a))) || !(term_name %in% names(coef(fit_b)))) {
    stop("Requested term was not found in one of the stratified fitted models: ", term_name)
  }

  estimate_a <- unname(coef(fit_a)[term_name])
  estimate_b <- unname(coef(fit_b)[term_name])
  se_a <- sqrt(diag(vcov(fit_a)))[term_name]
  se_b <- sqrt(diag(vcov(fit_b)))[term_name]
  var_a <- se_a^2
  var_b <- se_b^2
  std_error <- sqrt(var_a + var_b)

  df_a <- df.residual(fit_a)
  df_b <- df.residual(fit_b)
  df_resid <- (var_a + var_b)^2 / ((var_a^2 / df_a) + (var_b^2 / df_b))
  statistic <- (estimate_a - estimate_b) / std_error
  p_value <- 2 * pt(-abs(statistic), df = df_resid)
  crit <- qt(0.975, df = df_resid)

  tibble(
    contrast = contrast_name,
    estimate = estimate_a - estimate_b,
    std.error = std_error,
    statistic = statistic,
    p.value = p_value,
    conf.low = (estimate_a - estimate_b) - crit * std_error,
    conf.high = (estimate_a - estimate_b) + crit * std_error
  )
}

build_main_figure_panel_effects <- function(panel_score_long, panel_results_display) {
  contrast_lookup <- tibble(
    contrast = main_figure_forest_contrast_order,
    contrast_label = unname(main_figure_forest_contrast_labels[main_figure_forest_contrast_order]),
    contrast_order = seq_along(main_figure_forest_contrast_order)
  )

  map_dfr(main_figure_strict_panels, function(panel_name) {
    panel_df <- panel_score_long %>%
      filter(panel == panel_name)

    fit_lo <- lm(
      panel_score ~ Disease + Age + Sex + BMI,
      data = panel_df %>% filter(Age_group == "LO")
    )
    fit_eo <- lm(
      panel_score ~ Disease + Age + Sex + BMI,
      data = panel_df %>% filter(Age_group == "EO")
    )

    bind_rows(
      extract_lm_term(
        fit_lo,
        "DiseaseCRC",
        "LOCRC_vs_LOControl_ageadj"
      ),
      extract_lm_term(
        fit_eo,
        "DiseaseCRC",
        "EOCRC_vs_EOControl_ageadj"
      ),
      compare_independent_lm_terms(
        fit_eo,
        fit_lo,
        "DiseaseCRC",
        "EO_vs_LO_delta_beta_ageadj"
      )
    ) %>%
      mutate(
        panel = panel_name,
        panel_short = unname(main_figure_panel_short_labels[panel_name]),
        panel_label = unname(panel_labels[panel_name])
      )
  }) %>%
    mutate(
      adj.p.value = p.adjust(p.value, method = "BH")
    ) %>%
    left_join(contrast_lookup, by = "contrast") %>%
    left_join(
      panel_results_display %>%
        select(panel, contrast, limma_logFC = logFC, limma_P.Value = P.Value, limma_adj.P.Val = adj.P.Val),
      by = c("panel", "contrast")
    ) %>%
    mutate(
      contrast_label = factor(contrast_label, levels = rev(contrast_lookup$contrast_label)),
      panel_short = factor(
        panel_short,
        levels = unname(main_figure_panel_short_labels[main_figure_strict_panels])
      ),
      fdr_label = case_when(
        adj.p.value < 0.001 ~ "FDR < 0.001",
        adj.p.value < 0.01 ~ paste0("FDR = ", sprintf("%.3f", adj.p.value)),
        adj.p.value < 0.1 ~ paste0("FDR = ", sprintf("%.3f", adj.p.value)),
        TRUE ~ paste0("FDR = ", sprintf("%.2f", adj.p.value))
      )
    )
}

build_main_figure_canonical_effects <- function(canonical_targeted_pooled) {
  contrast_lookup <- tibble(
    contrast = main_figure_contrast_order,
    contrast_label = unname(main_figure_contrast_labels[main_figure_contrast_order])
  )

  effect_df <- canonical_targeted_pooled %>%
    filter(
      panel %in% main_figure_strict_panels,
      contrast %in% main_figure_contrast_order
    ) %>%
    inner_join(main_figure_metabolite_order, by = c("panel", "canonical_metabolite")) %>%
    left_join(contrast_lookup, by = "contrast") %>%
    mutate(
      panel_short = unname(main_figure_panel_short_labels[panel]),
      effect_label = paste0(
        sprintf("%.2f", mean_logFC),
        case_when(
          best_adj.P.Val < fdr_cutoff ~ "*",
          TRUE ~ ""
        )
      )
    )

  max_abs <- max(abs(effect_df$mean_logFC), na.rm = TRUE)
  effect_df %>%
    mutate(
      panel_short = factor(
        panel_short,
        levels = unname(main_figure_panel_short_labels[main_figure_strict_panels])
      ),
      contrast_label = factor(contrast_label, levels = unname(main_figure_contrast_labels[main_figure_contrast_order])),
      canonical_metabolite = factor(
        canonical_metabolite,
        levels = rev(main_figure_metabolite_order$canonical_metabolite[order(main_figure_metabolite_order$metabolite_order)])
      ),
      text_color = ifelse(abs(mean_logFC) > max_abs * 0.62, "white", "black"),
      max_abs = max_abs
    )
}

run_main_figure_leave_one_out <- function(feature_mat, panel_membership, eo_meta) {
  map_dfr(main_figure_strict_panels, function(panel_name) {
    ordered_metabolites <- main_figure_metabolite_order %>%
      filter(panel == panel_name) %>%
      arrange(metabolite_order) %>%
      pull(canonical_metabolite)

    scenario_tbl <- tibble(
      omitted_metabolite = c("Full panel", ordered_metabolites),
      omission_order = seq_along(c("Full panel", ordered_metabolites))
    )

    scenario_results <- map_dfr(scenario_tbl$omitted_metabolite, function(omitted_name) {
      membership_mod <- if (omitted_name == "Full panel") {
        panel_membership
      } else {
        panel_membership %>%
          filter(!(panel == panel_name & canonical_metabolite == omitted_name))
      }

      score_mat_mod <- build_panel_score_matrix(feature_mat, membership_mod)
      fit_df <- eo_meta %>%
        mutate(panel_score = as.numeric(score_mat_mod[panel_name, Sample_ID]))

      tidy(lm(panel_score ~ Disease + Age_centered + Sex + BMI, data = fit_df), conf.int = TRUE) %>%
        filter(term == "DiseaseCRC") %>%
        transmute(
          omitted_metabolite = omitted_name,
          estimate,
          conf.low,
          conf.high,
          p.value
        )
    }) %>%
      left_join(scenario_tbl, by = "omitted_metabolite") %>%
      mutate(
        panel = panel_name,
        panel_short = unname(main_figure_panel_short_labels[panel_name]),
        omission_label = ifelse(
          omitted_metabolite == "Full panel",
          "Full panel",
          paste("Drop", omitted_metabolite)
        )
      ) %>% left_join(panel_membership %>% distinct(panel, panel_label), by = 'panel') %>%
      arrange(omission_order) %>%
      mutate(
        full_panel_estimate = estimate[omitted_metabolite == "Full panel"][1],
        estimate_shift = estimate - full_panel_estimate,
        p_adj_within_panel = p.adjust(p.value, method = "BH"),
        sig_label = case_when(
          p_adj_within_panel < 0.001 ~ "***",
          p_adj_within_panel < 0.01 ~ "**",
          p_adj_within_panel < 0.05 ~ "*",
          TRUE ~ "ns"
        ),
        panel_omission = paste(panel_name, omission_label, sep = "||")
      )

    scenario_results
  }) %>%
    mutate(
      panel_omission = factor(panel_omission, levels = rev(unique(panel_omission))),
      panel_short = factor(
        panel_short,
        levels = unname(main_figure_panel_short_labels[main_figure_strict_panels])
      )
    )
}

plot_main_figure_panel_effects <- function(panel_effect_df) {
  text_df <- panel_effect_df %>%
    mutate(text_x = conf.high + 0.06)

  x_min <- min(panel_effect_df$conf.low, na.rm = TRUE) - 0.01
  x_max <- max(text_df$text_x, na.rm = TRUE) + 0.01

  ggplot(panel_effect_df, aes(x = estimate, y = contrast_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", 
               color = "grey55", linewidth = 0.6) +
    geom_errorbar(
      aes(xmin = conf.low, xmax = conf.high),
      orientation = "y",
      width = 0.18,
      linewidth = 0.6,
      color = "#4A4A4A"
    ) +
    geom_point(size = 2.8, color = "#B22222") +
    # geom_text(
    #   data = text_df,
    #   aes(x = text_x, y = contrast_label, label = fdr_label),
    #   inherit.aes = FALSE,
    #   hjust = 0,
    #   size = 5
    # ) +
    facet_wrap(~ panel_label, ncol = 1, scales = "free_y") +
    labs(
      # title = "A. Module effect sizes",
      x = "Adjusted panel-score difference",
      y = NULL
    ) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      axis.text.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      plot.margin = margin(5.5, 10, 5.5, 5.5),
      strip.text = element_text(size = 16),
      axis.title.x =  element_text(size = 16),
      
    ) +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off")
}

plot_main_figure_canonical_effects <- function(canonical_effect_df) {
  fill_limit <- max(canonical_effect_df$max_abs, na.rm = TRUE)

  ggplot(canonical_effect_df, aes(x = contrast_label, y = canonical_metabolite, fill = mean_logFC)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(
      aes(label = effect_label, color = text_color),
      size = 6,
      show.legend = FALSE
    ) +
    facet_wrap(~ panel_label, ncol = 1, scales = "free_y") +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-fill_limit, fill_limit)
    ) +
    scale_color_identity() +
    labs(
      # title = "B. Canonical metabolite effect map",
      # subtitle = "* FDR < 0.05",
      x = NULL,
      y = NULL,
      fill = "Mean log2 effect"
    ) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1,size = 16),
      axis.text.y = element_text(size = 16),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 16),
      legend.text = element_text(size = 12)
    )
}

plot_main_figure_leave_one_out <- function(loo_df) {
  full_ref <- loo_df %>%
    filter(omitted_metabolite == "Full panel") %>%
    select(panel_label, full_panel_estimate = estimate)

  ggplot(loo_df, aes(x = estimate, y = panel_omission)) +
    geom_vline(
      data = full_ref,
      aes(xintercept = full_panel_estimate),
      linetype = "dashed",
      color = "grey55",
      linewidth = 0.6
    ) +
    geom_errorbar(
      aes(xmin = conf.low, xmax = conf.high),
      orientation = "y",
      width = 0.18,
      linewidth = 0.6,
      color = "#4A4A4A"
    ) +
    geom_point(
      aes(shape = omitted_metabolite == "Full panel", 
          fill = p_adj_within_panel < fdr_cutoff),
      size = 2.8,
      color = "#B22222"
    ) +
    facet_wrap(~ panel_label, ncol = 1, scales = "free_y") +
    scale_shape_manual(values = c(`TRUE` = 23, `FALSE` = 21), guide = "none") +
    scale_fill_manual(
      values = c(`TRUE` = "#B22222", `FALSE` = "white"),
      guide = 'none'
      # labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "FDR >= 0.05"),
      # name = NULL
    ) +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) +
    labs(
      # title = "C. Leave-one-metabolite-out robustness",
      # subtitle = "Dashed line = full-panel EO effect",
      x = "EOCRC - EOControl adjusted effect",
      y = NULL
    ) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom",
      axis.text.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      strip.text = element_text(size = 16),
      axis.title.x = element_text(size = 16)
    )
}

build_species_support_strip_data <- function(species_candidates) {
  selected_species <- species_candidates %>%
    filter(
      panel == support_strip_target_panel,
      panel_fdr < support_strip_fdr_cutoff,
      disease_fdr < support_strip_fdr_cutoff,
      as.character(sign_panel_disease) %in% c("TRUE", "same_direction", "1")
    ) %>%
    mutate(
      panel_fdr = as.numeric(panel_fdr),
      disease_fdr = as.numeric(disease_fdr),
      lo_panel_fdr = as.numeric(lo_panel_fdr),
      lo_disease_fdr = as.numeric(lo_disease_fdr),
      eo_candidate_tier_support = case_when(
        panel_fdr < support_strip_fdr_cutoff &
          disease_fdr < support_strip_fdr_cutoff &
          as.character(sign_panel_disease) %in% c("TRUE", "same_direction", "1") &
          sensitivity_p < species_sensitivity_p_cutoff &
          as.character(sign_panel_sensitivity) %in% c("TRUE", "same_direction", "1") ~ "tier1_panel_disease_crc_support",
        panel_fdr < support_strip_fdr_cutoff &
          disease_fdr < support_strip_fdr_cutoff &
          as.character(sign_panel_disease) %in% c("TRUE", "same_direction", "1") ~ "tier2_panel_disease",
        panel_fdr < support_strip_fdr_cutoff ~ "tier3_panel_only",
        TRUE ~ "not_shortlisted"
      ),
      lo_candidate_tier_support = case_when(
        lo_panel_fdr < support_strip_fdr_cutoff &
          lo_disease_fdr < support_strip_fdr_cutoff &
          as.character(lo_sign_panel_disease) %in% c("TRUE", "same_direction", "1") &
          lo_sensitivity_p < species_sensitivity_p_cutoff &
          as.character(lo_sign_panel_sensitivity) %in% c("TRUE", "same_direction", "1") ~ "tier1_panel_disease_crc_support",
        lo_panel_fdr < support_strip_fdr_cutoff &
          lo_disease_fdr < support_strip_fdr_cutoff &
          as.character(lo_sign_panel_disease) %in% c("TRUE", "same_direction", "1") ~ "tier2_panel_disease",
        lo_panel_fdr < support_strip_fdr_cutoff ~ "tier3_panel_only",
        TRUE ~ "not_shortlisted"
      ),
      eo_shortlisted = eo_candidate_tier_support %in% support_strip_species_tiers,
      lo_shortlisted = lo_candidate_tier_support %in% support_strip_species_tiers,
      shortlist_origin = case_when(
        eo_shortlisted & lo_shortlisted ~ "supported_in_both_strata",
        eo_shortlisted ~ "EO_supported_only",
        lo_shortlisted ~ "LO_supported_only",
        TRUE ~ "neither"
      )
    ) %>%
    arrange(
      factor(shortlist_origin, levels = c("supported_in_both_strata", "EO_supported_only", "LO_supported_only", "neither")),
      pmin(panel_fdr, lo_panel_fdr, na.rm = TRUE),
      pmin(disease_fdr, lo_disease_fdr, na.rm = TRUE),
      desc(pmax(abs(as.numeric(panel_beta)), abs(as.numeric(lo_panel_beta)), na.rm = TRUE))
    ) %>%
    distinct(species_label, .keep_all = TRUE)

  support_df <- selected_species %>%
    transmute(
      item_label = gsub("_", " ", species_label),
      short_label = short_label,
      eo_candidate_tier = eo_candidate_tier_support,
      lo_candidate_tier = lo_candidate_tier_support,
      shortlist_origin,
      panel = panel,
      eo_prevalence = as.numeric(NA),
      contrast_key = "EO_disease_effect",
      effect = as.numeric(disease_beta),
      fdr = as.numeric(disease_fdr)
    ) %>%
    bind_rows(
      selected_species %>%
        transmute(
          item_label = gsub("_", " ", species_label),
          short_label = short_label,
          eo_candidate_tier = eo_candidate_tier_support,
          lo_candidate_tier = lo_candidate_tier_support,
          shortlist_origin,
          panel = panel,
          eo_prevalence = as.numeric(NA),
          contrast_key = "LO_disease_effect",
          effect = as.numeric(lo_disease_beta),
          fdr = as.numeric(lo_disease_fdr)
        )
    ) %>%
    bind_rows(
      selected_species %>%
        transmute(
          item_label = gsub("_", " ", species_label),
          short_label = short_label,
          eo_candidate_tier = eo_candidate_tier_support,
          lo_candidate_tier = lo_candidate_tier_support,
          shortlist_origin,
          panel = panel,
          eo_prevalence = as.numeric(NA),
          contrast_key = "EO_putrefactive_panel_association",
          effect = as.numeric(panel_beta),
          fdr = as.numeric(panel_fdr)
        )
    ) %>%
    bind_rows(
      selected_species %>%
        transmute(
          item_label = gsub("_", " ", species_label),
          short_label = short_label,
          eo_candidate_tier = eo_candidate_tier_support,
          lo_candidate_tier = lo_candidate_tier_support,
          shortlist_origin,
          panel = panel,
          eo_prevalence = as.numeric(NA),
          contrast_key = "LO_putrefactive_panel_association",
          effect = as.numeric(lo_panel_beta),
          fdr = as.numeric(lo_panel_fdr)
        )
    ) %>%
    mutate(
      contrast_label = recode(contrast_key, !!!support_strip_contrast_labels),
      sig_marker = case_when(
        is.na(fdr) ~ "",
        fdr < support_strip_fdr_cutoff ~ "*",
        TRUE ~ ""
      ),
      effect_label = if_else(is.na(effect), "NT", paste0(sprintf("%+.2f", effect), sig_marker))
    )

  max_abs <- max(abs(support_df$effect), na.rm = TRUE)
  order_levels <- selected_species %>%
    transmute(
      item_label = gsub("_", " ", species_label),
      eo_disease_beta = as.numeric(disease_beta),
      best_panel_fdr = pmin(as.numeric(panel_fdr), as.numeric(lo_panel_fdr), na.rm = TRUE),
      best_disease_fdr = pmin(as.numeric(disease_fdr), as.numeric(lo_disease_fdr), na.rm = TRUE),
      shortlist_origin = factor(shortlist_origin, levels = c("supported_in_both_strata", "EO_supported_only", "LO_supported_only", "neither"))
    ) %>%
    arrange(eo_disease_beta, best_disease_fdr, best_panel_fdr, shortlist_origin, item_label) %>%
    pull(item_label)

  support_df %>%
    mutate(
      item_label = factor(item_label, levels = rev(order_levels)),
      contrast_label = factor(contrast_label, levels = unname(support_strip_contrast_labels)),
      text_color = ifelse(abs(effect) > max_abs * 0.58, "white", "black"),
      text_color = ifelse(is.na(effect), "grey35", text_color),
      max_abs = max_abs,
      selection_rule = "Main-figure support strip restricted to EO-shortlisted putrefaction-core species: EO disease and EO panel association FDR < 0.05 with concordant direction; LO columns are contextual comparisons",
      shortlist_origin = factor(shortlist_origin, levels = c("supported_in_both_strata", "EO_supported_only", "LO_supported_only", "neither"))
    ) %>%
    arrange(desc(item_label), contrast_label)
}

plot_species_support_strip <- function(support_df) {
  fill_limit <- max(support_df$max_abs, na.rm = TRUE)

  ggplot(support_df, aes(x = contrast_label, y = item_label, fill = effect)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(
      aes(label = effect_label, color = text_color),
      size = 6,
      show.legend = FALSE
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-fill_limit, fill_limit)
    ) +
    scale_color_identity() +
    labs(
      # title = "D. Support strip",
      # subtitle = "EO shortlist selected by EO disease + EO spillover alignment; LO disease and LO spillover columns shown for age-stratified comparison; * FDR < 0.05",
      x = NULL,
      y = NULL,
      fill = "Effect"
    ) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      panel.grid = element_blank(),
      legend.text = element_text(size = 16)
    )
}
