
prepare_crossomic_matrix <- function(raw_df, feature_col, sample_map) {
  raw_df <- as.data.frame(raw_df)
  available_samples <- intersect(sample_map$Patient_ID, colnames(raw_df))

  if (length(available_samples) == 0) {
    stop("No overlapping cross-omics samples were found.")
  }

  rename_tbl <- sample_map %>%
    distinct(Patient_ID, Sample_ID)

  sample_ids <- rename_tbl$Sample_ID[match(available_samples, rename_tbl$Patient_ID)]

  prepared_df <- raw_df[, c(feature_col, available_samples), drop = FALSE]
  colnames(prepared_df) <- c("feature_id", sample_ids)

  feature_mat <- as.matrix(prepared_df[, -1, drop = FALSE])
  mode(feature_mat) <- "numeric"
  rownames(feature_mat) <- prepared_df$feature_id

  feature_mat
}

summarise_crossomic_qc <- function(
    feature_mat,
    sample_meta,
    subset_mask,
    group_var = "Disease",
    min_prevalence,
    min_group_prevalence,
    detect_threshold = 0) {
  subset_meta <- sample_meta[subset_mask, , drop = FALSE]
  subset_ids <- subset_meta$Sample_ID

  detect_mat <- feature_mat[, subset_ids, drop = FALSE] > detect_threshold
  overall_prevalence <- rowMeans(detect_mat, na.rm = TRUE)

  group_levels_local <- unique(as.character(subset_meta[[group_var]]))
  group_prev_df <- map_dfr(group_levels_local, function(group_name) {
    group_ids <- subset_meta$Sample_ID[as.character(subset_meta[[group_var]]) == group_name]

    tibble(
      feature_id = rownames(feature_mat),
      group_name = group_name,
      group_prevalence = rowMeans(feature_mat[, group_ids, drop = FALSE] > detect_threshold, na.rm = TRUE)
    )
  })

  group_summary <- group_prev_df %>%
    group_by(feature_id) %>%
    summarise(max_group_prevalence = max(group_prevalence, na.rm = TRUE), .groups = "drop")

  tibble(
    feature_id = rownames(feature_mat),
    overall_prevalence = overall_prevalence
  ) %>%
    left_join(group_summary, by = "feature_id") %>%
    mutate(
      keep_feature = overall_prevalence >= min_prevalence &
        max_group_prevalence >= min_group_prevalence
    ) %>%
    arrange(desc(keep_feature), desc(overall_prevalence), desc(max_group_prevalence))
}

clr_transform_matrix <- function(feature_mat, pseudocount = NULL) {
  feature_mat <- as.matrix(feature_mat)
  min_positive <- suppressWarnings(min(feature_mat[feature_mat > 0], na.rm = TRUE))

  if (!is.finite(min_positive)) {
    stop("CLR transform requires at least one positive value.")
  }

  if (is.null(pseudocount)) {
    pseudocount <- min_positive / 2
  }

  sample_mat <- t(feature_mat) + pseudocount
  sample_clr <- t(apply(sample_mat, 1, function(x) log(x) - mean(log(x))))
  colnames(sample_clr) <- rownames(feature_mat)
  rownames(sample_clr) <- colnames(feature_mat)

  t(sample_clr)
}

fit_single_coefficient_limma <- function(
    expr_mat,
    sample_meta,
    design_formula,
    coefficient_expression,
    contrast_name,
    model_label,
    test_type = c("ebayes", "treat"),
    lfc_threshold = 0) {
  fit_limma_contrasts(
    expr_mat = expr_mat,
    sample_meta = sample_meta,
    design_formula = design_formula,
    contrast_definitions = stats::setNames(coefficient_expression, contrast_name),
    model_label = model_label,
    test_type = match.arg(test_type),
    lfc_threshold = lfc_threshold
  )
}

extract_anchor_matrix <- function(feature_mat, anchor_ids, id_type = c("pathway", "ec")) {
  id_type <- match.arg(id_type)

  feature_key <- if (id_type == "pathway") {
    sub(":.*", "", rownames(feature_mat))
  } else {
    rownames(feature_mat)
  }

  keep_mask <- feature_key %in% anchor_ids
  anchor_mat <- feature_mat[keep_mask, , drop = FALSE]
  rownames(anchor_mat) <- feature_key[keep_mask]
  anchor_mat[!duplicated(rownames(anchor_mat)), , drop = FALSE]
}

extract_crossomic_candidates <- function(
    panel_assoc_results,
    disease_results,
    interaction_results,
    sensitivity_results,
    lo_sensitivity_results = NULL,
    lo_disease_results = NULL,
    lo_panel_assoc_results = NULL,
    feature_annotation) {
  candidate_seed <- bind_rows(
    panel_assoc_results %>%
      transmute(feature_id = metabolite_identification, panel = target_panel),
    if (is.null(lo_panel_assoc_results)) {
      tibble(feature_id = character(), panel = character())
    } else {
      lo_panel_assoc_results %>%
        transmute(feature_id = metabolite_identification, panel = target_panel)
    }
  ) %>%
    distinct(feature_id, panel)

  candidate_seed %>%
    left_join(
      panel_assoc_results %>%
        transmute(
          feature_id = metabolite_identification,
          panel = target_panel,
          panel_beta = logFC,
          panel_p = P.Value,
          panel_fdr = adj.P.Val
        ),
      by = c("feature_id", "panel")
    ) %>%
    left_join(
      disease_results %>%
        transmute(
          feature_id = metabolite_identification,
          disease_beta = logFC,
          disease_p = P.Value,
          disease_fdr = adj.P.Val
        ),
      by = "feature_id"
    ) %>%
    left_join(
      interaction_results %>%
        transmute(
          feature_id = metabolite_identification,
          interaction_beta = logFC,
          interaction_p = P.Value,
          interaction_fdr = adj.P.Val
        ),
      by = "feature_id"
    ) %>%
    left_join(
      sensitivity_results %>%
        transmute(
          feature_id = metabolite_identification,
          panel = target_panel,
          sensitivity_beta = logFC,
          sensitivity_p = P.Value,
          sensitivity_fdr = adj.P.Val
        ),
      by = c("feature_id", "panel")
    ) %>%
    left_join(
      if (is.null(lo_sensitivity_results)) {
        tibble(
          feature_id = character(),
          panel = character(),
          lo_sensitivity_beta = numeric(),
          lo_sensitivity_p = numeric(),
          lo_sensitivity_fdr = numeric()
        )
      } else {
        lo_sensitivity_results %>%
          transmute(
            feature_id = metabolite_identification,
            panel = target_panel,
            lo_sensitivity_beta = logFC,
            lo_sensitivity_p = P.Value,
            lo_sensitivity_fdr = adj.P.Val
          )
      },
      by = c("feature_id", "panel")
    ) %>%
    left_join(
      if (is.null(lo_disease_results)) {
        tibble(
          feature_id = character(),
          lo_disease_beta = numeric(),
          lo_disease_p = numeric(),
          lo_disease_fdr = numeric()
        )
      } else {
        lo_disease_results %>%
          transmute(
            feature_id = metabolite_identification,
            lo_disease_beta = logFC,
            lo_disease_p = P.Value,
            lo_disease_fdr = adj.P.Val
          )
      },
      by = "feature_id"
    ) %>%
    left_join(
      if (is.null(lo_panel_assoc_results)) {
        tibble(
          feature_id = character(),
          panel = character(),
          lo_panel_beta = numeric(),
          lo_panel_p = numeric(),
          lo_panel_fdr = numeric()
        )
      } else {
        lo_panel_assoc_results %>%
          transmute(
            feature_id = metabolite_identification,
            panel = target_panel,
            lo_panel_beta = logFC,
            lo_panel_p = P.Value,
            lo_panel_fdr = adj.P.Val
          )
      },
      by = c("feature_id", "panel")
    ) %>%
    left_join(feature_annotation, by = "feature_id") %>%
    mutate(
      sign_panel_disease = sign(panel_beta) == sign(disease_beta),
      sign_panel_sensitivity = sign(panel_beta) == sign(sensitivity_beta),
      lo_sign_panel_disease = sign(lo_panel_beta) == sign(lo_disease_beta),
      lo_sign_panel_sensitivity = sign(lo_panel_beta) == sign(lo_sensitivity_beta),
      candidate_tier = case_when(
        panel_fdr < species_candidate_fdr &
          disease_fdr < species_candidate_fdr &
          sign_panel_disease &
          sensitivity_p < species_sensitivity_p_cutoff &
          sign_panel_sensitivity ~ "tier1_panel_disease_crc_support",
        panel_fdr < species_candidate_fdr &
          disease_fdr < species_candidate_fdr &
          sign_panel_disease ~ "tier2_panel_disease",
        panel_fdr < species_candidate_fdr ~ "tier3_panel_only",
        TRUE ~ "not_shortlisted"
      ),
      lo_candidate_tier = case_when(
        lo_panel_fdr < species_candidate_fdr &
          lo_disease_fdr < species_candidate_fdr &
          lo_sign_panel_disease &
          lo_sensitivity_p < species_sensitivity_p_cutoff &
          lo_sign_panel_sensitivity ~ "tier1_panel_disease_crc_support",
        lo_panel_fdr < species_candidate_fdr &
          lo_disease_fdr < species_candidate_fdr &
          lo_sign_panel_disease ~ "tier2_panel_disease",
        lo_panel_fdr < species_candidate_fdr ~ "tier3_panel_only",
        TRUE ~ "not_shortlisted"
      ),
      shortlist_origin = case_when(
        candidate_tier %in% support_strip_species_tiers & lo_candidate_tier %in% support_strip_species_tiers ~ "supported_in_both_strata",
        candidate_tier %in% support_strip_species_tiers ~ "EO_supported_only",
        lo_candidate_tier %in% support_strip_species_tiers ~ "LO_supported_only",
        TRUE ~ "neither"
      )
    ) %>%
    arrange(
      panel,
      shortlist_origin,
      panel_fdr,
      disease_fdr,
      lo_panel_fdr,
      lo_disease_fdr,
      sensitivity_p,
      lo_sensitivity_p
    )
}

plot_candidate_species_heatmap <- function(candidate_tbl) {
  plot_tbl <- candidate_tbl %>%
    filter(candidate_tier != "not_shortlisted") %>%
    select(
      panel,
      target_panel_label,
      short_label,
      candidate_tier,
      disease_beta,
      panel_beta,
      interaction_beta
    ) %>%
    pivot_longer(
      cols = c(disease_beta, panel_beta, interaction_beta),
      names_to = "effect_type",
      values_to = "beta"
    ) %>%
    mutate(
      effect_type = factor(
        effect_type,
        levels = c("disease_beta", "panel_beta", "interaction_beta"),
        labels = c("EOCRC vs EOControl", "EO panel association", "Disease x age-group interaction")
      )
    )

  if (nrow(plot_tbl) == 0) {
    return(NULL)
  }

  row_order <- candidate_tbl %>%
    filter(candidate_tier != "not_shortlisted") %>%
    distinct(panel, short_label, candidate_tier, panel_fdr, disease_fdr) %>%
    arrange(panel, candidate_tier, panel_fdr, disease_fdr, short_label) %>%
    mutate(row_id = paste(panel, short_label, sep = "||"))

  plot_tbl <- plot_tbl %>%
    mutate(
      row_id = paste(panel, short_label, sep = "||"),
      row_id = factor(row_id, levels = rev(row_order$row_id))
    )

  ggplot(plot_tbl, aes(x = effect_type, y = row_id, fill = beta)) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(target_panel_label ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(
      low = "#1D4E89",
      mid = "white",
      high = "#B22222",
      midpoint = 0
    ) +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) +
    labs(
      title = "EO cross-omics candidate species",
      subtitle = "Rows are shortlisted by EO panel association, EO disease shift, and EOCRC sensitivity support",
      x = NULL,
      y = NULL,
      fill = "beta"
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1),
      strip.text.y = element_text(angle = 0)
    )
}

build_species_association_table <- function(
    panel_assoc_results,
    feature_annotation,
    qc_tbl) {
  assoc_tbl <- panel_assoc_results %>%
    rename(
      feature_id = metabolite_identification,
      panel = target_panel
    ) %>%
    left_join(feature_annotation, by = "feature_id") %>%
    left_join(
      qc_tbl %>%
        select(feature_id, overall_prevalence, max_group_prevalence),
      by = "feature_id"
    ) %>%
    mutate(
      target_panel_label = recode(panel, !!!panel_labels)
    )

  if ("t" %in% names(assoc_tbl)) {
    assoc_tbl <- assoc_tbl %>%
      mutate(
        se = if_else(!is.na(t) & abs(t) > 1e-8, abs(logFC / t), NA_real_),
        conf.low = logFC - 1.96 * se,
        conf.high = logFC + 1.96 * se
      )
  } else {
    assoc_tbl <- assoc_tbl %>%
      mutate(
        se = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_
      )
  }

  assoc_tbl %>%
    mutate(
      association_direction = case_when(
        adj.P.Val < fdr_cutoff & logFC > 0 ~ "positive",
        adj.P.Val < fdr_cutoff & logFC < 0 ~ "negative",
        TRUE ~ "nonsignificant"
      ),
      association_direction_fdr10 = "not_used"
    )
}

plot_top_associated_species <- function(association_tbl, fdr_threshold = fdr_cutoff, top_n = top_species_hits_per_direction) {
  plot_tbl <- association_tbl %>%
    filter(adj.P.Val < fdr_threshold) %>%
    mutate(
      direction_label = if_else(logFC > 0, "Positive association", "Negative association"),
      ci_low = if_else(is.finite(conf.low), conf.low, logFC),
      ci_high = if_else(is.finite(conf.high), conf.high, logFC)
    ) %>%
    group_by(panel, target_panel_label, direction_label) %>%
    arrange(adj.P.Val, desc(abs(logFC)), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup()

  if (nrow(plot_tbl) == 0) {
    return(NULL)
  }

  row_order <- plot_tbl %>%
    arrange(panel, direction_label, logFC, short_label) %>%
    mutate(row_id = paste(panel, direction_label, short_label, sep = "||")) %>%
    pull(row_id)

  plot_tbl <- plot_tbl %>%
    mutate(
      row_id = paste(panel, direction_label, short_label, sep = "||"),
      row_id = factor(row_id, levels = row_order)
    )

  ggplot(plot_tbl, aes(x = logFC, y = row_id, color = direction_label)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60", linewidth = 0.3) +
    geom_segment(aes(x = ci_low, xend = ci_high, yend = row_id), linewidth = 0.5, alpha = 0.8) +
    geom_point(aes(size = overall_prevalence), alpha = 0.9) +
    facet_grid(target_panel_label ~ direction_label, scales = "free_y", space = "free_y") +
    scale_color_manual(
      values = c(
        "Positive association" = "#B22222",
        "Negative association" = "#1D4E89"
      )
    ) +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|.*\\|\\|", "", x)) +
    labs(
      title = "Top EO panel-associated species",
      subtitle = "Effect sizes come from EO panel-association models adjusted for disease, age, sex, and BMI",
      x = "Adjusted beta for panel_score",
      y = NULL,
      size = "EO prevalence",
      color = NULL
    ) +
    theme(
      panel.grid.major.y = element_blank(),
      strip.text.y = element_text(angle = 0)
    )
}

build_taxonomic_cladogram_data <- function(association_tbl, fdr_threshold = fdr_cutoff) {
  sig_tbl <- association_tbl %>%
    filter(adj.P.Val < fdr_threshold) %>%
    mutate(direction_label = if_else(logFC > 0, "Positive", "Negative"))

  if (nrow(sig_tbl) == 0) {
    return(NULL)
  }

  map(unique(sig_tbl$panel), function(panel_name) {
    panel_tbl <- sig_tbl %>%
      filter(panel == panel_name) %>%
      arrange(phylum_label, family_label, genus_label, desc(logFC), short_label) %>%
      mutate(
        y = row_number(),
        x = 4,
        node_id = paste("species", feature_id, sep = "||"),
        parent_id = paste("genus", phylum_label, family_label, genus_label, sep = "||")
      )

    genus_tbl <- panel_tbl %>%
      group_by(panel, target_panel_label, phylum_label, family_label, genus_label) %>%
      summarise(
        y = mean(y),
        child_ymin = min(y),
        child_ymax = max(y),
        n_children = n(),
        .groups = "drop"
      ) %>%
      mutate(
        x = 3,
        node_id = paste("genus", phylum_label, family_label, genus_label, sep = "||"),
        parent_id = paste("family", phylum_label, family_label, sep = "||")
      )

    family_tbl <- genus_tbl %>%
      group_by(panel, target_panel_label, phylum_label, family_label) %>%
      summarise(
        y = mean(y),
        child_ymin = min(y),
        child_ymax = max(y),
        n_children = n(),
        .groups = "drop"
      ) %>%
      mutate(
        x = 2,
        node_id = paste("family", phylum_label, family_label, sep = "||"),
        parent_id = paste("phylum", phylum_label, sep = "||")
      )

    phylum_tbl <- family_tbl %>%
      group_by(panel, target_panel_label, phylum_label) %>%
      summarise(
        y = mean(y),
        child_ymin = min(y),
        child_ymax = max(y),
        n_children = n(),
        .groups = "drop"
      ) %>%
      mutate(
        x = 1,
        node_id = paste("phylum", phylum_label, sep = "||"),
        parent_id = NA_character_
      )

    vertical_segments <- bind_rows(
      phylum_tbl %>%
        transmute(panel, target_panel_label, x, xend = x, y = child_ymin, yend = child_ymax),
      family_tbl %>%
        transmute(panel, target_panel_label, x, xend = x, y = child_ymin, yend = child_ymax),
      genus_tbl %>%
        transmute(panel, target_panel_label, x, xend = x, y = child_ymin, yend = child_ymax)
    )

    horizontal_segments <- bind_rows(
      family_tbl %>%
        transmute(panel, target_panel_label, x = 1, xend = 2, y, yend = y),
      genus_tbl %>%
        transmute(panel, target_panel_label, x = 2, xend = 3, y, yend = y),
      panel_tbl %>%
        transmute(panel, target_panel_label, x = 3, xend = 4, y, yend = y)
    )

    list(
      vertical_segments = vertical_segments,
      horizontal_segments = horizontal_segments,
      species_points = panel_tbl,
      species_labels = panel_tbl %>%
        transmute(panel, target_panel_label, x = 4.10, y, label = short_label),
      internal_labels = bind_rows(
        phylum_tbl %>%
          transmute(panel, target_panel_label, x = 0.95, y, label = phylum_label, rank = "Phylum"),
        family_tbl %>%
          filter(n_children > 1) %>%
          transmute(panel, target_panel_label, x = 1.95, y, label = family_label, rank = "Family")
      )
    )
  }) %>%
    transpose() %>%
    map(bind_rows)
}

plot_taxonomic_cladogram <- function(cladogram_data, panel_name = NULL) {
  if (is.null(cladogram_data)) {
    return(NULL)
  }

  panel_names <- unique(as.character(cladogram_data$species_points$panel))
  panel_names <- panel_names[!is.na(panel_names)]
  if (length(panel_names) == 0) {
    return(NULL)
  }

  if (is.null(panel_name)) {
    panel_name <- panel_names[[1]]
  }

  species_df <- cladogram_data$species_points %>% filter(panel == panel_name)
  vertical_df <- cladogram_data$vertical_segments %>% filter(panel == panel_name)
  horizontal_df <- cladogram_data$horizontal_segments %>% filter(panel == panel_name)
  species_label_df <- cladogram_data$species_labels %>% filter(panel == panel_name)
  internal_label_df <- cladogram_data$internal_labels %>% filter(panel == panel_name)

  if (nrow(species_df) == 0) {
    return(NULL)
  }

  panel_title <- unique(species_df$target_panel_label)

  ggplot(species_df, aes(x = x, y = y)) +
    geom_segment(
      data = vertical_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 0.35,
      color = "grey55"
    ) +
    geom_segment(
      data = horizontal_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 0.35,
      color = "grey55"
    ) +
    geom_point(
      aes(fill = direction_label, size = overall_prevalence),
      shape = 21,
      color = "black",
      stroke = 0.25
    ) +
    geom_text(
      data = species_label_df,
      aes(x = x, y = y, label = label),
      hjust = 0,
      size = 2.8
    ) +
    geom_text(
      data = internal_label_df,
      aes(x = x, y = y, label = label),
      hjust = 1,
      size = 2.6,
      color = "grey20"
    ) +
    scale_fill_manual(values = c(Positive = "#B22222", Negative = "#1D4E89")) +
    scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("Phylum", "Family", "Genus", "Species"),
      limits = c(0.75, 4.75),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = panel_title,
      subtitle = "Species points mark significant positive or negative EO panel associations",
      x = NULL,
      y = NULL,
      fill = "Direction",
      size = "EO prevalence"
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

run_taxonomic_enrichment <- function(
    association_tbl,
    direction_var = "association_direction",
    significance_layer = "FDR<0.05",
    min_taxon_size = taxonomic_enrichment_min_taxon_size) {
  rank_cols <- c(family_label = "Family", genus_label = "Genus")

  map_dfr(names(rank_cols), function(rank_col) {
    map_dfr(unique(association_tbl$panel), function(panel_name) {
      panel_tbl <- association_tbl %>%
        filter(panel == panel_name, !is.na(.data[[rank_col]]), .data[[rank_col]] != "unclassified")

      if (nrow(panel_tbl) == 0) {
        return(tibble())
      }

      map_dfr(c("positive", "negative"), function(direction_name) {
        assoc_mask <- panel_tbl[[direction_var]] == direction_name
        total_assoc <- sum(assoc_mask, na.rm = TRUE)

        if (total_assoc == 0) {
          return(tibble())
        }

        taxon_counts <- panel_tbl %>%
          count(taxon_label = .data[[rank_col]], name = "taxon_total") %>%
          filter(taxon_total >= min_taxon_size)

        if (nrow(taxon_counts) == 0) {
          return(tibble())
        }

        map_dfr(taxon_counts$taxon_label, function(taxon_name) {
          a <- sum(panel_tbl[[rank_col]] == taxon_name & assoc_mask, na.rm = TRUE)
          b <- sum(panel_tbl[[rank_col]] == taxon_name & !assoc_mask, na.rm = TRUE)
          c <- sum(panel_tbl[[rank_col]] != taxon_name & assoc_mask, na.rm = TRUE)
          d <- sum(panel_tbl[[rank_col]] != taxon_name & !assoc_mask, na.rm = TRUE)

          fisher_fit <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")

          tibble(
            panel = panel_name,
            target_panel_label = unique(panel_tbl$target_panel_label),
            taxonomic_rank = rank_cols[[rank_col]],
            direction = direction_name,
            taxon_label = taxon_name,
            taxon_hits = a,
            taxon_total = a + b,
            panel_hits_total = total_assoc,
            background_total = nrow(panel_tbl),
            taxon_hit_fraction = if_else((a + b) > 0, a / (a + b), NA_real_),
            panel_hit_fraction = total_assoc / nrow(panel_tbl),
            odds_ratio = unname(fisher_fit$estimate),
            p.value = fisher_fit$p.value,
            significance_layer = significance_layer
          )
        })
      })
    })
  }) %>%
    group_by(panel, taxonomic_rank, direction, significance_layer) %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "BH")) %>%
    ungroup()
}

plot_taxonomic_enrichment <- function(enrichment_tbl, adj_threshold = fdr_cutoff) {
  plot_tbl <- enrichment_tbl %>%
    filter(adj.P.Val < adj_threshold) %>%
    group_by(panel, target_panel_label, taxonomic_rank, direction) %>%
    arrange(adj.P.Val, desc(taxon_hits), .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()

  if (nrow(plot_tbl) == 0) {
    return(NULL)
  }

  finite_or <- plot_tbl$odds_ratio[is.finite(plot_tbl$odds_ratio) & plot_tbl$odds_ratio > 0]
  max_or <- ifelse(length(finite_or) == 0, 4, max(finite_or))

  plot_tbl <- plot_tbl %>%
    mutate(
      odds_ratio_plot = case_when(
        is.infinite(odds_ratio) ~ max_or * 1.2,
        odds_ratio <= 0 ~ NA_real_,
        TRUE ~ odds_ratio
      ),
      log2_odds_ratio = log2(odds_ratio_plot),
      direction_label = if_else(direction == "positive", "Positive-associated taxa", "Negative-associated taxa"),
      row_id = paste(panel, taxonomic_rank, direction, taxon_label, sep = "||")
    )

  row_order <- plot_tbl %>%
    arrange(panel, taxonomic_rank, direction, log2_odds_ratio, taxon_label) %>%
    pull(row_id)

  plot_tbl <- plot_tbl %>%
    mutate(row_id = factor(row_id, levels = row_order))

  ggplot(plot_tbl, aes(x = log2_odds_ratio, y = row_id)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60", linewidth = 0.3) +
    geom_point(aes(size = taxon_hits, color = -log10(adj.P.Val))) +
    facet_grid(taxonomic_rank ~ direction_label + target_panel_label, scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|.*\\|\\|.*\\|\\|", "", x)) +
    scale_color_gradient(low = "#5B8E7D", high = "#B22222") +
    labs(
      title = "Taxonomic enrichment among EO panel-associated species",
      subtitle = "Fisher exact enrichment of family/genus among significant positive or negative panel-associated species",
      x = "log2 enrichment odds ratio",
      y = NULL,
      size = "Taxon hits",
      color = "-log10(FDR)"
    ) +
    theme(
      panel.grid.major.y = element_blank(),
      strip.text.y = element_text(angle = 0)
    )
}
