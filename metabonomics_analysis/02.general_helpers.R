
safe_file_tag <- function(x) {
  x %>%
    gsub("[^A-Za-z0-9]+", "_", .) %>%
    gsub("^_+|_+$", "", .)
}

extract_taxon_rank <- function(feature_ids, rank_code) {
  rank_pattern <- paste0(".*\\|", rank_code, "__([^|]+).*")
  ifelse(
    grepl(paste0("\\|", rank_code, "__"), feature_ids),
    sub(rank_pattern, "\\1", feature_ids),
    NA_character_
  )
}

build_species_taxonomy_annotation <- function(feature_ids) {
  tibble(feature_id = feature_ids) %>%
    mutate(
      phylum_label = extract_taxon_rank(feature_id, "p"),
      family_label = extract_taxon_rank(feature_id, "f"),
      genus_label = extract_taxon_rank(feature_id, "g"),
      species_label = extract_taxon_rank(feature_id, "s"),
      short_label = coalesce(species_label, feature_id)
    )
}

significance_label <- function(p_value) {
  case_when(
    is.na(p_value) ~ "NA",
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

save_gg <- function(plot_obj, filename, width = 10, height = 7) {
  output_file <- filename
  tryCatch(
    {
      ggsave(
        filename = output_file,
        plot = plot_obj,
        width = width,
        height = height,
        dpi = 300,
        bg = "white"
      )
    },
    error = function(e) {
      fallback_file <- sub("\\.pdf$", "_rerun.pdf", output_file, ignore.case = TRUE)
      message("Primary PDF write failed for ", output_file, "; retrying as ", fallback_file)
      ggsave(
        filename = fallback_file,
        plot = plot_obj,
        width = width,
        height = height,
        dpi = 300,
        bg = "white"
      )
    }
  )
}

run_betadisper_permutest <- function(dispersion_obj, permutations = 999) {
  if (!inherits(dispersion_obj, "betadisper")) {
    stop("run_betadisper_permutest() expects a vegan::betadisper object.")
  }

  # Call the vegan betadisper S3 method explicitly so the script does not depend
  # on the current interactive session's generic/method dispatch state.
  permutest_betadisper <- getS3method("permutest", "betadisper", envir = asNamespace("vegan"))
  permutest_betadisper(dispersion_obj, permutations = permutations)
}

pareto_scale_matrix <- function(mat) {
  mat <- as.matrix(mat)
  center_vec <- colMeans(mat, na.rm = TRUE)
  sd_vec <- apply(mat, 2, sd, na.rm = TRUE)
  scale_vec <- sqrt(sd_vec)
  scale_vec[is.na(scale_vec) | scale_vec == 0] <- 1
  scaled <- sweep(mat, 2, center_vec, "-")
  sweep(scaled, 2, scale_vec, "/")
}

z_score_features <- function(feature_by_sample_mat) {
  z_mat <- t(apply(feature_by_sample_mat, 1, function(x) {
    x <- as.numeric(x)
    if (all(is.na(x))) {
      return(rep(NA_real_, length(x)))
    }
    x_sd <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(x_sd) || x_sd == 0) {
      out <- rep(0, length(x))
      out[is.na(x)] <- NA_real_
      return(out)
    }
    (x - mean(x, na.rm = TRUE)) / x_sd
  }))
  rownames(z_mat) <- rownames(feature_by_sample_mat)
  colnames(z_mat) <- colnames(feature_by_sample_mat)
  z_mat
}

safe_colmeans_available <- function(mat) {
  mat <- as.matrix(mat)
  observed_n <- colSums(!is.na(mat))
  out <- colMeans(mat, na.rm = TRUE)
  out[observed_n == 0] <- NA_real_
  out
}

read_mode_file <- function(path, mode_name) {
  raw <- fread(path)
  sample_cols <- grep("^[oy]", names(raw), value = TRUE)
  raw %>%
    select(all_of(c("metabolite_identification", sample_cols))) %>%
    mutate(mode = mode_name)
}

score_mode_candidates <- function(lcms_df) {
  sample_cols <- setdiff(colnames(lcms_df), c("metabolite_identification", "mode"))

  lcms_df %>%
    rowwise() %>%
    mutate(
      n_detected = sum(!is.na(c_across(all_of(sample_cols))) & c_across(all_of(sample_cols)) > 0),
      missing_rate = mean(is.na(c_across(all_of(sample_cols)))),
      median_signal = median(c_across(all_of(sample_cols)), na.rm = TRUE),
      prefer_mode = case_when(
        grepl(
          "acid$|Acid$|butyr|propion|succin|lactate|acetate|citrate|fumar|malate",
          metabolite_identification,
          ignore.case = TRUE
        ) ~ "negative",
        grepl(
          "amine|pyridine|pyrid|ornithine|lysine|arginine|putrescine|cadaverine",
          metabolite_identification,
          ignore.case = TRUE
        ) ~ "positive",
        TRUE ~ NA_character_
      ),
      mode_rank = case_when(
        !is.na(prefer_mode) & mode == prefer_mode ~ 1,
        !is.na(prefer_mode) & mode != prefer_mode ~ 2,
        TRUE ~ 3
      )
    ) %>%
    ungroup() %>%
    group_by(metabolite_identification) %>%
    arrange(mode_rank, missing_rate, desc(median_signal), .by_group = TRUE) %>%
    mutate(
      selection_rank = row_number(),
      selected_mode = selection_rank == 1,
      selection_reason = case_when(
        selected_mode & !is.na(prefer_mode) & mode == prefer_mode ~ "preferred_ion_mode",
        selected_mode & !is.na(prefer_mode) & mode != prefer_mode ~ "fallback_nonpreferred_mode",
        selected_mode ~ "best_missingness_then_signal",
        TRUE ~ "not_selected"
      )
    ) %>%
    ungroup()
}

deduplicate_modes <- function(lcms_df) {
  scored_df <- if ("selection_rank" %in% names(lcms_df)) lcms_df else score_mode_candidates(lcms_df)
  scored_df %>%
    filter(selected_mode) %>%
    arrange(metabolite_identification)
}

extract_limma_results <- function(fit_obj, contrast_names, effect_cutoff, fdr_cutoff) {
  map_dfr(contrast_names, function(contrast_name) {
    topTable(fit_obj, coef = contrast_name, number = Inf, sort.by = "P") %>%
      rownames_to_column("metabolite_identification") %>%
      mutate(
        contrast = contrast_name,
        significant_fdr = adj.P.Val < fdr_cutoff,
        significant_fdr_fc = adj.P.Val < fdr_cutoff & abs(logFC) >= effect_cutoff,
        significant_fdr_0_10 = FALSE,
        significant_fdr_fc_0_10 = FALSE,
        direction = case_when(
          logFC > 0 ~ "higher",
          logFC < 0 ~ "lower",
          TRUE ~ "no_change"
        )
      )
  })
}

fit_limma_contrasts <- function(
    expr_mat,
    sample_meta,
    design_formula,
    contrast_definitions,
    model_label,
    test_type = c("treat", "ebayes"),
    lfc_threshold = 0) {
  test_type <- match.arg(test_type)

  design <- model.matrix(design_formula, data = sample_meta)
  colnames(design) <- make.names(colnames(design))

  contrast_matrix <- do.call(
    limma::makeContrasts,
    c(as.list(contrast_definitions), list(levels = design))
  )

  fit0 <- lmFit(expr_mat, design)
  fit1 <- contrasts.fit(fit0, contrast_matrix)

  if (test_type == "treat") {
    fit2 <- treat(fit1, lfc = lfc_threshold, robust = TRUE)
    result_tbl <- map_dfr(colnames(contrast_matrix), function(contrast_name) {
      topTreat(fit2, coef = contrast_name, number = Inf, sort.by = "P") %>%
        rownames_to_column("metabolite_identification") %>%
        mutate(
          contrast = contrast_name,
          model_label = model_label,
          test_type = test_type,
          lfc_threshold = lfc_threshold,
          significant_fdr = adj.P.Val < fdr_cutoff,
          significant_fdr_fc = adj.P.Val < fdr_cutoff,
          significant_fdr_0_10 = FALSE,
          significant_fdr_fc_0_10 = FALSE,
          direction = case_when(
            logFC > 0 ~ "higher",
            logFC < 0 ~ "lower",
            TRUE ~ "no_change"
          )
        )
    })
  } else {
    fit2 <- eBayes(fit1, robust = TRUE)
    result_tbl <- extract_limma_results(
      fit_obj = fit2,
      contrast_names = colnames(contrast_matrix),
      effect_cutoff = lfc_threshold,
      fdr_cutoff = fdr_cutoff
    ) %>%
      mutate(
        model_label = model_label,
        test_type = test_type,
        lfc_threshold = lfc_threshold
      )
  }

  list(
    results = result_tbl,
    fit = fit2,
    design = design,
    contrast_matrix = contrast_matrix
  )
}

select_volcano_labels <- function(df, targeted_membership, n_global = 6, n_targeted = 8) {
  top_global <- df %>%
    slice_min(order_by = adj.P.Val, n = n_global, with_ties = FALSE) %>%
    pull(metabolite_identification)

  top_targeted <- df %>%
    semi_join(targeted_membership, by = "metabolite_identification") %>%
    slice_min(order_by = adj.P.Val, n = n_targeted, with_ties = FALSE) %>%
    pull(metabolite_identification)

  unique(c(top_global, top_targeted))
}

plot_volcano <- function(df, contrast_name, targeted_membership) {
  label_features <- select_volcano_labels(df, targeted_membership)

  df <- df %>%
    mutate(
      neg_log10_fdr = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
      label = ifelse(metabolite_identification %in% label_features, metabolite_identification, NA_character_)
    )

  ggplot(df, aes(x = logFC, y = neg_log10_fdr)) +
    geom_point(
      aes(color = significant_fdr_fc),
      alpha = 0.75,
      size = 1.6
    ) +
    geom_vline(xintercept = c(-effect_cutoff, effect_cutoff), linetype = 2, color = "grey50") +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = 2, color = "grey50") +
    geom_text_repel(
      data = df %>% filter(!is.na(label)),
      aes(label = label),
      size = 3.1,
      max.overlaps = Inf,
      box.padding = 0.25,
      point.padding = 0.15,
      segment.color = "grey60"
    ) +
    scale_color_manual(values = c("FALSE" = "grey72", "TRUE" = "#C73E1D")) +
    labs(
      title = contrast_name,
      x = "log2 fold-change",
      y = "-log10(FDR)",
      color = paste0("FDR < ", fdr_cutoff, "\n& |log2FC| >= ", round(effect_cutoff, 3))
    ) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}
