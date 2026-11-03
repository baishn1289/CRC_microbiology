library(data.table)
library(sandwich)
library(tidyverse)

# Stratified disease models and direct CRC comparison.
fit_age_effects <- function(y, meta, panel_name) {
  d <- as.data.frame(meta)
  d$score <- as.numeric(y)
  d$Disease <- factor(d$Disease, levels = c("Control", "CRC"))
  d$Sex <- factor(d$Sex, levels = c("female", "male"))
  fits <- lapply(
    c("EO", "LO"),
    function(g) lm(score ~ Disease + Age + Sex + BMI, data = d[d$Age_group == g, ])
  )
  beta <- vapply(fits, function(f) unname(coef(f)["DiseaseCRC"]), numeric(1))
  estimates <- c(beta, beta[1] - beta[2])
  residual_df <- vapply(fits, df.residual, numeric(1))
  rows <- list()

  for (method in c("HC3_normal", "ordinary_lm")) {
    variances <- vapply(
      fits,
      function(f) {
        if (method == "HC3_normal")
          sandwich::vcovHC(f, type = "HC3")["DiseaseCRC", "DiseaseCRC"] else
          vcov(f)["DiseaseCRC", "DiseaseCRC"]
      },
      numeric(1)
    )
    se <- sqrt(c(variances, sum(variances)))
    df <- c(residual_df, sum(variances)^2 / sum(variances^2 / residual_df))
    critical <- if (method == "HC3_normal") rep(qnorm(0.975), 3) else qt(0.975, df)
    p <- if (method == "HC3_normal") 2 * pnorm(-abs(estimates / se)) else
      2 * pt(-abs(estimates / se), df)
    rows[[method]] <- data.table::data.table(
      panel = panel_name,
      contrast = c("EO", "LO", "EO_minus_LO"),
      inference = method,
      estimate = estimates,
      std_error = se,
      df = if (method == "HC3_normal") NA_real_ else df,
      ci_low = estimates - critical * se,
      ci_high = estimates + critical * se,
      p_value = p
    )
  }

  list(results = data.table::rbindlist(rows), fits = fits)
}

fit_crc_effect <- function(y, design, panel_name, model_name) {
  d <- as.data.frame(design)
  names(d) <- make.names(colnames(design))
  d$score <- as.numeric(y)
  fit <- lm(score ~ . - 1, data = d)
  beta <- unname(coef(fit)["Age_groupEO"])
  rows <- list()

  for (method in c("HC3_normal", "ordinary_lm")) {
    se <- sqrt(if (method == "HC3_normal")
      sandwich::vcovHC(fit, type = "HC3")["Age_groupEO", "Age_groupEO"] else
      vcov(fit)["Age_groupEO", "Age_groupEO"])
    critical <- if (method == "HC3_normal") qnorm(0.975) else qt(0.975, df.residual(fit))
    p <- if (method == "HC3_normal")
      2 * pnorm(-abs(beta / se)) else
      2 * pt(-abs(beta / se), df.residual(fit))
    rows[[method]] <- data.table::data.table(
      panel = panel_name,
      contrast = "EOCRC_minus_LOCRC",
      model = model_name,
      inference = method,
      n = nrow(design),
      n_EO = sum(design[, "Age_groupEO"] == 1),
      n_LO = sum(design[, "Age_groupEO"] == 0),
      estimate = beta,
      std_error = se,
      df = if (method == "HC3_normal") NA_real_ else df.residual(fit),
      ci_low = beta - critical * se,
      ci_high = beta + critical * se,
      p_value = p
    )
  }

  list(results = data.table::rbindlist(rows), fit = fit)
}

# Sample alignment, prevalence filtering and natural-log CLR.
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
  detect_threshold = 0
) {
  subset_meta <- sample_meta[subset_mask, , drop = FALSE]
  subset_ids <- subset_meta$Sample_ID
  detect_mat <- feature_mat[, subset_ids, drop = FALSE] > detect_threshold
  overall_prevalence <- rowMeans(detect_mat, na.rm = TRUE)
  group_levels_local <- unique(as.character(subset_meta[[group_var]]))
  group_prev_df <- map_dfr(
    group_levels_local,
    function(group_name) {
      group_ids <- subset_meta$Sample_ID[as.character(subset_meta[[group_var]]) == group_name]

      tibble(
        feature_id = rownames(feature_mat),
        group_name = group_name,
        group_prevalence = rowMeans(
          feature_mat[, group_ids, drop = FALSE] > detect_threshold,
          na.rm = TRUE
        )
      )
    }
  )

  group_summary <- group_prev_df %>%
    group_by(feature_id) %>%
    summarise(max_group_prevalence = max(group_prevalence, na.rm = TRUE), .groups = "drop")

  tibble(feature_id = rownames(feature_mat), overall_prevalence = overall_prevalence) %>%
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

# Retained score, coefficient and member-map plotting functions.
plot_panel_scores <- function(
  panel_score_long,
  panel_annotations,
  label_style = c("no_sig_only", "beta_sig"),
  ncol = 4,
  beta_offset = 0.08,
  sig_offset = 0.045,
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
      Group4 = factor(Group4, levels = c("LOControl", "LOCRC", "EOControl", "EOCRC"))
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

  ann$y_tick <- ann$y_position - bracket_depth * ann$y_span
  ann$y_beta <- ann$y_position + beta_offset * ann$y_span
  ann$y_sig <- ann$y_position - sig_offset * ann$y_span

  p <- ggplot(panel_score_long, aes(
    x = Group4,
    y = panel_score,
    fill = Group4
  )) +
    geom_jitter(
      aes(color = Group4),
      width = 0.15,
      alpha = 0.8,
      size = 0.55,
      show.legend = FALSE
    ) +
    geom_boxplot(
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
    facet_wrap(
      ~panel_label,
      scales = "free_y",
      ncol = ncol,
      labeller = as_labeller(function(x) {
        stringr::str_replace(x, "^Extended: ", "Extended panel\n")
      })
    ) +
    scale_fill_manual(
      values = group_palette,
      breaks = c("LOControl", "LOCRC", "EOControl", "EOCRC"),
      drop = FALSE
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
    labs(x = NULL, y = "Mean feature z-score") +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        hjust = 1,
        angle = 25,
        size = text_size
      ),
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
        aes(
          x = x_start,
          xend = x_end,
          y = y_position,
          yend = y_position
        ),
        inherit.aes = FALSE,
        linewidth = 0.45
      ) +
      geom_segment(
        data = ann,
        aes(
          x = x_start,
          xend = x_start,
          y = y_tick,
          yend = y_position
        ),
        inherit.aes = FALSE,
        linewidth = 0.45
      ) +
      geom_segment(
        data = ann,
        aes(
          x = x_end,
          xend = x_end,
          y = y_tick,
          yend = y_position
        ),
        inherit.aes = FALSE,
        linewidth = 0.45
      ) +
      geom_text(
        data = ann,
        aes(
          x = x_mid,
          y = y_beta,
          label = beta_label
        ),
        inherit.aes = FALSE,
        size = 4,
        lineheight = 0.92
      ) +
      geom_text(
        data = ann,
        aes(
          x = x_mid,
          y = y_sig,
          label = sig_label_plot
        ),
        inherit.aes = FALSE,
        size = 4,
        lineheight = 0.92
      )
  }

  p
}

plot_main_figure_panel_effects <- function(panel_effect_df) {
  text_df <- panel_effect_df %>%
    mutate(text_x = conf.high + 0.06)

  x_min <- min(panel_effect_df$conf.low, na.rm = TRUE) - 0.01
  x_max <- max(text_df$text_x, na.rm = TRUE) + 0.01

  ggplot(panel_effect_df, aes(x = estimate, y = contrast_label)) +
    geom_vline(
      xintercept = 0,
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
    geom_point(size = 2.8, color = "#B22222") +
    facet_wrap(~panel_label, ncol = 1, scales = "free_y") +
    labs(x = "Adjusted panel-score difference", y = NULL) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      axis.text.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      plot.margin = margin(5.5, 10, 5.5, 5.5),
      strip.text = element_text(size = 16),
      axis.title.x = element_text(size = 16),
    ) +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off")
}

plot_main_figure_canonical_effects <- function(canonical_effect_df) {
  fill_limit <- max(canonical_effect_df$max_abs, na.rm = TRUE)

  ggplot(
    canonical_effect_df,
    aes(
      x = contrast_label,
      y = canonical_metabolite,
      fill = mean_logFC
    )
  ) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(
      aes(label = effect_label, color = text_color),
      size = 6,
      show.legend = FALSE
    ) +
    facet_wrap(~panel_label, ncol = 1, scales = "free_y") +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-fill_limit, fill_limit)
    ) +
    scale_color_identity() +
    labs(
      x = NULL,
      y = NULL,
      fill = "Mean log2 effect"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 25,
        hjust = 1,
        size = 16
      ),
      axis.text.y = element_text(size = 16),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 16),
      legend.text = element_text(size = 12)
    )
}
