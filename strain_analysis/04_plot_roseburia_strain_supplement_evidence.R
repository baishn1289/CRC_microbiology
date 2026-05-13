#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(stringr)
  library(tibble)
})

env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

make_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

format_p <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  )
}

evidence_level <- function(p_value, q_value) {
  case_when(
    !is.na(q_value) & q_value < 0.10 ~ "FDR trend",
    !is.na(p_value) & p_value < 0.05 ~ "Nominal only",
    TRUE ~ "No support"
  )
}

evidence_palette <- c(
  "No support" = "#F1F0EA",
  "Nominal only" = "#E9B566",
  "FDR trend" = "#156B6F"
)

species_levels <- c("R. inulinivorans", "R. faecis")
species_labels <- c(
  "R. faecis" = expression(italic("R. faecis")),
  "R. inulinivorans" = expression(italic("R. inulinivorans"))
)

prepare_grid <- function(df, x_levels, target_x) {
  df %>%
    mutate(
      species_label = factor(species_label, levels = species_levels),
      stratum_label = factor(stratum_label, levels = x_levels),
      evidence = factor(
        evidence_level(p_value, q_bh_all_strain),
        levels = c("No support", "Nominal only", "FDR trend")
      ),
      is_target = stratum_label == target_x,
      is_fdr_trend = evidence == "FDR trend",
      is_nominal_only = evidence == "Nominal only",
      label = paste0(
        "p=", format_p(p_value),
        "\nq16=", format_p(q_bh_all_strain),
        "\nR2=", ifelse(is.na(r2), "NA", sprintf("%.4f", r2)),
        "\nn=", n_used_samples
      ),
      label_color = ifelse(evidence == "FDR trend", "white", "#202020")
    )
}

plot_evidence_heatmap <- function(df, title, subtitle, x_labels) {
  ggplot(df, aes(x = stratum_label, y = species_label)) +
    geom_tile(aes(fill = evidence), color = "#FFFFFF", linewidth = 1.8, width = 0.96, height = 0.90) +
    geom_tile(data = filter(df, is_target), fill = NA, color = "#9E2A2B", linewidth = 1.05, width = 0.96, height = 0.90) +
    geom_tile(data = filter(df, is_nominal_only), fill = NA, color = "#A05A00", linewidth = 0.95, linetype = "22", width = 0.96, height = 0.90) +
    geom_tile(data = filter(df, is_fdr_trend), fill = NA, color = "#0B1D2B", linewidth = 1.25, width = 0.96, height = 0.90) +
    geom_text(aes(label = label, color = label_color), size = 2.9, lineheight = 0.90) +
    scale_fill_manual(values = evidence_palette, drop = FALSE, name = "Evidence") +
    scale_color_identity() +
    scale_x_discrete(labels = x_labels, position = "top") +
    scale_y_discrete(labels = species_labels) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    coord_fixed(ratio = 0.56, clip = "off") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 12.5, color = "#202020", margin = margin(b = 3)),
      plot.subtitle = element_text(size = 9.2, color = "#4A4A4A", margin = margin(b = 6)),
      axis.text.x = element_text(face = "bold", size = 9.7, color = "#222222", margin = margin(b = 4), lineheight = 0.92),
      axis.text.y = element_text(size = 10.5, color = "#222222", margin = margin(r = 6)),
      panel.grid = element_blank(),
      legend.position = "none",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 10),
      plot.margin = margin(4, 10, 6, 10)
    )
}

make_card_plot <- function(cards) {
  ggplot(cards) +
    geom_rect(aes(xmin = x - 0.48, xmax = x + 0.48, ymin = 0.05, ymax = 0.95, fill = fill),
              color = "#D2D6D3", linewidth = 0.7) +
    geom_text(aes(x = x - 0.40, y = 0.72, label = title),
              hjust = 0, vjust = 1, size = 3.05, fontface = "bold", color = "#202020") +
    geom_text(aes(x = x - 0.40, y = 0.51, label = body),
              hjust = 0, vjust = 1, size = 2.35, lineheight = 0.96, color = "#303030") +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0.48, 3.52), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 10, 0, 10))
}

save_plot_dual <- function(plot_obj, output_dir, stem, width = 12.2, height = 7.7, dpi = 360) {
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")),
    plot_obj,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
  ggsave(
    file.path(output_dir, paste0(stem, ".png")),
    plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
}

project_root <- normalizePath(env_or_default("PROJECT_ROOT", getwd()), winslash = "/", mustWork = FALSE)
result_root <- env_or_default(
  "STRAINPHLAN_V5_RESULT_ROOT",
  file.path(project_root, "results/strain_analysis_v5")
)
summary_dir <- env_or_default(
  "STRAINPHLAN_ROSEBURIA_SUMMARY_DIR",
  file.path(result_root, "roseburia_conclusion_summary")
)
output_dir <- make_dir(env_or_default(
  "STRAINPHLAN_ROSEBURIA_STRAIN_SUPP_OUTPUT_DIR",
  file.path(summary_dir, "strain_supplement_evidence")
))

case_control <- read_tsv(file.path(summary_dir, "roseburia_strain_case_control_summary.tsv"), show_col_types = FALSE)
eolo <- read_tsv(file.path(summary_dir, "roseburia_strain_eolo_summary.tsv"), show_col_types = FALSE)

case_x <- c("Europe EO", "East Asia EO", "Europe LO", "East Asia LO")
eolo_x <- c("Europe CRC", "East Asia CRC", "Europe Control", "East Asia Control")

case_grid <- prepare_grid(case_control, case_x, "Europe EO")
eolo_grid <- prepare_grid(eolo, eolo_x, "Europe CRC")

rf_region <- read_tsv(
  file.path(result_root, "anpan_region_comparison/Europe_Asia_t__SGB4925/analysis_summary.tsv"),
  show_col_types = FALSE
) %>%
  filter(model_id == "01_phylogeny_region_overall") %>%
  transmute(species = "R. faecis", elpd_per_sample = pglmm_minus_base_elpd_per_sample)

ri_region <- read_tsv(
  file.path(result_root, "anpan_region_comparison_ri4940/Europe_Asia_t__SGB4940/analysis_summary.tsv"),
  show_col_types = FALSE
) %>%
  filter(model_id == "01_phylogeny_region_overall") %>%
  transmute(species = "R. inulinivorans", elpd_per_sample = pglmm_minus_base_elpd_per_sample)

region_card_text <- sprintf(
  "Clear in both species\ndELPD/sample %.3f / %.3f",
  rf_region$elpd_per_sample[[1]],
  ri_region$elpd_per_sample[[1]]
)

case_target_text <- case_grid %>%
  filter(stratum_label == "Europe EO") %>%
  arrange(species_label) %>%
  summarise(
    text = sprintf(
      "No support in both species\nq16 %.3f / %.3f",
      q_bh_all_strain[[1]],
      q_bh_all_strain[[2]]
    ),
    .groups = "drop"
  ) %>%
  pull(text)

exception_text <- paste(
  "East Asia-only weak signals",
  "q16=0.096; nominal q16=0.144",
  sep = "\n"
)

cards <- tibble(
  x = 1:3,
  fill = c("#EAF3F0", "#F7EDE8", "#F7F3DD"),
  title = c(
    "Europe-East Asia background",
    "Europe EO target tests",
    "Only weak positives elsewhere"
  ),
  body = c(
    region_card_text,
    case_target_text,
    exception_text
  )
)

case_plot <- plot_evidence_heatmap(
  case_grid,
  "A. Case-control strain structure within age-by-region strata",
  "Red outline = Europe EO target; black outline = q16 < 0.10; fill = evidence tier.",
  c(
    "Europe EO" = "Europe\nEO\n(target)",
    "East Asia EO" = "East Asia\nEO",
    "Europe LO" = "Europe\nLO",
    "East Asia LO" = "East Asia\nLO"
  )
)

eolo_plot <- plot_evidence_heatmap(
  eolo_grid,
  "B. EO vs LO strain structure within diagnosis-by-region strata",
  "Red outline = Europe CRC target contrast; dashed amber outline = nominal-only p < 0.05.",
  c(
    "Europe CRC" = "Europe\nCRC\n(target)",
    "East Asia CRC" = "East Asia\nCRC",
    "Europe Control" = "Europe\nControl",
    "East Asia Control" = "East Asia\nControl"
  )
)

card_plot <- make_card_plot(cards)

combined <- card_plot / case_plot / eolo_plot +
  plot_layout(heights = c(0.72, 1.85, 1.85)) +
  plot_annotation(
    title = "Roseburia strain follow-up: no Europe EOCRC-specific lineage replacement",
    subtitle = "Patristic-distance PERMANOVA on strain-retained samples with cohort-adjusted within-region models",
    caption = "Evidence tiers: teal = q16 < 0.10; amber = p < 0.05 but q16 >= 0.10; ivory = no support. q16 is BH-adjusted across all 16 strain-panel tiles. R2 is printed as effect size but not used as the main visual encoding."
  ) &
  theme(
    plot.title = element_text(face = "bold", size = 15.5, color = "#171717"),
    plot.subtitle = element_text(size = 10.2, color = "#444444"),
    plot.caption = element_text(size = 8.3, color = "#4B4B4B", hjust = 0)
  )

write_tsv(
  bind_rows(
    mutate(case_grid, panel = "case_control"),
    mutate(eolo_grid, panel = "eo_lo")
  ) %>%
    select(panel, species_label, stratum_label, evidence, p_value, q_bh_all_strain, r2, n_used_samples, is_target),
  file.path(output_dir, "roseburia_strain_supplement_evidence_source.tsv")
)

writeLines(
  c(
    "Roseburia strain supplement evidence figure",
    "",
    "Main conclusion:",
    "Europe EOCRC-specific dominant-lineage replacement is not supported by the current StrainPhlAn patristic-distance analyses.",
    "",
    "Visual encoding:",
    "- Cell fill is evidence tier, not R2. This avoids making small-sample, non-significant Europe EO R2 values look like strong support.",
    "- Red outlines mark the target tests relevant to the Europe EO abundance anomaly.",
    "- Black outline marks q16 < 0.10; dashed amber outline marks nominal-only p < 0.05."
  ),
  file.path(output_dir, "roseburia_strain_supplement_evidence_notes.txt")
)

save_plot_dual(combined, output_dir, "roseburia_strain_supplement_evidence_grid", width = 12.4, height = 7.8)

message("Wrote Roseburia strain supplement evidence figure to: ", output_dir)
