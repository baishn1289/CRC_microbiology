library(tidyverse)
library(data.table)
library(patchwork)
library(cowplot)

code_dir <- file.path(
  base_symbol,
  "downstream_analysis/src/github_code_availability/submission_version_261102/metabolomics_v35"
)
output_dir <- file.path(
  base_symbol,
  "downstream_analysis/output/submission_version_261102/metabolomics_v35"
)
analysis_dir <- file.path(base_symbol, "FUSCC_metabonomics_analysis")
v35_dir <- output_dir

source(file.path(code_dir, "00_helpers.R"))
figure_components <- readRDS(file.path(output_dir, "objects/figure_components.rds"))
main <- figure_components$main
supp_plots <- figure_components$supp_plots
supp <- readRDS(file.path(output_dir, "objects/supplementary_statistics.rds"))
species <- readRDS(file.path(output_dir, "objects/species_associations_v35.rds"))
theme_set(theme_bw(base_size = 6.5))
page_width <- 210 / 25.4
page_height <- 297 / 25.4
source_dir <- file.path(output_dir, "plots")

# Retain the selected species order and same-data CRC support flags.
selected <- copy(species$crc_support)
selected[, block := ifelse(disease_beta > 0, "EOCRC-enriched", "EOCRC-reduced")]
selected[, block_order := ifelse(disease_beta > 0, 1L, 2L)]
setorder(selected, block_order, -disease_beta)
selected[, row_position := .N - seq_len(.N) + 1 + ifelse(block_order == 1L, 2, 0)]
selected[, named_binomial := grepl("^[A-Z][a-z]+ [a-z]+$", species_label)]
selected[, label_expression := ifelse(
  named_binomial,
  paste0("italic(", encodeString(species_label, quote = '"'), ")"),
  encodeString(species_label, quote = '"')
)]
association <- merge(
  species$results[analysis_role == "panel_association" & panel == "D0"],
  selected[, .(feature_id, row_position, block, species_label_display = species_label)],
  by = "feature_id"
)
association[, label := ifelse(q_BH <= 0.1, "*", "")]
direction <- merge(
  CJ(feature_id = selected$feature_id, age_group = c("EO", "LO")),
  species$results[analysis_role == "disease", .(feature_id, age_group, estimate, q_BH)],
  by = c("feature_id", "age_group"),
  all.x = TRUE
)
direction <- merge(direction, selected[, .(feature_id, row_position)], by = "feature_id")
direction[, disease_direction := ifelse(
  is.na(estimate),
  NA_character_,
  ifelse(estimate > 0, "Higher in CRC", "Lower in CRC")
)]
direction[, supported := !is.na(q_BH) & q_BH <= 0.1]
nitrogen <- list(
  selected = selected,
  association = association,
  direction = direction
)

# Final A4 layout; ordinary R PDF devices only.
paper_plot <- function(
  p,
  text_mm = 2.1,
  point_mm = 0.8
) {
  p <- p + labs(caption = NULL, subtitle = NULL) +
    theme(
      axis.text = element_text(size = 6.3),
      axis.text.x = element_text(size = 6.3),
      axis.text.y = element_text(size = 6.3),
      axis.title = element_text(size = 6.4),
      axis.title.x = element_text(size = 6.4),
      axis.title.y = element_text(size = 6.4),
      strip.text = element_text(size = 7.0, face = "bold"),
      strip.text.x = element_text(size = 7.0, face = "bold"),
      strip.text.y = element_text(size = 7.0, face = "bold"),
      plot.title = element_text(size = 7.5, face = "bold"),
      plot.caption = element_blank(),
      legend.title = element_text(size = 6.3),
      legend.text = element_text(size = 6),
      legend.key.size = grid::unit(2, "mm"),
      plot.margin = margin(2, 2, 2, 2),
      panel.grid.minor = element_blank()
    )

  for (i in seq_along(p$layers)) {
    layer <- p$layers[[i]]
    if (inherits(layer$geom, "GeomText")) layer$aes_params$size <- text_mm
    if (inherits(layer$geom, "GeomPoint")) {
      layer$aes_params$size <- point_mm
      layer$aes_params$stroke <- 0.15
    }

    if (inherits(layer$geom, "GeomBoxplot") || inherits(layer$geom, "GeomErrorbar") ||
      inherits(layer$geom, "GeomSegment")) layer$aes_params$linewidth <- 0.25
  }

  p
}

score <- paper_plot(main$score_plot, point_mm = 0.35) + labs(y = "Mean member z-score")
score <- score + scale_x_discrete(labels = c(
  LOControl = "LO\nControl",
  LOCRC = "LO\nCRC",
  EOControl = "EO\nControl",
  EOCRC = "EO\nCRC"
)) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      size = 6.1
    )
  )
score$facet$params$labeller <- labeller(panel_label = c(
  "Arginine-ornithine-proline\nmetabolites" = "Arginine-ornithine-\nproline metabolites",
  "Aminovalerate-polyamine\nmetabolites" = "Aminovalerate-\npolyamine metabolites"
))
forest <- paper_plot(main$main_forest, point_mm = 0.9)

for (i in seq_along(forest$layers)) if (inherits(forest$layers[[i]]$geom, "GeomVline")) {
  forest$layers[[i]]$aes_params$linewidth <- 0.25
  forest$layers[[i]]$aes_params$colour <- "grey65"
  forest$layers[[i]]$aes_params$linetype <- "dashed"
}

forest$coordinates$limits$x <- c(min(forest$data$ci_low) - 0.05, max(forest$data$ci_high) + 0.36)
members <- paper_plot(main$canonical_plot, text_mm = 2.15) +
  guides(fill = guide_colourbar(barheight = grid::unit(18, "mm"), barwidth = grid::unit(1.5, "mm")))
top <- (score | forest | members) + plot_layout(widths = c(0.85, 1.25, 1.20)) +
  plot_annotation(tag_levels = list(c("B", "C", "D"))) &
  theme(plot.tag = element_text(size = 12, face = "bold"))

selected <- copy(nitrogen$selected)
association <- copy(nitrogen$association)
direction <- copy(nitrogen$direction)
y_limits <- c(0.4, max(selected$row_position) + 1.4)
boundary <- max(selected[block_order == 2L, row_position]) + 1
x_limits <- range(c(association$ci_low, association$ci_high)) + c(-0.15, 0.15)
row_theme <- theme_bw(base_size = 6.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_rect(linewidth = 0.25),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 6),
    axis.title.x = element_text(size = 6),
    plot.title = element_text(
      size = 7,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "none",
    plot.margin = margin(1.5, 1.5, 1.5, 1.5)
  )
names_plot <- ggplot(selected, aes(0, row_position)) +
  geom_text(
    aes(label = label_expression),
    parse = TRUE,
    hjust = 0,
    size = 2.15
  ) +
  annotate(
    "text",
    x = 0,
    y = max(selected$row_position) + 0.8,
    label = "EOCRC-enriched",
    hjust = 0,
    fontface = "bold",
    size = 2.0,
    color = "#A63A42"
  ) +
  annotate(
    "text",
    x = 0,
    y = max(selected[block_order == 2L, row_position]) + 0.8,
    label = "EOCRC-reduced",
    hjust = 0,
    fontface = "bold",
    size = 2.0,
    color = "#285D91"
  ) +
  scale_y_continuous(limits = y_limits, expand = expansion(mult = 0)) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  labs(title = "Species") +
  theme_void(base_size = 6.5) +
  theme(
    plot.title = element_text(
      size = 7,
      face = "bold",
      hjust = 0
    ),
    plot.margin = margin(1.5, 1.5, 1.5, 1.5)
  )
direction_plot <- ggplot(direction, aes(age_group, row_position)) +
  geom_hline(
    yintercept = boundary,
    color = "grey80",
    linewidth = 0.25
  ) +
  geom_point(
    data = direction[!is.na(estimate)],
    aes(
      shape = disease_direction,
      fill = disease_direction,
      alpha = supported
    ),
    size = 1.25,
    stroke = 0.15,
    color = "grey35"
  ) +
  geom_text(
    data = direction[is.na(estimate)],
    label = "QC",
    size = 1.8,
    color = "grey45"
  ) +
  scale_shape_manual(values = c("Higher in CRC" = 24, "Lower in CRC" = 25)) +
  scale_fill_manual(values = c("Higher in CRC" = "#C83E4D", "Lower in CRC" = "#3F5F8F")) +
  scale_alpha_manual(values = c("FALSE" = 0.3, "TRUE" = 1)) +
  scale_y_continuous(limits = y_limits, expand = expansion(mult = 0)) +
  labs(
    title = "Disease\ndirection",
    x = NULL,
    y = NULL
  ) +
  row_theme
age_forests <- list()

for (g in c("EO", "LO")) {
  d <- association[age_group == g]
  missing <- selected[!feature_id %in% d$feature_id]
  age_forests[[g]] <- ggplot(d, aes(estimate, row_position)) +
    geom_hline(
      yintercept = boundary,
      color = "grey80",
      linewidth = 0.25
    ) +
    geom_vline(
      xintercept = 0,
      linetype = 2,
      color = "grey65",
      linewidth = 0.25
    ) +
    geom_errorbar(
      aes(xmin = ci_low, xmax = ci_high),
      orientation = "y",
      width = 0.12,
      color = "grey45",
      linewidth = 0.25
    ) +
    geom_point(aes(color = significant_fdr), size = 1.1) +
    geom_text(aes(x = ci_high + 0.10, label = label), size = 2.25) +
    geom_text(
      data = missing,
      aes(x = 0, y = row_position),
      label = "QC fail",
      size = 1.9,
      color = "grey45"
    ) +
    scale_color_manual(values = c("FALSE" = "grey45", "TRUE" = "#C83E4D")) +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = 0)) +
    coord_cartesian(xlim = x_limits) +
    labs(
      title = paste0(g, " association"),
      x = "CLR change per mean\nmember z-score unit",
      y = NULL
    ) +
    row_theme
}

support_plot <- ggplot(selected, aes(0, row_position)) +
  geom_hline(
    yintercept = boundary,
    color = "grey80",
    linewidth = 0.25
  ) +
  geom_point(
    aes(shape = EOCRC_supported, color = EOCRC_supported),
    size = 1.6,
    stroke = 0.3
  ) +
  scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16)) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#237C50")) +
  scale_y_continuous(limits = y_limits, expand = expansion(mult = 0)) +
  coord_cartesian(xlim = c(-0.6, 0.6)) +
  labs(
    title = "EOCRC support\nq <= 0.1",
    x = NULL,
    y = NULL
  ) +
  row_theme +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
species_row <- plot_grid(
  names_plot,
  direction_plot,
  age_forests$EO,
  age_forests$LO,
  support_plot,
  nrow = 1,
  rel_widths = c(3.3, 0.85, 3.8, 3.8, 1.3),
  align = "h",
  axis = "tb"
)
main_page <- ggdraw() +
  draw_plot(
    top,
    x = 4 / 210,
    y = 111 / 297,
    width = 202 / 210,
    height = 91 / 297
  ) +
  draw_plot(
    species_row,
    x = 4 / 210,
    y = 6 / 297,
    width = 202 / 210,
    height = 100 / 297
  ) +
  draw_plot_label(
    "E",
    x = 1 / 210,
    y = 109 / 297,
    hjust = 0,
    vjust = 1,
    size = 12,
    fontface = "bold"
  )
grDevices::pdf(
  file.path(source_dir, "Fig4_B_E_A4.pdf"),
  width = page_width,
  height = page_height,
  useDingbats = FALSE
)
set.seed(20261034)
print(main_page)
grDevices::dev.off()

paper_heat <- paper_plot(supp_plots$heat) + theme(
  axis.text.y = element_text(size = 6.2),
  legend.position = "bottom"
) + guides(fill = guide_colourbar(
  barwidth = grid::unit(20, "mm"),
  barheight = grid::unit(1.6, "mm"),
  title.position = "left"
))
paper_membership <- paper_plot(supp_plots$membership_plot, point_mm = 0.7) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 6
    ),
    legend.position = "none"
  )
panel_key <- read_csv(
  file.path(v35_dir, "tables/supplementary_panel_display_key.csv"),
  show_col_types = FALSE
)
paper_key <- ggplot(
  panel_key,
  aes(
    0,
    rev(seq_len(nrow(panel_key))),
    label = paste0(display_key, "  ", str_wrap(panel_label, 26))
  )
) +
  geom_text(
    hjust = 0,
    size = 2.1,
    lineheight = 1.05
  ) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  labs(title = "Predefined panels") +
  theme_void(base_size = 6.5) +
  theme(
    plot.title = element_text(size = 7.5, face = "bold"),
    plot.margin = margin(2, 6, 2, 2)
  )
paper_A <- (paper_heat | paper_membership | paper_key) + plot_layout(widths = c(2.2, 1.1, 2.0))

score_long <- as_tibble(supp$score_long) %>%
  mutate(
    panel_label = factor(
      str_wrap(panel_label, 29),
      levels = str_wrap(panel_key$panel_label[panel_key$panel_tier == "extended_class_based"], 29)
    ),
    Group4 = factor(Group4, levels = c("LOControl", "LOCRC", "EOControl", "EOCRC"))
  )
annotations <- read_csv(
  file.path(v35_dir, "tables/supplementary_score_plot_annotations.csv"),
  show_col_types = FALSE
)
paper_B <- plot_panel_scores(
  score_long,
  annotations,
  label_style = "beta_sig",
  ncol = 3,
  beta_offset = 0.16,
  sig_offset = 0.07,
  text_size = 6.5
) + labs(y = "Mean member z-score")
paper_B <- paper_plot(paper_B, text_mm = 2.0, point_mm = 0.30) +
  scale_x_discrete(labels = c(
    LOControl = "LO\nControl",
    LOCRC = "LO\nCRC",
    EOControl = "EO\nControl",
    EOCRC = "EO\nCRC"
  )) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 6.5, face = "plain"),
    strip.text.x = element_text(size = 6.5, face = "plain"),
    axis.text.x = element_text(
      size = 5.9,
      angle = 0,
      hjust = 0.5
    )
  )

for (i in seq_along(paper_B$layers)) if (inherits(paper_B$layers[[i]]$geom, "GeomBoxplot")) {
  paper_B$layers[[i]]$aes_params$width <- 0.42
  paper_B$layers[[i]]$stat_params$width <- 0.42
}

effect_data <- as_tibble(supp_plots$effect_data) %>%
  left_join(panel_key %>%
    select(panel, display_key), by = "panel") %>%
  mutate(display_key = factor(display_key, levels = rev(paste0("E", seq_len(9)))))
span <- diff(range(c(effect_data$ci_low, effect_data$ci_high)))
paper_C <- ggplot(effect_data, aes(estimate, display_key)) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    color = "grey60",
    linewidth = 0.25
  ) +
  geom_errorbar(
    aes(xmin = ci_low, xmax = ci_high),
    orientation = "y",
    width = 0.12,
    color = "grey45",
    linewidth = 0.25
  ) +
  geom_point(color = "#C83E4D", size = 0.9) +
  geom_text(
    aes(x = max(ci_high) + span * 0.05, label = q_label),
    hjust = 0,
    size = 1.9
  ) +
  facet_wrap(~comparison, nrow = 1) +
  coord_cartesian(xlim = c(
    min(effect_data$ci_low) - span * 0.06,
    max(effect_data$ci_high) + span * 0.65
  )) +
  labs(
    x = "Mean member z-score difference",
    y = NULL,
    title = "Age-group effects"
  ) +
  theme_bw(base_size = 6.5) +
  theme(
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 6.4),
    strip.text = element_text(size = 6.5),
    plot.title = element_text(size = 7.5, face = "bold"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(2, 2, 2, 2)
  )
supplement_page <- ggdraw() +
  draw_plot(
    paper_A,
    x = 4 / 210,
    y = 153 / 297,
    width = 202 / 210,
    height = 139 / 297
  ) +
  draw_plot(
    paper_B,
    x = 4 / 210,
    y = 6 / 297,
    width = 117 / 210,
    height = 140 / 297
  ) +
  draw_plot(
    paper_C,
    x = 127 / 210,
    y = 12 / 297,
    width = 79 / 210,
    height = 130 / 297
  ) +
  draw_plot_label(
    c("A", "B", "C"),
    x = c(1, 1, 124) / 210,
    y = c(295, 150, 150) / 297,
    hjust = 0,
    vjust = 1,
    size = 12,
    fontface = "bold"
  )
grDevices::pdf(
  file.path(v35_dir, "plots/supplementary_figure_composite_v35.pdf"),
  width = page_width,
  height = page_height,
  useDingbats = FALSE
)
set.seed(20261034)
print(supplement_page)
grDevices::dev.off()
grDevices::png(
  file.path(v35_dir, "plots/supplementary_figure_composite_v35.png"),
  width = round(page_width * 300),
  height = round(page_height * 300),
  res = 300,
  type = "windows"
)
set.seed(20261034)
print(supplement_page)
grDevices::dev.off()
