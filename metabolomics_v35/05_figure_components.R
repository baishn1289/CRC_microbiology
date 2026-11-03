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
ctx <- readRDS(file.path(output_dir, "objects/fixed_scheme_context.rds"))
stats <- readRDS(file.path(output_dir, "objects/fixed_scheme_statistics.rds"))
supp <- readRDS(file.path(output_dir, "objects/supplementary_statistics.rds"))
theme_set(theme_bw(base_size = 12))

# Build the retained main-figure panels without exporting alternate layouts.
fdr_cutoff <- 0.1
panel_order <- c("U_P", "D0")
panel_labels <- c(
  U_P = "Arginine-ornithine-proline\nmetabolites",
  D0 = "Aminovalerate-polyamine\nmetabolites"
)
effect_order <- c("LO", "EO", "EO_minus_LO")
effect_labels <- c(
  LO = "LO effect",
  EO = "EO effect",
  EO_minus_LO = "EO - LO effect"
)
group_levels <- c("LOControl", "LOCRC", "EOControl", "EOCRC")
group_palette <- c(
  LOControl = "#B8C7DA",
  LOCRC = "#3F5F8F",
  EOControl = "#E8B8B0",
  EOCRC = "#C83E4D"
)
mean_df <- as_tibble(stats$mean_results) %>%
  mutate(
    panel_label = factor(unname(panel_labels[panel]), levels = unname(panel_labels)),
    contrast_label = factor(
      unname(effect_labels[contrast]),
      levels = rev(unname(effect_labels[effect_order]))
    ),
    conf.low = ci_low,
    conf.high = ci_high
  )
member_df <- as_tibble(stats$member_results)
score_long <- as_tibble(stats$score_long) %>%
  mutate(
    panel_label = factor(unname(panel_labels[panel]), levels = unname(panel_labels)),
    Group4 = factor(Group4, levels = group_levels)
  )
member_long <- as_tibble(stats$member_long) %>%
  mutate(
    panel_label = factor(unname(panel_labels[panel]), levels = unname(panel_labels)),
    Group4 = factor(Group4, levels = group_levels)
  )
membership_order <- unlist(ctx$members, use.names = FALSE)
main_forest <- plot_main_figure_panel_effects(mean_df) +
  geom_text(
    data = mean_df,
    aes(
      x = ci_high + 0.025,
      label = sprintf("q=%.3f%s", q_BH6, ifelse(significant_fdr, "*", ""))
    ),
    hjust = 0,
    size = 3.7,
    color = "grey15"
  ) +
  labs(
    x = "Adjusted mean member z-score difference",
    caption = "OLS + HC3, asymptotic 95% CI; * BH over 6 panel effects, q <= 0.1"
  ) +
  theme(
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    plot.caption = element_text(size = 9)
  )
main_forest$coordinates$limits$x <- c(min(mean_df$ci_low) - 0.04, max(mean_df$ci_high) + 0.25)
panel_ranges <- score_long %>%
  group_by(panel, panel_label) %>%
  summarise(
    y_max = max(panel_score),
    y_min = min(panel_score),
    .groups = "drop"
  ) %>%
  mutate(y_span = pmax(y_max - y_min, 0.35))
annotation <- bind_rows(
  mean_df %>%
    filter(contrast %in% c("EO", "LO")) %>%
    transmute(
      panel,
      contrast,
      estimate,
      q = q_BH6,
      x_start = if_else(contrast == "EO", 3, 1),
      x_end = if_else(contrast == "EO", 4, 2),
      comparison = "age-adjusted disease effect"
    ),
  as_tibble(stats$crc) %>%
    filter(model == "clinical_binary") %>%
    transmute(
      panel,
      contrast,
      estimate,
      q = q_BH2,
      x_start = 2,
      x_end = 4,
      comparison = "CRC-only Sex/BMI/Stage2/Site2 adjusted effect"
    )
) %>%
  left_join(panel_ranges, by = "panel") %>%
  mutate(
    x_mid = (x_start + x_end) / 2,
    y_position = y_max + y_span * if_else(contrast == "EOCRC_minus_LOCRC", 0.35, 0.15),
    beta_label = sprintf("beta=%+.2f", estimate),
    sig_label_plot = sprintf("q=%.3f%s", q, ifelse(q <= 0.1, "*", ""))
  )
set.seed(20261033)
score_plot <- plot_panel_scores(
  score_long,
  annotation,
  label_style = "no_sig_only",
  ncol = 1,
  text_size = 12
) +
  labs(y = "Mean member z-score")
canonical_df <- member_df %>%
  mutate(
    panel_label = factor(unname(panel_labels[panel]), levels = unname(panel_labels)),
    contrast_label = factor(
      unname(effect_labels[contrast]),
      levels = unname(effect_labels[effect_order])
    ),
    canonical_metabolite = factor(canonical_metabolite, levels = rev(membership_order)),
    mean_logFC = estimate,
    effect_label = paste0(sprintf("%.2f", estimate), ifelse(significant_fdr, "*", ""))
  )
canonical_df <- canonical_df %>%
  mutate(
    max_abs = max(abs(estimate)),
    text_color = if_else(abs(estimate) > max_abs * 0.62, "white", "black")
  )
canonical_plot <- plot_main_figure_canonical_effects(canonical_df) +
  facet_wrap(
    ~panel_label,
    ncol = 1,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_y_discrete(expand = expansion(mult = 0, add = 0)) +
  labs(
    fill = "Member z-score\neffect",
    caption = "Members: OLS + HC3; U_P BH15 (5 x 3); D0 BH9 (3 x 3); * q <= 0.1"
  ) +
  theme(
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12),
    plot.caption = element_text(size = 9)
  )
main <- list(
  main_forest = main_forest,
  score_plot = score_plot,
  canonical_plot = canonical_plot
)

# Supplementary heatmap shows the original member-z group means.
member_z <- as_tibble(t(supp$canonical_z), rownames = "Sample_ID") %>%
  pivot_longer(-Sample_ID, names_to = "canonical_metabolite", values_to = "member_z") %>%
  left_join(as_tibble(ctx$meta), by = "Sample_ID")
group_means <- member_z %>%
  group_by(canonical_metabolite, Group4) %>%
  summarise(
    mean_member_z = mean(member_z),
    n = n(),
    .groups = "drop"
  )
panel_key <- as_tibble(supp$lookup) %>%
  distinct(panel, panel_label, panel_tier) %>%
  mutate(display_key = c("U_P", "D0", paste0("E", seq_len(9))))
membership <- as_tibble(supp$lookup) %>%
  left_join(panel_key %>%
    select(panel, display_key), by = "panel")
row_order <- supp$canonical_order
membership <- membership %>%
  mutate(
    canonical_metabolite = factor(canonical_metabolite, levels = rev(row_order)),
    display_key = factor(display_key, levels = panel_key$display_key)
  )
group_means <- group_means %>%
  mutate(
    canonical_metabolite = factor(canonical_metabolite, levels = rev(row_order)),
    Group4 = factor(Group4, levels = group_levels)
  )
fill_limit <- max(abs(group_means$mean_member_z))
heat <- ggplot(group_means, aes(Group4, canonical_metabolite, fill = mean_member_z)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    limits = c(-fill_limit, fill_limit)
  ) +
  scale_y_discrete(expand = expansion(mult = 0)) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Mean member\nz-score",
    title = "33 unique predefined metabolites"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1),
    axis.text.y = element_text(size = 11),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 14)
  )
membership_plot <- ggplot(membership, aes(display_key, canonical_metabolite)) +
  geom_tile(
    fill = "grey98",
    color = "white",
    width = 0.94,
    height = 0.97
  ) +
  geom_point(size = 1.9, color = "grey25") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(expand = expansion(mult = 0)) +
  labs(
    x = NULL,
    y = NULL,
    title = "Panel membership"
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 10
    ),
    plot.title = element_text(size = 14, face = "bold"),
    panel.border = element_blank()
  )

score_long <- as_tibble(supp$score_long) %>%
  mutate(
    panel_label = factor(str_wrap(panel_label, 29), levels = str_wrap(unique(panel_label), 29)),
    Group4 = factor(Group4, levels = group_levels)
  )
panel_ranges <- score_long %>%
  group_by(panel, panel_label) %>%
  summarise(
    y_max = max(panel_score),
    y_min = min(panel_score),
    .groups = "drop"
  ) %>%
  mutate(y_span = pmax(y_max - y_min, 0.35))
annotations <- as_tibble(supp$effects) %>%
  filter(contrast %in% c("EO", "LO")) %>%
  transmute(
    panel,
    contrast,
    estimate,
    p_value,
    q = q_supp_panel_BH27,
    family = "supplementary panel BH27",
    x_start = if_else(contrast == "EO", 3, 1),
    x_end = if_else(contrast == "EO", 4, 2)
  ) %>%
  left_join(panel_ranges, by = "panel") %>%
  mutate(
    x_mid = (x_start + x_end) / 2,
    y_position = y_max + y_span * 0.13,
    beta_label = sprintf("beta=%+.2f", estimate),
    sig_label_plot = paste0(
      if_else(q < 0.001, "q<0.001", sprintf("q=%.3f", q)),
      ifelse(q <= 0.1, "*", "")
    )
  )

effect_data <- bind_rows(
  as_tibble(supp$effects) %>%
    filter(contrast == "EO_minus_LO") %>%
    transmute(
      panel,
      estimate,
      ci_low,
      ci_high,
      p_value,
      q = q_supp_panel_BH27,
      comparison = "Disease EO - LO\nBH27"
    ),
  as_tibble(supp$crc) %>%
    filter(model == "clinical_binary") %>%
    transmute(
      panel,
      estimate,
      ci_low,
      ci_high,
      p_value,
      q = q_supp_CRC_BH9,
      comparison = "CRC EO - LO\nBH9"
    )
) %>%
  left_join(panel_key %>%
    select(panel, panel_label), by = "panel") %>%
  mutate(
    panel_label = factor(
      str_wrap(panel_label, 23),
      levels = rev(str_wrap(supp$score_long$panel_label %>%
        unique(), 23))
    ),
    q_label = paste0(if_else(q < 0.001, "q<0.001", sprintf("q=%.3f", q)), ifelse(q <= 0.1, "*", ""))
  )

write_csv(group_means, file.path(output_dir, "tables/supplementary_member_group_means.csv"))
write_csv(membership, file.path(output_dir, "tables/supplementary_heatmap_membership.csv"))
write_csv(panel_key, file.path(output_dir, "tables/supplementary_panel_display_key.csv"))
write_csv(
  annotations,
  file.path(output_dir, "tables/supplementary_score_plot_annotations.csv")
)
supp_plots <- list(
  heat = heat,
  membership_plot = membership_plot,
  effect_data = effect_data
)
saveRDS(
  list(main = main, supp_plots = supp_plots),
  file.path(output_dir, "objects/figure_components.rds")
)
