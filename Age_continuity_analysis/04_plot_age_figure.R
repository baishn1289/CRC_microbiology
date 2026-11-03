library(tidyverse)
library(patchwork)
library(systemfonts)
library(ragg)

font_family <- "sans"
df_raw <- inputs$main
plot_df <- df_raw %>%
  mutate(
    row_id = rev(row_number()),
    display_label = gsub("_", " ", Species),
    HR = coef,
    CI_lower = ci_lo,
    CI_upper = ci_hi,
    Q = sprintf("%.3f", q_after_QC),
    point_fill = if_else(coef < 0, "Negative", "Positive")
  )
data_df <- plot_df
stripe_df <- plot_df %>%
  filter(row_id %% 2 == 0)
age_colors <- c(Negative = "#2563EB", Positive = "#D97706")
group_colors <- c(Control = "#007F7B", CRC = "#7B3F98")
p <- ggplot(plot_df, aes(y = row_id)) +
  geom_rect(
    data = stripe_df,
    aes(ymin = row_id - .5, ymax = row_id + .5),
    xmin = -3.35,
    xmax = 2.65,
    fill = "#F3F4F6",
    inherit.aes = FALSE
  ) +
  annotate(
    "segment",
    x = 0,
    xend = 0,
    y = .5,
    yend = 20.4,
    linetype = 2,
    linewidth = .35,
    color = "grey45"
  ) +
  geom_errorbar(
    data = data_df,
    aes(xmin = CI_lower, xmax = CI_upper),
    orientation = "y",
    width = .14,
    linewidth = .4,
    color = "#4B5563"
  ) +
  geom_point(
    data = data_df,
    aes(x = HR, fill = point_fill),
    shape = 22,
    size = 1.9,
    color = "black",
    stroke = .35
  ) +
  scale_fill_manual(values = age_colors, guide = "none") +
  geom_text(
    aes(x = -3.25, label = display_label),
    hjust = 0,
    fontface = "italic",
    family = font_family,
    size = 2.65
  ) +
  geom_text(aes(x = 1.35, label = Q), family = font_family, size = 2.6) +
  annotate(
    "text",
    x = -3.25,
    y = 20.9,
    label = "Species",
    hjust = 0,
    family = font_family,
    fontface = "bold",
    size = 3.2
  ) +
  annotate(
    "text",
    x = .02,
    y = 20.9,
    label = "beta / '10 years (95% CI)'",
    parse = TRUE,
    family = font_family,
    fontface = "bold",
    size = 3.15
  ) +
  annotate(
    "text",
    x = 1.35,
    y = 20.9,
    label = "q",
    family = font_family,
    fontface = "bold",
    size = 3.2
  ) +
  annotate(
    "text",
    x = 1.94,
    y = 20.9,
    label = "EO",
    family = font_family,
    fontface = "bold",
    size = 3.2
  ) +
  annotate(
    "text",
    x = 2.39,
    y = 20.9,
    label = "LO",
    family = font_family,
    fontface = "bold",
    size = 3.2
  ) +
  annotate(
    "segment",
    x = -3.35,
    xend = 2.65,
    y = 20.53,
    yend = 20.53,
    linewidth = .4
  ) +
  scale_x_continuous(
    limits = c(-3.35, 2.65),
    breaks = c(-.5, 0, .5, 1),
    labels = c("-0.5", "0", "0.5", "1.0"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(.35, 21.3), expand = c(0, 0)) +
  labs(
    x = NULL,
    y = NULL,
    tag = "A"
  ) +
  theme_void(base_size = 7.5, base_family = font_family) +
  theme(
    axis.text.x = element_text(size = 9, color = "black"),
    axis.ticks.x = element_line(linewidth = .25),
    axis.ticks.length = unit(.7, "mm"),
    plot.tag = element_text(size = 9.5, face = "bold"),
    plot.tag.position = "topleft",
    plot.margin = margin(2, 1, 2, 1, unit = "mm")
  )

crc_tiles <- list()

for (target in c("EO", "LO")) {
  available <- data_df[[paste0(target, "_available")]]
  disease_beta <- data_df[[paste0("original_", target, "_coef")]]
  tile_color <- ifelse(!available, "#EAECEF", ifelse(disease_beta > 0, "#F0D3AD", "#BED1E5"))
  tile_label <- ifelse(
    !available,
    "\u00d7",
    ifelse(data_df[[paste0(target, "_significant")]], "*", "")
  )
  position <- if (target == "EO") 1.94 else 2.39
  tile_data <- mutate(
    data_df,
    target = target,
    available = available,
    disease_beta = disease_beta,
    tile_x = position,
    tile_color = tile_color,
    tile_label = tile_label
  )
  p <- p + geom_tile(
    data = tile_data,
    aes(x = tile_x, y = row_id),
    width = .31,
    height = .60,
    fill = tile_color,
    color = "white"
  ) +
    geom_text(
      data = tile_data,
      aes(
        x = tile_x,
        y = row_id,
        label = tile_label
      ),
      family = font_family,
      size = 2.8,
      color = "#182535"
    )
  crc_tiles[[target]] <- select(
    tile_data,
    feature,
    Species,
    target,
    row_id,
    available,
    disease_beta,
    tile_x,
    tile_color,
    tile_label
  )
}

pair <- inputs$slopes %>%
  mutate(
    row_id = 21 - Figure_row,
    Group = factor(Group, levels = c("Control", "CRC")),
    y = row_id + if_else(Group == "Control", .17, -.17),
    point_fill = if_else(qc_pass, as.character(Group), "Flag")
  )
p_B <- ggplot(pair, aes(
  x = beta,
  y = y,
  color = Group
)) +
  geom_rect(
    data = stripe_df,
    aes(ymin = row_id - .5, ymax = row_id + .5),
    xmin = -1.05,
    xmax = 1.65,
    fill = "#F3F4F6",
    inherit.aes = FALSE
  ) +
  annotate(
    "segment",
    x = 0,
    xend = 0,
    y = .5,
    yend = 20.4,
    linetype = 2,
    color = "grey45",
    linewidth = .35
  ) +
  geom_errorbar(
    aes(xmin = ci_lo, xmax = ci_hi),
    orientation = "y",
    width = .12,
    linewidth = .4
  ) +
  geom_point(aes(fill = point_fill, shape = Group), size = 1.9, stroke = .4) +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_fill_manual(values = c(group_colors, Flag = "white"), guide = "none") +
  scale_shape_manual(values = c(Control = 22, CRC = 23), guide = "none") +
  annotate(
    "text",
    x = .30,
    y = 20.9,
    label = "beta / '10 years (95% CI)'",
    parse = TRUE,
    family = font_family,
    fontface = "bold",
    size = 9.5 / ggplot2::.pt
  ) +
  annotate(
    "segment",
    x = -1.05,
    xend = 1.65,
    y = 20.53,
    yend = 20.53,
    color = "black",
    linewidth = .4
  ) +
  scale_x_continuous(
    limits = c(-1.05, 1.65),
    breaks = c(-1, 0, 1),
    labels = c("-1", "0", "1"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(.35, 21.3), expand = c(0, 0)) +
  labs(
    x = NULL,
    y = NULL,
    tag = "B"
  ) +
  theme_void(base_size = 7.5, base_family = font_family) +
  theme(
    axis.text.x = element_text(size = 9, color = "black"),
    axis.ticks.x = element_line(linewidth = .25),
    axis.ticks.length = unit(.7, "mm"),
    plot.tag = element_text(size = 9.5, face = "bold"),
    plot.tag.position = "topleft",
    plot.margin = margin(2, 1, 2, 1, unit = "mm")
  )

legend_labels <- c("Filled: QC passed", "Open: QC not passed", "CRC", "Control")
legend_size_pt <- 8.1
legend_text_width_mm <- vapply(
  seq_along(legend_labels),
  function(j) {
    systemfonts::string_width(
      legend_labels[j],
      family = font_family,
      weight = if (j >= 3) "bold" else "normal",
      size = legend_size_pt,
      res = 720
    ) * 25.4 / 720
  },
  numeric(1)
)
legend_key_width_mm <- c(7.25, 7.25, 4.75, 4.75)
legend_item_width_mm <- legend_key_width_mm + legend_text_width_mm
legend_gap_mm <- 9
legend_total_width_mm <- sum(legend_item_width_mm) + 3 * legend_gap_mm
legend_item_x <- (190 - legend_total_width_mm) / 2 +
  c(0, cumsum(head(legend_item_width_mm + legend_gap_mm, -1)))
legend_text_x <- legend_item_x + legend_key_width_mm
legend_B <- ggplot() +
  annotate(
    "point",
    x = legend_item_x[1] + 1.25,
    y = .5,
    shape = 22,
    size = 2,
    stroke = .4,
    color = "#444444",
    fill = "#444444"
  ) +
  annotate(
    "point",
    x = legend_item_x[1] + 4.25,
    y = .5,
    shape = 23,
    size = 2,
    stroke = .4,
    color = "#444444",
    fill = "#444444"
  ) +
  annotate(
    "text",
    x = legend_text_x[1],
    y = .5,
    label = legend_labels[1],
    hjust = 0,
    family = font_family,
    color = "#333333",
    size = legend_size_pt / ggplot2::.pt
  ) +
  annotate(
    "point",
    x = legend_item_x[2] + 1.25,
    y = .5,
    shape = 22,
    size = 2,
    stroke = .4,
    color = "#444444",
    fill = "white"
  ) +
  annotate(
    "point",
    x = legend_item_x[2] + 4.25,
    y = .5,
    shape = 23,
    size = 2,
    stroke = .4,
    color = "#444444",
    fill = "white"
  ) +
  annotate(
    "text",
    x = legend_text_x[2],
    y = .5,
    label = legend_labels[2],
    hjust = 0,
    family = font_family,
    color = "#333333",
    size = legend_size_pt / ggplot2::.pt
  ) +
  annotate(
    "point",
    x = legend_item_x[3] + 1.25,
    y = .5,
    shape = 23,
    size = 2.1,
    stroke = .4,
    color = group_colors[["CRC"]],
    fill = group_colors[["CRC"]]
  ) +
  annotate(
    "text",
    x = legend_text_x[3],
    y = .5,
    label = legend_labels[3],
    hjust = 0,
    family = font_family,
    color = group_colors[["CRC"]],
    fontface = "bold",
    size = legend_size_pt / ggplot2::.pt
  ) +
  annotate(
    "point",
    x = legend_item_x[4] + 1.25,
    y = .5,
    shape = 22,
    size = 2.1,
    stroke = .4,
    color = group_colors[["Control"]],
    fill = group_colors[["Control"]]
  ) +
  annotate(
    "text",
    x = legend_text_x[4],
    y = .5,
    label = legend_labels[4],
    hjust = 0,
    family = font_family,
    color = group_colors[["Control"]],
    fontface = "bold",
    size = legend_size_pt / ggplot2::.pt
  ) +
  scale_x_continuous(limits = c(0, 190), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void(base_family = font_family) +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "mm"))
forest_row <- p + p_B + plot_layout(widths = c(.69, .31))
combined <- (forest_row / wrap_elements(full = legend_B)) + plot_layout(heights = c(1, .045)) +
  plot_annotation(theme = theme(plot.margin = margin(10, 1, 10, 1, unit = "mm")))
stem <- file.path(output_dir, "Figure_S_age_context")
ggsave(
  paste0(stem, ".pdf"),
  combined,
  width = 192,
  height = 232,
  units = "mm",
  device = grDevices::pdf,
  family = "Helvetica",
  useDingbats = FALSE
)
ggsave(
  paste0(stem, ".png"),
  combined,
  width = 192,
  height = 232,
  units = "mm",
  device = ragg::agg_png,
  dpi = 300,
  bg = "white"
)
