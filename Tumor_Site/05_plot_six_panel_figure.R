library(dplyr)
library(ggplot2)
library(grid)
library(ragg)

theme_set(theme_classic(base_size = 7, base_family = "sans") + theme(
  axis.line = element_line(linewidth = .3),
  axis.ticks = element_line(linewidth = .3),
  axis.text = element_text(colour = "black", size = 6.8),
  axis.title = element_text(size = 7),
  strip.background = element_blank(),
  strip.text = element_text(size = 7),
  panel.spacing = unit(1.5, "mm"),
  plot.margin = margin(0, 1, 0, 1, "mm")
))


matrix_rows <- coverage_matrix %>%
  mutate(
    x = .37 + (col_index - 1) * .155,
    y = 4 - row_index,
    cell_colour = ifelse(
      site_record == "Colon_unspecified",
      "#ECEDEF",
      ifelse(n > 0, "#E8F0F4", "white")
    ),
    label = ifelse(n > 0, as.character(n), "-")
  )
p_a <- ggplot(matrix_rows, aes(x, y)) +
  geom_tile(
    aes(fill = cell_colour),
    width = .145,
    height = .79,
    colour = "white",
    linewidth = .4
  ) +
  scale_fill_identity() +
  geom_text(aes(label = label), size = 2.7) +
  annotate(
    "text",
    x = .005,
    y = c(3, 2, 1),
    label = c("FUSCC-SHSD", "Yachida 2019", "Yang 2020"),
    hjust = 0,
    size = 2.35
  ) +
  annotate(
    "text",
    x = c(.37, .525, .68, .835),
    y = 4.1,
    label = c("Right\ncolon", "Left\ncolon", "Colon,\nunspecified", "Rectum"),
    size = 2.3
  ) +
  annotate(
    "text",
    x = .005,
    y = -.15,
    label = "Rectum - colon: 3 cohorts",
    hjust = 0,
    size = 2.35
  ) +
  annotate(
    "text",
    x = .005,
    y = -.95,
    label = "Sidedness and ternary pairs: 2 cohorts",
    hjust = 0,
    size = 2.35
  ) +
  annotate(
    "text",
    x = .005,
    y = -1.65,
    label = "Left side includes rectum",
    hjust = 0,
    size = 2.2
  ) +
  scale_x_continuous(limits = c(0, .94), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-2.1, 4.9), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

eo_order <- c(
  "Dialister pneumosintes",
  "Haemophilus parainfluenzae",
  "Schaalia odontolytica",
  "Veillonella rogosae",
  "Ruthenibacterium lactatiformans",
  "Dysosmobacter welbionis",
  "Prevotella stercorea"
)
b_rows <- eo_hits %>%
  arrange(match(taxon, eo_order)) %>%
  mutate(
    row_id = c(10.2, 9.2, 8.2, 5.8, 4.8, 3.8, 1.4),
    kind = "pooled",
    q_label = sprintf("%.4f", q_comparison)
  )
b_label <- ggplot(b_rows, aes(x = .01, y = row_id)) +
  geom_text(
    aes(label = taxon),
    hjust = 0,
    fontface = "italic",
    size = 2.35
  ) +
  annotate(
    "text",
    x = .01,
    y = c(11.35, 6.95, 2.55),
    label = c("Left side - right side", "Left colon - right colon", "Rectum - left colon"),
    hjust = 0,
    fontface = "bold",
    size = 2.35
  ) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(.3, 12.2), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))
b_core <- plotcase_forest(b_rows, c(-5.5, 6.2), c(-4, 0, 4), c(.3, 12.2))
b_q <- ggplot(b_rows, aes(x = .5, y = row_id)) +
  geom_text(aes(label = q_label), size = 2.35) +
  annotate(
    "text",
    x = .5,
    y = 11.35,
    label = "q_cmp",
    fontface = "bold",
    size = 2.4
  ) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(.3, 12.2), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

eo_taxa <- c("Dialister pneumosintes", "Haemophilus parainfluenzae", "Veillonella rogosae")
eo_core <- raw_panels <- list()
site_colours <- c(Right = "#D55E00", Left = "#0072B2")

for (taxon in eo_taxa) {
  d <- forest_rows %>%
    filter(taxon == .env$taxon) %>%
    mutate(row_id = 3:1)
  eo_core[[taxon]] <- plotcase_forest(d, c(-8, 2), c(-8, -4, 0, 2), c(.4, 3.8))
  levels_site <- if (taxon == "Veillonella rogosae")
    c("Right colon", "Left colon") else
    c("Right-sided", "Left-sided")
  raw <- raw_plot %>%
    filter(taxon == .env$taxon) %>%
    mutate(
      cohort_label = factor(cohort_label, levels = c("FUSCC", "Yachida 2019")),
      site_label = factor(site_label, levels = levels_site)
    )
  det <- detection_plot %>%
    filter(taxon == .env$taxon) %>%
    mutate(
      cohort_label = factor(cohort_label, levels = c("FUSCC", "Yachida 2019")),
      site_label = factor(site_label, levels = levels_site),
      label = paste0(n_detected, "/", n)
    )
  raw_panels[[taxon]] <- ggplot(raw, aes(site_label, display, colour = side)) +
    geom_boxplot(
      width = .45,
      outlier.shape = NA,
      linewidth = .35,
      fill = "white"
    ) +
    geom_point(
      position = position_jitter(
        width = .13,
        height = 0,
        seed = 20261016
      ),
      alpha = .6,
      size = .6
    ) +
    geom_hline(
      yintercept = log10(.001 + .0001),
      linetype = 3,
      linewidth = .25,
      colour = "grey65"
    ) +
    geom_text(
      data = det,
      aes(y = .85, label = label),
      colour = "black",
      size = 2.3
    ) +
    facet_wrap(~cohort_label, nrow = 1) +
    scale_colour_manual(values = site_colours, guide = "none") +
    scale_x_discrete(labels = c(
      "Right-sided" = "Right-\nsided",
      "Left-sided" = "Left-\nsided",
      "Right colon" = "Right\ncolon",
      "Left colon" = "Left\ncolon"
    )) +
    scale_y_continuous(
      limits = c(-4.15, 1.05),
      breaks = c(-4, -2, 0),
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = if (taxon == "Dialister pneumosintes") "log10(RA% + 0.0001)" else NULL) +
    theme(
      axis.text.x = element_text(size = 6.8),
      axis.title.y = element_text(size = 7),
      strip.text = element_text(size = 7)
    )
}

eo_labels <- label_column(eo_core[[1]]$data, c(.4, 3.8))
c_columns <- eo_hits %>%
  filter(taxon %in% eo_taxa) %>%
  arrange(match(taxon, eo_taxa)) %>%
  transmute(
    taxon,
    scope,
    comparison = c("Left side - right side", "Left side - right side", "Left colon - right colon"),
    n_samples = as.integer(n_samples),
    q_comparison
  )

lo_taxa <- c(
  "Veillonella parvula",
  "Veillonella dispar",
  "Veillonella atypica",
  "Clostridium sp AF20 17LB"
)
ld <- forest_rows %>%
  filter(panel == "f")
lo_limits <- c(floor(min(ld$ci_low) / 2) * 2, ceiling(max(ld$ci_high) / 2) * 2)
lo_cores <- list()

for (taxon in lo_taxa) {
  d <- ld %>%
    filter(taxon == .env$taxon) %>%
    mutate(row_id = 9:1)
  lo_cores[[taxon]] <- plotcase_forest(d, lo_limits, c(-8, -4, 0, 4, 8), c(.4, 9.8))
}

lo_labels <- label_column(lo_cores[[1]]$data, c(.4, 9.8), 2.4)
af20_labels <- label_column(lo_cores[[4]]$data, c(.4, 9.8), 2.4)

g_rows <- forest_rows %>%
  filter(panel == "g") %>%
  mutate(row_id = 6:1)
g_label <- label_column(g_rows, c(.4, 6.8), 2.4)
g_core <- plotcase_forest(g_rows, c(-5, 9.2), c(-4, 0, 4, 8), c(.4, 6.8))
figure_width <- 183
figure_height <- 245
plots <- list(
  a = p_a,
  b_label = b_label,
  b_core = b_core,
  b_q = b_q,
  eo_labels = eo_labels,
  lo_labels = lo_labels,
  af20_labels = af20_labels,
  g_label = g_label,
  g_core = g_core
)
positions <- data.frame(
  name = names(plots),
  x = c(.018, .461, .771, .945, .018, .018, .018, .525, .692),
  y = c(.785, .785, .785, .785, .642, .252, .028, .028, .028),
  width = c(.413, .303, .168, .043, .116, .156, .155, .160, .287),
  height = c(.185, .185, .185, .185, .068, .153, .155, .155, .155)
)

for (i in 1:3) {
  taxon <- eo_taxa[i]
  left <- c(.144, .424, .704)[i]
  nm <- paste0(c("eo_core_", "raw_"), i)
  plots[nm] <- list(eo_core[[taxon]], raw_panels[[taxon]])
  positions <- bind_rows(
    positions,
    data.frame(
      name = nm,
      x = c(left, left),
      y = c(.642, .480),
      width = c(.275, .275),
      height = c(.068, .155)
    )
  )
}

for (i in 1:4) {
  nm <- paste0("lo_core_", i)
  plots[[nm]] <- lo_cores[[i]]
  positions <- bind_rows(
    positions,
    data.frame(
      name = nm,
      x = if (i <= 3) c(.185, .455, .725)[i] else .186,
      y = if (i <= 3) .252 else .028,
      width = if (i <= 3) .254 else .294,
      height = if (i <= 3) .153 else .155
    )
  )
}

pdf(
  NULL,
  width = figure_width / 25.4,
  height = figure_height / 25.4,
  family = "Helvetica",
  useDingbats = FALSE
)
plot_grobs <- lapply(plots, ggplotGrob)

for (i in 2:3) {
  plot_grobs[[paste0("eo_core_", i)]]$heights <- plot_grobs$eo_core_1$heights
  plot_grobs[[paste0("eo_core_", i)]]$widths <- plot_grobs$eo_core_1$widths
}

for (pair in list(
  c("b_label", "b_core"),
  c("b_q", "b_core"),
  c("lo_labels", "lo_core_1"),
  c("af20_labels", "lo_core_4"),
  c("g_label", "g_core")
)) {
  plot_grobs[[pair[1]]]$heights <- plot_grobs[[pair[2]]]$heights
}

plot_grobs$eo_labels$heights <- plot_grobs$eo_core_1$heights
figure_parts <- vector("list", nrow(positions))
names(figure_parts) <- positions$name

for (i in seq_len(nrow(positions))) {
  box <- positions[i, ]
  grob <- plot_grobs[[box$name]]
  grob$name <- box$name
  grob$vp <- viewport(
    name = paste0("vp6_", box$name),
    x = box$x,
    y = box$y,
    width = box$width,
    height = box$height,
    just = c("left", "bottom")
  )
  figure_parts[[i]] <- grob
}

panel_letters <- LETTERS[1:6]
headers <- list()
headers <- c(
  headers,
  panel_title("A", "EOCRC sample coverage", .018, .985),
  panel_title("B", "EOCRC: exploratory associations", .461, .985),
  panel_title("C", "EOCRC", .018, .758, size = 7.5)
)

for (i in 1:3) {
  row <- c_columns[i, ]
  center <- c(.2815, .5615, .8415)[i]
  headers <- c(
    headers,
    list(
      textGrob(
        row$taxon,
        x = center,
        y = .758,
        just = c("center", "top"),
        gp = gpar(fontsize = 7.7, fontface = "italic")
      ),
      textGrob(
        row$comparison,
        x = center,
        y = .741,
        just = c("center", "top"),
        gp = gpar(fontsize = 7)
      ),
      textGrob(
        sprintf("n = %d; q_cmp = %.4f", row$n_samples, row$q_comparison),
        x = center,
        y = .724,
        just = c("center", "top"),
        gp = gpar(fontsize = 7)
      )
    )
  )
}

headers <- c(
  headers,
  list(
    textGrob(
      "D",
      x = .018,
      y = .469,
      just = c("left", "top"),
      gp = gpar(fontsize = 9, fontface = "bold")
    ),
    textGrob(
      "LOCRC: left side - right side (8 cohorts; 695 cases)",
      x = .018,
      y = .451,
      just = c("left", "top"),
      gp = gpar(fontsize = 7.3)
    )
  )
)

for (i in 1:3) {
  taxon <- lo_taxa[i]
  center <- c(.312, .582, .852)[i]
  q <- lo_hits$q_comparison[lo_hits$taxon == taxon]
  headers <- c(
    headers,
    list(
      textGrob(
        taxon,
        x = center,
        y = .431,
        just = c("center", "top"),
        gp = gpar(fontsize = 7.7, fontface = "italic")
      ),
      textGrob(
        sprintf("q_cmp = %.3g", q),
        x = center,
        y = .414,
        just = c("center", "top"),
        gp = gpar(fontsize = 7)
      )
    )
  )
}

af20 <- lo_hits[lo_hits$taxon == lo_taxa[4], ]
headers <- c(
  headers,
  panel_title("E", lo_taxa[4], .018, .228, TRUE),
  list(
    textGrob(
      "LOCRC: left side - right side",
      x = .039,
      y = .210,
      just = c("left", "top"),
      gp = gpar(fontsize = 7)
    ),
    textGrob(
      sprintf("q_cmp = %.4f; I2 = %.1f%%", af20$q_comparison, af20$I2),
      x = .039,
      y = .193,
      just = c("left", "top"),
      gp = gpar(fontsize = 7)
    )
  ),
  panel_title("F", "Colon cases: EO - LO", .525, .228),
  list(
    textGrob(
      "Clostridiales bacterium KLE1615",
      x = .546,
      y = .210,
      just = c("left", "top"),
      gp = gpar(fontsize = 7.5, fontface = "italic")
    ),
    textGrob(
      sprintf("5 cohorts; 408 cases; q_cmp = %.4f", age_hit$q_comparison),
      x = .546,
      y = .193,
      just = c("left", "top"),
      gp = gpar(fontsize = 7)
    )
  )
)
figure <- grobTree(children = do.call(gList, c(figure_parts, headers)))
invisible(dev.off())

stem <- file.path(output_dir, "Supplementary_Figure_Tumor_Site_six_panels")
pdf(
  paste0(stem, ".pdf"),
  width = figure_width / 25.4,
  height = figure_height / 25.4,
  family = "Helvetica",
  useDingbats = FALSE
)
grid.newpage()
grid.draw(figure)
dev.off()
ragg::agg_png(
  paste0(stem, ".png"),
  width = figure_width,
  height = figure_height,
  units = "mm",
  res = 600
)
grid.newpage()
grid.draw(figure)
dev.off()
