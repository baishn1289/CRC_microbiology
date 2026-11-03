library(tidyverse)
library(patchwork)
library(sysfonts)
library(showtext)
library(ragg)

# Retain the final clean figure layout and native PDF device.
figure_width <- 10
figure_height <- 12
font_add(
  "S1Arial",
  regular = "C:/Windows/Fonts/arial.ttf",
  bold = "C:/Windows/Fonts/arialbd.ttf",
  italic = "C:/Windows/Fonts/ariali.ttf",
  bolditalic = "C:/Windows/Fonts/arialbi.ttf"
)
theme_set(theme_classic(base_size = 12, base_family = "S1Arial"))
figure_theme <- theme(
  axis.line = element_line(linewidth = 0.35),
  axis.ticks = element_line(linewidth = 0.3),
  axis.text = element_text(size = 12, colour = "black"),
  axis.title = element_text(size = 13.2),
  plot.title = element_blank(),
  plot.subtitle = element_blank(),
  plot.caption = element_blank(),
  strip.background = element_blank(),
  strip.text = element_text(size = 12.2, face = "bold"),
  panel.spacing.x = unit(15, "pt"),
  plot.tag = element_text(size = 17, face = "bold"),
  plot.margin = margin(10, 12, 17, 9)
)

age_labels <- c(EO = "Early-onset (EO)", LO = "Late-onset (LO)")

alpha <- bind_rows(lapply(
  c("EO", "LO"),
  function(age) {
    fit <- ecology$alpha[[age]]
    bind_rows(fit$cohort_res %>%
      mutate(I2 = NA_real_), fit$meta_res) %>%
      mutate(
        Age_class = age,
        stars = case_when(
          pval < 0.001 ~ "***",
          pval < 0.01 ~ "**",
          pval < 0.05 ~ "*",
          TRUE ~ ""
        )
      )
  }
))
alpha <- alpha %>%
  group_by(Age_class) %>%
  arrange(Type != "Summary", coef, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(
    label = if_else(Type == "Summary", "Summary (REML)", Cohort),
    row_key = paste(Age_class, label, sep = "___")
  )
alpha$row_key <- factor(alpha$row_key, levels = alpha$row_key)
p_alpha <- ggplot(alpha, aes(x = coef, y = row_key)) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    linewidth = 0.35,
    colour = "grey55"
  ) +
  geom_errorbar(
    aes(xmin = lower, xmax = upper),
    orientation = "y",
    width = 0.14,
    linewidth = 0.4,
    colour = "grey30"
  ) +
  geom_point(
    data = filter(alpha, Type == "Individual"),
    shape = 21,
    fill = "grey35",
    colour = "black",
    size = 1.9,
    stroke = 0.35
  ) +
  geom_point(
    data = filter(alpha, Type == "Summary"),
    shape = 23,
    fill = "#C9252D",
    colour = "black",
    size = 3.5,
    stroke = 0.5
  ) +
  geom_text(
    data = filter(alpha, Type == "Individual", stars != ""),
    aes(x = upper, label = stars),
    nudge_x = 0.06,
    hjust = 0,
    size = 4.2,
    family = "S1Arial"
  ) +
  facet_wrap(
    ~Age_class,
    nrow = 1,
    scales = "free",
    labeller = as_labeller(age_labels)
  ) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    x = "Adjusted Shannon difference (CRC - Control)",
    y = NULL,
    tag = "A"
  ) +
  figure_theme +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 10.8)
  )

beta <- ecology$beta %>%
  select(-exploratory) %>%
  mutate(
    point_class = if_else(Type == "Summary", "Global", "Cohort"),
    label = if_else(Type == "Summary", "Global", Cohort)
  ) %>%
  group_by(Age_class) %>%
  arrange(Type != "Summary", coef, .by_group = TRUE) %>%
  ungroup() %>%
  mutate(row_key = paste(Age_class, label, sep = "___"))
beta$row_key <- factor(beta$row_key, levels = beta$row_key)
beta$point_class <- factor(beta$point_class, levels = c("Cohort", "Global"))
p_beta <- ggplot(beta, aes(x = coef, y = row_key)) +
  geom_hline(
    yintercept = 1.5,
    linetype = 2,
    colour = "grey80",
    linewidth = 0.3
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = coef,
      yend = row_key
    ),
    colour = "grey68",
    linewidth = 0.4
  ) +
  geom_point(
    aes(
      shape = point_class,
      fill = point_class,
      colour = point_class
    ),
    size = 2.6,
    stroke = 0.55
  ) +
  facet_wrap(
    ~Age_class,
    nrow = 1,
    scales = "free_y",
    labeller = as_labeller(age_labels)
  ) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_x_continuous(
    limits = c(0, ceiling((max(beta$coef) + 0.025) / 0.05) * 0.05),
    breaks = seq(0, 0.5, 0.05),
    labels = function(x) sprintf("%.0f%%", 100 * x),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_shape_manual(values = c(Cohort = 21, Global = 23), guide = "none") +
  scale_fill_manual(values = c(Cohort = "grey25", Global = "#C9252D"), guide = "none") +
  scale_colour_manual(values = c(Cohort = "black", Global = "black"), guide = "none") +
  labs(
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    x = "Disease-status R\u00b2",
    y = NULL,
    tag = "B"
  ) +
  figure_theme +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 10.8),
    legend.position = "none"
  )

disp <- ecology$dispersion %>%
  mutate(x_key = paste(Age_class, Group, sep = "___"))
ds <- ecology$dispersion_stats
x_labels <- setNames(
  rep(c("Control", "CRC"), length(ds$Age_class)),
  c(rbind(
    paste(ds$Age_class, "Control", sep = "___"),
    paste(ds$Age_class, "CRC", sep = "___")
  ))
)
disp$x_key <- factor(disp$x_key, levels = names(x_labels))
y_top <- max(disp$Distance) + 0.035
p_disp <- ggplot(disp, aes(
  x = x_key,
  y = Distance,
  fill = Group
)) +
  geom_boxplot(
    width = 0.42,
    outlier.shape = NA,
    alpha = 0.72,
    linewidth = 0.4
  ) +
  geom_point(
    position = position_jitter(
      width = 0.12,
      height = 0,
      seed = 123
    ),
    size = 0.5,
    alpha = 0.26,
    colour = "grey30"
  ) +
  facet_wrap(
    ~Age_class,
    nrow = 1,
    scales = "free_x",
    labeller = as_labeller(age_labels)
  ) +
  scale_x_discrete(labels = x_labels, expand = expansion(add = 0.55)) +
  scale_fill_manual(values = c(Control = "#9DB3D9", CRC = "#E97874"), guide = "none") +
  coord_cartesian(ylim = c(min(disp$Distance) - 0.015, y_top)) +
  labs(
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    x = NULL,
    y = "Distance to\nspatial median",
    tag = "C"
  ) +
  figure_theme +
  theme(axis.title.y = element_text(margin = margin(r = 3)))

figure_s1 <- p_alpha / p_beta / wrap_elements(full = p_disp) +
  plot_layout(heights = c(1.12, 1.18, 1.135)) +
  plot_annotation(
    theme = theme(
      text = element_text(family = "S1Arial"),
      plot.margin = margin(7, 8, 7, 8)
    )
  )

write.csv(beta, file.path(output_dir, "FigureS1B_plot_data.csv"), row.names = FALSE)
write.csv(alpha, file.path(output_dir, "FigureS1A_plot_data.csv"), row.names = FALSE)
write.csv(ds, file.path(output_dir, "FigureS1C_plot_data.csv"), row.names = FALSE)
pdf_path <- file.path(output_dir, "FigureS1.pdf")
grDevices::pdf(
  pdf_path,
  width = figure_width,
  height = figure_height,
  family = "Helvetica",
  useDingbats = FALSE
)
showtext::showtext_begin()
print(figure_s1)
showtext::showtext_end()
grDevices::dev.off()

for (dpi in c(450, 150)) {
  png_path <- file.path(
    output_dir,
    if (dpi == 450)
      "FigureS1.png" else "FigureS1_preview.png"
  )
  showtext::showtext_opts(dpi = dpi)
  ragg::agg_png(
    png_path,
    width = figure_width,
    height = figure_height,
    units = "in",
    res = dpi,
    background = "white"
  )
  showtext::showtext_begin()
  print(figure_s1)
  showtext::showtext_end()
  grDevices::dev.off()
}
