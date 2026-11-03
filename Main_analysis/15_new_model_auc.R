library(tidyverse)
library(patchwork)
library(grid)

out_dir <- file.path(base_dir, "downstream_analysis/output/picture")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!exists("plot_data_main")) {
  plot_data_main <- plot_data %>%
    filter(as.character(Feature) == "Microbiome Only")
}

method_levels <- c("Random Forest", "CatBoost", "TabPFN")
feature_levels <- c("Clinical Only", "Microbiome Only", "Combined")
group_levels <- c("LOCRC", "EOCRC")
summary_rows <- c("Mean", "Matched_EO_Mean")

group_cols <- c(LOCRC = "#3B5C8C", EOCRC = "#D0781F")
grey_text <- "#5A5A5A"
grey_point <- "#6F6F6F"
grey_line <- "#B8B8B8"

standardize_auc_data <- function(dat) {
  dat %>%
    as_tibble() %>%
    mutate(
      Cohort = as.character(Cohort),
      Method = factor(as.character(Method), levels = method_levels),
      Feature = factor(as.character(Feature), levels = feature_levels),
      Group = factor(as.character(Group), levels = group_levels),
      AUC = as.numeric(AUC)
    )
}

plot_data <- standardize_auc_data(plot_data)
plot_data_main <- standardize_auc_data(plot_data_main) %>%
  filter(Feature == "Microbiome Only")

main_points <- plot_data_main %>%
  filter(!Cohort %in% summary_rows)

cohort_mean <- main_points %>%
  group_by(Group, Cohort) %>%
  summarise(
    mean_auc = mean(AUC, na.rm = TRUE),
    n_non_missing = sum(!is.na(AUC)),
    .groups = "drop"
  ) %>%
  filter(n_non_missing > 0)

cohort_order_top <- cohort_mean %>%
  filter(Group == "LOCRC") %>%
  arrange(desc(mean_auc), Cohort) %>%
  pull(Cohort)

cohort_y_levels <- rev(cohort_order_top)
main_points <- main_points %>%
  mutate(Cohort = factor(Cohort, levels = cohort_y_levels))
cohort_mean <- cohort_mean %>%
  mutate(Cohort = factor(Cohort, levels = cohort_y_levels))

missing_eocrc <- setdiff(
  cohort_order_top,
  cohort_mean %>% filter(Group == "EOCRC") %>% pull(Cohort) %>% as.character()
)

missing_dash <- tibble(
  Cohort = factor(missing_eocrc, levels = cohort_y_levels),
  x = 0.508,
  xend = 0.520
)

n_locr <- cohort_mean %>%
  filter(Group == "LOCRC") %>%
  summarise(n = n_distinct(Cohort)) %>%
  pull(n)

n_eocr <- cohort_mean %>%
  filter(Group == "EOCRC") %>%
  summarise(n = n_distinct(Cohort)) %>%
  pull(n)

theme_panel <- function(base_size = 12, axis_lwd = 0.45) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_blank(),
      axis.line.x = element_line(linewidth = axis_lwd, colour = "black"),
      axis.line.y = element_line(linewidth = axis_lwd, colour = "black"),
      axis.ticks = element_line(linewidth = axis_lwd, colour = "black"),
      axis.text = element_text(colour = "black"),
      panel.grid.major.y = element_line(colour = "#E9E9E9", linewidth = 0.45),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = rel(1.25)),
      plot.margin = margin(5, 8, 5, 8)
    )
}

make_group_auc_plot <- function(group_name, show_y = TRUE,
                                show_legend = FALSE) {
  dat <- main_points %>%
    filter(Group == group_name)
  mean_dat <- cohort_mean %>%
    filter(Group == group_name)
  
  xmin <- (floor(min(dat$AUC)*10))/10
  
  p <- ggplot(dat, aes(x = AUC, y = Cohort)) +
    geom_vline(
      xintercept = 0.80,
      colour = grey_line,
      linewidth = 0.65,
      linetype = "dashed"
    ) +
    geom_point(
      aes(shape = Method),
      position = position_jitter(width = 0, height = 0.12, seed = 17),
      colour = grey_point,
      fill = "white",
      size = 2.45,
      stroke = 0.8,
      na.rm = TRUE
    ) +
    geom_point(
      data = mean_dat,
      aes(x = mean_auc, y = Cohort, shape = "Cohort mean"),
      colour = group_cols[[group_name]],
      fill = group_cols[[group_name]],
      size = 3.25,
      stroke = 0.6,
      inherit.aes = FALSE
    ) +
    scale_shape_manual(
      values = c(
        "Random Forest" = 21,
        "CatBoost" = 22,
        "TabPFN" = 24,
        "Cohort mean" = 23
      ),
      breaks = c(method_levels, "Cohort mean"),
      name = NULL
    ) +
    scale_x_continuous(
      breaks = seq(xmin, 1.0, by = 0.1),
      expand = expansion(mult = c(0, 0.005))
    ) +
    coord_cartesian(xlim = c(xmin, 1.0)) +
    scale_y_discrete(drop = FALSE) +
    labs(
      title = sprintf(
        "%s  (n = %d cohorts)",
        group_name,
        ifelse(group_name == "LOCRC", n_locr, n_eocr)
      ),
      x = "AUC",
      y = if (show_y) "Test cohort" else NULL
    ) +
    theme_panel(base_size = 12) +
    theme(
      plot.title = element_text(
        colour = group_cols[[group_name]],
        face = "bold",
        hjust = 0.5,
        size = 16,
        margin = margin(b = 8)
      ),
      legend.position = if (show_legend) c(0.70,0) else "none",
      legend.justification = c(0, 0),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = 9.5),
      legend.spacing.y = unit(0.05, "cm")
    ) +
    guides(
      shape = guide_legend(
        override.aes = list(
          colour = grey_point,
          fill = c("white", "white", "white", grey_point),
          size = c(2.8, 2.8, 2.8, 3.2),
          stroke = c(0.9, 0.9, 0.9, 0.8)
        )
      )
    )

  if (group_name == "EOCRC" && nrow(missing_dash) > 0) {
    p <- p +
      geom_segment(
        data = missing_dash,
        aes(x = x, xend = xend, y = Cohort, yend = Cohort),
        inherit.aes = FALSE,
        colour = "#C7C7C7",
        linewidth = 0.9,
        lineend = "butt"
      )
  }

  if (!show_y) {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()
      )
  }

  p
}

plot_a <- (make_group_auc_plot("LOCRC", show_y = TRUE, show_legend = FALSE) +
             plot_spacer()+
    make_group_auc_plot("EOCRC", show_y = FALSE, show_legend = TRUE)
  ) + plot_layout(widths = c(6,0.1, 5))

plot_a

algo_points <- main_points %>%
  filter(!is.na(AUC)) %>%
  mutate(Method = factor(Method, levels = method_levels))

algo_summary <- algo_points %>%
  group_by(Group, Method) %>%
  summarise(
    mean_auc = mean(AUC),
    sd_auc = sd(AUC),
    .groups = "drop"
  )

plot_b <- ggplot(algo_points, aes(x = Method, y = AUC, colour = Group)) +
  geom_hline(
    yintercept = 0.80,
    colour = grey_line,
    linewidth = 0.55,
    linetype = "dashed"
  ) +
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.08,
      jitter.height = 0,
      dodge.width = 0.46,
      seed = 17
    ),
    alpha = 0.23,
    size = 2.2,
    stroke = 0,
    show.legend = FALSE
  ) +
  geom_errorbar(
    data = algo_summary,
    aes(
      y = mean_auc,
      ymin = mean_auc - sd_auc,
      ymax = mean_auc + sd_auc,
      group = Group
    ),
    position = position_dodge(width = 0.46),
    width = 0.16,
    linewidth = 0.65
  ) +
  geom_point(
    data = algo_summary,
    aes(y = mean_auc, group = Group),
    position = position_dodge(width = 0.46),
    size = 3.2
  ) +
  scale_colour_manual(values = group_cols, name = NULL) +
  scale_x_discrete(labels = c("Random\nForest", "CatBoost", "TabPFN")) +
  scale_y_continuous(
    breaks = seq(0.5, 1.0, by = 0.1),
    expand = expansion(mult = c(0.02, 0))
  ) +
  coord_cartesian(ylim = c(0.48, 1.01), clip = "off")+
  labs(x = NULL, y = "AUC") +
  theme_panel(base_size = 12) +
  theme(
    # plot.title = element_text(face = "bold", hjust = 0, size = 15),
    axis.text.x = element_text(size = 10),
    legend.position = c(0.8, 0),
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key = element_blank(),  
    panel.grid.major.y = element_blank()
  )
plot_b

matched_diff <- cohort_mean %>%
  select(Group, Cohort, mean_auc) %>%
  mutate(Cohort = as.character(Cohort)) %>%
  pivot_wider(names_from = Group, values_from = mean_auc) %>%
  filter(!is.na(LOCRC), !is.na(EOCRC)) %>%
  mutate(
    delta = EOCRC - LOCRC,
    direction = if_else(delta >= 0, "EOCRC higher", "LOCRC higher"),
    cohort_label = Cohort
  ) %>%
  arrange(delta) %>%
  mutate(cohort_label = factor(cohort_label, levels = cohort_label))

mean_delta <- mean(matched_diff$delta)

supplement_pdf <- ggplot(matched_diff, aes(y = cohort_label)) +
  geom_vline(xintercept = 0, colour = "#5A5A5A", linewidth = 1.0) +
  geom_segment(
    aes(x = pmin(delta, 0), xend = pmax(delta, 0), yend = cohort_label),
    colour = grey_line,
    linewidth = 1.05
  ) +
  geom_point(
    aes(x = delta, colour = direction),
    size = 3.0,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = 0.19,
    y = 1.18,
    label = sprintf(
      "Mean \u0394 = %+0.3f\nmatched cohorts n = %d",
      mean_delta,
      nrow(matched_diff)
    ),
    hjust = 1,
    vjust = 0,
    colour = grey_text,
    size = 4.0
  ) +
  scale_colour_manual(
    values = c("EOCRC higher" = group_cols[["EOCRC"]], "LOCRC higher" = group_cols[["LOCRC"]])
  ) +
  scale_x_continuous(
    limits = c(-0.23, 0.23),
    breaks = seq(-0.2, 0.2, by = 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(title = "C  Matched-cohort difference", x = "EOCRC - LOCRC mean AUC", y = NULL) +
  theme_panel(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0, size = 15),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 11.5)
  )


feature_points <- plot_data %>%
  filter(!Cohort %in% summary_rows, !is.na(AUC)) %>%
  group_by(Group, Cohort, Feature) %>%
  summarise(mean_auc = mean(AUC), .groups = "drop") %>%
  mutate(
    Feature = factor(Feature, levels = feature_levels),
    feature_id = as.numeric(Feature),
    Group = factor(Group, levels = group_levels)
  )

feature_summary <- feature_points %>%
  group_by(Group, Feature) %>%
  summarise(
    mean_auc_summary = mean(mean_auc),
    sd_auc = sd(mean_auc),
    .groups = "drop"
  ) %>%
  rename(mean_auc = mean_auc_summary) %>%
  mutate(feature_id = as.numeric(Feature))

make_feature_plot <- function(group_name, show_y = TRUE) {
  dat <- feature_points %>%
    filter(Group == group_name)
  sum_dat <- feature_summary %>%
    filter(Group == group_name)

  ggplot(dat, aes(x = feature_id, y = mean_auc, group = Cohort)) +
    geom_hline(
      yintercept = 0.80,
      colour = grey_line,
      linewidth = 0.55,
      linetype = "dashed"
    ) +
    geom_line(
      colour = "#D5D5D5",
      linewidth = 0.65,
      alpha = 0.3
    ) +
    geom_point(
      colour = group_cols[[group_name]],
      fill = group_cols[[group_name]],
      alpha = 0.25,
      size = 2.4,
      stroke = 0
    ) +
    geom_segment(
      data = sum_dat,
      aes(
        x = feature_id,
        xend = feature_id,
        y = mean_auc - sd_auc,
        yend = mean_auc + sd_auc
      ),
      inherit.aes = FALSE,
      colour = group_cols[[group_name]],
      linewidth = 0.95
    ) +
    geom_segment(
      data = sum_dat,
      aes(
        x = feature_id - 0.06,
        xend = feature_id + 0.06,
        y = mean_auc - sd_auc,
        yend = mean_auc - sd_auc
      ),
      inherit.aes = FALSE,
      colour = group_cols[[group_name]],
      linewidth = 0.95
    ) +
    geom_segment(
      data = sum_dat,
      aes(
        x = feature_id - 0.06,
        xend = feature_id + 0.06,
        y = mean_auc + sd_auc,
        yend = mean_auc + sd_auc
      ),
      inherit.aes = FALSE,
      colour = group_cols[[group_name]],
      linewidth = 0.95
    ) +
    geom_point(
      data = sum_dat,
      aes(x = feature_id, y = mean_auc),
      inherit.aes = FALSE,
      shape = 23,
      colour = group_cols[[group_name]],
      fill = group_cols[[group_name]],
      size = 3.1,
      stroke = 0.6
    ) +
    scale_x_continuous(
      breaks = seq_along(feature_levels),
      labels = c("Clinical\nonly", "Microbiome\nonly", "Combined")
    ) +
    scale_y_continuous(
      breaks = seq(0.4, 1.0, by = 0.1),
      expand = expansion(mult = c(0.02, 0))
    ) +
    coord_cartesian(ylim = c(0.40, 1.02), clip = "off") +
    labs(
      title = group_name,
      x = NULL,
      y = if (show_y) "Cohort-mean AUC\nacross algorithms" else NULL
    ) +
    theme_panel(base_size = 12) +
    theme(
      plot.title = element_text(
        colour = group_cols[[group_name]],
        face = "bold",
        hjust = 0.5,
        size = 16,
        margin = margin(b = 9)
      ),
      axis.text.x = element_text(size = 10),
      panel.grid.major.y = element_blank()
    ) +
    {
      if (!show_y) {
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()
        )
      } else {
        theme()
      }
    }
}

plot_c <- (make_feature_plot("LOCRC", show_y = TRUE) +
             plot_spacer()+
    make_feature_plot("EOCRC", show_y = FALSE)
  ) +
  plot_layout(widths = c(1,0.1, 1)) +
  plot_annotation(
    theme = theme(
     plot.margin = margin(4, 8, 4, 8)
    )
  )
plot_c


bottom_row <- (free(plot_b, side = "l") | plot_spacer()|
    free(plot_c, side = "l")
) +
  plot_layout(widths = c(2,0.2,6))

bottom_row
main_pdf <- (plot_a / plot_spacer()/ 
               bottom_row
             )+
  plot_layout(heights = c(2.5, 0.05,1.15)) +
  plot_annotation(
    tag_levels = list(c("A", "","B", "C"))
  ) +
  plot_annotation(
    tag_levels = list(c("A", "", "B", "C"))
  ) &
  theme(
    plot.margin = margin(4,4, 4, 4),
    plot.tag = element_text(size = 18, face = "bold"),
    plot.tag.position = c(0.012, 0.98)
  )

main_pdf

pdf(file.path(out_dir, "17.1_model_auc_visual_scheme_main.pdf"), 
    width = 12, height = 9)
main_pdf
dev.off()

message("Saved: ", main_pdf_file)
message("Saved: ", supp_pdf_file)
