#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridGraphics)
  library(ggplotify)
  library(patchwork)
  library(readr)
  library(dplyr)
  library(tidyr)
})

env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

make_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

current_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hits <- grep(paste0("^", file_arg), args, value = TRUE)
  if (length(hits) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", hits[[1]]), winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

safe_ggsave <- function(filename, plot, ..., fallback_suffix = "_revised") {
  tryCatch(
    {
      ggsave(filename, plot, ...)
      filename
    },
    error = function(e) {
      fallback <- sub("(\\.[^.]+)$", paste0(fallback_suffix, "\\1"), filename)
      warning(
        sprintf("Could not write %s; writing %s instead. Original error: %s", filename, fallback, conditionMessage(e)),
        call. = FALSE
      )
      ggsave(fallback, plot, ...)
      fallback
    }
  )
}

source_with_env_flag <- function(script, env, flag_name, flag_value) {
  old_value <- Sys.getenv(flag_name, unset = NA_character_)
  on.exit({
    if (is.na(old_value)) {
      Sys.unsetenv(flag_name)
    } else {
      do.call(Sys.setenv, setNames(list(old_value), flag_name))
    }
  }, add = TRUE)
  do.call(Sys.setenv, setNames(list(flag_value), flag_name))
  sys.source(script, envir = env)
}

script_dir <- normalizePath(
  env_or_default(
    "STRAINPHLAN_V5_SCRIPT_DIR",
    current_script_dir()
  ),
  winslash = "/",
  mustWork = TRUE
)

project_root <- normalizePath(env_or_default("PROJECT_ROOT", dirname(script_dir)), winslash = "/", mustWork = FALSE)
result_root <- env_or_default(
  "STRAINPHLAN_V5_RESULT_ROOT",
  file.path(project_root, "results/strain_analysis_v5")
)

output_dir <- make_dir(env_or_default(
  "STRAINPHLAN_ROSEBURIA_FULL_SUPP_OUTPUT_DIR",
  file.path(result_root, "roseburia_full_supplement_figure")
))

region_script <- file.path(script_dir, "03_plot_roseburia_region_structure.R")
strain_evidence_script <- file.path(script_dir, "04_plot_roseburia_strain_supplement_evidence.R")

if (!file.exists(region_script)) {
  stop("Missing region structure script: ", region_script, call. = FALSE)
}
if (!file.exists(strain_evidence_script)) {
  stop("Missing strain evidence script: ", strain_evidence_script, call. = FALSE)
}

region_env <- new.env(parent = globalenv())
source_with_env_flag(
  region_script,
  region_env,
  "STRAINPHLAN_ROSEBURIA_REGION_RENDER_FILES",
  "false"
)

strain_env <- new.env(parent = globalenv())
sys.source(strain_evidence_script, envir = strain_env)

required_region_objects <- c( "overview_plot")
missing_region_objects <- setdiff(required_region_objects, ls(region_env))
if (length(missing_region_objects) > 0) {
  stop("Region script did not create required objects: ", paste(missing_region_objects, collapse = ", "), call. = FALSE)
}

required_region_data <- c("pcoa_tbl", "nn_summary", "region_summary")
missing_region_data <- setdiff(required_region_data, ls(region_env))
if (length(missing_region_data) > 0) {
  stop("Region script did not create required source data: ", paste(missing_region_data, collapse = ", "), call. = FALSE)
}

required_strain_objects <- c("case_plot", "eolo_plot")
missing_strain_objects <- setdiff(required_strain_objects, ls(strain_env))
if (length(missing_strain_objects) > 0) {
  stop("Strain evidence script did not create required objects: ", paste(missing_strain_objects, collapse = ", "), call. = FALSE)
}

required_strain_data <- c("case_grid", "eolo_grid")
missing_strain_data <- setdiff(required_strain_data, ls(strain_env))
if (length(missing_strain_data) > 0) {
  stop("Strain evidence script did not create required source data: ", paste(missing_strain_data, collapse = ", "), call. = FALSE)
}

# ==============================================================================
# Panels:
# A. Patristic-distance ordination
# B. Same-region nearest-neighbor rate
# C. Within-region strain-structure evidence
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Master colors
# ------------------------------------------------------------------------------
COL_OLD   <- "#2F3E7E"
COL_YOUNG <- "#F05A4A"
COL_TEXT  <- "#222222"

region_env$region_cols <- c(
  Europe    = COL_OLD,
  East_Asia = COL_YOUNG
)


make_compact_tree_grob <- function(data, panel_label = "") {
  local_data <- data
  local_panel_label <- panel_label
  
  grid::grid.grabExpr(
    gridGraphics::grid.echo(function() {
      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par), add = TRUE)
      
      par(
        mar = c(0.05, 0.05, 0.25, 0.05),
        oma = c(0, 0, 0, 0),
        xpd = NA
      )
      
      region_env$draw_region_ring_tree(
        local_data,
        panel_label = local_panel_label,
        title_line = -1.15
      )
    })
  )
}

make_tree_metric_plot <- function(species_label) {
  metric_row <- region_env$region_summary %>%
    filter(.data$species_label == !!species_label) %>%
    slice(1)

  if (nrow(metric_row) != 1) {
    stop("Cannot find one region-summary row for ", species_label, call. = FALSE)
  }

  metric_text <- sprintf(
    "dELPD/sample = %.3f   z = %.1f   n = %d",
    metric_row$pglmm_minus_base_elpd_per_sample,
    metric_row$pglmm_minus_base_z,
    metric_row$n_used_samples
  )

  ggplot() +
    annotate(
      "label",
      x = 0.5,
      y = 0.5,
      label = metric_text,
      linewidth = 0,
      fill = "white",
      color = COL_TEXT,
      size = 3.0,
      family = "",
      label.padding = unit(0.12, "lines")
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(0, 4, 1, 4))
}

make_region_legend_plot <- function() {
  legend_grob <- grid::grobTree(
    grid::textGrob(
      "Region ring",
      x = unit(0.02, "npc"),
      y = unit(0.58, "npc"),
      hjust = 0,
      gp = grid::gpar(fontface = "bold", fontsize = 10.5, col = COL_TEXT)
    ),
    grid::rectGrob(
      x = unit(0.06, "npc"),
      y = unit(0.50, "npc"),
      width = unit(4.2, "mm"),
      height = unit(4.2, "mm"),
      just = c("left", "center"),
      gp = grid::gpar(fill = COL_OLD, col = NA)
    ),
    grid::textGrob(
      "Europe",
      x = unit(0.28, "npc"),
      y = unit(0.50, "npc"),
      hjust = 0,
      gp = grid::gpar(fontsize = 9.5, col = COL_TEXT)
    ),
    grid::rectGrob(
      x = unit(0.06, "npc"),
      y = unit(0.455, "npc"),
      width = unit(4.2, "mm"),
      height = unit(4.2, "mm"),
      just = c("left", "center"),
      gp = grid::gpar(fill = COL_YOUNG, col = NA)
    ),
    grid::textGrob(
      "East Asia",
      x = unit(0.28, "npc"),
      y = unit(0.455, "npc"),
      hjust = 0,
      gp = grid::gpar(fontsize = 9.5, col = COL_TEXT)
    )
  )

  patchwork::wrap_elements(full = legend_grob)
}

rf_tree_plot <- ggplotify::as.ggplot(
  make_compact_tree_grob(region_env$species_data[[1]], "")
)

ri_tree_plot <- ggplotify::as.ggplot(
  make_compact_tree_grob(region_env$species_data[[2]], "")
)

rf_tree_block <- rf_tree_plot / make_tree_metric_plot("R. faecis") +
  plot_layout(heights = c(1, 0.085))

ri_tree_block <- ri_tree_plot / make_tree_metric_plot("R. inulinivorans") +
  plot_layout(heights = c(1, 0.085))

tree_row <- (
  rf_tree_block |
    ri_tree_block |
    make_region_legend_plot()
) +
  plot_layout(widths = c(1, 1, 0.22)) &
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )


tree_panel <- patchwork::wrap_elements(full = tree_row)


nn_tbl_revised <- region_env$nn_summary %>%
  select(species_label, observed_same_region, expected_random) %>%
  pivot_longer(
    c(observed_same_region, expected_random),
    names_to = "metric",
    values_to = "same_region_rate"
  ) %>%
  mutate(
    metric = recode(
      metric,
      expected_random = "Expected",
      observed_same_region = "Observed"
    ),
    metric = factor(metric, levels = c("Expected", "Observed")),
    label = paste0(round(100 * same_region_rate, 1), "%"),
    species_label = factor(
      species_label,
      levels = c("R. faecis", "R. inulinivorans")
    )
  )

nn_structure_plot <- ggplot(nn_tbl_revised, aes(metric, same_region_rate, fill = metric)) +
  geom_col(
    width = 0.60,
    color = "#333333",
    linewidth = 0.32
  ) +
  geom_text(
    aes(label = label),
    vjust = -0.35,
    size = 3.0,
    color = COL_TEXT
  ) +
  facet_wrap(~ species_label, nrow = 1) +
  scale_fill_manual(
    values = c(
      "Expected" = '#A7B4D2',
      "Observed" = '#F05A4A'
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    limits = c(0, 1.08),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Same-region nearest-neighbor rate",
    x = NULL,
    y = "Same-region rate"
  ) +
  theme_classic(base_size = 9.5) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold.italic", size = 9.5),
    axis.text.x = element_text(size = 8.6,angle = 30,vjust = 0.5,
                               hjust = 0.5, color = COL_TEXT),
    axis.text.y = element_text(size = 8.3, color = COL_TEXT),
    axis.title.y = element_text(size = 9.2),
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.30),
    plot.title = element_text(face = "bold", size = 11.2, margin = margin(b = 3)),
    plot.margin = margin(2, 6, 2, 6)
  )
nn_structure_plot

format_p2 <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    ifelse(x < 0.001, "<0.001", sprintf("%.3f", x))
  )
}

format_r2 <- function(x) {
  ifelse(
    is.na(x),
    "NA",
    sprintf("%.3f", x)
  )
}

classify_strain_evidence <- function(p_value, q_value) {
  case_when(
    # !is.na(q_value) & q_value < 0.10 ~ "q16 < 0.10",
    !is.na(p_value) & p_value < 0.05 ~ "Nominal p < 0.05",
    TRUE ~ "No evidence"
  )
}

make_heatmap_input <- function() {
  case_df <- strain_env$case_grid %>%
    mutate(
      evidence_panel = "Case-control strain structure",
      stratum_label = as.character(stratum_label)
    )
  
  eolo_df <- strain_env$eolo_grid %>%
    mutate(
      evidence_panel = "EO vs LO strain structure",
      stratum_label = as.character(stratum_label)
    )
  
  bind_rows(case_df, eolo_df) %>%
    mutate(
      evidence_panel = factor(
        evidence_panel,
        levels = c(
          "Case-control strain structure",
          "EO vs LO strain structure"
        )
      ),
      stratum_label = factor(
        stratum_label,
        levels = c(
          "Europe EO", "East Asia EO", "Europe LO", "East Asia LO",
          "Europe CRC", "East Asia CRC", "Europe Control", "East Asia Control"
        )
      ),
      species_label = factor(
        as.character(species_label),
        levels = c("R. faecis", "R. inulinivorans")
      ),
      evidence_status = classify_strain_evidence(p_value, q_bh_all_strain),
      evidence_status = factor(
        evidence_status,
        levels = c("No evidence", "Nominal p < 0.05"#, "q16 < 0.10"
                   )
      ),
      cell_label = paste0(
        "p=", format_p2(p_value),
        "\nq=", format_p2(q_bh_all_strain),
        "\nR²=", format_r2(r2),
        "\nn=", n_used_samples
      )
    )
}

make_heatmap_input_revised <- function() {
  make_heatmap_input() %>%
    mutate(
      cell_label = paste0(
        "p ", format_p2(p_value),
        "\nq16 ", format_p2(q_bh_all_strain),
        "\nR2 ", format_r2(r2)
      ),
      label_color =  COL_TEXT,
      label_face = if_else(evidence_status == "No evidence", "plain", "bold")
    )
}

heatmap_df <- make_heatmap_input_revised()

heatmap_x_labels <- c(
  "Europe EO" = "Europe\nEO",
  "East Asia EO" = "East Asia\nEO",
  "Europe LO" = "Europe\nLO",
  "East Asia LO" = "East Asia\nLO",
  "Europe CRC" = "Europe\nCRC",
  "East Asia CRC" = "East Asia\nCRC",
  "Europe Control" = "Europe\nControl",
  "East Asia Control" = "East Asia\nControl"
)

strain_evidence_heatmap <- ggplot(
  heatmap_df,
  aes(x = stratum_label, y = species_label)
) +
  geom_tile(
    aes(fill = evidence_status),
    color = "white",
    linewidth = 0.85,
    width = 0.98,
    height = 0.92
  ) +
  geom_tile(
    data = heatmap_df %>% filter(evidence_status != "No evidence"),
    fill = NA,
    color = "#4A342D",
    linewidth = 0.42,
    width = 0.98,
    height = 0.92,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = cell_label, color = label_color, fontface = label_face),
    size = 2.50,
    lineheight = 0.82,
    show.legend = FALSE
  ) +
  facet_grid(
    . ~ evidence_panel,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_x_discrete(
    labels = heatmap_x_labels,
    drop = TRUE,
    expand = expansion(add = 0)
  ) +
  scale_y_discrete(
    limits = rev(c("R. faecis", "R. inulinivorans")),
    labels = strain_env$species_labels,
    expand = expansion(add = 0)
  ) +
  scale_fill_manual(
    values = c(
      "No evidence" = "#F3F3F0",
      "Nominal p < 0.05" = "#F3C99B"#,
      # "q16 < 0.10" = "#A84E3A"
    ),
    name = NULL,
    drop = FALSE
  ) +
  scale_color_identity() +
  labs(
    title = "Within-region strain-structure evidence",
    x = NULL,
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid = element_blank(),
    panel.spacing.x = unit(0.60, "lines"),
    strip.background = element_blank(),
    strip.text.x = element_text(
      face = "bold",
      size = 10.7,
      color = COL_TEXT
    ),
    axis.text.x = element_text(
      face = "bold",
      size = 8.6,
      color = COL_TEXT,
      lineheight = 0.92,
      margin = margin(t = 4)
    ),
    axis.text.y = element_text(
      size = 9.6,
      color = COL_TEXT,
      margin = margin(r = 6)
    ),
    plot.title = element_text(
      face = "bold",
      size = 11.2,
      color = COL_TEXT,
      margin = margin(b = 5)
    ),
    legend.position = "top",
    legend.justification = "left",
    legend.direction = "horizontal",
    legend.key.width = unit(0.62, "lines"),
    legend.key.height = unit(0.62, "lines"),
    legend.text = element_text(size = 8.2, color = COL_TEXT),
    legend.margin = margin(-4, 0, -2, 0),
    plot.margin = margin(2, 6, 2, 6)
  )


bottom_row_inner <- nn_structure_plot +
  plot_spacer() +
  strain_evidence_heatmap +
  plot_layout(widths = c(0.50, 0.08, 1.00))


bottom_row <- plot_spacer() +
  bottom_row_inner +
  plot_spacer() +
  plot_layout(widths = c(0.006, 1, 0.006))

layout_design <- c(
    area(t = 1, l = 1,  b = 1, r = 10), 
    area(t = 2, l = 1,  b = 2, r = 1),   
    area(t = 3, l = 1, b = 3, r = 10) 
  )

full_supplement_plot <- (
  tree_panel /
    plot_spacer() /
    bottom_row
) +
  plot_layout(
    design = layout_design,
    heights = c(3, 0.15, 1)
  ) +
  plot_annotation(
    tag_levels = list(c("B", "C", "D"))
  ) &
  theme(
    plot.margin = margin(4, 4, 4, 4),
    plot.background = element_rect(fill = "white", color = NA),
    plot.tag = element_text(size = 16, face = "bold"),
    plot.tag.position = c(0.012, 0.98)
  )

# full_supplement_plot

# ------------------------------------------------------------------------------
# 5. Save
# ------------------------------------------------------------------------------

pdf_file <- safe_ggsave(
  file.path(output_dir, "roseburia_full_supplement_figure_revised.pdf"),
  full_supplement_plot,
  width = 10,
  height = 8,
  units = "in",
  bg = "white"
)
# '#7EA8CB'
# '#FF927F'
png_file <- safe_ggsave(
  file.path(output_dir, "roseburia_full_supplement_figure_revised.png"),
  full_supplement_plot,
  width = 14.0,
  height = 9.4,
  units = "in",
  dpi = 360,
  bg = "white"
)

saveRDS(
  list(
    full_supplement_plot = full_supplement_plot,
    tree_row = tree_row,
    nn_structure_plot = nn_structure_plot,
    strain_evidence_heatmap = strain_evidence_heatmap
  ),
  file.path(output_dir, "roseburia_full_supplement_figure_revised_components.rds")
)
