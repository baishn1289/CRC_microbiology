#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(ggplotify)
  library(grid)
  library(gridGraphics)
})

env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

env_or_default_num <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  parsed <- suppressWarnings(as.numeric(value))
  if (is.na(parsed)) {
    stop(sprintf("Environment variable %s must be numeric, got: %s", name, value), call. = FALSE)
  }
  parsed
}

env_or_default_bool <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    return(default)
  }
  value <- tolower(value)
  if (value %in% c("1", "true", "t", "yes", "y")) {
    return(TRUE)
  }
  if (value %in% c("0", "false", "f", "no", "n")) {
    return(FALSE)
  }
  stop(sprintf("Environment variable %s must be true/false, got: %s", name, value), call. = FALSE)
}

make_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

find_required_file <- function(directory, pattern, description) {
  files <- list.files(directory, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop(sprintf("Could not find %s in %s using pattern %s", description, directory, pattern), call. = FALSE)
  }
  files[[1]]
}

draw_annulus_sector <- function(theta_lower, theta_upper, r_inner, r_outer, col, border = NA, n = 6) {
  theta <- seq(theta_lower, theta_upper, length.out = max(2, n))
  polygon(
    x = c(r_outer * cos(theta), rev(r_inner * cos(theta))),
    y = c(r_outer * sin(theta), rev(r_inner * sin(theta))),
    col = col,
    border = border
  )
}

region_cols <- c(
  Europe = "#2266A5",
  East_Asia = "#D77A1F"
)

project_root <- normalizePath(env_or_default("PROJECT_ROOT", getwd()), winslash = "/", mustWork = FALSE)
result_root <- env_or_default("STRAINPHLAN_V5_RESULT_ROOT", file.path(project_root, "results/strain_analysis_v5"))

species_specs <- tribble(
  ~species_id, ~species_label, ~strain_dir, ~region_result_dir,
  "rf4925", "R. faecis",
  env_or_default("STRAINPHLAN_RF4925_STRAIN_DIR", file.path(project_root, "data/strainphlan/Europe_Asia_t__SGB4925")),
  env_or_default("STRAINPHLAN_RF4925_REGION_RESULT_DIR", file.path(result_root, "anpan_region_comparison/Europe_Asia_t__SGB4925")),
  "ri4940", "R. inulinivorans",
  env_or_default("STRAINPHLAN_RI4940_STRAIN_DIR", file.path(project_root, "data/strainphlan/Europe_Asia_t__SGB4940")),
  env_or_default("STRAINPHLAN_RI4940_REGION_RESULT_DIR", file.path(result_root, "anpan_region_comparison_ri4940/Europe_Asia_t__SGB4940"))
)

output_dir <- make_dir(env_or_default(
  "STRAINPHLAN_ROSEBURIA_REGION_STRUCTURE_OUTPUT_DIR",
  file.path(result_root, "roseburia_region_structure_paper")
))

single_tree_output_dir <- make_dir(env_or_default(
  "STRAINPHLAN_ROSEBURIA_REGION_SINGLE_TREE_OUTPUT_DIR",
  file.path(output_dir, "individual_region_trees")
))
single_tree_width <- env_or_default_num("STRAINPHLAN_ROSEBURIA_REGION_SINGLE_TREE_WIDTH", 5.4)
single_tree_height <- env_or_default_num("STRAINPHLAN_ROSEBURIA_REGION_SINGLE_TREE_HEIGHT", 5.4)
render_files <- env_or_default_bool("STRAINPHLAN_ROSEBURIA_REGION_RENDER_FILES", TRUE)

prepare_species_data <- function(spec) {
  tree_file <- find_required_file(spec$strain_dir, "^RAxML_bestTree\\..*\\.tre$", "RAxML best tree")
  meta <- read_tsv(file.path(spec$region_result_dir, "retained_metadata.tsv"), show_col_types = FALSE) %>%
    mutate(
      species_id = spec$species_id,
      species_label = spec$species_label,
      region = factor(region, levels = c("Europe", "East_Asia"))
    )

  tree_all <- read.tree(tree_file)
  sample_tree <- drop.tip(tree_all, tree_all$tip.label[str_detect(tree_all$tip.label, "^(GCF_|GCA_)")])
  tree_retained <- drop.tip(sample_tree, setdiff(sample_tree$tip.label, meta$tree_tip_label)) %>%
    ladderize()

  meta <- meta %>%
    filter(tree_tip_label %in% tree_retained$tip.label) %>%
    mutate(tree_tip_label = factor(tree_tip_label, levels = tree_retained$tip.label)) %>%
    arrange(tree_tip_label) %>%
    mutate(tree_tip_label = as.character(tree_tip_label))

  dist_mat <- cophenetic.phylo(tree_retained)
  dist_mat <- dist_mat[meta$tree_tip_label, meta$tree_tip_label, drop = FALSE]

  pcoa <- cmdscale(as.dist(dist_mat), eig = TRUE, k = 2)
  positive_eig <- pcoa$eig[pcoa$eig > 0]
  var_exp <- pcoa$eig[1:2] / sum(positive_eig) * 100
  pcoa_points <- as.data.frame(pcoa$points)
  colnames(pcoa_points)[seq_len(2)] <- c("PCo1", "PCo2")
  pcoa_tbl <- rownames_to_column(pcoa_points, var = "tree_tip_label") %>%
    as_tibble() %>%
    left_join(meta, by = "tree_tip_label") %>%
    mutate(
      species_label = spec$species_label,
      PCo1_var = var_exp[[1]],
      PCo2_var = var_exp[[2]]
    )

  nn_mat <- dist_mat
  diag(nn_mat) <- Inf
  nearest_idx <- apply(nn_mat, 1, which.min)
  nn_tbl <- tibble(
    tree_tip_label = rownames(nn_mat),
    nearest_tip_label = colnames(nn_mat)[nearest_idx],
    nearest_distance = nn_mat[cbind(seq_along(nearest_idx), nearest_idx)]
  ) %>%
    left_join(select(meta, tree_tip_label, region), by = "tree_tip_label") %>%
    left_join(
      select(meta, tree_tip_label, nearest_region = region),
      by = c("nearest_tip_label" = "tree_tip_label")
    ) %>%
    mutate(
      species_label = spec$species_label,
      same_region = region == nearest_region
    )

  region_counts <- table(meta$region)
  n_total <- sum(region_counts)
  expected_same_region <- sum(region_counts * (region_counts - 1)) / (n_total * (n_total - 1))
  nn_summary <- nn_tbl %>%
    summarise(
      species_id = spec$species_id,
      species_label = spec$species_label,
      n_samples = n(),
      observed_same_region = mean(same_region, na.rm = TRUE),
      expected_random = expected_same_region,
      enrichment = observed_same_region / expected_same_region,
      .groups = "drop"
    )

  idx <- which(upper.tri(dist_mat), arr.ind = TRUE)
  pair_tbl <- tibble(
    tree_tip_i = rownames(dist_mat)[idx[, 1]],
    tree_tip_j = colnames(dist_mat)[idx[, 2]],
    distance = dist_mat[idx]
  ) %>%
    left_join(select(meta, tree_tip_label, region_i = region), by = c("tree_tip_i" = "tree_tip_label")) %>%
    left_join(select(meta, tree_tip_label, region_j = region), by = c("tree_tip_j" = "tree_tip_label")) %>%
    mutate(
      species_label = spec$species_label,
      pair_group = case_when(
        region_i == "Europe" & region_j == "Europe" ~ "Within Europe",
        region_i == "East_Asia" & region_j == "East_Asia" ~ "Within East Asia",
        TRUE ~ "Between regions"
      ),
      pair_group = factor(pair_group, levels = c("Within Europe", "Within East Asia", "Between regions"))
    )

  region_summary <- read_tsv(file.path(spec$region_result_dir, "analysis_summary.tsv"), show_col_types = FALSE) %>%
    filter(model_id == "01_phylogeny_region_overall") %>%
    transmute(
      species_id = spec$species_id,
      species_label = spec$species_label,
      n_used_samples,
      pglmm_minus_base_elpd,
      pglmm_minus_base_elpd_per_sample,
      pglmm_minus_base_z,
      interpretation
    )

  list(
    spec = spec,
    tree = tree_retained,
    meta = meta,
    dist_mat = dist_mat,
    pcoa_tbl = pcoa_tbl,
    pair_tbl = pair_tbl,
    nn_tbl = nn_tbl,
    nn_summary = nn_summary,
    region_summary = region_summary
  )
}

species_data <- lapply(seq_len(nrow(species_specs)), function(i) prepare_species_data(species_specs[i, ]))

pcoa_tbl <- bind_rows(lapply(species_data, `[[`, "pcoa_tbl")) %>%
  mutate(
    species_label = factor(species_label, levels = species_specs$species_label),
    region_label = recode(as.character(region), Europe = "Europe", East_Asia = "East Asia"),
    region_label = factor(region_label, levels = c("Europe", "East Asia"))
  ) %>%
  arrange(region_label)

pair_tbl <- bind_rows(lapply(species_data, `[[`, "pair_tbl")) %>%
  mutate(
    species_label = factor(species_label, levels = species_specs$species_label),
    pair_group = factor(pair_group, levels = c("Within Europe", "Within East Asia", "Between regions"))
  )

nn_summary <- bind_rows(lapply(species_data, `[[`, "nn_summary")) %>%
  mutate(species_label = factor(species_label, levels = species_specs$species_label))

region_summary <- bind_rows(lapply(species_data, `[[`, "region_summary"))

axis_labels <- pcoa_tbl %>%
  group_by(species_label) %>%
  summarise(
    x_lab = sprintf("PCo1 (%.1f%%)", first(PCo1_var)),
    y_lab = sprintf("PCo2 (%.1f%%)", first(PCo2_var)),
    .groups = "drop"
  )

pcoa_plot <- ggplot(pcoa_tbl, aes(PCo1, PCo2, color = region_label)) +
  stat_ellipse(aes(group = region_label), linewidth = 0.55, alpha = 0.9, level = 0.75, show.legend = FALSE) +
  geom_point(size = 1.25, alpha = 0.58, stroke = 0) +
  facet_wrap(~ species_label, scales = "free") +
  scale_color_manual(
    values = c(Europe = region_cols[["Europe"]], `East Asia` = region_cols[["East_Asia"]]),
    breaks = c("Europe", "East Asia"),
    name = "Region"
  ) +
  labs(
    title = "A. Strain ordination from patristic distances",
    subtitle = "Each point is a strain-retained sample; ellipses summarize local regional concentration",
    x = "PCo1",
    y = "PCo2"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold.italic", size = 11),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 12.2),
    plot.subtitle = element_text(size = 9.0, color = "#555555")
  )

pair_summary <- pair_tbl %>%
  group_by(species_label, pair_group) %>%
  summarise(
    median_distance = median(distance, na.rm = TRUE),
    q25 = quantile(distance, 0.25, na.rm = TRUE),
    q75 = quantile(distance, 0.75, na.rm = TRUE),
    n_pairs = n(),
    .groups = "drop"
  )

pair_plot <- ggplot(pair_tbl, aes(pair_group, distance, fill = pair_group)) +
  geom_violin(width = 0.88, alpha = 0.85, color = NA, scale = "width") +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.96, linewidth = 0.35, color = "#1E1E1E") +
  facet_wrap(~ species_label, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "Within Europe" = "#D9E8F5",
      "Within East Asia" = "#F2D4B5",
      "Between regions" = "#4E5B61"
    ),
    guide = "none"
  ) +
  labs(
    title = "B. Pairwise patristic distances",
    subtitle = "Between-region distances are shown against within-region backgrounds",
    x = NULL,
    y = "Patristic distance"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold.italic", size = 11),
    axis.text.x = element_text(angle = 18, hjust = 1),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5, color = "#555555")
  )

nn_plot <- nn_summary %>%
  select(species_label, observed_same_region, expected_random, enrichment, n_samples) %>%
  pivot_longer(
    c(observed_same_region, expected_random),
    names_to = "metric",
    values_to = "same_region_rate"
  ) %>%
  mutate(
    metric = recode(metric, observed_same_region = "Observed\nnearest neighbor", expected_random = "Random\nexpectation"),
    metric = factor(metric, levels = c("Random\nexpectation", "Observed\nnearest neighbor")),
    label = paste0(round(100 * same_region_rate, 1), "%")
  ) %>%
  ggplot(aes(metric, same_region_rate, fill = metric)) +
  geom_col(width = 0.62, color = "#242424", linewidth = 0.35) +
  geom_text(aes(label = label), vjust = -0.35, size = 3.4) +
  facet_wrap(~ species_label) +
  scale_fill_manual(
    values = c("Random\nexpectation" = "#D8D4CB", "Observed\nnearest neighbor" = "#1F6F78"),
    guide = "none"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    limits = c(0, 1.08),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "B. Same-region nearest-neighbor enrichment",
    subtitle = "Observed local region clustering is compared with the sample-composition expectation",
    x = NULL,
    y = "Same-region nearest-neighbor rate"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold.italic", size = 11),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(face = "bold", size = 12.2),
    plot.subtitle = element_text(size = 9.0, color = "#555555")
  )

summary_label <- region_summary %>%
  mutate(
    line = sprintf(
      "%s: dELPD/sample %.3f; z %.1f; n=%d",
      species_label,
      pglmm_minus_base_elpd_per_sample,
      pglmm_minus_base_z,
      n_used_samples
    )
  ) %>%
  pull(line) %>%
  paste(collapse = "    |    ")

layout_design <- c(
  area(t = 1, l = 1,  b = 1, r = 4),   
  area(t = 2, l = 1,  b = 2, r = 3)
)

overview_plot <- (pcoa_plot / nn_plot) +
  plot_layout(    design = layout_design,heights = c(1, 1)#, guides = "collect"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9.0, color = "#444444"),
    plot.caption = element_text(size = 8.5, color = "#555555", hjust = 0)
  )
if (isTRUE(render_files)) {
  ggsave(
    file.path(output_dir, "roseburia_region_structure_overview.pdf"),
    overview_plot,
    width = 8,
    height = 7.6,
    units = "in",
    bg = "white"
  )
}


if (isTRUE(render_files)) {
  ggsave(
    file.path(output_dir, "roseburia_region_pairwise_distance_distributions.pdf"),
    pair_plot,
    width = 10.8,
    height = 5.2,
    units = "in",
    bg = "white"
  )
}

draw_region_ring_tree <- function(data, panel_label = "", title_line = -0.45,
                                  footer_line = -0.15, show_footer = FALSE) {
  tree <- data$tree
  meta <- data$meta
  species_label <- data$spec$species_label
  summary_row <- data$region_summary
  
  tip_depth <- max(ape::node.depth.edgelength(tree)[seq_along(tree$tip.label)], na.rm = TRUE)
  plot.phylo(
    tree,
    type = "fan",
    show.tip.label = FALSE,
    no.margin = TRUE,
    edge.color = adjustcolor("#202020", alpha.f = 0.70),
    edge.width = 0.34,
    x.lim = c(-1.15, 1.15) * tip_depth,
    y.lim = c(-1.15, 1.15) * tip_depth
  )
  plot_info <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  
  tip_tbl <- tibble(
    tree_tip_label = tree$tip.label,
    x = plot_info$xx[seq_along(tree$tip.label)],
    y = plot_info$yy[seq_along(tree$tip.label)]
  ) %>%
    mutate(
      angle = atan2(y, x),
      angle = if_else(angle < 0, angle + 2 * pi, angle)
    ) %>%
    left_join(select(meta, tree_tip_label, region), by = "tree_tip_label")
  
  order_idx <- order(tip_tbl$angle)
  angles_sorted <- tip_tbl$angle[order_idx]
  angle_ext <- c(tail(angles_sorted, 1) - 2 * pi, angles_sorted, head(angles_sorted, 1) + 2 * pi)
  theta_lower <- (angle_ext[1:length(angles_sorted)] + angle_ext[2:(length(angles_sorted) + 1)]) / 2
  theta_upper <- (angle_ext[2:(length(angles_sorted) + 1)] + angle_ext[3:(length(angles_sorted) + 2)]) / 2
  
  tip_tbl$theta_lower <- NA_real_
  tip_tbl$theta_upper <- NA_real_
  tip_tbl$theta_lower[order_idx] <- theta_lower
  tip_tbl$theta_upper[order_idx] <- theta_upper
  
  r_tip <- sqrt(tip_tbl$x^2 + tip_tbl$y^2)
  r_max <- max(r_tip, na.rm = TRUE)
  ring_inner <- r_max * 1.012
  ring_outer <- r_max * 1.068
  
  invisible(lapply(seq_len(nrow(tip_tbl)), function(i) {
    draw_annulus_sector(
      theta_lower = tip_tbl$theta_lower[[i]],
      theta_upper = tip_tbl$theta_upper[[i]],
      r_inner = ring_inner,
      r_outer = ring_outer,
      col = unname(region_cols[as.character(tip_tbl$region[[i]])])
    )
  }))
  
  title(
    main = str_squish(paste(panel_label, species_label)),
    line = title_line,
    cex.main = 0.95,
    font.main = 3
  )
  
  if (isTRUE(show_footer)) {
    mtext(
      sprintf(
        "Europe-East Asia structure: dELPD/sample %.3f; z %.1f; n=%d",
        summary_row$pglmm_minus_base_elpd_per_sample,
        summary_row$pglmm_minus_base_z,
        summary_row$n_used_samples
      ),
      side = 1,
      line = footer_line,
      cex = 0.72,
      col = "#333333"
    )
  }
}

draw_ring_trees_panel <- function(species_data, region_cols) {
  stopifnot(is.list(species_data), length(species_data) >= 2)
  stopifnot(all(c("Europe", "East_Asia") %in% names(region_cols)))
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  layout(matrix(c(1, 2, 3), nrow = 1), widths = c(1, 1, 0.26))
  par(mar = c(1.0, 0.2, 1.4, 0.2), oma = c(0.4, 0.4, 2.0, 0.4), xpd = NA)
  
  draw_region_ring_tree(species_data[[1]], panel_label = "A.", title_line = -0.35)
  draw_region_ring_tree(species_data[[2]], panel_label = "B.", title_line = -0.35)
  
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend(
    "center",
    legend = c("Europe", "East Asia"),
    fill = unname(region_cols[c("Europe", "East_Asia")]),
    border = NA,
    bty = "n",
    title = "Region ring",
    cex = 0.9,
    y.intersp = 1.25
  )
  
  mtext(
    "Region-colored Roseburia strain trees",
    side = 3,
    outer = TRUE,
    line = 0.45,
    cex = 1.10,
    font = 2
  )
}

make_ring_trees_grob <- function(species_data, region_cols) {
  local_species_data <- species_data
  local_region_cols <- region_cols
  grid::grid.grabExpr(
    gridGraphics::grid.echo(function() {
      draw_ring_trees_panel(local_species_data, local_region_cols)
    })
  )
}

make_single_ring_tree_grob <- function(data, panel_label = "", title_line = -0.45) {
  local_data <- data
  local_panel_label <- panel_label
  local_title_line <- title_line
  grid::grid.grabExpr(
    gridGraphics::grid.echo(function() {
      old_par <- par(no.readonly = TRUE)
      on.exit(par(old_par), add = TRUE)
      par(mar = c(1.0, 0.2, 1.6, 0.2), oma = c(0.2, 0.2, 0.2, 0.2), xpd = NA)
      draw_region_ring_tree(local_data, panel_label = local_panel_label, title_line = local_title_line)
    })
  )
}

ring_trees_grob <- make_ring_trees_grob(species_data, region_cols)
rf4925_region_tree_grob <- make_single_ring_tree_grob(species_data[[1]], panel_label = "A.")
ri4940_region_tree_grob <- make_single_ring_tree_grob(species_data[[2]], panel_label = "B.")
roseburia_region_tree_grobs <- list(
  rf4925 = rf4925_region_tree_grob,
  ri4940 = ri4940_region_tree_grob,
  combined = ring_trees_grob
)


ring_trees_plot <- ggplotify::as.ggplot(ring_trees_grob)
rf4925_region_tree_plot <- ggplotify::as.ggplot(rf4925_region_tree_grob)
ri4940_region_tree_plot <- ggplotify::as.ggplot(ri4940_region_tree_grob)
roseburia_region_tree_plots <- list(
  rf4925 = rf4925_region_tree_plot,
  ri4940 = ri4940_region_tree_plot,
  combined = ring_trees_plot
)

render_ring_trees <- function(filename, device_fun, width = 11.2, height = 5.9) {
  device_fun(filename, width = width, height = height)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = FALSE)

  draw_ring_trees_panel(species_data, region_cols)
}

render_single_ring_tree <- function(data, filename, device_fun, width = 5.4, height = 5.4) {
  device_fun(filename, width = width, height = height)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = FALSE)

  par(mar = c(1.0, 0.2, 1.6, 0.2), oma = c(0.2, 0.2, 0.2, 0.2), xpd = NA)
  draw_region_ring_tree(data, panel_label = "", title_line = -0.45)
}

if (isTRUE(render_files)) {
  render_ring_trees(
    file.path(output_dir, "roseburia_region_ring_fan_trees.pdf"),
    pdf,
    width = 11.2,
    height = 5.9
  )

  invisible(lapply(species_data, function(data) {
    species_id <- data$spec$species_id
    render_single_ring_tree(
      data,
      file.path(single_tree_output_dir, sprintf("roseburia_region_ring_tree_%s.pdf", species_id)),
      pdf,
      width = single_tree_width,
      height = single_tree_height
    )
    render_single_ring_tree(
      data,
      file.path(single_tree_output_dir, sprintf("roseburia_region_ring_tree_%s.png", species_id)),
      function(filename, width, height) png(filename, width = width, height = height, units = "in", res = 360, bg = "white"),
      width = single_tree_width,
      height = single_tree_height
    )
  }))

  render_ring_trees(
    file.path(output_dir, "roseburia_region_ring_fan_trees.png"),
    function(filename, width, height) png(filename, width = width, height = height, units = "in", res = 360, bg = "white"),
    width = 11.2,
    height = 5.9
  )
}

write_tsv(region_summary, file.path(output_dir, "roseburia_region_structure_anpan_summary.tsv"))
write_tsv(pair_summary, file.path(output_dir, "roseburia_region_pairwise_distance_summary.tsv"))
write_tsv(nn_summary, file.path(output_dir, "roseburia_region_nearest_neighbor_summary.tsv"))

writeLines(
  c(
    "Roseburia region structure paper visualization",
    "",
    "Primary use:",
    "Use roseburia_region_structure_overview as the statistical/visual summary of Europe-East Asia strain background differences.",
    "Use roseburia_region_ring_fan_trees only as a phylogenetic-context panel; it is descriptive and should not be used as disease-replacement evidence.",
    "Individual tree files are written to individual_region_trees by default. Override that path with STRAINPHLAN_ROSEBURIA_REGION_SINGLE_TREE_OUTPUT_DIR.",
    "",
    "Interpretation:",
    "Both species show clear Europe-East Asia strain background structure in Anpan LOO comparisons and PCoA/pairwise-distance summaries."
  ),
  file.path(output_dir, "roseburia_region_structure_notes.txt")
)

if (isTRUE(render_files)) {
  message("Wrote Roseburia region-structure paper figures to: ", output_dir)
} else {
  message("Created Roseburia region-structure R plot objects; file rendering was disabled.")
}
