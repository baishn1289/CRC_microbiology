suppressPackageStartupMessages({
  library(readxl)
  library(limma)
  library(tidyverse)
  library(data.table)
  library(broom)
  library(purrr)
  library(vegan)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
})

project_root <- normalizePath(Sys.getenv("EOCRC_PROJECT_ROOT", unset = getwd()), winslash = "/", mustWork = FALSE)
analysis_dir <- file.path(project_root, "FUSCC_metabonomics_analysis")
analysis_version <- Sys.getenv("EOCRC_METABO_VERSION", unset = "v28")

resolve_single_path <- function(paths, label) {
  hits <- unique(as.character(paths))
  hits <- hits[!is.na(hits) & nzchar(hits)]
  hits <- hits[file.exists(hits)]
  if (length(hits) != 1) {
    stop(label, ": expected exactly one existing file, found ", length(hits), ".")
  }
  normalizePath(hits[[1]], winslash = "/", mustWork = TRUE)
}

meta_path <- resolve_single_path(
  Sys.glob(file.path(project_root, "SRA_metadata*", "A30_FUSCC_Gut", "metadata_for_revision.xlsx")),
  "Metabolomics metadata workbook"
)

output_dir <- file.path(analysis_dir, paste0("extended_metabonomics_results_", analysis_version))
plot_dir <- file.path(output_dir, "plots")
table_dir <- file.path(output_dir, "tables")

version_output_path <- function(path) {
  if (!is.character(path) || length(path) != 1 || is.na(path) || !nzchar(path)) {
    return(path)
  }
  gsub("(_v)\\d+(?=\\.[A-Za-z0-9]+$)", paste0("_", analysis_version), path, perl = TRUE)
}

write_output_csv <- function(x, file, ...) {
  utils::write.csv(x = x, file = version_output_path(file), ...)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_bw(base_size = 12))

group_levels <- c("LOControl", "LOCRC", "EOControl", "EOCRC")
group_palette <- c(
  LOControl = "#4E79A7",
  LOCRC = "#D1495B",
  EOControl = "#76B7B2",
  EOCRC = "#EDAE49"
)

panel_levels <- c(
  "arginine_ornithine_routing_panel",
  "amino_acid_putrefaction_core_panel",
  "amino_acid_putrefaction_auxiliary_panel",
  "tryptophan_indole_branch_panel",
  "polyamine_lysine_putrefaction_panel",
  "biogenic_amine_panel",
  "indole_catabolism_panel",
  "aromatic_host_cometabolite_panel",
  "injurious_biogenic_amine_subpanel",
  "injurious_parent_phenol_subpanel",
  "injurious_host_processed_aromatic_toxin_subpanel"
)

panel_labels <- c(
  arginine_ornithine_routing_panel = "Arginine-ornithine routing axis",
  amino_acid_putrefaction_core_panel = "Amino-acid putrefaction core",
  amino_acid_putrefaction_auxiliary_panel = "Extended: amino-acid putrefaction auxiliary",
  tryptophan_indole_branch_panel = "Extended: tryptophan-indole branch",
  polyamine_lysine_putrefaction_panel = "Extended: polyamine / lysine putrefaction",
  biogenic_amine_panel = "Extended: biogenic amines",
  indole_catabolism_panel = "Extended: indole catabolism",
  aromatic_host_cometabolite_panel = "Extended: aromatic host-microbe co-metabolites",
  injurious_biogenic_amine_subpanel = "Extended: injurious biogenic amines",
  injurious_parent_phenol_subpanel = "Extended: parent phenolic toxins",
  injurious_host_processed_aromatic_toxin_subpanel = "Extended: host-processed aromatic toxins"
)

panel_tier_map <- c(
  arginine_ornithine_routing_panel = "strict_manuscript_aligned",
  amino_acid_putrefaction_core_panel = "strict_manuscript_aligned",
  amino_acid_putrefaction_auxiliary_panel = "extended_class_based",
  tryptophan_indole_branch_panel = "extended_class_based",
  polyamine_lysine_putrefaction_panel = "extended_class_based",
  biogenic_amine_panel = "extended_class_based",
  indole_catabolism_panel = "extended_class_based",
  aromatic_host_cometabolite_panel = "extended_class_based",
  injurious_biogenic_amine_subpanel = "extended_class_based",
  injurious_parent_phenol_subpanel = "extended_class_based",
  injurious_host_processed_aromatic_toxin_subpanel = "extended_class_based"
)

panel_tier_labels <- c(
  strict_manuscript_aligned = "Strict manuscript-aligned validation",
  extended_class_based = "Extended class-based support"
)

contrast_display_order <- c(
  "LOCRC_vs_LOControl_ageadj",
  "EOCRC_vs_EOControl_ageadj",
  "EOCRC_vs_LOCRC_crc_age_stage_site_adj"
)

main_figure_strict_panels <- c(
  "arginine_ornithine_routing_panel",
  "amino_acid_putrefaction_core_panel"
)

main_figure_panel_short_labels <- c(
  arginine_ornithine_routing_panel = "Routing axis",
  amino_acid_putrefaction_core_panel = "Putrefaction core"
)

main_figure_contrast_order <- c(
  "LOCRC_vs_LOControl_pooled_ageadj",
  "EOCRC_vs_EOControl_pooled_ageadj",
  "Disease_age_interaction_pooled_ageadj"
)

main_figure_contrast_labels <- c(
  LOCRC_vs_LOControl_pooled_ageadj = "LO effect",
  EOCRC_vs_EOControl_pooled_ageadj = "EO effect",
  Disease_age_interaction_pooled_ageadj = "EO > LO"
)

main_figure_forest_contrast_order <- c(
  "LOCRC_vs_LOControl_ageadj",
  "EOCRC_vs_EOControl_ageadj",
  "EO_vs_LO_delta_beta_ageadj"
)

main_figure_forest_contrast_labels <- c(
  LOCRC_vs_LOControl_ageadj = "LO effect",
  EOCRC_vs_EOControl_ageadj = "EO effect",
  EO_vs_LO_delta_beta_ageadj = "Delta beta (EO - LO)"
)

main_figure_metabolite_order <- tribble(
  ~panel, ~canonical_metabolite, ~metabolite_order, ~metabolite_group,
  "arginine_ornithine_routing_panel", "Arginine", 1, "Routing axis",
  "arginine_ornithine_routing_panel", "Citrulline", 2, "Routing axis",
  "arginine_ornithine_routing_panel", "Ornithine_related", 3, "Routing axis",
  "arginine_ornithine_routing_panel", "N-Acetylornithine", 4, "Routing axis",
  "amino_acid_putrefaction_core_panel", "5-Aminopentanoate", 5, "Putrefaction core",
  "amino_acid_putrefaction_core_panel", "N-Acetylputrescine", 6, "Putrefaction core",
  "amino_acid_putrefaction_core_panel", "N-Carbamoylputrescine", 7, "Putrefaction core"
)

support_strip_target_panel <- "amino_acid_putrefaction_core_panel"
support_strip_species_tiers <- c(
  "tier1_panel_disease_crc_support",
  "tier2_panel_disease"
)
support_strip_fdr_cutoff <- 0.05

support_strip_contrast_labels <- c(
  EO_disease_effect = "EO disease",
  EO_putrefactive_panel_association = "EO core assoc",
  LO_disease_effect = "LO disease",
  LO_putrefactive_panel_association = "LO core assoc"
)

fdr_cutoff <- 0.05
effect_cutoff <- log2(1.2)
min_present_fraction <- 0.70
max_missing_fraction <- 0.50
species_prevalence_min <- 0.10
species_group_prevalence_min <- 0.15
function_prevalence_min <- 0.10
function_group_prevalence_min <- 0.10
species_candidate_fdr <- 0.05
species_sensitivity_p_cutoff <- 0.05
top_species_hits_per_direction <- 10
taxonomic_enrichment_min_taxon_size <- 3
eo_target_panels <- c(
  "amino_acid_putrefaction_core_panel",
  "amino_acid_putrefaction_auxiliary_panel",
  "tryptophan_indole_branch_panel",
  "injurious_biogenic_amine_subpanel",
  "injurious_parent_phenol_subpanel",
  "injurious_host_processed_aromatic_toxin_subpanel"
)
manuscript_metacyc_anchor_ids <- c(
  "PWY-5004",
  "PWY-8187",
  "AST-PWY",
  "POLYAMSYN-PWY",
  "ARG+POLYAMINE-SYN",
  "MET-SAM-PWY",
  "HOMOSER-METSYN-PWY",
  "PWY-6328"
)
manuscript_ec_anchor_ids <- c(
  "3.5.3.1",
  "1.4.1.11",
  "1.4.1.4",
  "6.3.1.2",
  "4.1.99.1",
  "3.5.3.23",
  "3.5.3.6"
)

safe_file_tag <- function(x) {
  x %>%
    gsub("[^A-Za-z0-9]+", "_", .) %>%
    gsub("^_+|_+$", "", .)
}

extract_taxon_rank <- function(feature_ids, rank_code) {
  rank_pattern <- paste0(".*\\|", rank_code, "__([^|]+).*")
  ifelse(
    grepl(paste0("\\|", rank_code, "__"), feature_ids),
    sub(rank_pattern, "\\1", feature_ids),
    NA_character_
  )
}

build_species_taxonomy_annotation <- function(feature_ids) {
  tibble(feature_id = feature_ids) %>%
    mutate(
      phylum_label = extract_taxon_rank(feature_id, "p"),
      family_label = extract_taxon_rank(feature_id, "f"),
      genus_label = extract_taxon_rank(feature_id, "g"),
      species_label = extract_taxon_rank(feature_id, "s"),
      short_label = coalesce(species_label, feature_id)
    )
}

significance_label <- function(p_value) {
  case_when(
    is.na(p_value) ~ "NA",
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

save_gg <- function(plot_obj, filename, width = 10, height = 7) {
  output_file <- version_output_path(sub("\\.[^.]+$", ".pdf", filename))
  tryCatch(
    {
      ggsave(
        filename = output_file,
        plot = plot_obj,
        width = width,
        height = height,
        dpi = 300,
        bg = "white"
      )
    },
    error = function(e) {
      fallback_file <- sub("\\.pdf$", "_rerun.pdf", output_file, ignore.case = TRUE)
      message("Primary PDF write failed for ", output_file, "; retrying as ", fallback_file)
      ggsave(
        filename = fallback_file,
        plot = plot_obj,
        width = width,
        height = height,
        dpi = 300,
        bg = "white"
      )
    }
  )
}

run_betadisper_permutest <- function(dispersion_obj, permutations = 999) {
  if (!inherits(dispersion_obj, "betadisper")) {
    stop("run_betadisper_permutest() expects a vegan::betadisper object.")
  }

  # Call the vegan betadisper S3 method explicitly so the script does not depend
  # on the current interactive session's generic/method dispatch state.
  permutest_betadisper <- getS3method("permutest", "betadisper", envir = asNamespace("vegan"))
  permutest_betadisper(dispersion_obj, permutations = permutations)
}

pareto_scale_matrix <- function(mat) {
  mat <- as.matrix(mat)
  center_vec <- colMeans(mat, na.rm = TRUE)
  sd_vec <- apply(mat, 2, sd, na.rm = TRUE)
  scale_vec <- sqrt(sd_vec)
  scale_vec[is.na(scale_vec) | scale_vec == 0] <- 1
  scaled <- sweep(mat, 2, center_vec, "-")
  sweep(scaled, 2, scale_vec, "/")
}

z_score_features <- function(feature_by_sample_mat) {
  z_mat <- t(apply(feature_by_sample_mat, 1, function(x) {
    x <- as.numeric(x)
    if (all(is.na(x))) {
      return(rep(NA_real_, length(x)))
    }
    x_sd <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(x_sd) || x_sd == 0) {
      out <- rep(0, length(x))
      out[is.na(x)] <- NA_real_
      return(out)
    }
    (x - mean(x, na.rm = TRUE)) / x_sd
  }))
  rownames(z_mat) <- rownames(feature_by_sample_mat)
  colnames(z_mat) <- colnames(feature_by_sample_mat)
  z_mat
}

safe_colmeans_available <- function(mat) {
  mat <- as.matrix(mat)
  observed_n <- colSums(!is.na(mat))
  out <- colMeans(mat, na.rm = TRUE)
  out[observed_n == 0] <- NA_real_
  out
}

read_mode_file <- function(path, mode_name) {
  raw <- fread(path)
  sample_cols <- grep("^[oy]", names(raw), value = TRUE)
  raw %>%
    select(all_of(c("metabolite_identification", sample_cols))) %>%
    mutate(mode = mode_name)
}

score_mode_candidates <- function(lcms_df) {
  sample_cols <- setdiff(colnames(lcms_df), c("metabolite_identification", "mode"))

  lcms_df %>%
    rowwise() %>%
    mutate(
      n_detected = sum(!is.na(c_across(all_of(sample_cols))) & c_across(all_of(sample_cols)) > 0),
      missing_rate = mean(is.na(c_across(all_of(sample_cols)))),
      median_signal = median(c_across(all_of(sample_cols)), na.rm = TRUE),
      prefer_mode = case_when(
        grepl(
          "acid$|Acid$|butyr|propion|succin|lactate|acetate|citrate|fumar|malate",
          metabolite_identification,
          ignore.case = TRUE
        ) ~ "negative",
        grepl(
          "amine|pyridine|pyrid|ornithine|lysine|arginine|putrescine|cadaverine",
          metabolite_identification,
          ignore.case = TRUE
        ) ~ "positive",
        TRUE ~ NA_character_
      ),
      mode_rank = case_when(
        !is.na(prefer_mode) & mode == prefer_mode ~ 1,
        !is.na(prefer_mode) & mode != prefer_mode ~ 2,
        TRUE ~ 3
      )
    ) %>%
    ungroup() %>%
    group_by(metabolite_identification) %>%
    arrange(mode_rank, missing_rate, desc(median_signal), .by_group = TRUE) %>%
    mutate(
      selection_rank = row_number(),
      selected_mode = selection_rank == 1,
      selection_reason = case_when(
        selected_mode & !is.na(prefer_mode) & mode == prefer_mode ~ "preferred_ion_mode",
        selected_mode & !is.na(prefer_mode) & mode != prefer_mode ~ "fallback_nonpreferred_mode",
        selected_mode ~ "best_missingness_then_signal",
        TRUE ~ "not_selected"
      )
    ) %>%
    ungroup()
}

deduplicate_modes <- function(lcms_df) {
  scored_df <- if ("selection_rank" %in% names(lcms_df)) lcms_df else score_mode_candidates(lcms_df)
  scored_df %>%
    filter(selected_mode) %>%
    arrange(metabolite_identification)
}

extract_limma_results <- function(fit_obj, contrast_names, effect_cutoff, fdr_cutoff) {
  map_dfr(contrast_names, function(contrast_name) {
    topTable(fit_obj, coef = contrast_name, number = Inf, sort.by = "P") %>%
      rownames_to_column("metabolite_identification") %>%
      mutate(
        contrast = contrast_name,
        significant_fdr = adj.P.Val < fdr_cutoff,
        significant_fdr_fc = adj.P.Val < fdr_cutoff & abs(logFC) >= effect_cutoff,
        significant_fdr_0_10 = FALSE,
        significant_fdr_fc_0_10 = FALSE,
        direction = case_when(
          logFC > 0 ~ "higher",
          logFC < 0 ~ "lower",
          TRUE ~ "no_change"
        )
      )
  })
}

fit_limma_contrasts <- function(
    expr_mat,
    sample_meta,
    design_formula,
    contrast_definitions,
    model_label,
    test_type = c("treat", "ebayes"),
    lfc_threshold = 0) {
  test_type <- match.arg(test_type)

  design <- model.matrix(design_formula, data = sample_meta)
  colnames(design) <- make.names(colnames(design))

  contrast_matrix <- do.call(
    limma::makeContrasts,
    c(as.list(contrast_definitions), list(levels = design))
  )

  fit0 <- lmFit(expr_mat, design)
  fit1 <- contrasts.fit(fit0, contrast_matrix)

  if (test_type == "treat") {
    fit2 <- treat(fit1, lfc = lfc_threshold, robust = TRUE)
    result_tbl <- map_dfr(colnames(contrast_matrix), function(contrast_name) {
      topTreat(fit2, coef = contrast_name, number = Inf, sort.by = "P") %>%
        rownames_to_column("metabolite_identification") %>%
        mutate(
          contrast = contrast_name,
          model_label = model_label,
          test_type = test_type,
          lfc_threshold = lfc_threshold,
          significant_fdr = adj.P.Val < fdr_cutoff,
          significant_fdr_fc = adj.P.Val < fdr_cutoff,
          significant_fdr_0_10 = FALSE,
          significant_fdr_fc_0_10 = FALSE,
          direction = case_when(
            logFC > 0 ~ "higher",
            logFC < 0 ~ "lower",
            TRUE ~ "no_change"
          )
        )
    })
  } else {
    fit2 <- eBayes(fit1, robust = TRUE)
    result_tbl <- extract_limma_results(
      fit_obj = fit2,
      contrast_names = colnames(contrast_matrix),
      effect_cutoff = lfc_threshold,
      fdr_cutoff = fdr_cutoff
    ) %>%
      mutate(
        model_label = model_label,
        test_type = test_type,
        lfc_threshold = lfc_threshold
      )
  }

  list(
    results = result_tbl,
    fit = fit2,
    design = design,
    contrast_matrix = contrast_matrix
  )
}

select_volcano_labels <- function(df, targeted_membership, n_global = 6, n_targeted = 8) {
  top_global <- df %>%
    slice_min(order_by = adj.P.Val, n = n_global, with_ties = FALSE) %>%
    pull(metabolite_identification)

  top_targeted <- df %>%
    semi_join(targeted_membership, by = "metabolite_identification") %>%
    slice_min(order_by = adj.P.Val, n = n_targeted, with_ties = FALSE) %>%
    pull(metabolite_identification)

  unique(c(top_global, top_targeted))
}

plot_volcano <- function(df, contrast_name, targeted_membership) {
  label_features <- select_volcano_labels(df, targeted_membership)

  df <- df %>%
    mutate(
      neg_log10_fdr = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
      label = ifelse(metabolite_identification %in% label_features, metabolite_identification, NA_character_)
    )

  ggplot(df, aes(x = logFC, y = neg_log10_fdr)) +
    geom_point(
      aes(color = significant_fdr_fc),
      alpha = 0.75,
      size = 1.6
    ) +
    geom_vline(xintercept = c(-effect_cutoff, effect_cutoff), linetype = 2, color = "grey50") +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = 2, color = "grey50") +
    geom_text_repel(
      data = df %>% filter(!is.na(label)),
      aes(label = label),
      size = 3.1,
      max.overlaps = Inf,
      box.padding = 0.25,
      point.padding = 0.15,
      segment.color = "grey60"
    ) +
    scale_color_manual(values = c("FALSE" = "grey72", "TRUE" = "#C73E1D")) +
    labs(
      title = contrast_name,
      x = "log2 fold-change",
      y = "-log10(FDR)",
      color = paste0("FDR < ", fdr_cutoff, "\n& |log2FC| >= ", round(effect_cutoff, 3))
    ) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
}

# final canonicalized implementation
build_targeted_panel_membership <- function(feature_names) {
  feature_tbl <- tibble(metabolite_identification = feature_names)

  canonical_map <- tribble(
    ~metabolite_identification, ~canonical_metabolite,
    "Arginine", "Arginine",
    "Citrulline", "Citrulline",
    "L-Citrulline", "Citrulline",
    "D-Ornithine", "Ornithine_related",
    "DL-Ornithine", "Ornithine_related",
    "N-.alpha.-Acetyl-L-ornithine", "N-Acetylornithine",
    "N-Acetylornithine", "N-Acetylornithine",
    "N-Acetylputrescine", "N-Acetylputrescine",
    "5-Aminopentanoic acid", "5-Aminopentanoate",
    "5-Aminovaleric acid", "5-Aminopentanoate",
    "5-Aminopentanal", "5-Aminopentanal",
    "Indole-3-acetaldehyde", "Indole-3-acetaldehyde",
    "Indoleacetaldehyde", "Indole-3-acetaldehyde",
    "Indole-3-carboxyaldehyde", "Indole-3-carboxyaldehyde",
    "3-Formylindole", "Indole-3-carboxyaldehyde",
    "N-Carbamoylputrescine", "N-Carbamoylputrescine",
    "N-Acetylcadaverine", "N-Acetylcadaverine",
    "Lysine", "Lysine",
    "Tryptophan", "Tryptophan",
    "Urocanic acid", "Urocanic acid",
    "Gamma-Aminobutyric acid", "Gamma-Aminobutyric acid",
    "Histamine", "Histamine",
    "N-Acetylhistamine", "N-Acetylhistamine",
    "Tyramine", "Tyramine",
    "Tyramine-O-sulfate", "Tyramine-O-sulfate",
    "N-Methyltyramine", "N-Methyltyramine",
    "Indole", "Indole",
    "Indoleacetic acid", "Indole-3-acetic acid",
    "3-Indoleacetic acid", "Indole-3-acetic acid",
    "Indolelactic acid", "Indole-3-lactate",
    "Indole-3-acetamide", "Indole-3-acetamide",
    "Indole-3-carboxylic acid", "Indole-3-carboxylic acid",
    "1H-Indole-3-carboxylic acid", "Indole-3-carboxylic acid",
    "p-Cresol", "p-Cresol",
    "p-Cresol sulfate", "p-Cresol sulfate",
    "Phenol", "Phenol",
    "Phenol sulphate", "Phenol sulphate",
    "Indoxyl sulfate", "Indoxyl sulfate",
    "Phenylacetic acid", "Phenylacetic acid",
    "Phenylacetyl-L-glutamine", "Phenylacetyl-L-glutamine",
    "N-Phenylacetylglutamic acid", "N-Phenylacetylglutamic acid",
    "4-Hydroxyphenylacetylglutamine", "4-Hydroxyphenylacetylglutamine"
  )

  panel_catalog <- tribble(
    ~panel, ~canonical_metabolite,
    "arginine_ornithine_routing_panel", "Arginine",
    "arginine_ornithine_routing_panel", "Citrulline",
    "arginine_ornithine_routing_panel", "Ornithine_related",
    "arginine_ornithine_routing_panel", "N-Acetylornithine",
    "amino_acid_putrefaction_core_panel", "5-Aminopentanoate",
    "amino_acid_putrefaction_core_panel", "N-Acetylputrescine",
    "amino_acid_putrefaction_core_panel", "N-Carbamoylputrescine",
    "amino_acid_putrefaction_auxiliary_panel", "N-Acetylcadaverine",
    "amino_acid_putrefaction_auxiliary_panel", "5-Aminopentanal",
    "amino_acid_putrefaction_auxiliary_panel", "Gamma-Aminobutyric acid",
    "amino_acid_putrefaction_auxiliary_panel", "Urocanic acid",
    "tryptophan_indole_branch_panel", "Indole",
    "tryptophan_indole_branch_panel", "Indole-3-acetaldehyde",
    "tryptophan_indole_branch_panel", "Indole-3-carboxyaldehyde",
    "tryptophan_indole_branch_panel", "Indole-3-acetic acid",
    "tryptophan_indole_branch_panel", "Indole-3-lactate",
    "polyamine_lysine_putrefaction_panel", "N-Acetylputrescine",
    "polyamine_lysine_putrefaction_panel", "N-Carbamoylputrescine",
    "polyamine_lysine_putrefaction_panel", "N-Acetylcadaverine",
    "polyamine_lysine_putrefaction_panel", "5-Aminopentanoate",
    "biogenic_amine_panel", "Histamine",
    "biogenic_amine_panel", "N-Acetylhistamine",
    "biogenic_amine_panel", "Tyramine",
    "biogenic_amine_panel", "Tyramine-O-sulfate",
    "biogenic_amine_panel", "N-Methyltyramine",
    "indole_catabolism_panel", "Indole",
    "indole_catabolism_panel", "Indole-3-acetaldehyde",
    "indole_catabolism_panel", "Indole-3-acetamide",
    "indole_catabolism_panel", "Indole-3-carboxyaldehyde",
    "indole_catabolism_panel", "Indole-3-acetic acid",
    "indole_catabolism_panel", "Indole-3-lactate",
    "indole_catabolism_panel", "Indole-3-carboxylic acid",
    "aromatic_host_cometabolite_panel", "Phenylacetic acid",
    "aromatic_host_cometabolite_panel", "Phenylacetyl-L-glutamine",
    "aromatic_host_cometabolite_panel", "N-Phenylacetylglutamic acid",
    "aromatic_host_cometabolite_panel", "4-Hydroxyphenylacetylglutamine",
    "aromatic_host_cometabolite_panel", "Indoxyl sulfate",
    "injurious_biogenic_amine_subpanel", "Histamine",
    "injurious_biogenic_amine_subpanel", "Tyramine",
    "injurious_biogenic_amine_subpanel", "N-Methyltyramine",
    "injurious_parent_phenol_subpanel", "p-Cresol",
    "injurious_parent_phenol_subpanel", "Phenol",
    "injurious_host_processed_aromatic_toxin_subpanel", "p-Cresol sulfate",
    "injurious_host_processed_aromatic_toxin_subpanel", "Phenol sulphate",
    "injurious_host_processed_aromatic_toxin_subpanel", "Indoxyl sulfate"
  )

  feature_tbl %>%
    inner_join(canonical_map, by = "metabolite_identification") %>%
    inner_join(panel_catalog, by = "canonical_metabolite", relationship = "many-to-many") %>%
    distinct(metabolite_identification, canonical_metabolite, panel, .keep_all = TRUE) %>%
    mutate(
      panel = factor(panel, levels = panel_levels),
      panel_label = recode(panel, !!!panel_labels),
      panel_tier = recode(panel, !!!panel_tier_map),
      panel_tier_label = recode(panel_tier, !!!panel_tier_labels)
    ) %>%
    arrange(panel, canonical_metabolite, metabolite_identification)
}

build_focus_metabolite_review <- function(feature_names) {
  tibble(
    metabolite_identification = c(
      "5-Aminopentanoic acid",
      "5-Aminovaleric acid",
      "N-Acetylputrescine",
      "N-Carbamoylputrescine",
      "N-Acetylcadaverine",
      "5-Aminopentanal",
      "Gamma-Aminobutyric acid",
      "Urocanic acid",
      "Indole",
      "Indole-3-acetaldehyde",
      "Indole-3-carboxyaldehyde",
      "3-Formylindole",
      "Indoleacetic acid",
      "3-Indoleacetic acid",
      "Indolelactic acid",
      "Indoxyl sulfate",
      "4-Hydroxyphenylacetylglutamine",
      "Indole-3-carbinol",
      "5-Hydroxyindole-3-acetic acid",
      "L-Histidine",
      "L-Phenylalanine"
    ),
    final_decision = c(
      "keep_as_canonical_5_aminopentanoate",
      "collapse_into_5_aminopentanoate",
      "keep_putrefaction_core_tier1",
      "keep_putrefaction_core_tier1",
      "keep_putrefaction_auxiliary_tier2",
      "keep_putrefaction_auxiliary_tier2",
      "keep_putrefaction_auxiliary_tier2",
      "keep_putrefaction_auxiliary_tier2",
      "keep_indole_branch_supplementary",
      "keep_indole_branch_supplementary",
      "keep_indole_branch_supplementary",
      "collapse_into_indole_3_carboxyaldehyde",
      "keep_indole_branch_supplementary",
      "collapse_into_indole_3_acetic_acid",
      "keep_indole_branch_supplementary",
      "keep_extended_only",
      "keep_extended_only",
      "exclude_from_curated_panels",
      "exclude_from_curated_panels",
      "exclude_from_curated_panels",
      "exclude_from_curated_panels"
    ),
    representative_panel = c(
      "amino_acid_putrefaction_core_panel",
      "amino_acid_putrefaction_core_panel",
      "amino_acid_putrefaction_core_panel",
      "amino_acid_putrefaction_auxiliary_panel",
      "amino_acid_putrefaction_auxiliary_panel",
      "amino_acid_putrefaction_auxiliary_panel",
      "amino_acid_putrefaction_auxiliary_panel",
      "tryptophan_indole_branch_panel",
      "tryptophan_indole_branch_panel",
      "tryptophan_indole_branch_panel",
      "tryptophan_indole_branch_panel",
      "tryptophan_indole_branch_panel",
      "tryptophan_indole_branch_panel",
      "injurious_host_processed_aromatic_toxin_subpanel",
      "aromatic_host_cometabolite_panel",
      NA,
      NA,
      NA,
      NA,
      NA,
      NA
    ),
    biochemical_comment = c(
      "Primary lysine/cadaverine branch marker aligned with the amino-acid putrefaction mainline; retained as a first-tier core feature.",
      "Synonym-level duplicate of 5-aminopentanoic acid; retained in raw results but collapsed for panel scoring and heatmaps.",
      "Detectable putrescine-related spillover marker aligned with the arginine-ornithine-polyamine axis; retained as a first-tier core feature.",
      "Putrescine/agmatine-related intermediate retained as a same-axis support marker in the strict putrefaction core.",
      "Cadaverine itself is not detected, so this acetylated cadaverine signal is retained only as a second-tier auxiliary lysine-branch marker.",
      "Near-downstream lysine/cadaverine branch aldehyde retained only as a second-tier auxiliary marker.",
      "Broad amino-acid decarboxylation/putrefaction-adjacent readout retained only as a second-tier auxiliary marker because of lower pathway specificity.",
      "Direct product of histidine ammonia-lyase (EC 4.3.1.3); retained as an auxiliary histidine-derived ammonia-release marker rather than a core feature.",
      "Direct product of tryptophan indole-lyase (EC 4.1.99.1); retained only in the supplementary tryptophan-indole branch.",
      "Tryptophan-derived indole aldehyde marker retained only in the supplementary indole branch rather than the main amino-acid putrefaction core.",
      "Tryptophan-derived indole aldehyde marker retained only in the supplementary indole branch rather than the main amino-acid putrefaction core.",
      "Synonym-level duplicate of indole-3-carboxyaldehyde rather than an independent metabolite signal.",
      "Tryptophan-derived indole acid marker retained only in the supplementary indole branch.",
      "Synonym-level duplicate of indoleacetic acid and collapsed into indole-3-acetic acid.",
      "Tryptophan-derived indole lactate marker retained only in the supplementary indole branch.",
      "Strong gut microbiota-host co-metabolite toxin signal, but it extends beyond the exact manuscript EC/pathway anchors, so it is kept only in the extended panels.",
      "Aromatic host-microbe conjugate consistent with broadened aromatic fermentation support, but not specific enough for strict manuscript-aligned validation.",
      "Likely dominated by dietary cruciferous origin rather than a robust fecal putrefaction readout; excluded to avoid overclaiming.",
      "Predominantly a host serotonin catabolite rather than a clean microbial putrefaction marker; excluded.",
      "Represents broad amino-acid pool status rather than the arginine-ornithine precursor axis described in the manuscript; excluded from curated panels.",
      "Represents a broad amino-acid reservoir rather than a direct nitrogen-flux-diversion anchor; excluded from curated panels."
    )
  ) %>%
    mutate(
      present_in_dataset = metabolite_identification %in% feature_names
    )
}

build_panel_score_matrix <- function(feature_mat, panel_membership) {
  z_mat <- z_score_features(feature_mat)

  canonical_feature_map <- panel_membership %>%
    distinct(metabolite_identification, canonical_metabolite)

  canonical_feature_list <- split(
    canonical_feature_map$metabolite_identification,
    canonical_feature_map$canonical_metabolite
  )

  canonical_mat <- vapply(
    canonical_feature_list,
    function(mets) safe_colmeans_available(z_mat[unique(mets), , drop = FALSE]),
    numeric(ncol(z_mat))
  )
  canonical_mat <- t(canonical_mat)
  rownames(canonical_mat) <- names(canonical_feature_list)
  colnames(canonical_mat) <- colnames(z_mat)

  panel_list <- split(panel_membership$canonical_metabolite, panel_membership$panel)
  panel_list <- panel_list[lengths(panel_list) > 0]

  panel_mat <- vapply(
    panel_list,
    function(mets) safe_colmeans_available(canonical_mat[unique(mets), , drop = FALSE]),
    numeric(ncol(canonical_mat))
  )
  panel_mat <- t(panel_mat)
  rownames(panel_mat) <- names(panel_list)
  colnames(panel_mat) <- colnames(canonical_mat)
  panel_mat
}

build_targeted_canonical_long <- function(long_df, panel_membership) {
  canonical_map <- panel_membership %>%
    distinct(metabolite_identification, canonical_metabolite)

  long_df %>%
    inner_join(canonical_map, by = "metabolite_identification") %>%
    group_by(canonical_metabolite, Sample_ID, Group4, Disease, Age_group, Sex, BMI) %>%
    summarise(log2_intensity = mean(log2_intensity, na.rm = TRUE), .groups = "drop")
}

build_panel_score_coverage_table <- function(raw_long_df, panel_membership) {
  canonical_map <- panel_membership %>%
    distinct(panel, panel_label, panel_tier, panel_tier_label, metabolite_identification, canonical_metabolite)

  raw_long_df %>%
    inner_join(canonical_map, by = "metabolite_identification", relationship = "many-to-many") %>%
    mutate(has_value = !is.na(intensity) & intensity > 0) %>%
    group_by(Sample_ID, Group4, Disease, Age_group, Sex, BMI, panel, panel_label, panel_tier, panel_tier_label, canonical_metabolite) %>%
    summarise(has_value = any(has_value), .groups = "drop") %>%
    group_by(Sample_ID, Group4, Disease, Age_group, Sex, BMI, panel, panel_label, panel_tier, panel_tier_label) %>%
    summarise(
      n_observed = sum(has_value),
      n_total = dplyr::n(),
      observed_fraction = n_observed / n_total,
      .groups = "drop"
    )
}

summarise_panel_score_coverage <- function(panel_coverage_tbl) {
  panel_coverage_tbl %>%
    group_by(panel, panel_label, panel_tier, panel_tier_label, Group4, Disease, Age_group) %>%
    summarise(
      n_samples = dplyr::n(),
      mean_observed_fraction = mean(observed_fraction, na.rm = TRUE),
      median_observed_fraction = median(observed_fraction, na.rm = TRUE),
      min_observed_fraction = min(observed_fraction, na.rm = TRUE),
      complete_coverage_rate = mean(observed_fraction >= 1, na.rm = TRUE),
      .groups = "drop"
    )
}

# final canonicalized heatmap implementation
plot_targeted_heatmap <- function(long_df, significant_targeted_panel_df, id_col = "canonical_metabolite") {
  if (nrow(significant_targeted_panel_df) == 0) {
    return(NULL)
  }

  order_tbl <- significant_targeted_panel_df %>%
    distinct(panel, panel_label, panel_tier, panel_tier_label, .data[[id_col]]) %>%
    arrange(panel, .data[[id_col]]) %>%
    mutate(panel_feature = paste(panel, .data[[id_col]], sep = "||"))

  heatmap_df <- long_df %>%
    filter(.data[[id_col]] %in% order_tbl[[id_col]]) %>%
    group_by(.data[[id_col]], Group4) %>%
    summarise(mean_log2 = mean(log2_intensity, na.rm = TRUE), .groups = "drop") %>%
    group_by(.data[[id_col]]) %>%
    mutate(group_mean_z = as.numeric(scale(mean_log2))) %>%
    ungroup() %>%
    inner_join(order_tbl, by = id_col, relationship = "many-to-many") %>%
    mutate(
      Group4 = factor(Group4, levels = group_levels),
      panel_feature = factor(panel_feature, levels = rev(order_tbl$panel_feature))
    )

  ggplot(heatmap_df, aes(x = Group4, y = panel_feature, fill = group_mean_z)) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(panel_label ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(
      low = "#1D4E89",
      mid = "white",
      high = "#B22222",
      midpoint = 0
    ) +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) +
    labs(
      title = "Group-mean z-scores for significant curated metabolites",
      x = NULL,
      y = NULL,
      fill = "z-score"
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      strip.text.y = element_text(angle = 0)
    )
}

plot_panel_member_heatmap <- function(long_df, panel_lookup, feature_results, id_col = "canonical_metabolite") {
  order_tbl <- panel_lookup %>%
    distinct(panel, panel_label, panel_tier, panel_tier_label, .data[[id_col]]) %>%
    arrange(panel, .data[[id_col]]) %>%
    mutate(panel_feature = paste(panel, .data[[id_col]], sep = "||"))

  if (nrow(order_tbl) == 0) {
    return(NULL)
  }

  sig_col <- if ("significant_fdr" %in% names(feature_results)) {
    "significant_fdr"
  } else if ("any_significant_fdr" %in% names(feature_results)) {
    "any_significant_fdr"
  } else {
    stop("feature_results must contain either significant_fdr or any_significant_fdr.")
  }

  feature_sig <- feature_results %>%
    group_by(.data[[id_col]]) %>%
    summarise(any_fdr_sig = any(.data[[sig_col]], na.rm = TRUE), .groups = "drop")

  heatmap_df <- long_df %>%
    filter(.data[[id_col]] %in% order_tbl[[id_col]]) %>%
    group_by(.data[[id_col]], Group4) %>%
    summarise(mean_log2 = mean(log2_intensity, na.rm = TRUE), .groups = "drop") %>%
    group_by(.data[[id_col]]) %>%
    mutate(group_mean_z = as.numeric(scale(mean_log2))) %>%
    ungroup() %>%
    inner_join(order_tbl, by = id_col, relationship = "many-to-many") %>%
    left_join(feature_sig, by = id_col) %>%
    mutate(
      Group4 = factor(Group4, levels = group_levels),
      panel_feature = factor(panel_feature, levels = rev(order_tbl$panel_feature)),
      sig_flag = ifelse(any_fdr_sig, "FDR<0.05 in >=1 contrast", "Not FDR-significant")
    )

  ggplot(heatmap_df, aes(x = Group4, y = panel_feature, fill = group_mean_z#, alpha = sig_flag
                         )) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(
       panel_label~.,
      scales = "free_y",
      space = "free_y",
      labeller = as_labeller(function(x) {
        stringr::str_replace(x, "^Extended: ", "Extended panel\n")
      })
    )+
    # facet_grid(panel_label ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(
      low = "#1D4E89",
      mid = "white",
      high = "#B22222",
      midpoint = 0,
      guide = guide_colorbar(
        barwidth = unit(6, "cm"),
        barheight = unit(0.45, "cm"),
        title.position = "left"
        # ticks.colour = "black",
        # frame.colour = "black"
      )
    ) +
    # scale_alpha_manual(values = c("FDR<0.05 in >=1 contrast" = 1, "Not FDR-significant" = 0.55)) +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) +
    labs(
      title = "Group-mean z-scores for all curated metabolites",
      subtitle = "Canonicalized metabolite labels remove synonym-level double counting",
      x = NULL,
      y = NULL,
      fill = "z-score"#,
      # alpha = NULL
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size =10,angle = 30, hjust = 1),
      axis.text.y = element_text(size =10),
      strip.text.y = element_text(angle = 0,size =10),
      title = element_text(size =12),
      plot.subtitle = element_text(size =10),
      # strip.background = element_blank(),
      legend.position = "bottom"
    )
}

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
    detect_threshold = 0) {
  subset_meta <- sample_meta[subset_mask, , drop = FALSE]
  subset_ids <- subset_meta$Sample_ID

  detect_mat <- feature_mat[, subset_ids, drop = FALSE] > detect_threshold
  overall_prevalence <- rowMeans(detect_mat, na.rm = TRUE)

  group_levels_local <- unique(as.character(subset_meta[[group_var]]))
  group_prev_df <- map_dfr(group_levels_local, function(group_name) {
    group_ids <- subset_meta$Sample_ID[as.character(subset_meta[[group_var]]) == group_name]

    tibble(
      feature_id = rownames(feature_mat),
      group_name = group_name,
      group_prevalence = rowMeans(feature_mat[, group_ids, drop = FALSE] > detect_threshold, na.rm = TRUE)
    )
  })

  group_summary <- group_prev_df %>%
    group_by(feature_id) %>%
    summarise(max_group_prevalence = max(group_prevalence, na.rm = TRUE), .groups = "drop")

  tibble(
    feature_id = rownames(feature_mat),
    overall_prevalence = overall_prevalence
  ) %>%
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

fit_single_coefficient_limma <- function(
    expr_mat,
    sample_meta,
    design_formula,
    coefficient_expression,
    contrast_name,
    model_label,
    test_type = c("ebayes", "treat"),
    lfc_threshold = 0) {
  fit_limma_contrasts(
    expr_mat = expr_mat,
    sample_meta = sample_meta,
    design_formula = design_formula,
    contrast_definitions = stats::setNames(coefficient_expression, contrast_name),
    model_label = model_label,
    test_type = match.arg(test_type),
    lfc_threshold = lfc_threshold
  )
}

extract_anchor_matrix <- function(feature_mat, anchor_ids, id_type = c("pathway", "ec")) {
  id_type <- match.arg(id_type)

  feature_key <- if (id_type == "pathway") {
    sub(":.*", "", rownames(feature_mat))
  } else {
    rownames(feature_mat)
  }

  keep_mask <- feature_key %in% anchor_ids
  anchor_mat <- feature_mat[keep_mask, , drop = FALSE]
  rownames(anchor_mat) <- feature_key[keep_mask]
  anchor_mat[!duplicated(rownames(anchor_mat)), , drop = FALSE]
}

extract_crossomic_candidates <- function(
    panel_assoc_results,
    disease_results,
    interaction_results,
    sensitivity_results,
    lo_sensitivity_results = NULL,
    lo_disease_results = NULL,
    lo_panel_assoc_results = NULL,
    feature_annotation) {
  panel_assoc_results %>%
    transmute(
      feature_id = metabolite_identification,
      panel = target_panel,
      panel_beta = logFC,
      panel_p = P.Value,
      panel_fdr = adj.P.Val
    ) %>%
    left_join(
      disease_results %>%
        transmute(
          feature_id = metabolite_identification,
          disease_beta = logFC,
          disease_p = P.Value,
          disease_fdr = adj.P.Val
        ),
      by = "feature_id"
    ) %>%
    left_join(
      interaction_results %>%
        transmute(
          feature_id = metabolite_identification,
          interaction_beta = logFC,
          interaction_p = P.Value,
          interaction_fdr = adj.P.Val
        ),
      by = "feature_id"
    ) %>%
    left_join(
      sensitivity_results %>%
        transmute(
          feature_id = metabolite_identification,
          panel = target_panel,
          sensitivity_beta = logFC,
          sensitivity_p = P.Value,
          sensitivity_fdr = adj.P.Val
        ),
      by = c("feature_id", "panel")
    ) %>%
    left_join(
      if (is.null(lo_sensitivity_results)) {
        tibble(
          feature_id = character(),
          panel = character(),
          lo_sensitivity_beta = numeric(),
          lo_sensitivity_p = numeric(),
          lo_sensitivity_fdr = numeric()
        )
      } else {
        lo_sensitivity_results %>%
          transmute(
            feature_id = metabolite_identification,
            panel = target_panel,
            lo_sensitivity_beta = logFC,
            lo_sensitivity_p = P.Value,
            lo_sensitivity_fdr = adj.P.Val
          )
      },
      by = c("feature_id", "panel")
    ) %>%
    left_join(
      if (is.null(lo_disease_results)) {
        tibble(
          feature_id = character(),
          lo_disease_beta = numeric(),
          lo_disease_p = numeric(),
          lo_disease_fdr = numeric()
        )
      } else {
        lo_disease_results %>%
          transmute(
            feature_id = metabolite_identification,
            lo_disease_beta = logFC,
            lo_disease_p = P.Value,
            lo_disease_fdr = adj.P.Val
          )
      },
      by = "feature_id"
    ) %>%
    left_join(
      if (is.null(lo_panel_assoc_results)) {
        tibble(
          feature_id = character(),
          panel = character(),
          lo_panel_beta = numeric(),
          lo_panel_p = numeric(),
          lo_panel_fdr = numeric()
        )
      } else {
        lo_panel_assoc_results %>%
          transmute(
            feature_id = metabolite_identification,
            panel = target_panel,
            lo_panel_beta = logFC,
            lo_panel_p = P.Value,
            lo_panel_fdr = adj.P.Val
          )
      },
      by = c("feature_id", "panel")
    ) %>%
    left_join(feature_annotation, by = "feature_id") %>%
    mutate(
      sign_panel_disease = sign(panel_beta) == sign(disease_beta),
      sign_panel_sensitivity = sign(panel_beta) == sign(sensitivity_beta),
      lo_sign_panel_disease = sign(lo_panel_beta) == sign(lo_disease_beta),
      lo_sign_panel_sensitivity = sign(lo_panel_beta) == sign(lo_sensitivity_beta),
      candidate_tier = case_when(
        panel_fdr < species_candidate_fdr &
          disease_fdr < species_candidate_fdr &
          sign_panel_disease &
          sensitivity_p < species_sensitivity_p_cutoff &
          sign_panel_sensitivity ~ "tier1_panel_disease_crc_support",
        panel_fdr < species_candidate_fdr &
          disease_fdr < species_candidate_fdr &
          sign_panel_disease ~ "tier2_panel_disease",
        panel_fdr < species_candidate_fdr ~ "tier3_panel_only",
        TRUE ~ "not_shortlisted"
      ),
      lo_candidate_tier = case_when(
        lo_panel_fdr < species_candidate_fdr &
          lo_disease_fdr < species_candidate_fdr &
          lo_sign_panel_disease &
          lo_sensitivity_p < species_sensitivity_p_cutoff &
          lo_sign_panel_sensitivity ~ "tier1_panel_disease_crc_support",
        lo_panel_fdr < species_candidate_fdr &
          lo_disease_fdr < species_candidate_fdr &
          lo_sign_panel_disease ~ "tier2_panel_disease",
        lo_panel_fdr < species_candidate_fdr ~ "tier3_panel_only",
        TRUE ~ "not_shortlisted"
      ),
      shortlist_origin = case_when(
        candidate_tier %in% support_strip_species_tiers & lo_candidate_tier %in% support_strip_species_tiers ~ "EO_and_LO",
        candidate_tier %in% support_strip_species_tiers ~ "EO_only",
        lo_candidate_tier %in% support_strip_species_tiers ~ "LO_only",
        TRUE ~ "neither"
      )
    ) %>%
    arrange(
      panel,
      shortlist_origin,
      panel_fdr,
      disease_fdr,
      lo_panel_fdr,
      lo_disease_fdr,
      sensitivity_p,
      lo_sensitivity_p
    )
}

plot_candidate_species_heatmap <- function(candidate_tbl) {
  plot_tbl <- candidate_tbl %>%
    filter(candidate_tier != "not_shortlisted") %>%
    select(
      panel,
      target_panel_label,
      short_label,
      candidate_tier,
      disease_beta,
      panel_beta,
      interaction_beta
    ) %>%
    pivot_longer(
      cols = c(disease_beta, panel_beta, interaction_beta),
      names_to = "effect_type",
      values_to = "beta"
    ) %>%
    mutate(
      effect_type = factor(
        effect_type,
        levels = c("disease_beta", "panel_beta", "interaction_beta"),
        labels = c("EOCRC vs EOControl", "EO panel association", "Disease x age-group interaction")
      )
    )

  if (nrow(plot_tbl) == 0) {
    return(NULL)
  }

  row_order <- candidate_tbl %>%
    filter(candidate_tier != "not_shortlisted") %>%
    distinct(panel, short_label, candidate_tier, panel_fdr, disease_fdr) %>%
    arrange(panel, candidate_tier, panel_fdr, disease_fdr, short_label) %>%
    mutate(row_id = paste(panel, short_label, sep = "||"))

  plot_tbl <- plot_tbl %>%
    mutate(
      row_id = paste(panel, short_label, sep = "||"),
      row_id = factor(row_id, levels = rev(row_order$row_id))
    )

  ggplot(plot_tbl, aes(x = effect_type, y = row_id, fill = beta)) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(target_panel_label ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient2(
      low = "#1D4E89",
      mid = "white",
      high = "#B22222",
      midpoint = 0
    ) +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) +
    labs(
      title = "EO cross-omics candidate species",
      subtitle = "Rows are shortlisted by EO panel association, EO disease shift, and EOCRC sensitivity support",
      x = NULL,
      y = NULL,
      fill = "beta"
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 25, hjust = 1),
      strip.text.y = element_text(angle = 0)
    )
}

build_species_association_table <- function(
    panel_assoc_results,
    feature_annotation,
    qc_tbl) {
  assoc_tbl <- panel_assoc_results %>%
    rename(
      feature_id = metabolite_identification,
      panel = target_panel
    ) %>%
    left_join(feature_annotation, by = "feature_id") %>%
    left_join(
      qc_tbl %>%
        select(feature_id, overall_prevalence, max_group_prevalence),
      by = "feature_id"
    ) %>%
    mutate(
      target_panel_label = recode(panel, !!!panel_labels)
    )

  if ("t" %in% names(assoc_tbl)) {
    assoc_tbl <- assoc_tbl %>%
      mutate(
        se = if_else(!is.na(t) & abs(t) > 1e-8, abs(logFC / t), NA_real_),
        conf.low = logFC - 1.96 * se,
        conf.high = logFC + 1.96 * se
      )
  } else {
    assoc_tbl <- assoc_tbl %>%
      mutate(
        se = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_
      )
  }

  assoc_tbl %>%
    mutate(
      association_direction = case_when(
        adj.P.Val < fdr_cutoff & logFC > 0 ~ "positive",
        adj.P.Val < fdr_cutoff & logFC < 0 ~ "negative",
        TRUE ~ "nonsignificant"
      ),
      association_direction_fdr10 = "not_used"
    )
}

plot_top_associated_species <- function(association_tbl, fdr_threshold = fdr_cutoff, top_n = top_species_hits_per_direction) {
  plot_tbl <- association_tbl %>%
    filter(adj.P.Val < fdr_threshold) %>%
    mutate(
      direction_label = if_else(logFC > 0, "Positive association", "Negative association"),
      ci_low = if_else(is.finite(conf.low), conf.low, logFC),
      ci_high = if_else(is.finite(conf.high), conf.high, logFC)
    ) %>%
    group_by(panel, target_panel_label, direction_label) %>%
    arrange(adj.P.Val, desc(abs(logFC)), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup()

  if (nrow(plot_tbl) == 0) {
    return(NULL)
  }

  row_order <- plot_tbl %>%
    arrange(panel, direction_label, logFC, short_label) %>%
    mutate(row_id = paste(panel, direction_label, short_label, sep = "||")) %>%
    pull(row_id)

  plot_tbl <- plot_tbl %>%
    mutate(
      row_id = paste(panel, direction_label, short_label, sep = "||"),
      row_id = factor(row_id, levels = row_order)
    )

  ggplot(plot_tbl, aes(x = logFC, y = row_id, color = direction_label)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60", linewidth = 0.3) +
    geom_segment(aes(x = ci_low, xend = ci_high, yend = row_id), linewidth = 0.5, alpha = 0.8) +
    geom_point(aes(size = overall_prevalence), alpha = 0.9) +
    facet_grid(target_panel_label ~ direction_label, scales = "free_y", space = "free_y") +
    scale_color_manual(
      values = c(
        "Positive association" = "#B22222",
        "Negative association" = "#1D4E89"
      )
    ) +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|.*\\|\\|", "", x)) +
    labs(
      title = "Top EO panel-associated species",
      subtitle = "Effect sizes come from EO panel-association models adjusted for disease, age, sex, and BMI",
      x = "Adjusted beta for panel_score",
      y = NULL,
      size = "EO prevalence",
      color = NULL
    ) +
    theme(
      panel.grid.major.y = element_blank(),
      strip.text.y = element_text(angle = 0)
    )
}

build_taxonomic_cladogram_data <- function(association_tbl, fdr_threshold = fdr_cutoff) {
  sig_tbl <- association_tbl %>%
    filter(adj.P.Val < fdr_threshold) %>%
    mutate(direction_label = if_else(logFC > 0, "Positive", "Negative"))

  if (nrow(sig_tbl) == 0) {
    return(NULL)
  }

  map(unique(sig_tbl$panel), function(panel_name) {
    panel_tbl <- sig_tbl %>%
      filter(panel == panel_name) %>%
      arrange(phylum_label, family_label, genus_label, desc(logFC), short_label) %>%
      mutate(
        y = row_number(),
        x = 4,
        node_id = paste("species", feature_id, sep = "||"),
        parent_id = paste("genus", phylum_label, family_label, genus_label, sep = "||")
      )

    genus_tbl <- panel_tbl %>%
      group_by(panel, target_panel_label, phylum_label, family_label, genus_label) %>%
      summarise(
        y = mean(y),
        child_ymin = min(y),
        child_ymax = max(y),
        n_children = n(),
        .groups = "drop"
      ) %>%
      mutate(
        x = 3,
        node_id = paste("genus", phylum_label, family_label, genus_label, sep = "||"),
        parent_id = paste("family", phylum_label, family_label, sep = "||")
      )

    family_tbl <- genus_tbl %>%
      group_by(panel, target_panel_label, phylum_label, family_label) %>%
      summarise(
        y = mean(y),
        child_ymin = min(y),
        child_ymax = max(y),
        n_children = n(),
        .groups = "drop"
      ) %>%
      mutate(
        x = 2,
        node_id = paste("family", phylum_label, family_label, sep = "||"),
        parent_id = paste("phylum", phylum_label, sep = "||")
      )

    phylum_tbl <- family_tbl %>%
      group_by(panel, target_panel_label, phylum_label) %>%
      summarise(
        y = mean(y),
        child_ymin = min(y),
        child_ymax = max(y),
        n_children = n(),
        .groups = "drop"
      ) %>%
      mutate(
        x = 1,
        node_id = paste("phylum", phylum_label, sep = "||"),
        parent_id = NA_character_
      )

    vertical_segments <- bind_rows(
      phylum_tbl %>%
        transmute(panel, target_panel_label, x, xend = x, y = child_ymin, yend = child_ymax),
      family_tbl %>%
        transmute(panel, target_panel_label, x, xend = x, y = child_ymin, yend = child_ymax),
      genus_tbl %>%
        transmute(panel, target_panel_label, x, xend = x, y = child_ymin, yend = child_ymax)
    )

    horizontal_segments <- bind_rows(
      family_tbl %>%
        transmute(panel, target_panel_label, x = 1, xend = 2, y, yend = y),
      genus_tbl %>%
        transmute(panel, target_panel_label, x = 2, xend = 3, y, yend = y),
      panel_tbl %>%
        transmute(panel, target_panel_label, x = 3, xend = 4, y, yend = y)
    )

    list(
      vertical_segments = vertical_segments,
      horizontal_segments = horizontal_segments,
      species_points = panel_tbl,
      species_labels = panel_tbl %>%
        transmute(panel, target_panel_label, x = 4.10, y, label = short_label),
      internal_labels = bind_rows(
        phylum_tbl %>%
          transmute(panel, target_panel_label, x = 0.95, y, label = phylum_label, rank = "Phylum"),
        family_tbl %>%
          filter(n_children > 1) %>%
          transmute(panel, target_panel_label, x = 1.95, y, label = family_label, rank = "Family")
      )
    )
  }) %>%
    transpose() %>%
    map(bind_rows)
}

plot_taxonomic_cladogram <- function(cladogram_data, panel_name = NULL) {
  if (is.null(cladogram_data)) {
    return(NULL)
  }

  panel_names <- unique(as.character(cladogram_data$species_points$panel))
  panel_names <- panel_names[!is.na(panel_names)]
  if (length(panel_names) == 0) {
    return(NULL)
  }

  if (is.null(panel_name)) {
    panel_name <- panel_names[[1]]
  }

  species_df <- cladogram_data$species_points %>% filter(panel == panel_name)
  vertical_df <- cladogram_data$vertical_segments %>% filter(panel == panel_name)
  horizontal_df <- cladogram_data$horizontal_segments %>% filter(panel == panel_name)
  species_label_df <- cladogram_data$species_labels %>% filter(panel == panel_name)
  internal_label_df <- cladogram_data$internal_labels %>% filter(panel == panel_name)

  if (nrow(species_df) == 0) {
    return(NULL)
  }

  panel_title <- unique(species_df$target_panel_label)

  ggplot(species_df, aes(x = x, y = y)) +
    geom_segment(
      data = vertical_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 0.35,
      color = "grey55"
    ) +
    geom_segment(
      data = horizontal_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      linewidth = 0.35,
      color = "grey55"
    ) +
    geom_point(
      aes(fill = direction_label, size = overall_prevalence),
      shape = 21,
      color = "black",
      stroke = 0.25
    ) +
    geom_text(
      data = species_label_df,
      aes(x = x, y = y, label = label),
      hjust = 0,
      size = 2.8
    ) +
    geom_text(
      data = internal_label_df,
      aes(x = x, y = y, label = label),
      hjust = 1,
      size = 2.6,
      color = "grey20"
    ) +
    scale_fill_manual(values = c(Positive = "#B22222", Negative = "#1D4E89")) +
    scale_x_continuous(
      breaks = c(1, 2, 3, 4),
      labels = c("Phylum", "Family", "Genus", "Species"),
      limits = c(0.75, 4.75),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = panel_title,
      subtitle = "Species points mark significant positive or negative EO panel associations",
      x = NULL,
      y = NULL,
      fill = "Direction",
      size = "EO prevalence"
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
}

run_taxonomic_enrichment <- function(
    association_tbl,
    direction_var = "association_direction",
    significance_layer = "FDR<0.05",
    min_taxon_size = taxonomic_enrichment_min_taxon_size) {
  rank_cols <- c(family_label = "Family", genus_label = "Genus")

  map_dfr(names(rank_cols), function(rank_col) {
    map_dfr(unique(association_tbl$panel), function(panel_name) {
      panel_tbl <- association_tbl %>%
        filter(panel == panel_name, !is.na(.data[[rank_col]]), .data[[rank_col]] != "unclassified")

      if (nrow(panel_tbl) == 0) {
        return(tibble())
      }

      map_dfr(c("positive", "negative"), function(direction_name) {
        assoc_mask <- panel_tbl[[direction_var]] == direction_name
        total_assoc <- sum(assoc_mask, na.rm = TRUE)

        if (total_assoc == 0) {
          return(tibble())
        }

        taxon_counts <- panel_tbl %>%
          count(taxon_label = .data[[rank_col]], name = "taxon_total") %>%
          filter(taxon_total >= min_taxon_size)

        if (nrow(taxon_counts) == 0) {
          return(tibble())
        }

        map_dfr(taxon_counts$taxon_label, function(taxon_name) {
          a <- sum(panel_tbl[[rank_col]] == taxon_name & assoc_mask, na.rm = TRUE)
          b <- sum(panel_tbl[[rank_col]] == taxon_name & !assoc_mask, na.rm = TRUE)
          c <- sum(panel_tbl[[rank_col]] != taxon_name & assoc_mask, na.rm = TRUE)
          d <- sum(panel_tbl[[rank_col]] != taxon_name & !assoc_mask, na.rm = TRUE)

          fisher_fit <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")

          tibble(
            panel = panel_name,
            target_panel_label = unique(panel_tbl$target_panel_label),
            taxonomic_rank = rank_cols[[rank_col]],
            direction = direction_name,
            taxon_label = taxon_name,
            taxon_hits = a,
            taxon_total = a + b,
            panel_hits_total = total_assoc,
            background_total = nrow(panel_tbl),
            taxon_hit_fraction = if_else((a + b) > 0, a / (a + b), NA_real_),
            panel_hit_fraction = total_assoc / nrow(panel_tbl),
            odds_ratio = unname(fisher_fit$estimate),
            p.value = fisher_fit$p.value,
            significance_layer = significance_layer
          )
        })
      })
    })
  }) %>%
    group_by(panel, taxonomic_rank, direction, significance_layer) %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "BH")) %>%
    ungroup()
}

plot_taxonomic_enrichment <- function(enrichment_tbl, adj_threshold = fdr_cutoff) {
  plot_tbl <- enrichment_tbl %>%
    filter(adj.P.Val < adj_threshold) %>%
    group_by(panel, target_panel_label, taxonomic_rank, direction) %>%
    arrange(adj.P.Val, desc(taxon_hits), .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()

  if (nrow(plot_tbl) == 0) {
    return(NULL)
  }

  finite_or <- plot_tbl$odds_ratio[is.finite(plot_tbl$odds_ratio) & plot_tbl$odds_ratio > 0]
  max_or <- ifelse(length(finite_or) == 0, 4, max(finite_or))

  plot_tbl <- plot_tbl %>%
    mutate(
      odds_ratio_plot = case_when(
        is.infinite(odds_ratio) ~ max_or * 1.2,
        odds_ratio <= 0 ~ NA_real_,
        TRUE ~ odds_ratio
      ),
      log2_odds_ratio = log2(odds_ratio_plot),
      direction_label = if_else(direction == "positive", "Positive-associated taxa", "Negative-associated taxa"),
      row_id = paste(panel, taxonomic_rank, direction, taxon_label, sep = "||")
    )

  row_order <- plot_tbl %>%
    arrange(panel, taxonomic_rank, direction, log2_odds_ratio, taxon_label) %>%
    pull(row_id)

  plot_tbl <- plot_tbl %>%
    mutate(row_id = factor(row_id, levels = row_order))

  ggplot(plot_tbl, aes(x = log2_odds_ratio, y = row_id)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60", linewidth = 0.3) +
    geom_point(aes(size = taxon_hits, color = -log10(adj.P.Val))) +
    facet_grid(taxonomic_rank ~ direction_label + target_panel_label, scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|.*\\|\\|.*\\|\\|", "", x)) +
    scale_color_gradient(low = "#5B8E7D", high = "#B22222") +
    labs(
      title = "Taxonomic enrichment among EO panel-associated species",
      subtitle = "Fisher exact enrichment of family/genus among significant positive or negative panel-associated species",
      x = "log2 enrichment odds ratio",
      y = NULL,
      size = "Taxon hits",
      color = "-log10(FDR)"
    ) +
    theme(
      panel.grid.major.y = element_blank(),
      strip.text.y = element_text(angle = 0)
    )
}

run_panel_association_models <- function(
    expr_mat,
    sample_meta_base,
    panel_score_long,
    target_panels,
    design_formula,
    coefficient_expression,
    model_label) {
  map_dfr(target_panels, function(panel_name) {
    panel_meta <- sample_meta_base %>%
      left_join(
        panel_score_long %>%
          filter(panel == panel_name) %>%
          select(Sample_ID, panel_score),
        by = "Sample_ID"
      )

    fit_single_coefficient_limma(
      expr_mat = expr_mat[, panel_meta$Sample_ID, drop = FALSE],
      sample_meta = panel_meta,
      design_formula = design_formula,
      coefficient_expression = coefficient_expression,
      contrast_name = paste0(panel_name, "_association"),
      model_label = model_label,
      test_type = "ebayes",
      lfc_threshold = 0
    )$results %>%
      mutate(
        target_panel = panel_name,
        target_panel_label = unname(panel_labels[panel_name])
      )
  })
}

build_panel_score_annotations <- function(panel_score_long, panel_results) {
  contrast_positions <- tibble(
    contrast = contrast_display_order,
    x_start = c(1, 3, 2),
    x_end = c(2, 4, 4)
  )

  panel_ranges <- panel_score_long %>%
    group_by(panel, panel_label) %>%
    summarise(
      y_max = max(panel_score, na.rm = TRUE),
      y_min = min(panel_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(y_span = pmax(y_max - y_min, 0.35))

  panel_results %>%
    filter(contrast %in% contrast_display_order) %>%
    left_join(contrast_positions, by = "contrast") %>%
    left_join(panel_ranges, by = c("panel", "panel_label")) %>%
    mutate(
      contrast = factor(contrast, levels = contrast_display_order),
      rank_within_panel = as.numeric(contrast),
      y_position = y_max + y_span * (0.12 + 0.12 * rank_within_panel),
      x_mid = (x_start + x_end) / 2,
      sig_label = significance_label(adj.P.Val),
      beta_label = if_else(is.na(logFC), "NA", sprintf("%+.2f", logFC)),
      annot_label = paste0("\u03b2=", beta_label, "\n", sig_label)
    )
}

plot_panel_scores <- function(panel_score_long,
                              panel_annotations,
                              label_style = c("no_sig_only", "beta_sig"),
                              ncol = 4,
                              beta_offset = 0.08,
                              sig_offset  = 0.045,
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
      Group4 = factor(
        Group4,
        levels = c("LOControl", "LOCRC", "EOControl", "EOCRC")
      )
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
  ann$y_sig  <- ann$y_position - sig_offset  * ann$y_span
  
  p <- ggplot(panel_score_long, aes(x = Group4, y = panel_score, fill = Group4)) +
    geom_jitter(
      aes(color = Group4),
      width = 0.15,
      alpha = 0.8,
      size = 0.55,
      show.legend = FALSE
    ) + geom_boxplot(
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
    # facet_wrap(~ panel_label, scales = "free_y", ncol = ncol) +
    facet_wrap(
      ~ panel_label,
      scales = "free_y",
      ncol = ncol,
      labeller = as_labeller(function(x) {
        stringr::str_replace(x, "^Extended: ", "Extended panel\n")
      })
    )+
    scale_fill_manual(
      values = group_palette,
      breaks = c("LOControl", "LOCRC", "EOControl", "EOCRC"),
      drop = FALSE
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
    labs(
      x = NULL,
      y = "Mean feature z-score"
    ) +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(hjust = 1,angle = 25, size = text_size),
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
        aes(x = x_start, xend = x_end, y = y_position, yend = y_position),
        inherit.aes = FALSE,
        linewidth = 0.45
      ) +
      geom_segment(
        data = ann,
        aes(x = x_start, xend = x_start, y = y_tick, yend = y_position),
        inherit.aes = FALSE,
        linewidth = 0.45
      ) +
      geom_segment(
        data = ann,
        aes(x = x_end, xend = x_end, y = y_tick, yend = y_position),
        inherit.aes = FALSE,
        linewidth = 0.45
      ) +
      geom_text(
        data = ann,
        aes(x = x_mid, y = y_beta, label = beta_label),
        inherit.aes = FALSE,
        size = 4,
        lineheight = 0.92
      ) +
      geom_text(
        data = ann,
        aes(x = x_mid, y = y_sig, label = sig_label_plot),
        inherit.aes = FALSE,
        size = 4,
        lineheight = 0.92
      )
  }
  
  p
}

plot_focus_boxplots <- function(long_df, feature_order) {
  if (length(feature_order) == 0) {
    return(NULL)
  }

  plot_df <- long_df %>%
    filter(metabolite_identification %in% feature_order) %>%
    mutate(
      Group4 = factor(Group4, levels = group_levels),
      metabolite_identification = factor(metabolite_identification, levels = feature_order)
    )

  ggplot(plot_df, aes(x = Group4, y = log2_intensity, fill = Group4)) +
    geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.85) +
    geom_jitter(width = 0.15, alpha = 0.3, size = 0.6, color = "black") +
    facet_wrap(~ metabolite_identification, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = group_palette) +
    labs(
      title = "Boxplots for selected nitrogen / putrefaction-related metabolites",
      x = NULL,
      y = "log2 intensity"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

extract_lm_term <- function(fit, term_name, contrast_name) {
  if (!(term_name %in% names(coef(fit)))) {
    stop("Requested term was not found in the fitted model: ", term_name)
  }

  estimate <- unname(coef(fit)[term_name])
  std_error <- sqrt(diag(vcov(fit)))[term_name]
  df_resid <- df.residual(fit)
  statistic <- estimate / std_error
  p_value <- 2 * pt(-abs(statistic), df = df_resid)
  crit <- qt(0.975, df = df_resid)

  tibble(
    contrast = contrast_name,
    estimate = estimate,
    std.error = std_error,
    statistic = statistic,
    p.value = p_value,
    conf.low = estimate - crit * std_error,
    conf.high = estimate + crit * std_error
  )
}

compare_independent_lm_terms <- function(fit_a, fit_b, term_name, contrast_name) {
  if (!(term_name %in% names(coef(fit_a))) || !(term_name %in% names(coef(fit_b)))) {
    stop("Requested term was not found in one of the stratified fitted models: ", term_name)
  }

  estimate_a <- unname(coef(fit_a)[term_name])
  estimate_b <- unname(coef(fit_b)[term_name])
  se_a <- sqrt(diag(vcov(fit_a)))[term_name]
  se_b <- sqrt(diag(vcov(fit_b)))[term_name]
  var_a <- se_a^2
  var_b <- se_b^2
  std_error <- sqrt(var_a + var_b)

  df_a <- df.residual(fit_a)
  df_b <- df.residual(fit_b)
  df_resid <- (var_a + var_b)^2 / ((var_a^2 / df_a) + (var_b^2 / df_b))
  statistic <- (estimate_a - estimate_b) / std_error
  p_value <- 2 * pt(-abs(statistic), df = df_resid)
  crit <- qt(0.975, df = df_resid)

  tibble(
    contrast = contrast_name,
    estimate = estimate_a - estimate_b,
    std.error = std_error,
    statistic = statistic,
    p.value = p_value,
    conf.low = (estimate_a - estimate_b) - crit * std_error,
    conf.high = (estimate_a - estimate_b) + crit * std_error
  )
}

build_main_figure_panel_effects <- function(panel_score_long, panel_results_display) {
  contrast_lookup <- tibble(
    contrast = main_figure_forest_contrast_order,
    contrast_label = unname(main_figure_forest_contrast_labels[main_figure_forest_contrast_order]),
    contrast_order = seq_along(main_figure_forest_contrast_order)
  )

  map_dfr(main_figure_strict_panels, function(panel_name) {
    panel_df <- panel_score_long %>%
      filter(panel == panel_name)

    fit_lo <- lm(
      panel_score ~ Disease + Age + Sex + BMI,
      data = panel_df %>% filter(Age_group == "LO")
    )
    fit_eo <- lm(
      panel_score ~ Disease + Age + Sex + BMI,
      data = panel_df %>% filter(Age_group == "EO")
    )

    bind_rows(
      extract_lm_term(
        fit_lo,
        "DiseaseCRC",
        "LOCRC_vs_LOControl_ageadj"
      ),
      extract_lm_term(
        fit_eo,
        "DiseaseCRC",
        "EOCRC_vs_EOControl_ageadj"
      ),
      compare_independent_lm_terms(
        fit_eo,
        fit_lo,
        "DiseaseCRC",
        "EO_vs_LO_delta_beta_ageadj"
      )
    ) %>%
      mutate(
        panel = panel_name,
        panel_short = unname(main_figure_panel_short_labels[panel_name]),
        panel_label = unname(panel_labels[panel_name])
      )
  }) %>%
    mutate(
      adj.p.value = p.adjust(p.value, method = "BH")
    ) %>%
    left_join(contrast_lookup, by = "contrast") %>%
    left_join(
      panel_results_display %>%
        select(panel, contrast, limma_logFC = logFC, limma_P.Value = P.Value, limma_adj.P.Val = adj.P.Val),
      by = c("panel", "contrast")
    ) %>%
    mutate(
      contrast_label = factor(contrast_label, levels = rev(contrast_lookup$contrast_label)),
      panel_short = factor(
        panel_short,
        levels = unname(main_figure_panel_short_labels[main_figure_strict_panels])
      ),
      fdr_label = case_when(
        adj.p.value < 0.001 ~ "FDR < 0.001",
        adj.p.value < 0.01 ~ paste0("FDR = ", sprintf("%.3f", adj.p.value)),
        adj.p.value < 0.1 ~ paste0("FDR = ", sprintf("%.3f", adj.p.value)),
        TRUE ~ paste0("FDR = ", sprintf("%.2f", adj.p.value))
      )
    )
}

build_main_figure_canonical_effects <- function(canonical_targeted_pooled) {
  contrast_lookup <- tibble(
    contrast = main_figure_contrast_order,
    contrast_label = unname(main_figure_contrast_labels[main_figure_contrast_order])
  )

  effect_df <- canonical_targeted_pooled %>%
    filter(
      panel %in% main_figure_strict_panels,
      contrast %in% main_figure_contrast_order
    ) %>%
    inner_join(main_figure_metabolite_order, by = c("panel", "canonical_metabolite")) %>%
    left_join(contrast_lookup, by = "contrast") %>%
    mutate(
      panel_short = unname(main_figure_panel_short_labels[panel]),
      effect_label = paste0(
        sprintf("%.2f", mean_logFC),
        case_when(
          best_adj.P.Val < fdr_cutoff ~ "*",
          TRUE ~ ""
        )
      )
    )

  max_abs <- max(abs(effect_df$mean_logFC), na.rm = TRUE)
  effect_df %>%
    mutate(
      panel_short = factor(
        panel_short,
        levels = unname(main_figure_panel_short_labels[main_figure_strict_panels])
      ),
      contrast_label = factor(contrast_label, levels = unname(main_figure_contrast_labels[main_figure_contrast_order])),
      canonical_metabolite = factor(
        canonical_metabolite,
        levels = rev(main_figure_metabolite_order$canonical_metabolite[order(main_figure_metabolite_order$metabolite_order)])
      ),
      text_color = ifelse(abs(mean_logFC) > max_abs * 0.62, "white", "black"),
      max_abs = max_abs
    )
}

run_main_figure_leave_one_out <- function(feature_mat, panel_membership, eo_meta) {
  map_dfr(main_figure_strict_panels, function(panel_name) {
    ordered_metabolites <- main_figure_metabolite_order %>%
      filter(panel == panel_name) %>%
      arrange(metabolite_order) %>%
      pull(canonical_metabolite)

    scenario_tbl <- tibble(
      omitted_metabolite = c("Full panel", ordered_metabolites),
      omission_order = seq_along(c("Full panel", ordered_metabolites))
    )

    scenario_results <- map_dfr(scenario_tbl$omitted_metabolite, function(omitted_name) {
      membership_mod <- if (omitted_name == "Full panel") {
        panel_membership
      } else {
        panel_membership %>%
          filter(!(panel == panel_name & canonical_metabolite == omitted_name))
      }

      score_mat_mod <- build_panel_score_matrix(feature_mat, membership_mod)
      fit_df <- eo_meta %>%
        mutate(panel_score = as.numeric(score_mat_mod[panel_name, Sample_ID]))

      tidy(lm(panel_score ~ Disease + Age_centered + Sex + BMI, data = fit_df), conf.int = TRUE) %>%
        filter(term == "DiseaseCRC") %>%
        transmute(
          omitted_metabolite = omitted_name,
          estimate,
          conf.low,
          conf.high,
          p.value
        )
    }) %>%
      left_join(scenario_tbl, by = "omitted_metabolite") %>%
      mutate(
        panel = panel_name,
        panel_short = unname(main_figure_panel_short_labels[panel_name]),
        omission_label = ifelse(
          omitted_metabolite == "Full panel",
          "Full panel",
          paste("Drop", omitted_metabolite)
        )
      ) %>% left_join(panel_membership %>% distinct(panel, panel_label), by = 'panel') %>%
      arrange(omission_order) %>%
      mutate(
        full_panel_estimate = estimate[omitted_metabolite == "Full panel"][1],
        estimate_shift = estimate - full_panel_estimate,
        p_adj_within_panel = p.adjust(p.value, method = "BH"),
        sig_label = case_when(
          p_adj_within_panel < 0.001 ~ "***",
          p_adj_within_panel < 0.01 ~ "**",
          p_adj_within_panel < 0.05 ~ "*",
          TRUE ~ "ns"
        ),
        panel_omission = paste(panel_name, omission_label, sep = "||")
      )

    scenario_results
  }) %>%
    mutate(
      panel_omission = factor(panel_omission, levels = rev(unique(panel_omission))),
      panel_short = factor(
        panel_short,
        levels = unname(main_figure_panel_short_labels[main_figure_strict_panels])
      )
    )
}

plot_main_figure_panel_effects <- function(panel_effect_df) {
  text_df <- panel_effect_df %>%
    mutate(text_x = conf.high + 0.06)

  x_min <- min(panel_effect_df$conf.low, na.rm = TRUE) - 0.01
  x_max <- max(text_df$text_x, na.rm = TRUE) + 0.01

  ggplot(panel_effect_df, aes(x = estimate, y = contrast_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", 
               color = "grey55", linewidth = 0.6) +
    geom_errorbar(
      aes(xmin = conf.low, xmax = conf.high),
      orientation = "y",
      width = 0.18,
      linewidth = 0.6,
      color = "#4A4A4A"
    ) +
    geom_point(size = 2.8, color = "#B22222") +
    # geom_text(
    #   data = text_df,
    #   aes(x = text_x, y = contrast_label, label = fdr_label),
    #   inherit.aes = FALSE,
    #   hjust = 0,
    #   size = 5
    # ) +
    facet_wrap(~ panel_label, ncol = 1, scales = "free_y") +
    labs(
      # title = "A. Module effect sizes",
      x = "Adjusted panel-score difference",
      y = NULL
    ) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      axis.text.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      plot.margin = margin(5.5, 10, 5.5, 5.5),
      strip.text = element_text(size = 16),
      axis.title.x =  element_text(size = 16),
      
    ) +
    coord_cartesian(xlim = c(x_min, x_max), clip = "off")
}

plot_main_figure_canonical_effects <- function(canonical_effect_df) {
  fill_limit <- max(canonical_effect_df$max_abs, na.rm = TRUE)

  ggplot(canonical_effect_df, aes(x = contrast_label, y = canonical_metabolite, fill = mean_logFC)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(
      aes(label = effect_label, color = text_color),
      size = 6,
      show.legend = FALSE
    ) +
    facet_wrap(~ panel_label, ncol = 1, scales = "free_y") +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-fill_limit, fill_limit)
    ) +
    scale_color_identity() +
    labs(
      # title = "B. Canonical metabolite effect map",
      # subtitle = "* FDR < 0.05",
      x = NULL,
      y = NULL,
      fill = "Mean log2 effect"
    ) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1,size = 16),
      axis.text.y = element_text(size = 16),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 16),
      legend.text = element_text(size = 12)
    )
}

plot_main_figure_leave_one_out <- function(loo_df) {
  full_ref <- loo_df %>%
    filter(omitted_metabolite == "Full panel") %>%
    select(panel_label, full_panel_estimate = estimate)

  ggplot(loo_df, aes(x = estimate, y = panel_omission)) +
    geom_vline(
      data = full_ref,
      aes(xintercept = full_panel_estimate),
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
    geom_point(
      aes(shape = omitted_metabolite == "Full panel", 
          fill = p_adj_within_panel < fdr_cutoff),
      size = 2.8,
      color = "#B22222"
    ) +
    facet_wrap(~ panel_label, ncol = 1, scales = "free_y") +
    scale_shape_manual(values = c(`TRUE` = 23, `FALSE` = 21), guide = "none") +
    scale_fill_manual(
      values = c(`TRUE` = "#B22222", `FALSE` = "white"),
      guide = 'none'
      # labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "FDR >= 0.05"),
      # name = NULL
    ) +
    scale_y_discrete(labels = function(x) sub("^.*\\|\\|", "", x)) +
    labs(
      # title = "C. Leave-one-metabolite-out robustness",
      # subtitle = "Dashed line = full-panel EO effect",
      x = "EOCRC - EOControl adjusted effect",
      y = NULL
    ) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom",
      axis.text.y = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      strip.text = element_text(size = 16),
      axis.title.x = element_text(size = 16)
    )
}

build_species_support_strip_data <- function(species_candidates) {
  selected_species <- species_candidates %>%
    filter(
      panel == support_strip_target_panel,
      (
        (
          panel_fdr < support_strip_fdr_cutoff &
            disease_fdr < support_strip_fdr_cutoff &
            as.character(sign_panel_disease) %in% c("TRUE", "same_direction", "1")
        ) |
          (
            lo_panel_fdr < support_strip_fdr_cutoff &
              lo_disease_fdr < support_strip_fdr_cutoff &
              as.character(lo_sign_panel_disease) %in% c("TRUE", "same_direction", "1")
          )
      )
    ) %>%
    mutate(
      panel_fdr = as.numeric(panel_fdr),
      disease_fdr = as.numeric(disease_fdr),
      lo_panel_fdr = as.numeric(lo_panel_fdr),
      lo_disease_fdr = as.numeric(lo_disease_fdr),
      eo_candidate_tier_support = case_when(
        panel_fdr < support_strip_fdr_cutoff &
          disease_fdr < support_strip_fdr_cutoff &
          as.character(sign_panel_disease) %in% c("TRUE", "same_direction", "1") &
          sensitivity_p < species_sensitivity_p_cutoff &
          as.character(sign_panel_sensitivity) %in% c("TRUE", "same_direction", "1") ~ "tier1_panel_disease_crc_support",
        panel_fdr < support_strip_fdr_cutoff &
          disease_fdr < support_strip_fdr_cutoff &
          as.character(sign_panel_disease) %in% c("TRUE", "same_direction", "1") ~ "tier2_panel_disease",
        panel_fdr < support_strip_fdr_cutoff ~ "tier3_panel_only",
        TRUE ~ "not_shortlisted"
      ),
      lo_candidate_tier_support = case_when(
        lo_panel_fdr < support_strip_fdr_cutoff &
          lo_disease_fdr < support_strip_fdr_cutoff &
          as.character(lo_sign_panel_disease) %in% c("TRUE", "same_direction", "1") &
          lo_sensitivity_p < species_sensitivity_p_cutoff &
          as.character(lo_sign_panel_sensitivity) %in% c("TRUE", "same_direction", "1") ~ "tier1_panel_disease_crc_support",
        lo_panel_fdr < support_strip_fdr_cutoff &
          lo_disease_fdr < support_strip_fdr_cutoff &
          as.character(lo_sign_panel_disease) %in% c("TRUE", "same_direction", "1") ~ "tier2_panel_disease",
        lo_panel_fdr < support_strip_fdr_cutoff ~ "tier3_panel_only",
        TRUE ~ "not_shortlisted"
      ),
      eo_shortlisted = eo_candidate_tier_support %in% support_strip_species_tiers,
      lo_shortlisted = lo_candidate_tier_support %in% support_strip_species_tiers,
      shortlist_origin = case_when(
        eo_shortlisted & lo_shortlisted ~ "EO_and_LO",
        eo_shortlisted ~ "EO_only",
        lo_shortlisted ~ "LO_only",
        TRUE ~ "neither"
      )
    ) %>%
    arrange(
      factor(shortlist_origin, levels = c("EO_and_LO", "EO_only", "LO_only", "neither")),
      pmin(panel_fdr, lo_panel_fdr, na.rm = TRUE),
      pmin(disease_fdr, lo_disease_fdr, na.rm = TRUE),
      desc(pmax(abs(as.numeric(panel_beta)), abs(as.numeric(lo_panel_beta)), na.rm = TRUE))
    ) %>%
    distinct(species_label, .keep_all = TRUE)

  support_df <- selected_species %>%
    transmute(
      item_label = species_label,
      short_label = short_label,
      eo_candidate_tier = eo_candidate_tier_support,
      lo_candidate_tier = lo_candidate_tier_support,
      shortlist_origin,
      panel = panel,
      eo_prevalence = as.numeric(NA),
      contrast_key = "EO_disease_effect",
      effect = as.numeric(disease_beta),
      fdr = as.numeric(disease_fdr)
    ) %>%
    bind_rows(
      selected_species %>%
        transmute(
          item_label = species_label,
          short_label = short_label,
          eo_candidate_tier = eo_candidate_tier_support,
          lo_candidate_tier = lo_candidate_tier_support,
          shortlist_origin,
          panel = panel,
          eo_prevalence = as.numeric(NA),
          contrast_key = "LO_disease_effect",
          effect = as.numeric(lo_disease_beta),
          fdr = as.numeric(lo_disease_fdr)
        )
    ) %>%
    bind_rows(
      selected_species %>%
        transmute(
          item_label = species_label,
          short_label = short_label,
          eo_candidate_tier = eo_candidate_tier_support,
          lo_candidate_tier = lo_candidate_tier_support,
          shortlist_origin,
          panel = panel,
          eo_prevalence = as.numeric(NA),
          contrast_key = "EO_putrefactive_panel_association",
          effect = as.numeric(panel_beta),
          fdr = as.numeric(panel_fdr)
        )
    ) %>%
    bind_rows(
      selected_species %>%
        transmute(
          item_label = species_label,
          short_label = short_label,
          eo_candidate_tier = eo_candidate_tier_support,
          lo_candidate_tier = lo_candidate_tier_support,
          shortlist_origin,
          panel = panel,
          eo_prevalence = as.numeric(NA),
          contrast_key = "LO_putrefactive_panel_association",
          effect = as.numeric(lo_panel_beta),
          fdr = as.numeric(lo_panel_fdr)
        )
    ) %>%
    mutate(
      contrast_label = recode(contrast_key, !!!support_strip_contrast_labels),
      sig_marker = case_when(
        is.na(fdr) ~ "",
        fdr < support_strip_fdr_cutoff ~ "*",
        TRUE ~ ""
      ),
      effect_label = if_else(is.na(effect), "NA", paste0(sprintf("%+.2f", effect), sig_marker))
    )

  max_abs <- max(abs(support_df$effect), na.rm = TRUE)
  order_levels <- selected_species %>%
    transmute(
      item_label = species_label,
      eo_disease_beta = as.numeric(disease_beta),
      best_panel_fdr = pmin(as.numeric(panel_fdr), as.numeric(lo_panel_fdr), na.rm = TRUE),
      best_disease_fdr = pmin(as.numeric(disease_fdr), as.numeric(lo_disease_fdr), na.rm = TRUE),
      shortlist_origin = factor(shortlist_origin, levels = c("EO_and_LO", "EO_only", "LO_only", "neither"))
    ) %>%
    arrange(eo_disease_beta, best_disease_fdr, best_panel_fdr, shortlist_origin, item_label) %>%
    pull(item_label)

  support_df %>%
    mutate(
      item_label = factor(item_label, levels = rev(order_levels)),
      contrast_label = factor(contrast_label, levels = unname(support_strip_contrast_labels)),
      text_color = ifelse(abs(effect) > max_abs * 0.58, "white", "black"),
      max_abs = max_abs,
      selection_rule = "Union of EO and LO putrefaction-core shortlist species using support-strip significance threshold FDR < 0.05; EO and LO disease/panel columns shown side by side",
      shortlist_origin = factor(shortlist_origin, levels = c("EO_and_LO", "EO_only", "LO_only", "neither"))
    ) %>%
    arrange(desc(item_label), contrast_label)
}

plot_species_support_strip <- function(support_df) {
  fill_limit <- max(support_df$max_abs, na.rm = TRUE)

  ggplot(support_df, aes(x = contrast_label, y = item_label, fill = effect)) +
    geom_tile(color = "white", linewidth = 0.35) +
    geom_text(
      aes(label = effect_label, color = text_color),
      size = 6,
      show.legend = FALSE
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-fill_limit, fill_limit)
    ) +
    scale_color_identity() +
    labs(
      # title = "D. Support strip",
      # subtitle = "EO shortlist selected by EO disease + EO spillover alignment; LO disease and LO spillover columns shown for age-stratified comparison; * FDR < 0.05",
      x = NULL,
      y = NULL,
      fill = "Effect"
    ) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, size = 16),
      axis.text.y = element_text(size = 16),
      panel.grid = element_blank(),
      legend.text = element_text(size = 16)
    )
}

A30_metagenomic_id <- read_xlsx(meta_path, sheet = 2)
full_meta <- read_xlsx(meta_path, sheet = 1)

sample_id_col <- names(full_meta)[10]
batch_col <- names(full_meta)[11]
stage_col <- names(full_meta)[20]
tumor_site_col <- names(full_meta)[21]

FUSCC_cl2 <- full_meta %>%
  filter(ID %in% A30_metagenomic_id$`*ID`) %>%
  transmute(
    Patient_ID = ID,
    Sex = factor(gender),
    Age = as.numeric(age),
    BMI = as.numeric(BMI),
    Group = Group,
    Sample_ID = .data[[sample_id_col]],
    Batch = factor(.data[[batch_col]]),
    Stage_raw = .data[[stage_col]],
    Tumor_Site_raw = .data[[tumor_site_col]],
    Cohort = "A30_FUSCC",
    Country = "China",
    Diagnosis = Group
  ) %>%
  mutate(
    Group = ifelse(Group == 'CTRL',"Control",Group ),
    Disease = factor(Group, levels = c("Control", "CRC")),
    Age_group = factor(ifelse(Age < 50, "EO", "LO"), levels = c("LO", "EO")),
    Group4 = factor(paste0(Age_group, Disease), levels = group_levels),
    Age_within_group = Age - ave(Age, Age_group, FUN = function(x) mean(x, na.rm = TRUE)),
    Tumor_Site3 = factor(case_when(
      Tumor_Site_raw == 0 ~ "right_colon",
      Tumor_Site_raw == 1 ~ "left_colon",
      Tumor_Site_raw == 2 ~ "rectum",
      TRUE ~ "missing"
    ), levels = c("right_colon", "left_colon", "rectum", "missing")),
    Tumor_Site2 = factor(case_when(
      Tumor_Site3 %in% c("right_colon", "left_colon") ~ "colon",
      Tumor_Site3 == "rectum" ~ "rectum",
      TRUE ~ "missing"
    ), levels = c("colon", "rectum", "missing")),
    Stage_chr = str_trim(as.character(Stage_raw)),
    Stage_num = suppressWarnings(readr::parse_number(Stage_chr)),
    Stage_collapsed = factor(case_when(
      !is.na(Stage_num) & Stage_num <= 2 ~ "Stage0_2",
      !is.na(Stage_num) & Stage_num >= 3 ~ "Stage3_4",
      TRUE ~ "Missing_or_pCR"
    ), levels = c("Stage0_2", "Stage3_4", "Missing_or_pCR"))
  )

stopifnot(
  nrow(FUSCC_cl2) > 0,
  all(!is.na(FUSCC_cl2$Sample_ID)),
  !any(duplicated(FUSCC_cl2$Sample_ID))
)

write_output_csv(
  FUSCC_cl2 %>%
    group_by(Group4) %>%
    summarise(
      n = n(),
      age_mean = mean(Age, na.rm = TRUE),
      age_sd = sd(Age, na.rm = TRUE),
      bmi_mean = mean(BMI, na.rm = TRUE),
      bmi_sd = sd(BMI, na.rm = TRUE),
      female_n = sum(Sex == "female", na.rm = TRUE),
      male_n = sum(Sex == "male", na.rm = TRUE),
      .groups = "drop"
    ),
  file.path(table_dir, "clinical_group_summary_v28.csv"),
  row.names = FALSE
)

write_output_csv(
  FUSCC_cl2 %>%
    filter(Disease == "CRC") %>%
    count(Age_group, Stage_collapsed) %>%
    tidyr::pivot_wider(names_from = Stage_collapsed, values_from = n, values_fill = 0),
  file.path(table_dir, "crc_stage_by_age_group_v28.csv"),
  row.names = FALSE
)

write_output_csv(
  FUSCC_cl2 %>%
    filter(Disease == "CRC") %>%
    count(Age_group, Tumor_Site3) %>%
    tidyr::pivot_wider(names_from = Tumor_Site3, values_from = n, values_fill = 0),
  file.path(table_dir, "crc_tumor_site_by_age_group_v28.csv"),
  row.names = FALSE
)

LC_MS_negative <- read_mode_file(file.path(analysis_dir, "OMIX006518-09.csv"), "negative")
LC_MS_positive <- read_mode_file(file.path(analysis_dir, "OMIX006518-08.csv"), "positive")
mode_dedup_scored <- score_mode_candidates(bind_rows(LC_MS_negative, LC_MS_positive))
LC_MS_selected <- deduplicate_modes(mode_dedup_scored)

write_output_csv(
  mode_dedup_scored %>% arrange(metabolite_identification, selection_rank),
  file.path(table_dir, "targeted_mode_dedup_decisions.csv"),
  row.names = FALSE
)

sample_cols <- grep("^[oy]", names(LC_MS_selected), value = TRUE)

LC_long <- LC_MS_selected %>%
  select(all_of(c("metabolite_identification", sample_cols))) %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "Sample_ID",
    values_to = "intensity"
  ) %>%
  left_join(FUSCC_cl2, by = "Sample_ID")

metab_qc <- LC_long %>%
  group_by(metabolite_identification) %>%
  summarise(
    overall_missing_rate = mean(is.na(intensity)),
    overall_detect_rate = mean(!is.na(intensity) & intensity > 0),
    .groups = "drop"
  ) %>%
  left_join(
    LC_long %>%
      group_by(metabolite_identification, Group4) %>%
      summarise(group_detect_rate = mean(!is.na(intensity) & intensity > 0), .groups = "drop") %>%
      group_by(metabolite_identification) %>%
      summarise(max_group_detect_rate = max(group_detect_rate), .groups = "drop"),
    by = "metabolite_identification"
  ) %>%
  mutate(
    keep_feature = overall_missing_rate <= max_missing_fraction &
      max_group_detect_rate >= min_present_fraction
  )

write_output_csv(metab_qc, file.path(table_dir, "metabolite_qc_summary.csv"), row.names = FALSE)

keep_metabs <- metab_qc %>%
  filter(keep_feature) %>%
  pull(metabolite_identification)

LC_long2 <- LC_long %>%
  filter(metabolite_identification %in% keep_metabs) %>%
  group_by(metabolite_identification) %>%
  mutate(
    min_pos = suppressWarnings(min(intensity[intensity > 0], na.rm = TRUE)),
    min_pos = ifelse(is.infinite(min_pos), NA_real_, min_pos),
    intensity_imp = case_when(
      !is.na(intensity) & intensity > 0 ~ intensity,
      is.na(min_pos) ~ NA_real_,
      TRUE ~ min_pos / 2
    ),
    log2_intensity = log2(intensity_imp)
  ) %>%
  ungroup()

feature_mat <- LC_long2 %>%
  select(metabolite_identification, Sample_ID, log2_intensity) %>%
  distinct() %>%
  pivot_wider(names_from = Sample_ID, values_from = log2_intensity) %>%
  as.data.frame()

rownames(feature_mat) <- feature_mat$metabolite_identification
feature_mat$metabolite_identification <- NULL
feature_mat <- as.matrix(feature_mat)
feature_mat <- feature_mat[, FUSCC_cl2$Sample_ID, drop = FALSE]

sample_mat_log2 <- t(feature_mat)
sample_mat_pareto <- pareto_scale_matrix(sample_mat_log2)

pca_obj <- prcomp(sample_mat_pareto, center = FALSE, scale. = FALSE)
pca_var <- round(100 * summary(pca_obj)$importance[2, 1:2], 1)

pca_df <- as.data.frame(pca_obj$x[, 1:5]) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(FUSCC_cl2, by = "Sample_ID")

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group4, shape = Disease)) +
  geom_point(size = 2.0, alpha = 0.75) +
  stat_ellipse(aes(group = Group4), linewidth = 0.6, alpha = 0.5) +
  scale_color_manual(values = group_palette) +
  labs(
    title = "Pareto-scaled PCA",
    x = paste0("PC1 (", pca_var[1], "%)"),
    y = paste0("PC2 (", pca_var[2], "%)")
  ) +
  theme(panel.grid.minor = element_blank())

save_gg(pca_plot, file.path(plot_dir, "pca_group4.pdf"), width = 8.5, height = 6.5)

adonis_group4_df <- as.data.frame(
  vegan::adonis2(
    sample_mat_pareto ~ Group4 + Age_within_group + Sex + BMI,
    data = FUSCC_cl2,
    method = "euclidean",
    by = "margin",
    permutations = 999
  )
) %>%
  rownames_to_column("term")

adonis_interaction_df <- as.data.frame(
  vegan::adonis2(
    sample_mat_pareto ~ Disease * Age_group + Age_within_group + Sex + BMI,
    data = FUSCC_cl2,
    method = "euclidean",
    by = "terms",
    permutations = 999
  )
) %>%
  rownames_to_column("term")

write_output_csv(adonis_group4_df, file.path(table_dir, "adonis2_group4_overall.csv"), row.names = FALSE)
write_output_csv(
  adonis_interaction_df,
  file.path(table_dir, "adonis2_disease_age_interaction_terms.csv"),
  row.names = FALSE
)

dispersion_obj <- vegan::betadisper(dist(sample_mat_pareto), group = FUSCC_cl2$Group4)
dispersion_anova <- as.data.frame(anova(dispersion_obj)) %>%
  rownames_to_column("term")
dispersion_perm <- as.data.frame(run_betadisper_permutest(dispersion_obj, permutations = 999)$tab) %>%
  rownames_to_column("term")

write_output_csv(dispersion_anova, file.path(table_dir, "betadisper_anova.csv"), row.names = FALSE)
write_output_csv(dispersion_perm, file.path(table_dir, "betadisper_permutest.csv"), row.names = FALSE)

dispersion_df <- FUSCC_cl2 %>%
  mutate(distance_to_centroid = dispersion_obj$distances)

dispersion_plot <- ggplot(dispersion_df, aes(x = Group4, y = distance_to_centroid, fill = Group4)) +
  geom_boxplot(outlier.shape = NA, width = 0.65, alpha = 0.85) +
  geom_jitter(width = 0.15, alpha = 0.35, size = 0.8, color = "black") +
  scale_fill_manual(values = group_palette) +
  labs(
    title = "Multivariate dispersion by four-group status",
    x = NULL,
    y = "Distance to centroid"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank()
  )

save_gg(dispersion_plot, file.path(plot_dir, "group_dispersion_boxplot.pdf"), width = 8, height = 6)

panel_membership <- build_targeted_panel_membership(rownames(feature_mat))
write_output_csv(panel_membership, file.path(table_dir, "targeted_panel_membership.csv"), row.names = FALSE)
write_output_csv(
  panel_membership %>%
    distinct(panel, panel_label, panel_tier, panel_tier_label, canonical_metabolite) %>%
    arrange(panel, canonical_metabolite),
  file.path(table_dir, "targeted_panel_catalog_v28.csv"),
  row.names = FALSE
)
write_output_csv(
  build_focus_metabolite_review(rownames(feature_mat)),
  file.path(table_dir, "panel_design_review_v28.csv"),
  row.names = FALSE
)
write_output_csv(
  mode_dedup_scored %>%
    filter(selected_mode) %>%
    semi_join(
      panel_membership %>%
        filter(panel_tier == "strict_manuscript_aligned") %>%
        distinct(metabolite_identification),
      by = "metabolite_identification"
    ) %>%
    arrange(metabolite_identification, mode) %>%
    select(
      metabolite_identification,
      mode,
      prefer_mode,
      mode_rank,
      selection_rank,
      selected_mode,
      selection_reason,
      n_detected,
      missing_rate,
      median_signal
    ),
  file.path(table_dir, "targeted_mode_dedup_strict_metabolites.csv"),
  row.names = FALSE
)

panel_score_coverage_by_sample <- build_panel_score_coverage_table(
  raw_long_df = LC_long %>% filter(metabolite_identification %in% rownames(feature_mat)),
  panel_membership = panel_membership
)
write_output_csv(
  panel_score_coverage_by_sample,
  file.path(table_dir, "panel_score_coverage_by_sample.csv"),
  row.names = FALSE
)
write_output_csv(
  summarise_panel_score_coverage(panel_score_coverage_by_sample),
  file.path(table_dir, "panel_score_coverage_summary_by_group.csv"),
  row.names = FALSE
)

eo_meta <- FUSCC_cl2 %>%
  filter(Age_group == "EO") %>%
  mutate(Age_centered = Age - mean(Age, na.rm = TRUE))

lo_meta <- FUSCC_cl2 %>%
  filter(Age_group == "LO") %>%
  mutate(Age_centered = Age - mean(Age, na.rm = TRUE))

crc_meta <- FUSCC_cl2 %>%
  filter(Disease == "CRC") %>%
  group_by(Age_group) %>%
  mutate(Age_within_crc_group = Age - mean(Age, na.rm = TRUE)) %>%
  ungroup()

write_output_csv(
  tibble(
    model_label = c(
      "EO_stratified_main",
      "LO_stratified_main",
      "pooled_interaction_ageadj",
      "CRC_only_sensitivity"
    ),
    n_samples = c(nrow(eo_meta), nrow(lo_meta), nrow(FUSCC_cl2), nrow(crc_meta))
  ),
  file.path(table_dir, "model_sample_sizes_v28.csv"),
  row.names = FALSE
)

eo_feature_fit <- fit_limma_contrasts(
  expr_mat = feature_mat[, eo_meta$Sample_ID, drop = FALSE],
  sample_meta = eo_meta,
  design_formula = ~ Disease + Age_centered + Sex + BMI,
  contrast_definitions = c(EOCRC_vs_EOControl_ageadj = "DiseaseCRC"),
  model_label = "EO_stratified_main",
  test_type = "treat",
  lfc_threshold = effect_cutoff
)

lo_feature_fit <- fit_limma_contrasts(
  expr_mat = feature_mat[, lo_meta$Sample_ID, drop = FALSE],
  sample_meta = lo_meta,
  design_formula = ~ Disease + Age_centered + Sex + BMI,
  contrast_definitions = c(LOCRC_vs_LOControl_ageadj = "DiseaseCRC"),
  model_label = "LO_stratified_main",
  test_type = "treat",
  lfc_threshold = effect_cutoff
)

pooled_feature_fit <- fit_limma_contrasts(
  expr_mat = feature_mat,
  sample_meta = FUSCC_cl2,
  design_formula = ~ Disease * Age_group + Age_within_group + Sex + BMI,
  contrast_definitions = c(
    LOCRC_vs_LOControl_pooled_ageadj = "DiseaseCRC",
    EOCRC_vs_EOControl_pooled_ageadj = "DiseaseCRC + DiseaseCRC.Age_groupEO",
    Disease_age_interaction_pooled_ageadj = "DiseaseCRC.Age_groupEO"
  ),
  model_label = "pooled_interaction_ageadj",
  test_type = "treat",
  lfc_threshold = effect_cutoff
)

crc_feature_fit <- fit_limma_contrasts(
  expr_mat = feature_mat[, crc_meta$Sample_ID, drop = FALSE],
  sample_meta = crc_meta,
  design_formula = ~ Age_group + Age_within_crc_group + Sex + BMI + Tumor_Site3 + Stage_collapsed,
  contrast_definitions = c(EOCRC_vs_LOCRC_crc_age_stage_site_adj = "Age_groupEO"),
  model_label = "CRC_only_sensitivity",
  test_type = "treat",
  lfc_threshold = effect_cutoff
)

feature_results_display <- bind_rows(
  eo_feature_fit$results,
  lo_feature_fit$results,
  crc_feature_fit$results
)

feature_results_pooled <- pooled_feature_fit$results
feature_results_all <- bind_rows(feature_results_display, feature_results_pooled)

write_output_csv(
  feature_results_display,
  file.path(table_dir, "feature_results_display_models_treat.csv"),
  row.names = FALSE
)
write_output_csv(
  feature_results_pooled,
  file.path(table_dir, "feature_results_pooled_interaction_treat.csv"),
  row.names = FALSE
)
write_output_csv(
  feature_results_all,
  file.path(table_dir, "feature_results_all_models_treat.csv"),
  row.names = FALSE
)

write_output_csv(
  feature_results_display %>% filter(significant_fdr),
  file.path(table_dir, "feature_results_display_models_treat_fdr_lt_0_05.csv"),
  row.names = FALSE
)

targeted_results_display <- feature_results_display %>%
  inner_join(
    panel_membership,
    by = "metabolite_identification",
    relationship = "many-to-many"
  ) %>%
  arrange(contrast, adj.P.Val, canonical_metabolite, metabolite_identification)

targeted_results_pooled <- feature_results_pooled %>%
  inner_join(
    panel_membership,
    by = "metabolite_identification",
    relationship = "many-to-many"
  ) %>%
  arrange(contrast, adj.P.Val, canonical_metabolite, metabolite_identification)

targeted_results_all <- bind_rows(targeted_results_display, targeted_results_pooled)

write_output_csv(
  targeted_results_display,
  file.path(table_dir, "targeted_panel_feature_results_display_models.csv"),
  row.names = FALSE
)
write_output_csv(
  targeted_results_pooled,
  file.path(table_dir, "targeted_panel_feature_results_pooled_interaction.csv"),
  row.names = FALSE
)
write_output_csv(
  targeted_results_all,
  file.path(table_dir, "targeted_panel_feature_results_all_models.csv"),
  row.names = FALSE
)

significant_targeted_display <- targeted_results_display %>%
  filter(significant_fdr)

canonical_targeted_display <- targeted_results_display %>%
  group_by(model_label, contrast, panel, panel_label, panel_tier, panel_tier_label, canonical_metabolite) %>%
  summarise(
    best_adj.P.Val = min(adj.P.Val, na.rm = TRUE),
    best_P.Value = min(P.Value, na.rm = TRUE),
    any_significant_fdr = any(significant_fdr),
    any_significant_fdr_fc = any(significant_fdr_fc),
    mean_logFC = mean(logFC, na.rm = TRUE),
    n_raw_features = n_distinct(metabolite_identification),
    raw_features = paste(sort(unique(metabolite_identification)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(direction = ifelse(mean_logFC >= 0, "higher", "lower"))

canonical_targeted_pooled <- targeted_results_pooled %>%
  group_by(model_label, contrast, panel, panel_label, panel_tier, panel_tier_label, canonical_metabolite) %>%
  summarise(
    best_adj.P.Val = min(adj.P.Val, na.rm = TRUE),
    best_P.Value = min(P.Value, na.rm = TRUE),
    any_significant_fdr = any(significant_fdr),
    any_significant_fdr_fc = any(significant_fdr_fc),
    mean_logFC = mean(logFC, na.rm = TRUE),
    n_raw_features = n_distinct(metabolite_identification),
    raw_features = paste(sort(unique(metabolite_identification)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(direction = ifelse(mean_logFC >= 0, "higher", "lower"))

significant_targeted_display_canonical <- canonical_targeted_display %>%
  filter(any_significant_fdr)

write_output_csv(
  significant_targeted_display,
  file.path(table_dir, "targeted_panel_feature_results_display_models_fdr_lt_0_05.csv"),
  row.names = FALSE
)
write_output_csv(
  canonical_targeted_display,
  file.path(table_dir, "targeted_panel_feature_results_display_models_canonical_summary.csv"),
  row.names = FALSE
)
write_output_csv(
  significant_targeted_display_canonical,
  file.path(table_dir, "targeted_panel_feature_results_display_models_canonical_fdr_lt_0_05.csv"),
  row.names = FALSE
)
write_output_csv(
  canonical_targeted_display %>% filter(any_significant_fdr),
  file.path(table_dir, "targeted_panel_feature_results_display_models_canonical_fdr_lt_0_05_only.csv"),
  row.names = FALSE
)
write_output_csv(
  canonical_targeted_pooled,
  file.path(table_dir, "targeted_panel_feature_results_pooled_interaction_canonical_summary.csv"),
  row.names = FALSE
)

panel_score_mat <- build_panel_score_matrix(feature_mat, panel_membership)

eo_panel_fit <- fit_limma_contrasts(
  expr_mat = panel_score_mat[, eo_meta$Sample_ID, drop = FALSE],
  sample_meta = eo_meta,
  design_formula = ~ Disease + Age_centered + Sex + BMI,
  contrast_definitions = c(EOCRC_vs_EOControl_ageadj = "DiseaseCRC"),
  model_label = "EO_stratified_main",
  test_type = "ebayes",
  lfc_threshold = 0
)

lo_panel_fit <- fit_limma_contrasts(
  expr_mat = panel_score_mat[, lo_meta$Sample_ID, drop = FALSE],
  sample_meta = lo_meta,
  design_formula = ~ Disease + Age_centered + Sex + BMI,
  contrast_definitions = c(LOCRC_vs_LOControl_ageadj = "DiseaseCRC"),
  model_label = "LO_stratified_main",
  test_type = "ebayes",
  lfc_threshold = 0
)

pooled_panel_fit <- fit_limma_contrasts(
  expr_mat = panel_score_mat,
  sample_meta = FUSCC_cl2,
  design_formula = ~ Disease * Age_group + Age_within_group + Sex + BMI,
  contrast_definitions = c(
    LOCRC_vs_LOControl_pooled_ageadj = "DiseaseCRC",
    EOCRC_vs_EOControl_pooled_ageadj = "DiseaseCRC + DiseaseCRC.Age_groupEO",
    Disease_age_interaction_pooled_ageadj = "DiseaseCRC.Age_groupEO"
  ),
  model_label = "pooled_interaction_ageadj",
  test_type = "ebayes",
  lfc_threshold = 0
)

crc_panel_fit <- fit_limma_contrasts(
  expr_mat = panel_score_mat[, crc_meta$Sample_ID, drop = FALSE],
  sample_meta = crc_meta,
  design_formula = ~ Age_group + Age_within_crc_group + Sex + BMI + Tumor_Site3 + Stage_collapsed,
  contrast_definitions = c(EOCRC_vs_LOCRC_crc_age_stage_site_adj = "Age_groupEO"),
  model_label = "CRC_only_sensitivity",
  test_type = "ebayes",
  lfc_threshold = 0
)

panel_results_display <- bind_rows(
  eo_panel_fit$results,
  lo_panel_fit$results,
  crc_panel_fit$results
) %>%
  rename(panel = metabolite_identification) %>%
  mutate(
    panel_label = recode(panel, !!!panel_labels),
    panel_tier = recode(panel, !!!panel_tier_map),
    panel_tier_label = recode(panel_tier, !!!panel_tier_labels)
  )

panel_results_pooled <- pooled_panel_fit$results %>%
  rename(panel = metabolite_identification) %>%
  mutate(
    panel_label = recode(panel, !!!panel_labels),
    panel_tier = recode(panel, !!!panel_tier_map),
    panel_tier_label = recode(panel_tier, !!!panel_tier_labels)
  )

panel_results_all <- bind_rows(panel_results_display, panel_results_pooled)

write_output_csv(
  panel_results_display,
  file.path(table_dir, "targeted_panel_score_results_display_models.csv"),
  row.names = FALSE
)
write_output_csv(
  panel_results_pooled,
  file.path(table_dir, "targeted_panel_score_results_pooled_interaction.csv"),
  row.names = FALSE
)
write_output_csv(
  panel_results_all,
  file.path(table_dir, "targeted_panel_score_results_all_models.csv"),
  row.names = FALSE
)

panel_score_long <- as.data.frame(t(panel_score_mat)) %>%
  rownames_to_column("Sample_ID") %>%
  pivot_longer(
    cols = -Sample_ID,
    names_to = "panel",
    values_to = "panel_score"
  ) %>%
  mutate(
    panel = factor(panel, levels = panel_levels),
    panel_label = recode(panel, !!!panel_labels),
    panel_tier = recode(panel, !!!panel_tier_map),
    panel_tier_label = recode(panel_tier, !!!panel_tier_labels)
  ) %>%
  left_join(FUSCC_cl2, by = "Sample_ID")

write_output_csv(panel_score_long, file.path(table_dir, "targeted_panel_scores_by_sample.csv"), row.names = FALSE)

panel_score_annotations <- build_panel_score_annotations(panel_score_long, panel_results_display)
write_output_csv(
  panel_score_annotations,
  file.path(table_dir, "targeted_panel_score_annotations_v28.csv"),
  row.names = FALSE
)

panel_score_plot <- plot_panel_scores(panel_score_long, 
                                      panel_score_annotations, 
                                      label_style = "beta_sig",
                                      ncol = 3,eolocr_lift = 0.15,
                                      text_size = 12)
save_gg(panel_score_plot, file.path(plot_dir, "targeted_panel_scores.pdf"), 
        width = 15, height = 13)

strict_panel_score_plot <- plot_panel_scores(
  panel_score_long %>% filter(panel_tier == "strict_manuscript_aligned"),
  panel_score_annotations %>% filter(panel_tier == "strict_manuscript_aligned"),
  label_style = "no_sig_only",beta_offset = 0.03,sig_offset = 0.04,
  ncol = 1
)
save_gg(
  strict_panel_score_plot,
  file.path(plot_dir, "targeted_panel_scores_strict.pdf"),
  width = 4.5,
  height = 6.5
)

extended_panel_score_plot <- plot_panel_scores(
  panel_score_long %>% filter(panel_tier == "extended_class_based"),
  panel_score_annotations %>% filter(panel_tier == "extended_class_based"),
  label_style = "beta_sig",ncol = 5,
  beta_offset = 0.08,sig_offset = 0.08,eolocr_lift = 0.2,
  text_size = 10
)
save_gg(
  extended_panel_score_plot,
  file.path(plot_dir, "targeted_panel_scores_extended.pdf"),
  width = 12,
  height = 7
)

volcano_plots <- map(
  contrast_display_order,
  ~ plot_volcano(
    df = feature_results_display %>% filter(contrast == .x),
    contrast_name = .x,
    targeted_membership = panel_membership
  )
)

volcano_panel <- wrap_plots(volcano_plots, ncol = 2) +
  plot_annotation(title = "Differential metabolite volcano plots for display models")

save_gg(volcano_panel, file.path(plot_dir, "volcano_plots_display_models.pdf"), width = 14, height = 10)

targeted_canonical_long <- build_targeted_canonical_long(LC_long2, panel_membership)
write_output_csv(
  targeted_canonical_long,
  file.path(table_dir, "targeted_panel_canonical_scores_by_sample.csv"),
  row.names = FALSE
)

targeted_heatmap_plot <- plot_targeted_heatmap(
  long_df = targeted_canonical_long,
  significant_targeted_panel_df = significant_targeted_display_canonical
)

if (!is.null(targeted_heatmap_plot)) {
  save_gg(
    targeted_heatmap_plot,
    file.path(plot_dir, "significant_targeted_metabolite_heatmap.pdf"),
    width = 8.5,
    height = 10
  )
}

all_panel_member_heatmap <- plot_panel_member_heatmap(
  long_df = targeted_canonical_long,
  panel_lookup = panel_membership,
  feature_results = canonical_targeted_display
)

if (!is.null(all_panel_member_heatmap)) {
  save_gg(
    all_panel_member_heatmap,
    file.path(plot_dir, "all_panel_members_heatmap.pdf"),
    width = 11,
    height = 14
  )
}

strict_targeted_heatmap_plot <- plot_targeted_heatmap(
  long_df = targeted_canonical_long,
  significant_targeted_panel_df = significant_targeted_display_canonical %>%
    filter(panel_tier == "strict_manuscript_aligned")
)

if (!is.null(strict_targeted_heatmap_plot)) {
  save_gg(
    strict_targeted_heatmap_plot,
    file.path(plot_dir, "significant_targeted_metabolite_heatmap_strict.pdf"),
    width = 8.2,
    height = 5.8
  )
}

extended_targeted_heatmap_plot <- plot_targeted_heatmap(
  long_df = targeted_canonical_long,
  significant_targeted_panel_df = significant_targeted_display_canonical %>%
    filter(panel_tier == "extended_class_based")
)

if (!is.null(extended_targeted_heatmap_plot)) {
  save_gg(
    extended_targeted_heatmap_plot,
    file.path(plot_dir, "significant_targeted_metabolite_heatmap_extended.pdf"),
    width = 8.8,
    height = 10.5
  )
}

strict_all_panel_member_heatmap <- plot_panel_member_heatmap(
  long_df = targeted_canonical_long,
  panel_lookup = panel_membership %>% filter(panel_tier == "strict_manuscript_aligned"),
  feature_results = canonical_targeted_display %>% filter(panel_tier == "strict_manuscript_aligned")
)

if (!is.null(strict_all_panel_member_heatmap)) {
  save_gg(
    strict_all_panel_member_heatmap,
    file.path(plot_dir, "all_panel_members_heatmap_strict.pdf"),
    width = 8.4,
    height = 6.2
  )
}

extended_all_panel_member_heatmap <- plot_panel_member_heatmap(
  long_df = targeted_canonical_long,
  panel_lookup = panel_membership %>% filter(panel_tier == "extended_class_based"),
  feature_results = canonical_targeted_display %>% filter(panel_tier == "extended_class_based")
)

if (!is.null(extended_all_panel_member_heatmap)) {
  save_gg(
    extended_all_panel_member_heatmap,
    file.path(plot_dir, "all_panel_members_heatmap_extended.pdf"),
    width = 9.2,
    height = 11
  )
}

predefined_focus <- c(
  "5-Aminopentanoic acid",
  "N-Acetylputrescine",
  "Citrulline",
  "D-Ornithine",
  "N-Acetylornithine",
  "Indole-3-acetaldehyde",
  "Indole-3-carboxyaldehyde",
  "Indoxyl sulfate",
  "4-Hydroxyphenylacetylglutamine",
  "Phenylacetic acid",
  "N-Acetylputrescine",
  "N-Acetylcadaverine",
  "Histamine",
  "Tyramine",
  "p-Cresol",
  "p-Cresol sulfate",
  "Phenol",
  "Phenol sulphate",
  "Indole",
  "Indole-3-acetamide",
  "Arginine"
)

representative_feature_map <- panel_membership %>%
  group_by(canonical_metabolite) %>%
  summarise(
    representative_feature = sort(unique(metabolite_identification))[1],
    .groups = "drop"
  )

top_targeted_features <- significant_targeted_display_canonical %>%
  group_by(canonical_metabolite) %>%
  summarise(best_fdr = min(best_adj.P.Val), .groups = "drop") %>%
  arrange(best_fdr) %>%
  slice_head(n = 12) %>%
  left_join(representative_feature_map, by = "canonical_metabolite") %>%
  pull(representative_feature)

focus_features <- unique(c(predefined_focus, top_targeted_features))
focus_features <- focus_features[focus_features %in% rownames(feature_mat)]
focus_features <- focus_features[seq_len(min(length(focus_features), 15))]

focus_boxplot <- plot_focus_boxplots(LC_long2, focus_features)

if (!is.null(focus_boxplot)) {
  save_gg(
    focus_boxplot,
    file.path(plot_dir, "focus_metabolite_boxplots.pdf"),
    width = 13,
    height = 13
  )
}

main_figure_panel_effects <- build_main_figure_panel_effects(
  panel_score_long = panel_score_long,
  panel_results_display = panel_results_display
)

main_figure_canonical_effects <- build_main_figure_canonical_effects(
  canonical_targeted_pooled = canonical_targeted_pooled
)

main_figure_leave_one_out <- run_main_figure_leave_one_out(
  feature_mat = feature_mat,
  panel_membership = panel_membership,
  eo_meta = eo_meta
)

write_output_csv(
  main_figure_panel_effects,
  file.path(table_dir, "main_figure_panel_effects_v28.csv"),
  row.names = FALSE
)
write_output_csv(
  main_figure_canonical_effects %>%
    select(-text_color, -max_abs),
  file.path(table_dir, "main_figure_canonical_effects_v28.csv"),
  row.names = FALSE
)
write_output_csv(
  main_figure_leave_one_out,
  file.path(table_dir, "main_figure_leave_one_out_v28.csv"),
  row.names = FALSE
)

main_panel_effect_plot <- plot_main_figure_panel_effects(main_figure_panel_effects %>% 
                                                           mutate(panel_label =  factor(panel_label,
                                                                                        levels = c('Arginine-ornithine routing axis',
                                                                                                   'Amino-acid putrefaction core')  )))
main_canonical_effect_plot <- plot_main_figure_canonical_effects(main_figure_canonical_effects)
main_leave_one_out_plot <- plot_main_figure_leave_one_out(main_figure_leave_one_out)

save_gg(
  main_panel_effect_plot,
  file.path(plot_dir, "main_figure_panel_effects_v28.pdf"),
  width = 8,
  height = 8
)
save_gg(
  main_canonical_effect_plot,
  file.path(plot_dir, "main_figure_canonical_effects_v28.pdf"),
  width = 10,
  height = 8
)
save_gg(
  main_leave_one_out_plot,
  file.path(plot_dir, "main_figure_leave_one_out_v28.pdf"),
  width = 11,
  height = 8
)

main_figure_composite <- ((main_panel_effect_plot / main_leave_one_out_plot) | main_canonical_effect_plot) +
  plot_layout(widths = c(1.05, 1.2), guides = "collect") +
  plot_annotation(
    title = "EOCRC nitrogen diversion axis: effect size, metabolite order, and robustness",
    subtitle = paste(
      "LO effect = LOCRC - LOControl; EO effect = EOCRC - EOControl; EO > LO = disease-by-age interaction.",
      "Leave-one-out sensitivity uses the EO age-adjusted model."
    )
  ) &
  theme(legend.position = "bottom")

contrast_summary <- feature_results_display %>%
  group_by(model_label, contrast) %>%
  summarise(
    n_fdr_sig = sum(significant_fdr),
    n_fdr_fc_sig = sum(significant_fdr_fc),
    .groups = "drop"
  )

pooled_interaction_summary <- feature_results_pooled %>%
  group_by(model_label, contrast) %>%
  summarise(
    n_fdr_sig = sum(significant_fdr),
    n_fdr_fc_sig = sum(significant_fdr_fc),
    .groups = "drop"
  )

targeted_contrast_summary <- targeted_results_display %>%
  group_by(model_label, contrast, panel, panel_label, panel_tier, panel_tier_label) %>%
  summarise(
    n_fdr_sig = n_distinct(canonical_metabolite[significant_fdr]),
    n_fdr_fc_sig = n_distinct(canonical_metabolite[significant_fdr_fc]),
    n_panel_members = n_distinct(canonical_metabolite),
    .groups = "drop"
  )

write_output_csv(contrast_summary, file.path(table_dir, "contrast_level_hit_counts_display_models.csv"), row.names = FALSE)
write_output_csv(
  pooled_interaction_summary,
  file.path(table_dir, "contrast_level_hit_counts_pooled_interaction.csv"),
  row.names = FALSE
)
write_output_csv(
  targeted_contrast_summary,
  file.path(table_dir, "targeted_panel_hit_counts_by_contrast.csv"),
  row.names = FALSE
)

writeLines(
  c(
    paste0("Formal cross-omics framework (", analysis_version, ")"),
    "1. Primary discovery is restricted to the EO stratum (EOCRC vs EOControl), adjusted for continuous age, sex, and BMI.",
    "2. Candidate species must first associate with the target metabolite panel after covariate adjustment, rather than being selected only by crude abundance differences.",
    "3. EO disease association is evaluated separately so that panel-linked species can be distinguished from background panel correlates.",
    "4. EOCRC-only sensitivity models add tumor site and stage to test whether panel-species links persist within EO cancers.",
    "5. The former injurious_luminal_panel is split into biogenic amine, parent phenol, and host-processed aromatic toxin subpanels to avoid chemistry-driven signal cancellation.",
    "6. Pooled Disease x Age_group interaction is retained as secondary EO-specificity support, not as the primary discovery screen.",
    "7. HUMAnN MetaCyc/EC anchors are evaluated at the community level as mechanistic support, not as proof that a specific species encodes the function.",
    "8. Formal significance is defined only by FDR < 0.05 in this version; the former FDR < 0.10 exploratory layer is disabled."
  ),
  file.path(output_dir, "crossomics_analysis_framework_v28.txt")
)

load(file.path(analysis_dir, "FUSCC_metabonomics.Rdata"))

omics_sample_map <- FUSCC_cl2 %>%
  select(Patient_ID, Sample_ID) %>%
  distinct()

eo_crc_only_meta <- FUSCC_cl2 %>%
  filter(Age_group == "EO", Disease == "CRC") %>%
  mutate(Age_centered_crc = Age - mean(Age, na.rm = TRUE))

lo_crc_only_meta <- FUSCC_cl2 %>%
  filter(Age_group == "LO", Disease == "CRC") %>%
  mutate(Age_centered_crc = Age - mean(Age, na.rm = TRUE))

species_raw_mat <- prepare_crossomic_matrix(
  raw_df = meta_list_species[["FUSCC-SHSD"]],
  feature_col = "clade_name",
  sample_map = omics_sample_map
)
species_raw_mat <- species_raw_mat[, FUSCC_cl2$Sample_ID, drop = FALSE]

species_qc <- summarise_crossomic_qc(
  feature_mat = species_raw_mat,
  sample_meta = FUSCC_cl2,
  subset_mask = FUSCC_cl2$Age_group == "EO",
  group_var = "Disease",
  min_prevalence = species_prevalence_min,
  min_group_prevalence = species_group_prevalence_min
)
write_output_csv(species_qc, file.path(table_dir, "species_qc_summary_v28.csv"), row.names = FALSE)

keep_species <- species_qc %>%
  filter(keep_feature) %>%
  pull(feature_id)

species_clr_mat <- clr_transform_matrix(species_raw_mat[keep_species, , drop = FALSE])
species_annotation <- build_species_taxonomy_annotation(rownames(species_clr_mat))

species_eo_disease_fit <- fit_single_coefficient_limma(
  expr_mat = species_clr_mat[, eo_meta$Sample_ID, drop = FALSE],
  sample_meta = eo_meta,
  design_formula = ~ Disease + Age_centered + Sex + BMI,
  coefficient_expression = "DiseaseCRC",
  contrast_name = "EOCRC_vs_EOControl_ageadj",
  model_label = "species_EO_stratified_main",
  test_type = "ebayes",
  lfc_threshold = 0
)

species_lo_disease_fit <- fit_single_coefficient_limma(
  expr_mat = species_clr_mat[, lo_meta$Sample_ID, drop = FALSE],
  sample_meta = lo_meta,
  design_formula = ~ Disease + Age_centered + Sex + BMI,
  coefficient_expression = "DiseaseCRC",
  contrast_name = "LOCRC_vs_LOControl_ageadj",
  model_label = "species_LO_stratified_main",
  test_type = "ebayes",
  lfc_threshold = 0
)

species_pooled_interaction_fit <- fit_single_coefficient_limma(
  expr_mat = species_clr_mat,
  sample_meta = FUSCC_cl2,
  design_formula = ~ Disease * Age_group + Age_within_group + Sex + BMI,
  coefficient_expression = "DiseaseCRC.Age_groupEO",
  contrast_name = "Disease_age_interaction_pooled_ageadj",
  model_label = "species_pooled_interaction_ageadj",
  test_type = "ebayes",
  lfc_threshold = 0
)

species_panel_assoc_results <- run_panel_association_models(
  expr_mat = species_clr_mat,
  sample_meta_base = eo_meta,
  panel_score_long = panel_score_long,
  target_panels = eo_target_panels,
  design_formula = ~ panel_score + Disease + Age_centered + Sex + BMI,
  coefficient_expression = "panel_score",
  model_label = "species_EO_panel_association"
)

species_lo_panel_assoc_results <- run_panel_association_models(
  expr_mat = species_clr_mat,
  sample_meta_base = lo_meta,
  panel_score_long = panel_score_long,
  target_panels = eo_target_panels,
  design_formula = ~ panel_score + Disease + Age_centered + Sex + BMI,
  coefficient_expression = "panel_score",
  model_label = "species_LO_panel_association"
)

species_panel_sensitivity_results <- run_panel_association_models(
  expr_mat = species_clr_mat,
  sample_meta_base = eo_crc_only_meta,
  panel_score_long = panel_score_long,
  target_panels = eo_target_panels,
  design_formula = ~ panel_score + Age_centered_crc + Sex + BMI + Tumor_Site3 + Stage_collapsed,
  coefficient_expression = "panel_score",
  model_label = "species_EOCRC_panel_sensitivity"
)

species_lo_panel_sensitivity_results <- run_panel_association_models(
  expr_mat = species_clr_mat,
  sample_meta_base = lo_crc_only_meta,
  panel_score_long = panel_score_long,
  target_panels = eo_target_panels,
  design_formula = ~ panel_score + Age_centered_crc + Sex + BMI + Tumor_Site3 + Stage_collapsed,
  coefficient_expression = "panel_score",
  model_label = "species_LOCRC_panel_sensitivity"
)

write_output_csv(
  species_eo_disease_fit$results %>%
    left_join(species_annotation, by = c("metabolite_identification" = "feature_id")),
  file.path(table_dir, "species_eo_disease_results_v28.csv"),
  row.names = FALSE
)
write_output_csv(
  species_lo_disease_fit$results %>%
    left_join(species_annotation, by = c("metabolite_identification" = "feature_id")),
  file.path(table_dir, "species_lo_disease_results_v28.csv"),
  row.names = FALSE
)
write_output_csv(
  species_pooled_interaction_fit$results %>%
    left_join(species_annotation, by = c("metabolite_identification" = "feature_id")),
  file.path(table_dir, "species_pooled_interaction_results_v28.csv"),
  row.names = FALSE
)
write_output_csv(
  species_panel_assoc_results %>%
    left_join(species_annotation, by = c("metabolite_identification" = "feature_id")),
  file.path(table_dir, "species_panel_association_results_v28.csv"),
  row.names = FALSE
)
write_output_csv(
  species_lo_panel_assoc_results %>%
    left_join(species_annotation, by = c("metabolite_identification" = "feature_id")),
  file.path(table_dir, "species_lo_panel_association_results_v28.csv"),
  row.names = FALSE
)
write_output_csv(
  species_panel_sensitivity_results %>%
    left_join(species_annotation, by = c("metabolite_identification" = "feature_id")),
  file.path(table_dir, "species_panel_sensitivity_results_v28.csv"),
  row.names = FALSE
)
write_output_csv(
  species_lo_panel_sensitivity_results %>%
    left_join(species_annotation, by = c("metabolite_identification" = "feature_id")),
  file.path(table_dir, "species_lo_panel_sensitivity_results_v28.csv"),
  row.names = FALSE
)

species_panel_assoc_annotated <- build_species_association_table(
  panel_assoc_results = species_panel_assoc_results,
  feature_annotation = species_annotation,
  qc_tbl = species_qc
)

write_output_csv(
  species_panel_assoc_annotated,
  file.path(table_dir, "species_panel_association_annotated_v28.csv"),
  row.names = FALSE
)

write_output_csv(
  species_panel_assoc_annotated %>%
    filter(adj.P.Val < fdr_cutoff) %>%
    arrange(panel, adj.P.Val, desc(abs(logFC))),
  file.path(table_dir, "species_panel_top_hits_primary_v28.csv"),
  row.names = FALSE
)

top_species_plot <- plot_top_associated_species(
  species_panel_assoc_annotated,
  fdr_threshold = fdr_cutoff,
  top_n = top_species_hits_per_direction
)
if (!is.null(top_species_plot)) {
  save_gg(
    top_species_plot,
    file.path(plot_dir, "species_top_associated_forest_v28.pdf"),
    width = 12.5,
    height = 9.5
  )
}

species_cladogram_primary_data <- build_taxonomic_cladogram_data(
  species_panel_assoc_annotated,
  fdr_threshold = fdr_cutoff
)
if (!is.null(species_cladogram_primary_data)) {
  primary_cladogram_panels <- unique(as.character(species_cladogram_primary_data$species_points$panel))
  primary_cladogram_panels <- primary_cladogram_panels[!is.na(primary_cladogram_panels)]

  if (length(primary_cladogram_panels) == 1) {
    primary_plot <- plot_taxonomic_cladogram(
      species_cladogram_primary_data,
      panel_name = primary_cladogram_panels[[1]]
    )
    if (!is.null(primary_plot)) {
      save_gg(primary_plot, file.path(plot_dir, "species_taxonomic_cladogram_primary_v28.pdf"), width = 12, height = 9)
    }
  } else {
    walk(primary_cladogram_panels, function(panel_name) {
      panel_plot <- plot_taxonomic_cladogram(species_cladogram_primary_data, panel_name = panel_name)
      if (!is.null(panel_plot)) {
        output_name <- paste0(
          "species_taxonomic_cladogram_primary_",
          safe_file_tag(panel_name),
          "_v28.pdf"
        )
        save_gg(panel_plot, file.path(plot_dir, output_name), width = 12, height = 9)
      }
    })
  }
}

species_taxonomic_enrichment_primary <- run_taxonomic_enrichment(
  association_tbl = species_panel_assoc_annotated,
  direction_var = "association_direction",
  significance_layer = "FDR<0.05"
)

write_output_csv(
  species_taxonomic_enrichment_primary,
  file.path(table_dir, "species_taxonomic_enrichment_primary_v28.csv"),
  row.names = FALSE
)

species_taxonomic_enrichment_primary_plot <- plot_taxonomic_enrichment(
  species_taxonomic_enrichment_primary,
  adj_threshold = fdr_cutoff
)
if (!is.null(species_taxonomic_enrichment_primary_plot)) {
  save_gg(
    species_taxonomic_enrichment_primary_plot,
    file.path(plot_dir, "species_taxonomic_enrichment_primary_v28.pdf"),
    width = 14,
    height = 8
  )
}

species_candidates <- extract_crossomic_candidates(
  panel_assoc_results = species_panel_assoc_results,
  disease_results = species_eo_disease_fit$results,
  interaction_results = species_pooled_interaction_fit$results,
  sensitivity_results = species_panel_sensitivity_results,
  lo_sensitivity_results = species_lo_panel_sensitivity_results,
  lo_disease_results = species_lo_disease_fit$results,
  lo_panel_assoc_results = species_lo_panel_assoc_results,
  feature_annotation = species_annotation
) %>%
  mutate(
    target_panel_label = recode(panel, !!!panel_labels),
    eo_specificity_supported = interaction_fdr < fdr_cutoff & sign(interaction_beta) == sign(disease_beta)
  )

write_output_csv(species_candidates, file.path(table_dir, "species_candidate_summary_v28.csv"), row.names = FALSE)
write_output_csv(
  species_candidates %>%
    count(panel, target_panel_label, candidate_tier, eo_specificity_supported),
  file.path(table_dir, "species_candidate_counts_v28.csv"),
  row.names = FALSE
)
write_output_csv(
  species_candidates %>%
    filter(candidate_tier != "not_shortlisted"),
  file.path(table_dir, "species_candidate_shortlist_v28.csv"),
  row.names = FALSE
)

candidate_species_heatmap <- plot_candidate_species_heatmap(species_candidates)
if (!is.null(candidate_species_heatmap)) {
  save_gg(
    candidate_species_heatmap,
    file.path(plot_dir, "species_candidate_heatmap_v28.pdf"),
    width = 11.5,
    height = 9.5
  )
}

metacyc_raw_mat <- prepare_crossomic_matrix(
  raw_df = metacyc_list[["A30"]],
  feature_col = names(metacyc_list[["A30"]])[1],
  sample_map = omics_sample_map
)
metacyc_raw_mat <- metacyc_raw_mat[, FUSCC_cl2$Sample_ID, drop = FALSE]
metacyc_anchor_raw <- extract_anchor_matrix(
  feature_mat = metacyc_raw_mat,
  anchor_ids = manuscript_metacyc_anchor_ids,
  id_type = "pathway"
)

ec_raw_mat <- prepare_crossomic_matrix(
  raw_df = ec_list[["A30"]],
  feature_col = names(ec_list[["A30"]])[1],
  sample_map = omics_sample_map
)
ec_raw_mat <- ec_raw_mat[, FUSCC_cl2$Sample_ID, drop = FALSE]
ec_anchor_raw <- extract_anchor_matrix(
  feature_mat = ec_raw_mat,
  anchor_ids = manuscript_ec_anchor_ids,
  id_type = "ec"
)

metacyc_anchor_qc <- summarise_crossomic_qc(
  feature_mat = metacyc_anchor_raw,
  sample_meta = FUSCC_cl2,
  subset_mask = FUSCC_cl2$Age_group == "EO",
  group_var = "Disease",
  min_prevalence = function_prevalence_min,
  min_group_prevalence = function_group_prevalence_min
)
ec_anchor_qc <- summarise_crossomic_qc(
  feature_mat = ec_anchor_raw,
  sample_meta = FUSCC_cl2,
  subset_mask = FUSCC_cl2$Age_group == "EO",
  group_var = "Disease",
  min_prevalence = function_prevalence_min,
  min_group_prevalence = function_group_prevalence_min
)

write_output_csv(metacyc_anchor_qc, file.path(table_dir, "metacyc_anchor_qc_v28.csv"), row.names = FALSE)
write_output_csv(ec_anchor_qc, file.path(table_dir, "ec_anchor_qc_v28.csv"), row.names = FALSE)

metacyc_anchor_clr <- clr_transform_matrix(
  metacyc_anchor_raw[metacyc_anchor_qc$feature_id[metacyc_anchor_qc$keep_feature], , drop = FALSE]
)
ec_anchor_clr <- clr_transform_matrix(
  ec_anchor_raw[ec_anchor_qc$feature_id[ec_anchor_qc$keep_feature], , drop = FALSE]
)

metacyc_anchor_disease_fit <- fit_single_coefficient_limma(
  expr_mat = metacyc_anchor_clr[, eo_meta$Sample_ID, drop = FALSE],
  sample_meta = eo_meta,
  design_formula = ~ Disease + Age_centered + Sex + BMI,
  coefficient_expression = "DiseaseCRC",
  contrast_name = "EOCRC_vs_EOControl_ageadj",
  model_label = "metacyc_anchor_EO_stratified_main",
  test_type = "ebayes",
  lfc_threshold = 0
)
ec_anchor_disease_fit <- fit_single_coefficient_limma(
  expr_mat = ec_anchor_clr[, eo_meta$Sample_ID, drop = FALSE],
  sample_meta = eo_meta,
  design_formula = ~ Disease + Age_centered + Sex + BMI,
  coefficient_expression = "DiseaseCRC",
  contrast_name = "EOCRC_vs_EOControl_ageadj",
  model_label = "ec_anchor_EO_stratified_main",
  test_type = "ebayes",
  lfc_threshold = 0
)

metacyc_anchor_panel_results <- run_panel_association_models(
  expr_mat = metacyc_anchor_clr,
  sample_meta_base = eo_meta,
  panel_score_long = panel_score_long,
  target_panels = eo_target_panels,
  design_formula = ~ panel_score + Disease + Age_centered + Sex + BMI,
  coefficient_expression = "panel_score",
  model_label = "metacyc_anchor_EO_panel_association"
)
ec_anchor_panel_results <- run_panel_association_models(
  expr_mat = ec_anchor_clr,
  sample_meta_base = eo_meta,
  panel_score_long = panel_score_long,
  target_panels = eo_target_panels,
  design_formula = ~ panel_score + Disease + Age_centered + Sex + BMI,
  coefficient_expression = "panel_score",
  model_label = "ec_anchor_EO_panel_association"
)

metacyc_anchor_results_combined <- bind_rows(
  metacyc_anchor_disease_fit$results %>% mutate(anchor_type = "MetaCyc_pathway"),
  metacyc_anchor_panel_results %>% mutate(anchor_type = "MetaCyc_pathway")
)

ec_anchor_results_combined <- bind_rows(
  ec_anchor_disease_fit$results %>% mutate(anchor_type = "EC"),
  ec_anchor_panel_results %>% mutate(anchor_type = "EC")
)

write_output_csv(
  metacyc_anchor_results_combined,
  file.path(table_dir, "metacyc_anchor_results_v28.csv"),
  row.names = FALSE
)
write_output_csv(
  ec_anchor_results_combined,
  file.path(table_dir, "ec_anchor_results_v28.csv"),
  row.names = FALSE
)

candidate_species_ids <- species_candidates %>%
  filter(candidate_tier != "not_shortlisted") %>%
  pull(feature_id) %>%
  unique()

if (length(candidate_species_ids) > 0) {
  candidate_anchor_support <- bind_rows(
    map_dfr(eo_target_panels, function(panel_name) {
      panel_species_ids <- species_candidates %>%
        filter(panel == panel_name, candidate_tier != "not_shortlisted") %>%
        pull(feature_id) %>%
        unique()

      if (length(panel_species_ids) == 0 || nrow(metacyc_anchor_clr) == 0) {
        return(tibble())
      }

      map_dfr(panel_species_ids, function(species_id) {
        map_dfr(rownames(metacyc_anchor_clr), function(anchor_id) {
          lm_df <- eo_meta %>%
            mutate(
              species_value = species_clr_mat[species_id, Sample_ID],
              anchor_value = metacyc_anchor_clr[anchor_id, Sample_ID]
            )

          tidy(lm(anchor_value ~ species_value + Disease + Age_centered + Sex + BMI, data = lm_df)) %>%
            filter(term == "species_value") %>%
            transmute(
              panel = panel_name,
              target_panel_label = unname(panel_labels[panel_name]),
              anchor_type = "MetaCyc_pathway",
              anchor_id = anchor_id,
              feature_id = species_id,
              estimate,
              statistic,
              p.value
            )
        })
      })
    }),
    map_dfr(eo_target_panels, function(panel_name) {
      panel_species_ids <- species_candidates %>%
        filter(panel == panel_name, candidate_tier != "not_shortlisted") %>%
        pull(feature_id) %>%
        unique()

      if (length(panel_species_ids) == 0 || nrow(ec_anchor_clr) == 0) {
        return(tibble())
      }

      map_dfr(panel_species_ids, function(species_id) {
        map_dfr(rownames(ec_anchor_clr), function(anchor_id) {
          lm_df <- eo_meta %>%
            mutate(
              species_value = species_clr_mat[species_id, Sample_ID],
              anchor_value = ec_anchor_clr[anchor_id, Sample_ID]
            )

          tidy(lm(anchor_value ~ species_value + Disease + Age_centered + Sex + BMI, data = lm_df)) %>%
            filter(term == "species_value") %>%
            transmute(
              panel = panel_name,
              target_panel_label = unname(panel_labels[panel_name]),
              anchor_type = "EC",
              anchor_id = anchor_id,
              feature_id = species_id,
              estimate,
              statistic,
              p.value
            )
        })
      })
    })
  ) %>%
    group_by(panel, anchor_type) %>%
    mutate(adj.P.Val = p.adjust(p.value, method = "BH")) %>%
    ungroup() %>%
    left_join(species_annotation, by = "feature_id")

  write_output_csv(
    candidate_anchor_support,
    file.path(table_dir, "candidate_species_anchor_support_v28.csv"),
    row.names = FALSE
  )
} else {
  write_output_csv(
    tibble(),
    file.path(table_dir, "candidate_species_anchor_support_v28.csv"),
    row.names = FALSE
  )
}

main_figure_support_strip <- build_species_support_strip_data(species_candidates = species_candidates)

write_output_csv(
  main_figure_support_strip %>% select(-text_color, -max_abs),
  file.path(table_dir, "main_figure_support_strip_species_v28.csv"),
  row.names = FALSE
)

support_strip_plot <- plot_species_support_strip(main_figure_support_strip)
save_gg(
  support_strip_plot,
  file.path(plot_dir, "main_figure_support_strip_species_v28.pdf"),
  width =8,
  height = 5
)


# =========================
# =========================
left_margin   <- 0.02
right_margin  <- 0.02
top_margin    <- 0.02
bottom_margin <- 0.02

col_gap <- 0.015
row_gap <- 0.04

# =========================
# =========================
top_ratios <- c(3, 4, 6)

usable_width_top <- 1 - left_margin - right_margin - 2 * col_gap
unit_top <- usable_width_top / sum(top_ratios)

w_B <- top_ratios[1] * unit_top
w_C <- top_ratios[2] * unit_top
w_D <- top_ratios[3] * unit_top

x_B <- left_margin
x_C <- x_B + w_B + col_gap
x_D <- x_C + w_C + col_gap

# =========================
# =========================
row_ratios <- c(1, 1.8)
usable_height <- 1 - top_margin - bottom_margin - row_gap

h_top <- usable_height * row_ratios[1] / sum(row_ratios)
h_bottom <- usable_height * row_ratios[2] / sum(row_ratios)

y_bottom <- bottom_margin
y_top <- y_bottom + h_bottom + row_gap

# =========================
# =========================
w_E <- usable_width_top * 10 / 14
x_E <- left_margin

# x_E <- left_margin + (usable_width_top - w_E) / 2

# =========================
# =========================
main_figure_composite <- ggdraw() +
  draw_plot(strict_panel_score_plot,     x = x_B, y = y_top,    width = w_B, height = h_top) +
  draw_plot(main_panel_effect_plot,      x = x_C, y = y_top,    width = w_C, height = h_top) +
  draw_plot(main_canonical_effect_plot,  x = x_D, y = y_top,    width = w_D, height = h_top) +
  draw_plot(support_strip_plot,          x = x_E, y = y_bottom, width = w_E, height = h_bottom) +
  draw_plot_label(
    label = c("B", "C", "D", "E"),
    x = c(x_B, x_C, x_D, x_E),
    y = c(y_top + h_top, y_top + h_top, y_top + h_top, y_bottom + h_bottom),
    hjust = -0.05,
    vjust = 1.1,
    fontface = "bold",
    size = 26
  )

save_gg(
  main_figure_composite,
  file.path(plot_dir, "main_figure_nitrogen_diversion_v28.pdf"),
  width =22,
  height = 18
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

layout_design <- c(
  area(t = 1, l = 1,  b = 1, r = 2.5),  
  area(t = 2, l = 1,  b = 2, r = 6) 
)
supplementary_figure_composite <- all_panel_member_heatmap/free(extended_panel_score_plot, side = "l")+
  plot_layout(
    design = layout_design,
    heights = c(2,1)
  )+ plot_annotation(
    tag_levels = list(c("A", "B"))
  ) &
  theme(
    plot.tag = element_text(size = 26, face = "bold"),
    plot.tag.position = c(0.02, 1)
  )


save_gg(
  supplementary_figure_composite,
  file.path(plot_dir, "supplementary_figure_composite_v28.pdf"),
  width =15,
  height = 22
)

