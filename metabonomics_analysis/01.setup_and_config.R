
suppressPackageStartupMessages({
  library(readxl)
  library(limma)
  library(edgeR)
  library(tidyverse)
  library(data.table)
  library(broom)
  library(purrr)
  library(emmeans)
  library(vegan)
  library(ggrepel)
  library(patchwork)
})

disk_symbol <- "J:/"
base_symbol <- "/home/global_EOCRC_microbiology/"
project_root <- base_symbol
analysis_version <- "v29"
downstream_dir <- file.path(project_root, "downstream_analysis")
analysis_dir <- file.path(project_root, "FUSCC_metabonomics_analysis")

resolve_single_path <- function(paths, label) {
  hits <- unique(as.character(paths))
  hits <- hits[!is.na(hits) & nzchar(hits)]
  hits <- hits[file.exists(hits)]
  if (length(hits) != 1) {
    stop(label, ": expected exactly one existing file, found ", length(hits), ".")
  }
  normalizePath(hits[[1]], winslash = "/", mustWork = TRUE)
}

resolve_cohort_file <- function(cohort_dir, file_spec, label, optional = FALSE) {
  if (is.na(file_spec) || !nzchar(file_spec)) {
    if (optional) {
      return(NA_character_)
    }
    stop(label, ": missing file specification.")
  }

  candidate_paths <- if (grepl("[*?]", file_spec)) {
    Sys.glob(file.path(cohort_dir, file_spec))
  } else {
    file.path(cohort_dir, file_spec)
  }

  hits <- unique(candidate_paths[file.exists(candidate_paths)])
  if (optional && length(hits) == 0) {
    return(NA_character_)
  }
  if (length(hits) != 1) {
    stop(label, ": expected exactly one matching file under ", cohort_dir, ", found ", length(hits), ".")
  }
  normalizePath(hits[[1]], winslash = "/", mustWork = TRUE)
}

meta_path <- resolve_single_path(
  Sys.glob(file.path(project_root, "SRA_metadata*", "A30_FUSCC_Gut", "metadata_for_revision.xlsx")),
  "Metabolomics metadata workbook"
)

output_dir <- file.path(downstream_dir, "output", "metabonomics_analysis", paste0("extended_metabonomics_results_", analysis_version))
plot_dir <- file.path(output_dir, "plots")
table_dir <- file.path(output_dir, "tables")

rna_root_dir <- file.path(project_root, "FUSCC_RNA_seq")
rna_cohort_config <- tribble(
  ~cohort_id, ~cohort_label, ~cohort_role, ~subdir, ~metadata_spec, ~qc_spec, ~qc_mode,
  "microbiome", "Discovery: matched microbiome-RNA cohort", "discovery", "Microbiome_RNA_Seq", "62*.xlsx", NA_character_, "count_only",
  "cancer_cell", "Validation: cancer-cell RNA cohort", "validation", "Cancer_cell_RNA_Seq", "meta.csv", "rna_info.csv", "count_only",
  "multi_omics", "Validation: multi-omics RNA cohort", "validation", "Multi_omics_RNA_Seq", "*20240813.csv", NA_character_, "strict_qc_label"
) %>%
  mutate(
    cohort_dir = file.path(rna_root_dir, subdir),
    annotation_path = map_chr(cohort_dir, ~ resolve_cohort_file(.x, "1_genes_fpkm_expression.txt", paste("RNA annotation for", basename(.x)))),
    count_path = map_chr(cohort_dir, ~ resolve_cohort_file(.x, "gene_count_matrix.txt", paste("RNA count matrix for", basename(.x)))),
    metadata_path = map2_chr(cohort_dir, metadata_spec, ~ resolve_cohort_file(.x, .y, paste("RNA metadata for", basename(.x)))),
    qc_path = map2_chr(cohort_dir, qc_spec, ~ resolve_cohort_file(.x, .y, paste("RNA QC file for", basename(.x)), optional = TRUE))
  )
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

bulk_rna_host_signature_definitions <- tribble(
  ~module, ~gene_symbol, ~module_label, ~module_category, ~module_order,
  "host_cationic_aa_uptake_signature", "SLC7A1", "Host cationic AA uptake", "host_signature", 1,
  "host_cationic_aa_uptake_signature", "SLC6A14", "Host cationic AA uptake", "host_signature", 1,
  "host_cationic_aa_uptake_signature", "SLC38A2", "Host cationic AA uptake", "host_signature", 1,
  "host_arginine_ornithine_diversion_signature", "ARG2", "Host arginine-to-ornithine diversion", "host_signature", 2,
  "host_arginine_ornithine_diversion_signature", "OAT", "Host arginine-to-ornithine diversion", "host_signature", 2,
  "host_arginine_ornithine_diversion_signature", "SLC25A15", "Host arginine-to-ornithine diversion", "host_signature", 2,
  "host_polyamine_synthesis_signature", "ODC1", "Host polyamine synthesis", "host_signature", 3,
  "host_polyamine_synthesis_signature", "AMD1", "Host polyamine synthesis", "host_signature", 3,
  "host_polyamine_synthesis_signature", "AZIN1", "Host polyamine synthesis", "host_signature", 3,
  "host_polyamine_synthesis_signature", "SRM", "Host polyamine synthesis", "host_signature", 3,
  "host_polyamine_catabolic_stress_signature", "SAT1", "Host polyamine catabolic stress", "host_signature", 4,
  "host_polyamine_catabolic_stress_signature", "SMOX", "Host polyamine catabolic stress", "host_signature", 4,
  "host_polyamine_catabolic_stress_signature", "PAOX", "Host polyamine catabolic stress", "host_signature", 4
)

bulk_rna_discovery_target_panels <- c(
  "arginine_ornithine_routing_panel",
  "amino_acid_putrefaction_core_panel"
)

bulk_rna_cohort_levels <- c(
  "microbiome",
  "cancer_cell",
  "multi_omics"
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
