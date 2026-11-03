library(data.table)
library(sandwich)

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
prepared <- readRDS(file.path(output_dir, "objects/prepared_data.rds"))
ctx <- prepared$ctx
prep <- prepared$prep
fixed <- readRDS(file.path(output_dir, "objects/fixed_scheme_statistics.rds"))
original <- fread(file.path(
  analysis_dir,
  "extended_metabonomics_results_v29/tables/targeted_panel_membership.csv"
))

# Keep the nine exploratory definitions and their separate correction families.
extended_lookup <- unique(original[
  panel_tier == "extended_class_based",
  .(panel, panel_label, panel_tier, panel_tier_label, canonical_metabolite)
])
extended_panels <- unique(extended_lookup$panel)
stopifnot(length(extended_panels) == 9L)
display_labels <- c(
  amino_acid_putrefaction_auxiliary_panel = "Amino-acid catabolites",
  tryptophan_indole_branch_panel = "Tryptophan-indole metabolites",
  polyamine_lysine_putrefaction_panel = "Polyamine-lysine metabolites",
  biogenic_amine_panel = "Biogenic amines",
  indole_catabolism_panel = "Indole metabolites",
  aromatic_host_cometabolite_panel = "Aromatic co-metabolites",
  injurious_biogenic_amine_subpanel = "Selected biogenic amines",
  injurious_parent_phenol_subpanel = "Parent phenolic compounds",
  injurious_host_processed_aromatic_toxin_subpanel = "Conjugated aromatic metabolites"
)
extended_lookup[, panel_label := unname(display_labels[panel])]
stopifnot(!anyNA(extended_lookup$panel_label))
fixed_lookup <- rbindlist(lapply(
  names(ctx$members),
  function(p) data.table(
    panel = p,
    panel_label = if (p == "U_P")
      "Arginine-ornithine-proline metabolites" else
      "Aminovalerate-polyamine metabolites",
    panel_tier = "fixed_core",
    panel_tier_label = "Fixed representative core",
    canonical_metabolite = ctx$members[[p]]
  )
))
lookup <- rbindlist(list(fixed_lookup, extended_lookup))
name_map <- unique(original[, .(canonical_metabolite, metabolite_identification)])
mapping <- merge(
  name_map,
  prep$feature_qc,
  by.x = "metabolite_identification",
  by.y = "v29_label",
  all.x = TRUE,
  allow.cartesian = TRUE
)
stopifnot(!anyNA(mapping$peak_id), all(mapping$peak_id %in% rownames(prep$log_mats$PQN_log2)))
mapping <- mapping[
  !canonical_metabolite %in% ctx$mapping$canonical_metabolite,
  .(
    canonical_metabolite,
    peak_id,
    source_name = metabolite_identification,
    mode,
    qc_cv_percent,
    detection_fraction,
    keep_primary,
    mapping_rule = "v29_selected_synonym_peaks_equal_average"
  )
]
core_mapping <- ctx$mapping[, .(
  canonical_metabolite,
  peak_id,
  source_name = metabolite_identification,
  mode,
  qc_cv_percent,
  detection_fraction = sample_detection_fraction,
  keep_primary = pass_v29_detection,
  mapping_rule = "frozen_representative_core"
)]
mapping <- unique(rbindlist(list(core_mapping, mapping)))
canonical_order <- unique(lookup$canonical_metabolite)
stopifnot(setequal(canonical_order, mapping$canonical_metabolite), all(mapping$keep_primary))
log_mat <- prep$log_mats$PQN_log2
peak_z <- t(scale(t(log_mat)))
canonical_z <- canonical_log <- matrix(
  NA_real_,
  length(canonical_order),
  nrow(ctx$meta),
  dimnames = list(canonical_order, ctx$meta$Sample_ID)
)

for (met in canonical_order) {
  ids <- unique(mapping[canonical_metabolite == met, peak_id])
  canonical_z[met, ] <- colMeans(peak_z[ids, ctx$meta$Sample_ID, drop = FALSE])
  canonical_log[met, ] <- colMeans(log_mat[ids, ctx$meta$Sample_ID, drop = FALSE])
}

stopifnot(
  max(abs(canonical_z[rownames(ctx$canonical_z), ] - ctx$canonical_z)) < 1e-12,
  max(abs(canonical_log[rownames(ctx$member_log2), ] - ctx$member_log2)) < 1e-12
)
members <- split(extended_lookup$canonical_metabolite, extended_lookup$panel)
scores <- t(vapply(
  members,
  function(m) colMeans(canonical_z[unique(m), , drop = FALSE]),
  numeric(nrow(ctx$meta))
))
colnames(scores) <- ctx$meta$Sample_ID
fwrite(mapping, file.path(v35_dir, "tables/supplementary_peak_mapping.csv"), bom = TRUE)
fwrite(lookup, file.path(v35_dir, "tables/supplementary_panel_membership.csv"), bom = TRUE)
catalog <- extended_lookup[, .(
  panel_label = unique(panel_label),
  n_members = uniqueN(canonical_metabolite),
  members = paste(unique(canonical_metabolite), collapse = ";")
), by = panel]
fwrite(catalog, file.path(v35_dir, "tables/supplementary_panel_catalog.csv"), bom = TRUE)

age_fits <- lapply(names(members), function(p) fit_age_effects(scores[p, ], ctx$meta, p))
names(age_fits) <- names(members)
effects <- rbindlist(lapply(age_fits, function(x) x$results[inference == "HC3_normal"]))
effects[, c("preprocessing", "analysis_role") := list("PQN_log2", "supplementary_exploratory")]
effects[, q_supp_panel_BH27 := p.adjust(p_value, "BH")]
effects[, significant_fdr := q_supp_panel_BH27 <= 0.1]
effects[, correction_family := "9 extended panels x 3 disease effects = 27"]
fwrite(
  effects,
  file.path(v35_dir, "tables/supplementary_panel_effects_HC3_BH27.csv"),
  bom = TRUE
)
designs <- list(
  clinical_multicategory = fixed$crc_design_multicategory,
  clinical_binary = fixed$crc_design_binary
)
crc_fits <- list()

for (model_name in names(designs)) {
  for (p in names(members)) {
    crc_fits[[paste(model_name, p)]] <- fit_crc_effect(
      scores[p, rownames(designs[[model_name]])],
      designs[[model_name]],
      p,
      model_name
    )
  }
}

crc <- rbindlist(lapply(crc_fits, function(x) x$results[inference == "HC3_normal"]))
crc[, preprocessing := "PQN_log2"]
crc[, analysis_role := ifelse(
  model == "clinical_binary",
  "supplementary_primary_encoding",
  "supplementary_encoding_sensitivity"
)]
crc[, q_supp_CRC_BH9 := p.adjust(p_value, "BH"), by = model]
crc[, significant_fdr := q_supp_CRC_BH9 <= 0.1]
fwrite(crc, file.path(v35_dir, "tables/supplementary_CRC_effects_HC3_BH9.csv"), bom = TRUE)

member_fits <- lapply(canonical_order, function(m) fit_age_effects(canonical_z[m, ], ctx$meta, m))
names(member_fits) <- canonical_order
member_effects <- rbindlist(lapply(member_fits, function(x) x$results[inference == "HC3_normal"]))
setnames(member_effects, "panel", "canonical_metabolite")
member_effects[, c(
  "preprocessing",
  "analysis_role"
) := list(
  "PQN_log2",
  "supplementary_exploratory"
)]
member_effects[, q_supp_member_BH := p.adjust(p_value, "BH")]
member_effects[, significant_fdr := q_supp_member_BH <= 0.1]
member_effects[, correction_family := paste0(
  length(canonical_order),
  " unique members x 3 effects = ",
  .N
)]
fwrite(
  member_effects,
  file.path(v35_dir, "tables/supplementary_member_effects_HC3_BH.csv"),
  bom = TRUE
)

score_long <- melt(
  data.table(Sample_ID = ctx$meta$Sample_ID, as.data.table(t(scores))),
  id.vars = "Sample_ID",
  variable.name = "panel",
  value.name = "panel_score",
  variable.factor = FALSE
)
score_long <- merge(
  merge(score_long, unique(extended_lookup[, .(panel, panel_label)]), by = "panel"),
  ctx$meta,
  by = "Sample_ID",
  sort = FALSE
)
member_long <- merge(
  melt(
    data.table(Sample_ID = ctx$meta$Sample_ID, as.data.table(t(canonical_log))),
    id.vars = "Sample_ID",
    variable.name = "canonical_metabolite",
    value.name = "log2_intensity",
    variable.factor = FALSE
  ),
  ctx$meta,
  by = "Sample_ID",
  sort = FALSE
)
raw <- readRDS(file.path(analysis_dir, "v30/peak_audit_inputs_v30.rds"))$raw_all
raw <- as.data.table(raw)
coverage <- rbindlist(lapply(
  names(members),
  function(p) {
    detected <- vapply(
      members[[p]],
      function(m) {
        ids <- mapping[canonical_metabolite == m, peak_id]
        x <- as.matrix(raw[match(ids, peak_id), ctx$meta$Sample_ID, with = FALSE])
        colSums(is.finite(x) & x > 0) > 0
      },
      logical(nrow(ctx$meta))
    )
    data.table(
      Sample_ID = ctx$meta$Sample_ID,
      panel = p,
      observed_members = rowSums(detected),
      total_members = ncol(detected),
      observed_fraction = rowMeans(detected)
    )
  }
))
fwrite(coverage, file.path(v35_dir, "tables/supplementary_score_coverage.csv"), bom = TRUE)
fwrite(
  score_long,
  file.path(v35_dir, "tables/supplementary_panel_scores_long.csv"),
  bom = TRUE
)
fwrite(member_long, file.path(v35_dir, "tables/supplementary_members_long.csv"), bom = TRUE)
saveRDS(
  list(
    lookup = lookup,
    extended_lookup = extended_lookup,
    mapping = mapping,
    canonical_order = canonical_order,
    canonical_z = canonical_z,
    canonical_log = canonical_log,
    members = members,
    scores = scores,
    effects = effects,
    crc = crc,
    member_effects = member_effects,
    score_long = score_long,
    member_long = member_long,
    coverage = coverage,
    age_fits = age_fits,
    crc_fits = crc_fits,
    member_fits = member_fits
  ),
  file.path(v35_dir, "objects/supplementary_statistics.rds")
)
