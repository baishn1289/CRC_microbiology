library(data.table)
library(sandwich)
library(limma)

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
panels <- c("U_P", "D0")
contrast_order <- c("EO", "LO", "EO_minus_LO")
fdr_cutoff <- 0.1

# Panel-score effects: six tests in one reporting family.
age_fits <- lapply(panels, function(p) fit_age_effects(ctx$scores[p, ], ctx$meta, p))
names(age_fits) <- panels
all_methods <- rbindlist(lapply(age_fits, function(x) x$results))
all_methods[, q_BH6 := p.adjust(p_value, "BH"), by = inference]
all_methods[, significant_fdr := q_BH6 <= fdr_cutoff]
mean_results <- copy(all_methods[inference == "HC3_normal"])
mean_results[, c("outcome_id", "preprocessing", "variant", "correction_family", "analysis_role") :=
  list(
    paste0("representative_core__", panel),
    "PQN_log2",
    "representative_core",
    "2 fixed panels x 3 disease contrasts = 6",
    "primary"
  )]
mean_results <- mean_results[order(match(panel, panels), match(contrast, contrast_order))]
fwrite(mean_results, file.path(output_dir, "tables/panel_mean_effects_BH6.csv"))

# Clinical comparisons use the same 302 complete cases and two encodings.
designs <- prepared$crc_designs
crc_fits <- list()

for (model_name in names(designs)) {
  design <- designs[[model_name]]

  for (p in panels) {
    crc_fits[[paste(model_name, p)]] <- fit_crc_effect(
      ctx$scores[p, rownames(design)],
      design,
      p,
      model_name
    )
  }
}

crc_all_methods <- rbindlist(lapply(crc_fits, function(x) x$results))
crc_all_methods[, q_BH2 := p.adjust(p_value, "BH"), by = .(model, inference)]
crc_all_methods[, significant_fdr := q_BH2 <= fdr_cutoff]
crc <- copy(crc_all_methods[inference == "HC3_normal"])
crc[, formula := unname(prepared$crc_formulas[model])]
crc[, preprocessing := "PQN_log2"]
crc[, variant := "representative_core"]
crc[, analysis_role := ifelse(model == "clinical_binary", "primary", "sensitivity")]
fwrite(crc, file.path(output_dir, "tables/CRC_direct_comparison_BH2.csv"))

# Members: all three effects remain together within U_P BH15 and D0 BH9.
member_order <- unlist(ctx$members, use.names = FALSE)
member_fits <- lapply(member_order, function(m) fit_age_effects(ctx$canonical_z[m, ], ctx$meta, m))
names(member_fits) <- member_order
member_results <- rbindlist(lapply(member_fits, function(x) x$results[inference == "HC3_normal"]))
setnames(member_results, "panel", "canonical_metabolite")
member_results[, panel := vapply(
  canonical_metabolite,
  function(m) names(ctx$members)[vapply(ctx$members, function(ms) m %in% ms, logical(1))],
  character(1)
)]
member_results <- member_results[order(
  match(panel, panels),
  match(canonical_metabolite, member_order),
  match(contrast, contrast_order)
)]
member_results[, member_family_n := .N, by = panel]
member_results[, q_member_within_panel := p.adjust(p_value, "BH"), by = panel]
member_results[, significant_fdr := q_member_within_panel <= fdr_cutoff]
member_results[, c("preprocessing", "variant", "analysis_role") :=
  list("PQN_log2", "representative_core", "primary_member")]
member_results[, correction_family := ifelse(
  panel == "U_P",
  "U_P: 5 members x 3 effects = 15 (BH15)",
  "D0: 3 members x 3 effects = 9 (BH9)"
)]
stopifnot(
  all(member_results[panel == "U_P", member_family_n] == 15L),
  all(member_results[panel == "D0", member_family_n] == 9L)
)
stopifnot(
  nrow(member_results) == 24L,
  !anyDuplicated(member_results[, .(canonical_metabolite, contrast)])
)
fwrite(
  member_results,
  file.path(v35_dir, "tables/member_effects_within_panel_BH15_BH9.csv"),
  bom = TRUE
)

# Secondary set tests retain their original seeds and rotation count.
set_rows <- settings <- list()

for (j in seq_along(contrast_order)) {
  contrast_name <- contrast_order[j]
  seed <- 20261011L + j
  set.seed(seed)
  tab <- limma::mroast(
    ctx$background_z,
    index = ctx$set_ids,
    design = ctx$design,
    contrast = ctx$contrasts[, contrast_name],
    set.statistic = "mean",
    nrot = 19999L,
    approx.zscore = TRUE,
    robust = TRUE,
    trend = FALSE,
    midp = FALSE,
    sort = "none"
  )
  result <- as.data.table(tab, keep.rownames = "panel")[, .(
    panel,
    NGenes,
    PropDown,
    PropUp,
    Direction,
    PValue,
    PValue.Mixed
  )]
  result[, contrast := contrast_name]
  set_rows[[contrast_name]] <- result
  settings[[contrast_name]] <- data.table(
    contrast = contrast_name,
    seed = seed,
    nrot = 19999L,
    n_features = nrow(ctx$background_z),
    n_samples = ncol(ctx$background_z),
    n_sets = 2L,
    set_statistic = "mean",
    midp = FALSE,
    robust = TRUE,
    trend = FALSE,
    approx_zscore = TRUE
  )
}

sets <- list(results = rbindlist(set_rows))
sets$results[, q_Mixed_BH6 := p.adjust(PValue.Mixed, "BH")]
sets$results[, q_directional_BH6 := p.adjust(PValue, "BH")]
sets$results[, significant_Mixed := q_Mixed_BH6 <= fdr_cutoff]
sets$results[, significant_directional := q_directional_BH6 <= fdr_cutoff]
sets$results[, MCSE_mixed_approx := sqrt(PValue.Mixed * (1 - PValue.Mixed) / 20000)]
sets$results[, c("preprocessing", "variant", "analysis_role") :=
  list("PQN_log2", "representative_core", "secondary_member_set")]
sets$settings <- rbindlist(settings)
fwrite(sets$results, file.path(v35_dir, "tables/member_set_results_BH6.csv"), bom = TRUE)
fwrite(sets$settings, file.path(v35_dir, "tables/member_set_test_settings.csv"), bom = TRUE)

# P2 leave-one-member-out diagnostics do not redefine the panel.
loo_fits <- list()
loo_rows <- list()

for (omitted in ctx$members$D0) {
  remaining <- setdiff(ctx$members$D0, omitted)
  score <- colMeans(ctx$canonical_z[remaining, , drop = FALSE])
  fit_bundle <- fit_age_effects(score, ctx$meta, "D0")
  loo_fits[[omitted]] <- fit_bundle
  loo_rows[[omitted]] <- copy(fit_bundle$results[inference == "HC3_normal"])[
    ,
    c(
      "omitted_member",
      "remaining_members",
      "analysis_role"
    ) := list(
      omitted,
      length(remaining),
      "leave_one_out_diagnostic"
    )
  ]
}

loo <- rbindlist(loo_rows)
loo[, q_diagnostic_BH3 := p.adjust(p_value, "BH"), by = contrast]
loo[, correction_family := "3 omissions per contrast; diagnostic only"]
full_reference <- copy(mean_results[panel == "D0"])[
  ,
  c(
    "omitted_member",
    "remaining_members",
    "analysis_role"
  ) := list(
    "none",
    3L,
    "full_panel_primary_reference"
  )
]
loo <- rbindlist(list(full_reference, loo), fill = TRUE)
fwrite(loo, file.path(v35_dir, "tables/D0_leave_one_out_HC3_BH3.csv"), bom = TRUE)

score_long <- rbindlist(lapply(
  panels,
  function(p) data.table(
    Sample_ID = ctx$meta$Sample_ID,
    panel = p,
    panel_score = ctx$scores[p, ]
  )
))
score_long <- merge(
  score_long,
  ctx$meta,
  by = "Sample_ID",
  sort = FALSE
)
z_long <- melt(
  data.table(Sample_ID = ctx$meta$Sample_ID, as.data.table(t(ctx$canonical_z))),
  id.vars = "Sample_ID",
  variable.name = "canonical_metabolite",
  value.name = "member_z",
  variable.factor = FALSE
)
log_long <- melt(
  data.table(Sample_ID = ctx$meta$Sample_ID, as.data.table(t(ctx$member_log2))),
  id.vars = "Sample_ID",
  variable.name = "canonical_metabolite",
  value.name = "log2_intensity",
  variable.factor = FALSE
)
member_long <- merge(
  merge(z_long, log_long, by = c("Sample_ID", "canonical_metabolite")),
  ctx$meta,
  by = "Sample_ID",
  sort = FALSE
)
member_long[, panel := ifelse(canonical_metabolite %in% ctx$members$U_P, "U_P", "D0")]

fwrite(score_long, file.path(output_dir, "tables/fixed_panel_scores_long.csv"))
fwrite(member_long, file.path(output_dir, "tables/fixed_members_long.csv"))
saveRDS(
  list(
    mean_results = mean_results,
    member_results = member_results,
    crc = crc,
    leave_one_out = loo,
    score_long = score_long,
    member_long = member_long,
    crc_design_multicategory = designs$clinical_multicategory,
    crc_design_binary = designs$clinical_binary
  ),
  file.path(output_dir, "objects/fixed_scheme_statistics.rds")
)
