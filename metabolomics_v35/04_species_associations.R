library(data.table)
library(tidyverse)
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
ctx <- readRDS(file.path(output_dir, "objects/fixed_scheme_context.rds"))
fixed <- readRDS(file.path(output_dir, "objects/fixed_scheme_statistics.rds"))
meta <- as.data.frame(ctx$meta)
meta$Disease <- factor(meta$Disease, levels = c("Control", "CRC"))
meta$Sex <- factor(meta$Sex, levels = c("female", "male"))

# Species prevalence is filtered separately by age stratum before union CLR.
raw_env <- new.env()
load(file.path(analysis_dir, "FUSCC_metabonomics.Rdata"), envir = raw_env)
raw <- prepare_crossomic_matrix(
  raw_env$meta_list_species[["FUSCC-SHSD"]],
  "clade_name",
  as_tibble(meta[, c("Patient_ID", "Sample_ID")])
)
stopifnot(setequal(colnames(raw), meta$Sample_ID), !anyDuplicated(rownames(raw)))
raw <- raw[, meta$Sample_ID, drop = FALSE]
stopifnot(all(is.finite(raw)), all(raw >= 0))

qc <- lapply(
  c("EO", "LO"),
  function(g) as.data.table(summarise_crossomic_qc(
    raw,
    meta,
    meta$Age_group == g,
    "Disease",
    min_prevalence = 0.10,
    min_group_prevalence = 0.15
  ))
)
names(qc) <- c("EO", "LO")
eligible <- lapply(qc, function(d) d[keep_feature == TRUE, feature_id])
union_ids <- union(eligible$EO, eligible$LO)
both_ids <- intersect(eligible$EO, eligible$LO)
clr <- clr_transform_matrix(raw[union_ids, , drop = FALSE])
pseudocount <- min(raw[union_ids, , drop = FALSE][raw[union_ids, , drop = FALSE] > 0]) / 2
stopifnot(
  length(eligible$EO) == 368L,
  length(eligible$LO) == 362L,
  length(union_ids) == 392L,
  length(both_ids) == 338L,
  max(abs(colMeans(clr))) < 1e-10
)
fwrite(
  data.table(
    species_set = c("EO", "LO", "CLR_union", "both_strata"),
    n_species = c(368L, 362L, 392L, 338L),
    pseudocount = pseudocount
  ),
  file.path(v35_dir, "tables/species_QC_and_CLR_settings_v35.csv"),
  bom = TRUE
)

fit_species_matrix <- function(y, x, term, age_group, panel, analysis_role) {
  stopifnot(identical(colnames(y), rownames(x)), qr(x)$rank == ncol(x))
  fit <- lm.fit(x, t(y))
  bread <- solve(crossprod(x))
  leverage <- rowSums((x %*% bread) * x)
  influence <- drop(x %*% bread[, term])
  variance <- colSums((fit$residuals / (1 - leverage))^2 * influence^2)
  estimate <- as.numeric(fit$coefficients[term, ])
  se <- sqrt(variance)

  list(
    results = data.table(
      feature_id = rownames(y),
      age_group = age_group,
      panel = panel,
      analysis_role = analysis_role,
      n = ncol(y),
      estimate = estimate,
      std_error = se,
      ci_low = estimate - qnorm(0.975) * se,
      ci_high = estimate + qnorm(0.975) * se,
      p_value = 2 * pnorm(-abs(estimate / se)),
      inference = "OLS_HC3_normal"
    ),
    design = x,
    term = term,
    feature_ids = rownames(y)
  )
}

bundles <- list()

for (g in c("EO", "LO")) {
  d <- meta[meta$Age_group == g, , drop = FALSE]
  x <- model.matrix(~ Disease + Age + Sex + BMI, d)
  rownames(x) <- d$Sample_ID
  y <- clr[eligible[[g]], d$Sample_ID, drop = FALSE]
  bundles[[paste0(
    "disease_",
    g
  )]] <- fit_species_matrix(
    y,
    x,
    "DiseaseCRC",
    g,
    NA_character_,
    "disease"
  )

  for (p in c("U_P", "D0")) {
    xa <- cbind(x, panel_score = ctx$scores[p, d$Sample_ID])
    bundles[[paste("association", g, p, sep = "_")]] <-
      fit_species_matrix(y, xa, "panel_score", g, p, "panel_association")
  }

  clinical_ids <- intersect(d$Sample_ID, rownames(fixed$crc_design_binary))
  xc <- fixed$crc_design_binary[clinical_ids, colnames(fixed$crc_design_binary) != "Age_groupEO", drop = FALSE]
  xc <- cbind(xc, Age = meta$Age[match(clinical_ids, meta$Sample_ID)])

  for (p in c("U_P", "D0")) {
    xa <- cbind(xc, panel_score = ctx$scores[p, clinical_ids])
    bundles[[paste("CRC", g, p, sep = "_")]] <-
      fit_species_matrix(
        clr[eligible[[g]], clinical_ids, drop = FALSE],
        xa,
        "panel_score",
        g,
        p,
        "CRC_sensitivity"
      )
  }
}

all_results <- rbindlist(lapply(bundles, function(b) b$results))
all_results[, q_BH := p.adjust(p_value, "BH"), by = analysis_role]
all_results[, family_n := .N, by = analysis_role]
all_results[, significant_fdr := q_BH <= 0.1]
all_results[, species_label := gsub("_", " ", sub("^.*\\|s__", "", feature_id))]
stopifnot(
  all_results[analysis_role == "disease", .N] == 730L,
  all_results[analysis_role == "panel_association", .N] == 1460L,
  all_results[analysis_role == "CRC_sensitivity", .N] == 1460L
)
all_results[, preprocessing := "PQN_log2_panel_score; natural_log_CLR_species"]
fwrite(
  all_results,
  file.path(v35_dir, "tables/species_all_results_HC3_v35.csv"),
  bom = TRUE
)

eo <- all_results[analysis_role == "disease" & age_group == "EO" & feature_id %in% both_ids]
lo <- all_results[analysis_role == "disease" & age_group == "LO" & feature_id %in% both_ids]
delta <- merge(
  eo,
  lo,
  by = c("feature_id", "species_label"),
  suffixes = c("_EO", "_LO")
)
delta[, estimate := estimate_EO - estimate_LO]
delta[, std_error := sqrt(std_error_EO^2 + std_error_LO^2)]
delta[, c("contrast", "p_value", "ci_low", "ci_high") := list(
  "EO_minus_LO",
  2 * pnorm(-abs(estimate / std_error)),
  estimate - qnorm(0.975) * std_error,
  estimate + qnorm(0.975) * std_error
)]
delta[, q_BH338 := p.adjust(p_value, "BH")]
fwrite(
  delta[, .(
    feature_id,
    species_label,
    contrast,
    estimate,
    std_error,
    ci_low,
    ci_high,
    p_value,
    q_BH338
  )],
  file.path(v35_dir, "tables/species_disease_Delta_HC3_BH338_v35.csv"),
  bom = TRUE
)

eo_disease <- all_results[
  analysis_role == "disease" & age_group == "EO",
  .(
    feature_id,
    species_label,
    disease_beta = estimate,
    disease_p = p_value,
    disease_q = q_BH
  )
]
eo_assoc <- all_results[
  analysis_role == "panel_association" & age_group == "EO" & panel == "D0",
  .(
    feature_id,
    association_beta = estimate,
    association_p = p_value,
    association_q = q_BH
  )
]
candidates <- merge(eo_disease, eo_assoc, by = "feature_id")
candidates[, selected_for_strip := disease_q <= 0.1 & association_q <= 0.1 &
  sign(disease_beta) == sign(association_beta)]
setorder(candidates, disease_beta)
selected_ids <- candidates[selected_for_strip == TRUE, feature_id]
stopifnot(length(selected_ids) > 0)
strip <- all_results[feature_id %in% selected_ids &
  (analysis_role == "disease" | (analysis_role == "panel_association" & panel == "D0"))]
strip[, column := paste(age_group, ifelse(analysis_role == "disease", "disease", "D0 assoc"))]
strip[, correction_family := paste0(analysis_role, " BH", family_n)]
strip[, effect_label := paste0(sprintf("%+.2f", estimate), ifelse(significant_fdr, "*", ""))]
strip[, selection_rule := "EO D0 association and EO disease q<=0.1 with concordant signs; post-selection display"]
fwrite(
  candidates,
  file.path(v35_dir, "tables/species_EO_D0_candidates_v35.csv"),
  bom = TRUE
)
fwrite(
  strip,
  file.path(v35_dir, "tables/main_figure_support_strip_species_v35.csv"),
  bom = TRUE
)
crc_support <- merge(
  candidates[selected_for_strip == TRUE],
  all_results[
    analysis_role == "CRC_sensitivity" & age_group == "EO" & panel == "D0",
    .(
      feature_id,
      EOCRC_beta = estimate,
      EOCRC_ci_low = ci_low,
      EOCRC_ci_high = ci_high,
      EOCRC_p = p_value,
      EOCRC_q = q_BH,
      EOCRC_n = n
    )
  ],
  by = "feature_id"
)
crc_support[, EOCRC_concordant := sign(EOCRC_beta) == sign(association_beta)]
crc_support[, EOCRC_supported := EOCRC_concordant & EOCRC_q <= 0.1]
crc_support[, EOCRC_nominal := EOCRC_concordant & EOCRC_p < 0.05]
crc_support[, selection_rule := "Original EO disease/association double-FDR selection; CRC support is post-selection sensitivity"]
fwrite(
  crc_support,
  file.path(v35_dir, "tables/species_selected_CRC_support_v35.csv"),
  bom = TRUE
)
saveRDS(
  list(
    raw = raw,
    clr = clr,
    qc = qc,
    eligible = eligible,
    both_ids = both_ids,
    pseudocount = pseudocount,
    results = all_results,
    delta = delta,
    candidates = candidates,
    strip = strip,
    crc_support = crc_support,
    selected_ids = selected_ids,
    bundles = bundles
  ),
  file.path(v35_dir, "objects/species_associations_v35.rds")
)
