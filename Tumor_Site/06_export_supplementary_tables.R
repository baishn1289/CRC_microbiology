library(dplyr)
library(tidyr)

# Extend only the supplementary-table scope to both two-cohort interactions.
table_results <- all_results
interaction_scopes <- c(
  "ternary__interaction__joint__Rectum_vs_Left_Colon",
  "ternary__interaction__joint__Rectum_vs_Right_Colon"
)

for (sc in interaction_scopes) {
  rows <- which(table_results$scope == sc)
  a <- table_results[rows, ]
  valid <- a$status %in% c("ok", "descriptive_k2") & a$k == 2L &
    is.finite(a$estimate) & is.finite(a$se) & a$se > 0 &
    is.finite(a$p_value) & a$p_value >= 0 & a$p_value <= 1
  members <- rows[valid]
  table_results$q_comparison[members] <- p.adjust(table_results$p_value[members], "BH")
  table_results$n_comparison[rows] <- length(members)
  table_results$in_BH_family[members] <- table_results$reportable[members] <- TRUE
  table_results$policy_tier[rows] <- "two_cohort_interaction_exploratory"
  table_results$full_meta_QC_applied[rows] <- FALSE
  table_results$BH_family_policy[rows] <- "estimable_two_cohort_tests"
  table_results$analysis_status[rows] <- "excluded_exploratory_unestimable"
  table_results$analysis_status[members] <- "eligible_exploratory_estimable"
  table_results$P_for_BH[members] <- table_results$p_value[members]
  table_results$comparison_q_lt010[members] <- table_results$q_comparison[members] < .1
}

comparison_summary <- table_results %>%
  group_by(scope, policy_tier) %>%
  summarise(
    target_cohorts = first(target_cohorts),
    target_n_samples = first(target_n_samples),
    n_BH = sum(in_BH_family),
    n_hits = if (any(in_BH_family)) sum(q_comparison < .1, na.rm = TRUE) else NA_integer_,
    min_q = if (any(in_BH_family)) min(q_comparison, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  arrange(scope) %>%
  mutate(comparison_id = sprintf("TS%02d", row_number()))
scope_parts <- strsplit(comparison_summary$scope, "__", fixed = TRUE)
comparison_summary$classification <- vapply(scope_parts, `[`, character(1), 1)
comparison_summary$mode <- vapply(scope_parts, `[`, character(1), 2)
comparison_summary$age <- vapply(scope_parts, `[`, character(1), 3)
comparison_summary$contrast <- vapply(scope_parts, `[`, character(1), 4)
site_names <- c(
  Colon = "colon",
  Rectum = "rectum",
  Left_Colon = "left colon",
  Right_Colon = "right colon",
  left_sided = "left side",
  right_sided = "right side"
)

for (i in seq_len(nrow(comparison_summary))) {
  anatomical_sites <- unname(site_names[strsplit(
    comparison_summary$contrast[i],
    "_vs_",
    fixed = TRUE
  )[[1]]])
  anatomical_sites[1] <- paste0(
    toupper(substr(anatomical_sites[1], 1, 1)),
    substring(anatomical_sites[1], 2)
  )
  contrast <- paste(anatomical_sites, collapse = " - ")
  if (comparison_summary$mode[i] == "age") contrast <- paste0(contrast, " cases: EOCRC - LOCRC")
  if (comparison_summary$mode[i] == "interaction")
    contrast <- paste0("(", contrast, ") in EOCRC - (", contrast, ") in LOCRC")
  comparison_summary$contrast[i] <- contrast
}

comparison_summary$classification <- c(
  binary = "Binary",
  ternary = "Ternary",
  sidedness = "Sidedness"
)[comparison_summary$classification]
comparison_summary$evidence_tier <- ifelse(
  comparison_summary$policy_tier == "multi_ge3_fullQC50",
  "Multi-cohort full QC",
  "Two-cohort exploratory"
)
candidates <- table_results %>%
  filter(reportable, in_BH_family, q_comparison < .1) %>%
  left_join(
    select(comparison_summary, scope, comparison_id, table_evidence_tier = evidence_tier),
    by = "scope"
  ) %>%
  arrange(comparison_id, q_comparison)
overview_table <- comparison_summary %>%
  transmute(
    `Comparison ID` = comparison_id,
    `Analysis type` = ifelse(
      mode == "site",
      paste0(age, "CRC site association"),
      ifelse(mode == "age", "Within-site age association", "Age-by-site interaction")
    ),
    Classification = classification,
    Contrast = contrast,
    `Evidence tier` = evidence_tier,
    `Eligible cohorts` = target_cohorts,
    `Eligible cases` = target_n_samples,
    `BH tests` = n_BH,
    `Candidates (q < 0.10)` = n_hits,
    `Minimum comparison-wise q` = min_q,
    `Analysis status` = "Included"
  )
candidate_table <- candidates %>%
  left_join(
    select(overview_table, `Comparison ID`, `Analysis type`, Classification, Contrast),
    by = c("comparison_id" = "Comparison ID")
  ) %>%
  transmute(
    `Comparison ID` = comparison_id,
    `Analysis type`,
    Classification,
    Contrast,
    Taxon = taxon,
    `Evidence tier` = table_evidence_tier,
    `Contributing cohorts` = k,
    `Analyzed cases` = n_samples,
    `Pooled beta` = estimate,
    `95% CI` = sprintf("[%.3f, %.3f]", ci_low, ci_high),
    `Raw P` = p_value,
    `Comparison-wise q` = q_comparison,
    `I-squared (%)` = I2,
    `Maximum cohort weight (%)` = max_weight_percent,
    `Direction agreement (%)` = 100 * qc_consistency_score,
    `BH tests` = n_comparison
  )
write.csv(
  overview_table,
  file.path(output_dir, "Supplementary_Table_Tumor_Site_Analysis_Overview.csv"),
  row.names = FALSE,
  na = "NA"
)
write.csv(
  candidate_table,
  file.path(output_dir, "Supplementary_Table_Tumor_Site_Candidate_Associations.csv"),
  row.names = FALSE,
  na = "NA"
)
write.csv(
  filter(table_results, scope %in% interaction_scopes, in_BH_family),
  file.path(output_dir, "Two_cohort_interaction_results.csv"),
  row.names = FALSE,
  na = "NA"
)
