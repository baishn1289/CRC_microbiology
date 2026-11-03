library(dplyr)
library(tibble)

# Apply the final 50% meta-QC, then comparison-wise BH in estimable families.
stopifnot(is.function(process_meta_with_rescue))
duplicate_scopes <- c(
  "ternary__age__joint__Rectum" = "binary__age__joint__Rectum",
  "sidedness__age__joint__right_sided" = "ternary__age__joint__Right_Colon"
)
meta_qc <- direct_meta %>%
  mutate(unique_scope = ifelse(
    scope %in% names(duplicate_scopes),
    unname(duplicate_scopes[scope]),
    scope
  ))
qc50 <- apply_meta_qc(meta_qc, direct_local, process_meta_with_rescue, weight_limit = 50)
d <- meta_qc %>%
  left_join(
    qc50 %>%
      select(
        scope,
        feature,
        qc_pass_50 = qc_pass,
        qc_failure_reasons_50 = qc_failure_reasons,
        qc_nonzero_cohorts,
        qc_consistency_score
      ),
    by = c("scope", "feature")
  ) %>%
  left_join(target_coverage, by = "scope")
all_results <- estimable_bh(d)
all_results$full_meta_QC_applied <- all_results$policy_tier == "multi_ge3_fullQC50"
all_results$analysis_status <- ifelse(
  all_results$policy_tier == "excluded_two_cohort_interaction",
  "outside_policy",
  ifelse(all_results$reportable, "eligible", "excluded")
)
all_results$BH_family_policy <- ifelse(
  all_results$policy_tier == "multi_ge3_fullQC50",
  "full_QC50_passed",
  ifelse(
    all_results$policy_tier == "EO_two_cohort_exploratory",
    "estimable_two_cohort_tests",
    "outside_policy"
  )
)
all_results$comparison_q_lt010 <- ifelse(
  all_results$in_BH_family,
  all_results$q_comparison < .1,
  NA
)
all_results <- all_results[, setdiff(names(all_results), c("unique_scope"))]
write.csv(
  all_results,
  file.path(output_dir, "all_comparisons.csv"),
  row.names = FALSE,
  na = "NA"
)
