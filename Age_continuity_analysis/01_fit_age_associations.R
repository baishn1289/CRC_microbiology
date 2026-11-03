library(dplyr)
library(readr)
library(metafor)

code_dir <- file.path(
  base_symbol,
  "downstream_analysis/src/github_code_availability/submission_version_261102/Age_continuity_analysis"
)
output_dir <- file.path(
  base_symbol,
  "downstream_analysis/output/submission_version_261102/Age_continuity_analysis"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(code_dir, "00_helpers.R"))

# Construct the original percentage-abundance inputs in the current R session.
meta <- as.data.frame(metadata_allcohorts)
tables <- meta_list_metadata_intersect_samples
stopifnot(!anyDuplicated(meta$Sample_ID), !anyDuplicated(paste(meta$Cohort, meta$Patient_ID)))
meta$Age <- as.numeric(meta$Age)
meta$BMI <- as.numeric(meta$BMI)
meta$Sex <- factor(as.character(meta$Sex), levels = c("female", "male"))
meta$D <- as.integer(meta$Group == "CRC")
meta$Age10 <- (meta$Age - 50) / 10
meta$BMI_c <- meta$BMI - 25
cohort_data <- list()

for (key in names(tables)) {
  tab <- as.data.frame(tables[[key]])
  ids <- intersect(colnames(tab)[-1], meta$Sample_ID)
  mm <- meta[match(ids, meta$Sample_ID), , drop = FALSE]
  stopifnot(length(unique(mm$Cohort)) == 1L)
  species_rows <- grepl("\\|s__[^|]+$", tab$clade_name)
  xx <- as.matrix(tab[species_rows, ids, drop = FALSE])
  storage.mode(xx) <- "double"
  rownames(xx) <- tab$clade_name[species_rows]
  stopifnot(!anyDuplicated(rownames(xx)), all(is.finite(xx)), all(xx >= 0))
  mm$complete <- complete.cases(mm[, c("Age", "Sex", "BMI")]) &
    mm$Age >= 18 & mm$Age <= 100 & mm$Group %in% c("Control", "CRC")
  cohort_data[[unique(mm$Cohort)]] <- list(metadata = mm, raw = xx)
}

targets <- c("Common_age", "CRC_age", "Control_age", "Interaction_v7")
all_features <- sort(unique(unlist(lapply(cohort_data, function(dat) rownames(dat$raw)))))

# Define target-specific technical families before fitting.
designs <- list()
votes <- list()

for (cohort in names(cohort_data)) {
  dat <- cohort_data[[cohort]]

  for (target in targets[1:3]) {
    key <- paste(cohort, target, sep = "::")
    obj <- prepare_design(dat, target)
    designs[[key]] <- obj
    xx <- dat$raw[, obj$keep, drop = FALSE]
    detected <- rowSums(xx > .001)
    votes[[key]] <- data.frame(
      feature = rownames(xx),
      Cohort = cohort,
      Target = target,
      n = ncol(xx),
      detected = detected,
      prevalence = detected / ncol(xx),
      cohort_eligible = obj$audit$eligible,
      nominated = obj$audit$eligible & detected >= 3 & detected / ncol(xx) >= .05
    )
  }
}

votes <- bind_rows(votes)
nomination <- votes %>%
  group_by(Target, feature) %>%
  summarise(nominated_cohorts = sum(nominated), .groups = "drop")
screens <- lapply(
  targets[1:3],
  function(target) sort(nomination$feature[nomination$Target == target &
    nomination$nominated_cohorts >= 3])
)
names(screens) <- targets[1:3]
screens$Interaction_v7 <- sort(union(screens$CRC_age, screens$Control_age))
tested_features <- sort(unique(unlist(screens)))
pseudocount <- setNames(rep(Inf, length(tested_features)), tested_features)

for (dat in cohort_data) {
  xx <- dat$raw[, dat$metadata$complete, drop = FALSE]

  for (feature in intersect(tested_features, rownames(xx))) {
    positive <- xx[feature, ][xx[feature, ] > 0]
    if (length(positive)) pseudocount[feature] <- min(pseudocount[feature], min(positive) / 2)
  }
}

stopifnot(all(is.finite(pseudocount)), all(pseudocount > 0))
local_rows <- list()

for (cohort in names(cohort_data)) {
  for (target in targets[1:3]) {
    key <- paste(cohort, target, sep = "::")
    obj <- fit_target(cohort_data[[cohort]], designs[[key]], target, screens[[target]], pseudocount)
    local_rows[[key]] <- obj$local
  }
}

local <- bind_rows(local_rows)
stopifnot(all(local$residual_df >= 10))
write_csv(local, file.path(output_dir, "local_age_slopes.csv"))

pair <- inner_join(
  filter(local, Target == "CRC_age"),
  filter(local, Target == "Control_age"),
  by = c("feature", "Cohort"),
  suffix = c("_CRC", "_Control")
)
interaction <- transmute(
  pair,
  feature,
  Cohort,
  Target = "Interaction_v7",
  beta = beta_CRC - beta_Control,
  variance = variance_CRC + variance_Control,
  n = n_CRC + n_Control,
  n_control = n_Control,
  n_crc = n_CRC,
  residual_df = pmin(residual_df_CRC, residual_df_Control),
  residual_df_crc = residual_df_CRC,
  residual_df_control = residual_df_Control,
  formula = "within_cohort_CRC_minus_Control; independent_HC3_variances_added"
)
write_csv(interaction, file.path(output_dir, "local_interaction_differences.csv"))
all_local <- bind_rows(local, interaction)
pooled <- list()
by_target_feature <- split(all_local, paste(all_local$Target, all_local$feature, sep = "::"))

for (target in targets) {
  for (feature in screens[[target]]) {
    dd <- by_target_feature[[paste(target, feature, sep = "::")]]
    if (is.null(dd)) dd <- all_local[FALSE, ]
    obj <- pool_target(dd, feature, target)
    pooled[[length(pooled) + 1L]] <- obj$row
  }
}

result <- bind_rows(pooled)
result$q_after_QC <- NA_real_
result$qc_family_n <- result$screen_family_n <- NA_integer_

for (target in targets) {
  idx <- which(result$Target == target)
  passed <- idx[result$qc_pass[idx]]
  result$qc_family_n[idx] <- length(passed)
  result$screen_family_n[idx] <- length(idx)
  result$q_after_QC[passed] <- p.adjust(result$p[passed], "BH")
}

result$candidate_q010 <- result$qc_pass & !is.na(result$q_after_QC) & result$q_after_QC < .1
write_csv(result, file.path(output_dir, "all_targets_all_species.csv"))

for (target in targets) {
  rr <- result[result$Target == target, ]
  write_csv(rr, file.path(output_dir, paste0(target, "_all_species.csv")))
  write_csv(filter(rr, qc_pass), file.path(output_dir, paste0(target, "_QC_passed.csv")))
  write_csv(
    filter(rr, candidate_q010),
    file.path(output_dir, paste0(target, "_candidates_q010.csv"))
  )
}

summary <- result %>%
  group_by(Target) %>%
  summarise(
    screen_n = n(),
    estimable_n = sum(status == "ok"),
    QC_passed_n = sum(qc_pass),
    nominal_P005_QC = sum(qc_pass & p < .05, na.rm = TRUE),
    q010_QC = sum(candidate_q010),
    q005_QC = sum(q_after_QC < .05, na.rm = TRUE),
    .groups = "drop"
  )
summary <- summary[match(targets, summary$Target), ]
write_csv(summary, file.path(output_dir, "discovery_summary.csv"))
common_candidates <- result$feature[result$Target == "Common_age" & result$candidate_q010]
write_csv(
  filter(result, feature %in% common_candidates),
  file.path(output_dir, "Common_candidates_all_four_targets.csv")
)
write_csv(
  filter(result, Target == "Interaction_v7", feature %in% common_candidates, qc_pass, p < .05),
  file.path(output_dir, "Common_candidates_nominal_v7_interaction.csv")
)


print(summary, width = Inf)
