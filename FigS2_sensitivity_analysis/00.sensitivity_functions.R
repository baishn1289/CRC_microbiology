library(tidyverse)
library(MMUPHin)
library(metafor)

# Main-analysis functions and objects remain in the shared R session.
sensitivity_qc <- function(fit) {
  z <- process_meta_with_rescue(
    fit,
    i2_tier1 = 50,
    i2_tier2 = 80,
    weight_limit = 40,
    base_consistency = 0.6,
    rescue_consistency = 0.75
  )
  z$qval.fdr <- NA_real_
  keep <- which(z$pass_qc & is.finite(z$pval))
  z$qval.fdr[keep] <- p.adjust(z$pval[keep], "BH")
  z %>%
    mutate(feature = Pathway, significant = pass_qc & !is.na(qval.fdr) & qval.fdr < 0.1)
}

# The same model and masking rules apply to each sensitivity refit.
sensitivity_fit <- function(metadata, abundance, output_path) {
  md <- metadata %>%
    mutate(
      Group = factor(as.character(Group), levels = c("Control", "CRC")),
      Cohort = factor(Cohort),
      Age = as.numeric(as.character(Age)),
      Sex = factor(Sex),
      BMI = as.numeric(as.character(BMI))
    ) %>%
    column_to_rownames("Sample_ID")
  mat <- abundance %>%
    column_to_rownames("clade_name") %>%
    as.matrix()
  mat <- apply_step2_masking(mat, md, "Cohort")
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  MMUPHin::lm_meta(
    feature_abd = as.data.frame(mat[, rownames(md), drop = FALSE]),
    exposure = "Group",
    batch = "Cohort",
    covariates = c("Age", "Sex", "BMI"),
    data = md,
    control = list(
      rma_method = "HS",
      transform = "LOG",
      normalization = "NONE",
      min_prevalence = 0,
      min_abundance = 0.001,
      forest_plot = NULL,
      verbose = FALSE,
      output = output_path
    )
  )
}

# Retain age-by-cohort strata with both disease groups represented.
sensitivity_metadata <- function(metadata) {
  metadata <- metadata %>%
    filter(Group %in% c("Control", "CRC"), Age_class %in% c("EO", "LO"))
  eligible <- metadata %>%
    count(Cohort, Age_class, Group) %>%
    group_by(Cohort, Age_class) %>%
    summarise(eligible = all(c("Control", "CRC") %in% Group) && all(n >= 3), .groups = "drop") %>%
    filter(eligible)
  semi_join(metadata, eligible, by = c("Cohort", "Age_class"))
}
