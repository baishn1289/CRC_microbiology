library(dplyr)
library(tidyr)
library(metafor)
library(Maaslin2)
library(ggplot2)
library(grid)

code_dir <- file.path(
  base_symbol,
  "downstream_analysis/src/github_code_availability/submission_version_261102/Tumor_Site"
)
output_dir <- file.path(
  base_symbol,
  "downstream_analysis/output/submission_version_261102/Tumor_Site"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(code_dir, "00_helpers.R"))

# Use the original prepared matrix and clinical eligibility, never old fitted results.
if (!exists("abundance_full") ||
  !exists("metadata_all") ||
  !"complete_primary" %in% names(metadata_all)) {
  prepared_input <- new.env()
  load(
    file.path(base_symbol, "Tumor_Site/v5/inputs/V4_prepared.Rdata"),
    envir = prepared_input
  )
  abundance_full <- prepared_input$abundance_full
  metadata_all <- prepared_input$metadata_all
}

stopifnot(
  !anyDuplicated(metadata_all$Sample_ID),
  !anyDuplicated(colnames(abundance_full)),
  setequal(metadata_all$Sample_ID, colnames(abundance_full)),
  nrow(abundance_full) == 3603L
)
abundance_full <- abundance_full[, metadata_all$Sample_ID, drop = FALSE]
stopifnot(max(abs(colSums(abundance_full) - 100)) < 1e-10)
metadata_all$Age_class <- factor(as.character(metadata_all$Age_class), levels = c("LO", "EO"))
site_metadata <- read.csv(
  file.path(base_symbol, "Tumor_Site/metadata_allcohorts_manual_sort_261016.csv"),
  colClasses = "character",
  na.strings = character(),
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
site_index <- match(metadata_all$Sample_ID, site_metadata$Sample_ID)
stopifnot(!anyDuplicated(site_metadata$Sample_ID), !anyNA(site_index))

for (field in paste0("Tumor_Site_", c("binary", "ternary", "sidedness"))) {
  metadata_all[[field]] <- site_metadata[[field]][site_index]
}

levels_site <- list(
  binary = c("Colon", "Rectum"),
  ternary = c("Right_Colon", "Left_Colon", "Rectum"),
  sidedness = c("right_sided", "left_sided")
)
prepared_sets <- build_site_sets(metadata_all, levels_site)
sets <- lapply(
  prepared_sets$sets,
  prepare_site_set,
  abundance = abundance_full,
  metadata = metadata_all
)
target_coverage <- bind_rows(lapply(
  sets,
  function(s) data.frame(
    scope = s$scope,
    target_cohorts = n_distinct(s$metadata$Cohort),
    target_n_samples = nrow(s$metadata)
  )
))
