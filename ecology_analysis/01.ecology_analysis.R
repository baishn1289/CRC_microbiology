library(tidyverse)
library(vegan)
library(MMUPHin)
library(broom)
library(metafor)
library(permute)

# Continue in the R session used for the main microbiome analysis.
code_dir <- file.path(
  base_symbol,
  "downstream_analysis/src/github_code_availability/submission_version_261102/ecology_analysis"
)
output_dir <- file.path(
  base_symbol,
  "downstream_analysis/output/submission_version_261102/ecology_analysis"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Retain cohorts with at least three participants per disease group.
world_meta <- list()

for (age in c("EO", "LO")) {
  meta_age <- metadata_allcohorts %>%
    filter(Age_class == age)
  eligible <- meta_age %>%
    count(Cohort, Group) %>%
    group_by(Cohort) %>%
    summarise(keep = all(c("CRC", "Control") %in% Group) && all(n >= 3), .groups = "drop") %>%
    filter(keep) %>%
    pull(Cohort)
  world_meta[[age]] <- meta_age %>%
    filter(Cohort %in% eligible) %>%
    mutate(Group = factor(Group, levels = c("Control", "CRC")))
}

# Reuse the main consensus function with the executed ecological thresholds.
abundance <- run_consensus_EO_LO_combine(
  meta_list = meta_list_metadata_intersect_samples,
  meta_EO = world_meta$EO,
  meta_LO = world_meta$LO,
  min_cohorts = 3,
  MIN_STEP1_COUNT = 3,
  min_prevalence = 0.05,
  abund_threshold = 0.001
)
primary <- list()

for (age in c("EO", "LO")) {
  feat <- abundance[[age]] %>%
    column_to_rownames("clade_name") %>%
    as.matrix()
  mask_meta <- world_meta[[age]] %>%
    column_to_rownames("Sample_ID")
  feat <- apply_step2_masking(feat, mask_meta, "Cohort")
  result <- run_beta_meta_final(
    feat,
    world_meta[[age]],
    file.path(output_dir, age),
    keep_details = TRUE
  )
  primary[[age]] <- result
}

beta_stats <- bind_rows(
  primary$EO$stats %>%
    mutate(Age_class = "EO"),
  primary$LO$stats %>%
    mutate(Age_class = "LO")
)

# Use the global PERMANOVA samples and distances for spatial-median PERMDISP.
dispersion <- list()
dispersion_stats <- list()

for (age in c("EO", "LO")) {
  m <- primary[[age]]$metadata
  d <- primary[[age]]$distance
  stopifnot(identical(attr(d, "Labels"), rownames(m)))
  bd <- vegan::betadisper(
    d,
    m$Group,
    type = "median",
    bias.adjust = FALSE
  )
  set.seed(123)
  perm <- vegan::permutest(
    bd,
    permutations = permute::how(nperm = 999, blocks = m$Cohort)
  )
  distances <- data.frame(
    Sample_ID = rownames(m),
    Age_class = age,
    Cohort = m$Cohort,
    Group = m$Group,
    Distance = as.numeric(bd$distances)
  )
  means <- tapply(distances$Distance, distances$Group, mean)
  dispersion[[age]] <- distances
  dispersion_stats[[age]] <- data.frame(
    Age_class = age,
    n_control = sum(m$Group == "Control"),
    n_crc = sum(m$Group == "CRC"),
    mean_control = unname(means["Control"]),
    mean_crc = unname(means["CRC"]),
    mean_diff = unname(means["CRC"] - means["Control"]),
    permdisp_p = perm$tab["Groups", "Pr(>F)"],
    center = "spatial_median",
    permutations = 999,
    restriction = "within_Cohort",
    distance_source = "same_as_primary_global_PERMANOVA"
  )
}

dispersion <- bind_rows(dispersion)
dispersion_stats <- bind_rows(dispersion_stats)

# Compute species-level Shannon and combine adjusted cohort effects by REML.
alpha_all_species <- calc_alpha_from_list(meta_list_metadata_intersect_samples)
alpha_fit <- list()

for (age in c("EO", "LO")) {
  alpha_data <- alpha_all_species %>%
    select(-Cohort) %>%
    inner_join(world_meta[[age]], by = "Sample_ID")
  alpha_fit[[age]] <- run_alpha_meta_standard(alpha_data, metric = "Shannon")
}

alpha_stats <- bind_rows(lapply(
  c("EO", "LO"),
  function(age)
    bind_rows(alpha_fit[[age]]$cohort_res, alpha_fit[[age]]$meta_res) %>%
      mutate(Age_class = age)
))

ecology <- list(
  beta = beta_stats,
  primary = primary,
  dispersion = dispersion,
  dispersion_stats = dispersion_stats,
  alpha = alpha_fit
)

# Export model-specific counts and full-precision statistical results.
beta_rows <- list()

for (age in c("EO", "LO")) {
  p <- primary[[age]]

  for (i in seq_len(nrow(p$stats))) {
    row <- p$stats[i, ]
    global <- row$Type == "Summary"
    test <- if (global) p$tests$global else p$tests$cohort[[row$Cohort]]
    beta_rows[[length(beta_rows) + 1]] <- data.frame(
      Age_stratum = age,
      Estimate_type = if (global) "Global sample-pooled" else "Cohort",
      Cohort = if (global) "Global" else row$Cohort,
      N_control = row$n_control,
      N_CRC = row$n_crc,
      N_total = row$n,
      Group_df = test["Group", "Df"],
      Residual_df = test["Residual", "Df"],
      Group_R2 = row$coef,
      Pseudo_F = test["Group", "F"],
      P_nominal = row$pval,
      Covariates = gsub(";", ", ", row$covariates, fixed = TRUE),
      Input = if (global) "TSS + MMUPHin" else "TSS",
      Permutations = nrow(attr(test, "F.perm"))
    )
  }
}

b <- bind_rows(beta_rows)
a <- alpha_stats %>%
  transmute(
    Age_stratum = Age_class,
    Estimate_type = if_else(Type == "Summary", "REML summary", "Cohort"),
    Cohort = if_else(Type == "Summary", "Overall", Cohort),
    N_total = n,
    Mean_difference = coef,
    SE = stderr,
    CI95_lower = lower,
    CI95_upper = upper,
    P_nominal = pval,
    Residual_df = residual_df
  )
c_table <- dispersion_stats %>%
  transmute(
    Age_stratum = Age_class,
    N_control = n_control,
    N_CRC = n_crc,
    N_total = n_control + n_crc,
    Mean_control = mean_control,
    Mean_CRC = mean_crc,
    CRC_minus_control = mean_diff,
    P_PERMDISP = permdisp_p,
    Permutations = permutations
  )

tables <- list(
  B_PERMANOVA = b,
  A_Shannon = a,
  C_PERMDISP = c_table
)

for (name in names(tables)) write.csv(
  tables[[name]],
  file.path(output_dir, paste0(name, ".csv")),
  row.names = FALSE,
  na = ""
)
write.csv(
  alpha_all_species,
  file.path(output_dir, "alpha_species_per_sample.csv"),
  row.names = FALSE
)
write.csv(dispersion, file.path(output_dir, "PERMDISP_distances.csv"), row.names = FALSE)
