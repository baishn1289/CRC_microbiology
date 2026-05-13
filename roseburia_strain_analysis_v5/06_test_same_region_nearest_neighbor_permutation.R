#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

env_or_default_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed)) stop(sprintf("%s must be an integer.", name), call. = FALSE)
  parsed
}

current_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hits <- grep(paste0("^", file_arg), args, value = TRUE)
  if (length(hits) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", hits[[1]]), winslash = "/", mustWork = TRUE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

source_with_env_flag <- function(script, env, flag_name, flag_value) {
  old_value <- Sys.getenv(flag_name, unset = NA_character_)
  on.exit({
    if (is.na(old_value)) {
      Sys.unsetenv(flag_name)
    } else {
      do.call(Sys.setenv, setNames(list(old_value), flag_name))
    }
  }, add = TRUE)
  do.call(Sys.setenv, setNames(list(flag_value), flag_name))
  sys.source(script, envir = env)
}

same_region_permutation_test <- function(species_data, n_perm, seed) {
  set.seed(seed)
  nn <- species_data$nn_tbl
  observed <- mean(nn$same_region, na.rm = TRUE)
  regions <- as.character(nn$region)
  nearest <- match(nn$nearest_tip_label, nn$tree_tip_label)

  perm_rates <- replicate(n_perm, {
    perm_region <- sample(regions, length(regions), replace = FALSE)
    mean(perm_region == perm_region[nearest], na.rm = TRUE)
  })

  tibble(
    species_id = species_data$spec$species_id,
    species_label = species_data$spec$species_label,
    n_samples = nrow(nn),
    n_permutations = n_perm,
    observed_same_region = observed,
    expected_random = species_data$nn_summary$expected_random,
    permuted_mean = mean(perm_rates),
    permuted_sd = sd(perm_rates),
    permuted_max = max(perm_rates),
    empirical_p_greater = (sum(perm_rates >= observed) + 1) / (n_perm + 1)
  )
}

script_dir <- normalizePath(env_or_default("STRAINPHLAN_V5_SCRIPT_DIR", current_script_dir()), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(env_or_default("PROJECT_ROOT", dirname(script_dir)), winslash = "/", mustWork = FALSE)
result_root <- env_or_default("STRAINPHLAN_V5_RESULT_ROOT", file.path(project_root, "results/strain_analysis_v5"))
output_dir <- env_or_default("STRAINPHLAN_ROSEBURIA_REGION_STRUCTURE_OUTPUT_DIR", file.path(result_root, "roseburia_region_structure_paper"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

region_script <- file.path(script_dir, "03_plot_roseburia_region_structure.R")
region_env <- new.env(parent = globalenv())
source_with_env_flag(region_script, region_env, "STRAINPHLAN_ROSEBURIA_REGION_RENDER_FILES", "false")

n_perm <- env_or_default_int("STRAINPHLAN_NN_N_PERM", 10000)
seed <- env_or_default_int("STRAINPHLAN_NN_SEED", 20260508)

permutation_summary <- bind_rows(lapply(region_env$species_data, same_region_permutation_test, n_perm = n_perm, seed = seed))
write_tsv(permutation_summary, file.path(output_dir, "roseburia_region_nearest_neighbor_permutation_test.tsv"))
print(permutation_summary)
