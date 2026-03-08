library(tidyverse)
library(metafor)

extract_cohort_stats <- function(fits_list, age_group_label) {
  map_dfr(names(fits_list), function(cohort) {
    df <- fits_list[[cohort]]
    req_cols <- c("metadata", "value", "feature", "coef", "stderr")
    if(!all(req_cols %in% colnames(df))) {
      warning(sprintf("Cohort %s: missing required cols; skipping.", cohort))
      return(NULL)
    }
    df %>%
      filter(metadata == "Group", value == "CRC") %>%
      transmute(
        feature = as.character(feature),
        coef = as.numeric(coef),
        stderr = as.numeric(stderr),
        Cohort = cohort,
        AgeGroup = age_group_label,
        Variance = stderr^2
      )
  })
}

calc_interaction_p_robust <- function(data, target_id, min_cohort_each_group = 2, min_total_rows = 5) {
  sub_df <- data %>% filter(Pathway_ID == target_id)
  n_EO <- sub_df %>% filter(AgeGroup == "EO") %>% pull(Cohort) %>% unique() %>% length()
  n_LO <- sub_df %>% filter(AgeGroup == "LO") %>% pull(Cohort) %>% unique() %>% length()
  n_total_rows <- nrow(sub_df)
  fail_return <- tibble(
    Pathway_ID = target_id,
    P_interaction = NA_real_,
    Coef_diff_EO_vs_LO = NA_real_,
    CI_lower = NA_real_,
    CI_upper = NA_real_,
    n_cohort_EO = n_EO,
    n_cohort_LO = n_LO,
    n_rows = n_total_rows,
    I2 = NA_real_
  )
  
  if (n_EO < min_cohort_each_group || n_LO < min_cohort_each_group || n_total_rows < min_total_rows) {
    return(fail_return)
  }

  sub_df <- sub_df %>% mutate(AgeGroup = factor(AgeGroup, levels = c("LO", "EO")))

  out <- tryCatch({
    res <- rma(yi = coef, vi = Variance, mods = ~ AgeGroup, data = sub_df, method = "REML")
    tab <- coef(summary(res))
    if (nrow(tab) < 2) {
      return(fail_return)
    }
    moderator_row <- tab[2, , drop = FALSE] 
    est <- as.numeric(moderator_row[ , "estimate"])
    pval <- as.numeric(moderator_row[ , "pval"])
    ci_lb <- as.numeric(moderator_row[ , "ci.lb"])
    ci_ub <- as.numeric(moderator_row[ , "ci.ub"])
    # 估计 I2（总体异质性）
    I2_val <- tryCatch({
      res$I2
    }, error = function(e) NA_real_)
    
    tibble(
      Pathway_ID = target_id,
      P_interaction = pval,
      Coef_diff_EO_vs_LO = est,
      CI_lower = ci_lb,
      CI_upper = ci_ub,
      n_cohort_EO = n_EO,
      n_cohort_LO = n_LO,
      n_rows = n_total_rows,
      I2 = I2_val
    )
  }, error = function(e) {
    warning(sprintf("rma failed for %s : %s", target_id, e$message))
    return(fail_return)
  })
  return(out)
}

species_raw_EO <- extract_cohort_stats(fit_meta.crc_EO$maaslin_fits, "EO")
species_raw_LO <- extract_cohort_stats(fit_meta.crc_LO$maaslin_fits, "LO")

combined_cohort_data <- bind_rows(species_raw_EO, species_raw_LO) %>%
  mutate(Pathway_ID = str_trim(str_extract(feature, "^[^:]+")))

target_pathways_interaction <- unique(combined_cohort_data$Pathway_ID)
cat("Total data points:", nrow(combined_cohort_data), "\n")

species_interaction_list <- map_dfr(target_pathways_interaction, 
                                    ~ calc_interaction_p_robust(combined_cohort_data, .x))

species_interaction_results <- bind_rows(species_interaction_list)

species_interaction_results_final <- species_interaction_results %>%
  filter(!is.na(P_interaction)) %>%
  arrange(P_interaction) %>%
  mutate(
    Nominal_Significant = ifelse(P_interaction < 0.05, "Yes (*)", "No"),
    Species = str_remove(str_extract(Pathway_ID, 's__[A-Za-z0-9_]+'), '^s__'),
    `95%_CI_Diff` = sprintf("[%.3f, %.3f]", CI_lower, CI_upper)
  ) %>% 
  select(Species, Pathway = Pathway_ID, Coef_diff_EO_vs_LO, `95%_CI_Diff`, 
         P_interaction, Nominal_Significant#, everything()
         )
