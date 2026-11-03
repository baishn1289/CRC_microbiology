library(tidyverse)
library(MMUPHin)
library(metafor)

read_tsv <- readr::read_tsv
write_tsv <- readr::write_tsv

# Select discovery taxa using the shared prevalence and closure rules.
run_discovery_consensus <- function(
  meta_list,
  discovery_eo,
  discovery_lo,
  min_prevalence = 0.05,
  min_cohorts = 3,
  min_detected_samples = 3,
  abundance_threshold = 0.001
) {
  out <- run_consensus_EO_LO_combine(
    meta_list,
    discovery_eo,
    discovery_lo,
    min_prevalence = min_prevalence,
    min_cohorts = min_cohorts,
    MIN_STEP1_COUNT = min_detected_samples,
    abund_threshold = abundance_threshold
  )
  votes <- lapply(
    list(discovery_eo, discovery_lo),
    function(md) {
      nominated <- unlist(
        lapply(
          meta_list,
          function(tab) {
            z <- tab %>%
              filter(str_detect(clade_name, "\\|s__[^|]+$")) %>%
              column_to_rownames("clade_name") %>%
              select(any_of(md$Sample_ID)) %>%
              mutate(across(everything(), as.numeric))
            if (!ncol(z))
              return(character())
            n <- rowSums(z > abundance_threshold)
            rownames(z)[n >= min_detected_samples & n / ncol(z) >= min_prevalence]
          }
        ),
        use.names = FALSE
      )
      counts <- table(nominated)
      names(counts)[counts >= min_cohorts]
    }
  )

  list(
    merged = out$merged,
    keep_eo = votes[[1]],
    keep_lo = votes[[2]],
    keep_final = unique(c(votes[[1]], votes[[2]]))
  )
}

# Reuse main meta-QC and preserve discovery-specific boundary handling.
process_discovery_qc <- function(fit_object) {
  z <- process_meta_with_rescue(
    fit_object,
    i2_tier1 = 50,
    i2_tier2 = 80,
    weight_limit = 40,
    base_consistency = 0.6,
    rescue_consistency = 0.75
  )
  missing_support <- is.na(z$n_total_cohorts)
  z$n_total_cohorts[missing_support] <- 0L
  z$n_same_dir[missing_support] <- 0L
  z$n_opp_dir[missing_support] <- 0L
  z$pass_min_cohort[missing_support] <- FALSE
  z$pass_qc[missing_support] <- FALSE
  z$qc_status[missing_support & z$pass_weight] <- "Fail_Too_Few_Cohorts"
  z$is_rescue_zone <- z$is_in_rescue_zone
  z$qval_fdr <- NA_real_
  keep <- which(z$pass_qc & is.finite(z$pval))
  z$qval_fdr[keep] <- p.adjust(z$pval[keep], "BH")
  z
}

# Require at least three participants in both disease groups.
get_eligible_cohorts <- function(metadata) {
  cohort_counts <- metadata %>%
    count(Cohort, Age_class, Group) %>%
    pivot_wider(
      names_from = c(Age_class, Group),
      values_from = n,
      values_fill = 0,
      names_glue = "{Age_class}{Group}"
    )

  for (column_name in c("EOControl", "EOCRC", "LOControl", "LOCRC")) {
    if (!column_name %in% colnames(cohort_counts)) {
      cohort_counts[[column_name]] <- 0L
    }
  }

  list(
    eo = cohort_counts %>%
      filter(EOControl >= 3, EOCRC >= 3) %>%
      pull(Cohort) %>%
      as.character(),
    lo = cohort_counts %>%
      filter(LOControl >= 3, LOCRC >= 3) %>%
      pull(Cohort) %>%
      as.character(),
    counts = cohort_counts
  )
}

# Set the reference levels and sample alignment used by the models.
format_model_metadata <- function(metadata) {
  formatted <- metadata %>%
    filter(Group %in% c("CRC", "Control")) %>%
    mutate(
      Group = factor(as.character(Group), levels = c("Control", "CRC")),
      Cohort = as.factor(Cohort),
      Age = suppressWarnings(as.numeric(as.character(Age))),
      Sex = as.factor(Sex),
      BMI = suppressWarnings(as.numeric(as.character(BMI)))
    )
  if (anyDuplicated(formatted$Sample_ID)) {
    stop("模型 metadata 中存在重复 Sample_ID。")
  }

  formatted %>%
    column_to_rownames("Sample_ID")
}

# Fit an age-stratum meta-analysis with the original main-model settings.
fit_mmuphin_meta <- function(metadata, abundance_matrix, output_label) {
  model_metadata <- format_model_metadata(metadata)
  sample_ids <- rownames(model_metadata)
  missing_samples <- setdiff(sample_ids, colnames(abundance_matrix))
  if (length(missing_samples) > 0) {
    stop("丰度矩阵缺少样本：", paste(missing_samples, collapse = ", "))
  }

  model_abundance <- abundance_matrix[, sample_ids, drop = FALSE] %>%
    apply_step2_masking(meta_data = model_metadata) %>%
    as.data.frame()
  output_dir <- file.path(mmuphin_output_dir, output_label)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  suppressMessages(MMUPHin::lm_meta(
    feature_abd = model_abundance,
    exposure = "Group",
    batch = "Cohort",
    covariates = c("Age", "Sex", "BMI"),
    data = model_metadata,
    control = list(
      rma_method = "HS",
      transform = "LOG",
      normalization = "NONE",
      min_prevalence = 0,
      min_abundance = 0.001,
      forest_plot = NULL,
      verbose = FALSE,
      output = output_dir
    )
  ))
}

# Freeze shared significant taxa using discovery results only.
select_discovery_species <- function(processed_eo, processed_lo) {
  significant_eo <- processed_eo %>%
    filter(pass_qc, qval_fdr < 0.1) %>%
    transmute(
      feature = feature,
      Species_EO = str_remove(
        str_extract(feature, "s__[A-Za-z0-9_]+"),
        "^s__"
      ),
      discovery_beta_EO = coef,
      discovery_se_EO = stderr,
      discovery_q_EO = qval_fdr,
      discovery_k_EO = k,
      discovery_I2_EO = I2,
      discovery_max_weight_EO = max_weight,
      discovery_qc_EO = qc_status
    )
  significant_lo <- processed_lo %>%
    filter(pass_qc, qval_fdr < 0.1) %>%
    transmute(
      feature = feature,
      Species_LO = str_remove(
        str_extract(feature, "s__[A-Za-z0-9_]+"),
        "^s__"
      ),
      discovery_beta_LO = coef,
      discovery_se_LO = stderr,
      discovery_q_LO = qval_fdr,
      discovery_k_LO = k,
      discovery_I2_LO = I2,
      discovery_max_weight_LO = max_weight,
      discovery_qc_LO = qc_status
    )
  inner_join(significant_eo, significant_lo, by = "feature") %>%
    mutate(
      Species = coalesce(Species_EO, Species_LO),
      discovery_direction = sign(discovery_beta_EO),
      discovery_direction_label = if_else(
        discovery_direction >
          0,
        "CRC-enriched",
        "CRC-depleted"
      ),
      discovery_direction_concordant = sign(discovery_beta_EO) == sign(discovery_beta_LO)
    ) %>%
    filter(discovery_direction_concordant) %>%
    select(-Species_EO, -Species_LO) %>%
    arrange(Species)
}

# Construct the three nested discovery selection tiers.
build_discovery_tier_species <- function(processed_eo, processed_lo) {
  qc_eo <- processed_eo %>%
    filter(pass_qc) %>%
    transmute(
      feature,
      Species_EO = str_remove(
        str_extract(feature, "s__[A-Za-z0-9_]+"),
        "^s__"
      ),
      discovery_beta_EO = coef,
      discovery_se_EO = stderr,
      discovery_q_EO = qval_fdr,
      discovery_k_EO = k,
      discovery_I2_EO = I2,
      discovery_max_weight_EO = max_weight,
      discovery_qc_EO = qc_status
    )
  qc_lo <- processed_lo %>%
    filter(pass_qc) %>%
    transmute(
      feature,
      Species_LO = str_remove(
        str_extract(feature, "s__[A-Za-z0-9_]+"),
        "^s__"
      ),
      discovery_beta_LO = coef,
      discovery_se_LO = stderr,
      discovery_q_LO = qval_fdr,
      discovery_k_LO = k,
      discovery_I2_LO = I2,
      discovery_max_weight_LO = max_weight,
      discovery_qc_LO = qc_status
    )
  inner_join(qc_eo, qc_lo, by = "feature") %>%
    mutate(
      Species = coalesce(Species_EO, Species_LO),
      discovery_direction = sign(discovery_beta_EO),
      discovery_direction_label = if_else(discovery_direction > 0, "CRC-enriched", "CRC-depleted"),
      discovery_direction_concordant = sign(discovery_beta_EO) ==
        sign(discovery_beta_LO),
      tier_all_directional_qc_pass = discovery_direction_concordant,
      tier_one_side_significant = discovery_direction_concordant &
        (discovery_q_EO < 0.1 | discovery_q_LO < 0.1),
      tier_strict_core = discovery_direction_concordant &
        discovery_q_EO < 0.1 & discovery_q_LO < 0.1
    ) %>%
    filter(tier_all_directional_qc_pass) %>%
    select(-Species_EO, -Species_LO) %>%
    arrange(Species)
}

# Extract held-out coefficients without selecting on validation results.
extract_validation_meta <- function(fit_object, suffix) {
  extracted <- fit_object$meta_fits %>%
    as_tibble() %>%
    rowwise() %>%
    mutate(max_weight = ifelse(
      all(is.na(c_across(starts_with("weight_")))),
      NA_real_,
      max(c_across(starts_with("weight_")), na.rm = TRUE)
    )) %>%
    ungroup() %>%
    select(
      feature,
      coef,
      stderr,
      pval,
      k,
      tau2,
      I2,
      max_weight
    )
  names(extracted)[-1] <- paste0(names(extracted)[-1], "_", suffix)
  extracted
}

# Retain cohort estimates for the held-out sensitivity analyses.
extract_local_validation <- function(fit_object, age_stratum) {
  bind_rows(fit_object$maaslin_fits) %>%
    filter(metadata == "Group", value == "CRC") %>%
    transmute(
      age_stratum = age_stratum,
      Cohort = batch,
      feature,
      local_beta = coef,
      local_se = stderr,
      local_p = pval,
      n_model = N,
      n_nonzero = N.not.zero
    )
}

# Summarise finite differences with the original test definitions.
summarise_test_statistics <- function(values) {
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(tibble(
      n = 0L,
      median = NA_real_,
      mean = NA_real_,
      n_positive = 0L,
      proportion_positive = NA_real_,
      wilcoxon_p_one_sided = NA_real_,
      wilcoxon_method = NA_character_,
      sign_test_p_one_sided = NA_real_
    ))
  }

  nonzero_values <- values[values != 0]
  exact_possible <- length(nonzero_values) > 0 &&
    length(nonzero_values) <= 50 &&
    !anyDuplicated(abs(nonzero_values))
  wilcoxon_result <- tryCatch(
    wilcox.test(
      values,
      mu = 0,
      alternative = "greater",
      exact = exact_possible,
      correct = FALSE
    ),
    error = function(error) NULL
  )
  sign_p <- if (length(nonzero_values) > 0) {
    binom.test(
      sum(nonzero_values > 0),
      length(nonzero_values),
      p = 0.5,
      alternative = "greater"
    )$p.value
  } else {
    1
  }

  tibble(
    n = length(values),
    median = median(values),
    mean = mean(values),
    n_positive = sum(values > 0),
    proportion_positive = mean(values >
      0),
    wilcoxon_p_one_sided = if (is.null(wilcoxon_result)) {
      NA_real_
    } else {
      wilcoxon_result$p.value
    },
    wilcoxon_method = if (exact_possible)
      "exact"
    else "asymptotic_no_cc",
    sign_test_p_one_sided = sign_p
  )
}

# Resample taxa for the original conditional median interval.
bootstrap_median_ci <- function(values, repetitions) {
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(c(NA_real_, NA_real_))
  }

  if (length(values) == 1) {
    return(rep(values, 2))
  }

  bootstrap_medians <- replicate(
    repetitions,
    median(sample(values, length(values), replace = TRUE))
  )
  as.numeric(quantile(
    bootstrap_medians,
    probs = c(0.025, 0.975),
    names = FALSE,
    type = 6
  ))
}
