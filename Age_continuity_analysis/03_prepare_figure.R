library(dplyr)
library(readr)
library(tibble)

# Reuse the main disease-QC function for the original EO/LO background only.
stopifnot(is.function(process_meta_with_rescue))
main <- result %>%
  filter(Target == "Common_age", qc_pass, q_after_QC < .1) %>%
  arrange(q_after_QC, p, Species, feature) %>%
  mutate(Figure_row = row_number(), coef = beta)

for (stratum in c("EO", "LO")) {
  disease <- process_meta_with_rescue(
    get(paste0("fit_meta.crc_", stratum)),
    i2_tier1 = 50,
    i2_tier2 = 80,
    weight_limit = 40,
    base_consistency = .60,
    rescue_consistency = .75
  )
  published <- disease %>%
    filter(grepl("^Pass", qc_status)) %>%
    mutate(q_publication = p.adjust(pval, "fdr"))
  raw_index <- match(main$feature, disease$Pathway)
  qc_index <- match(main$feature, published$Pathway)
  main[[paste0("original_", stratum, "_coef")]] <- disease$coef[raw_index]
  main[[paste0("original_", stratum, "_pval")]] <- disease$pval[raw_index]
  main[[paste0(
    "original_",
    stratum,
    "_q_original_publication_rule"
  )]] <- published$q_publication[qc_index]
  main[[paste0(stratum, "_available")]] <- !is.na(qc_index)
  main[[paste0(stratum, "_significant")]] <- !is.na(qc_index) &
    published$q_publication[qc_index] < .1
}

pair <- result %>%
  filter(Target %in% c("CRC_age", "Control_age"), feature %in% main$feature) %>%
  mutate(Group = sub("_age$", "", Target), Figure_row = match(feature, main$feature)) %>%
  arrange(Figure_row, Group)
stopifnot(nrow(main) == 20L, nrow(pair) == 40L)
inputs <- list(main = main, slopes = pair)
