library(dplyr)
library(tidyr)

frozen <- all_results
meta <- direct_meta
local_results <- direct_local
eo_scopes <- c(
  "binary__site__EO__Rectum_vs_Colon",
  "sidedness__site__EO__left_sided_vs_right_sided",
  "ternary__site__EO__Left_Colon_vs_Right_Colon",
  "ternary__site__EO__Rectum_vs_Left_Colon",
  "ternary__site__EO__Rectum_vs_Right_Colon"
)
details <- bind_rows(lapply(sets[eo_scopes], function(s) s$metadata %>%
  mutate(scope = s$scope))) %>%
  count(scope, Cohort, Site, name = "n")
raw_blocks <- list()
raw_taxa <- c("Dialister pneumosintes", "Haemophilus parainfluenzae", "Veillonella rogosae")

for (taxon in raw_taxa) {
  candidate <- frozen %>%
    filter(
      reportable,
      in_BH_family,
      q_comparison < .1,
      .data$taxon == .env$taxon,
      grepl("__site__EO__", scope, fixed = TRUE)
    )
  stopifnot(nrow(candidate) == 1L)
  s <- sets[[candidate$scope]]
  raw <- abundance_full[candidate$feature, s$metadata$Sample_ID]
  raw_blocks[[taxon]] <- s$metadata %>%
    transmute(
      scope = candidate$scope,
      feature = candidate$feature,
      taxon = .env$taxon,
      Cohort,
      cohort_label = ifelse(Cohort == "FUSCC-SHSD", "FUSCC", "Yachida 2019"),
      Sample_ID,
      Site = as.character(Site),
      side = ifelse(grepl("[Rr]ight", Site), "Right", "Left"),
      site_label = if (grepl("ternary", candidate$scope))
        paste(side, "colon") else
        paste0(side, "-sided"),
      abundance_percent = as.numeric(raw),
      detected = raw > .001,
      display = log10(raw + .0001)
    )
}

raw_plot <- bind_rows(raw_blocks)
detection_plot <- raw_plot %>%
  group_by(scope, feature, taxon, Cohort, cohort_label, side, site_label) %>%
  summarise(
    n = n(),
    n_detected = sum(detected),
    n_zero = sum(abundance_percent == 0),
    .groups = "drop"
  )
j <- match(paste(frozen$scope, frozen$feature), paste(meta$scope, meta$feature))
stopifnot(
  !anyNA(j),
  identical(frozen$estimate, meta$estimate[j]),
  identical(frozen$p_value, meta$p_value[j])
)
plot_results <- frozen %>%
  mutate(tau2 = meta$tau2[j], Q_p = meta$Q_p[j])
all_hits <- plot_results %>%
  filter(reportable, in_BH_family, !is.na(q_comparison), q_comparison < .1)
eo_scopes <- c(
  "binary__site__EO__Rectum_vs_Colon",
  "sidedness__site__EO__left_sided_vs_right_sided",
  "ternary__site__EO__Left_Colon_vs_Right_Colon",
  "ternary__site__EO__Rectum_vs_Left_Colon",
  "ternary__site__EO__Rectum_vs_Right_Colon"
)
eo_hits <- all_hits %>%
  filter(scope %in% eo_scopes)
lo_hits <- all_hits %>%
  filter(scope == "sidedness__site__LO__left_sided_vs_right_sided")
age_hit <- all_hits %>%
  filter(scope == "binary__age__joint__Colon")
scope_summary <- frozen %>%
  group_by(scope, policy_tier) %>%
  summarise(
    target_cohorts = first(target_cohorts),
    target_n_samples = first(target_n_samples),
    n_original = n(),
    n_BH = sum(in_BH_family),
    n_hits = if (any(in_BH_family)) sum(comparison_q_lt010, na.rm = TRUE) else NA_integer_,
    min_q = if (any(in_BH_family)) min(q_comparison, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

binary <- details %>%
  filter(scope == eo_scopes[1]) %>%
  select(Cohort, Site, n) %>%
  pivot_wider(names_from = Site, values_from = n)
ternary <- details %>%
  filter(scope == eo_scopes[3]) %>%
  select(Cohort, Site, n) %>%
  pivot_wider(names_from = Site, values_from = n)
coverage_matrix <- binary %>%
  left_join(ternary, by = "Cohort") %>%
  mutate(
    across(c(Right_Colon, Left_Colon), ~ replace_na(.x, 0L)),
    Colon_unspecified = Colon - Right_Colon - Left_Colon
  ) %>%
  select(Cohort, Right_Colon, Left_Colon, Colon_unspecified, Rectum) %>%
  pivot_longer(-Cohort, names_to = "site_record", values_to = "n") %>%
  mutate(
    row_index = match(Cohort, c("FUSCC-SHSD", "YachidaS_2019", "Yang_2020")),
    col_index = match(site_record, c("Right_Colon", "Left_Colon", "Colon_unspecified", "Rectum"))
  ) %>%
  arrange(row_index, col_index)
coverage_by_comparison <- details %>%
  group_by(scope) %>%
  summarise(
    n = sum(n),
    k = n_distinct(Cohort),
    .groups = "drop"
  ) %>%
  left_join(scope_summary %>%
    select(scope, n_BH, n_hits), by = "scope")

cohort_order <- c(
  "FUSCC-SHSD",
  "YachidaS_2019",
  "Wirbel_2019",
  "Feng_2015",
  "Vogtmann_2016",
  "Yu_2017",
  "Zeller_2014",
  "Piccinno_2025_IIGM_IT",
  "Yang_2020"
)
cohort_names <- setNames(
  c(
    "FUSCC",
    "Yachida 2019",
    "Wirbel 2019",
    "Feng 2015",
    "Vogtmann 2016",
    "Yu 2017",
    "Zeller 2014",
    "Piccinno Italy",
    "Yang 2020"
  ),
  cohort_order
)
selected <- bind_rows(
  eo_hits %>%
    filter(taxon %in% c(
      "Dialister pneumosintes",
      "Haemophilus parainfluenzae",
      "Veillonella rogosae"
    )) %>%
    mutate(panel = case_when(
      taxon == "Dialister pneumosintes" ~ "c",
      taxon == "Haemophilus parainfluenzae" ~ "d",
      TRUE ~ "e"
    )),
  lo_hits %>%
    mutate(panel = "f"),
  age_hit %>%
    mutate(panel = "g")
)
forest_blocks <- list()

for (i in seq_len(nrow(selected))) {
  row <- selected[i, ]
  z <- local_results %>%
    filter(scope == row$scope, feature == row$feature, status == "ok") %>%
    arrange(match(Cohort, cohort_order))
  w <- 1 / (z$se^2 + row$tau2)
  forest_blocks[[i]] <- bind_rows(
    z %>%
      transmute(
        scope,
        feature,
        taxon,
        panel = row$panel,
        Cohort,
        Level = unname(cohort_names[Cohort]),
        n_model,
        n_detected,
        estimate,
        se,
        ci_low,
        ci_high,
        p_value,
        q_comparison = NA_real_,
        weight = 100 * w / sum(w),
        kind = "cohort"
      ),
    row %>%
      transmute(
        scope,
        feature,
        taxon,
        panel,
        Cohort = "HS pooled",
        Level = "HS pooled",
        n_model = n_samples,
        n_detected = NA_integer_,
        estimate,
        se,
        ci_low,
        ci_high,
        p_value,
        q_comparison,
        weight = NA_real_,
        kind = "pooled"
      )
  )
}

forest_rows <- bind_rows(forest_blocks)
