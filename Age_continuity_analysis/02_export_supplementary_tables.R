library(dplyr)
library(readr)

# Project current unified estimates; never recalculate q within displayed rows.
age_table <- result %>%
  mutate(target_order = match(
    Target,
    c("Common_age", "CRC_age", "Control_age", "Interaction_v7")
  )) %>%
  arrange(target_order, q_after_QC, p, Species, feature) %>%
  transmute(
    Species,
    Analysis = case_when(
      Target == "Common_age" ~ "Common age",
      Target == "Interaction_v7" ~ "Age-by-CRC interaction",
      TRUE ~ "Independent age"
    ),
    Population = case_when(
      Target == "Common_age" ~ "CRC and Control",
      Target == "Interaction_v7" ~ "CRC minus Control",
      TRUE ~ sub("_age$", "", Target)
    ),
    `Inference method` = method,
    `Cohorts (k)` = k,
    N = n,
    `Beta (log2 / 10 years)` = beta,
    SE = se,
    `95% CI lower` = ci_lo,
    `95% CI upper` = ci_hi,
    P = p,
    q = ifelse(qc_pass, q_after_QC, NA_real_),
    `Target QC` = ifelse(qc_pass, "Pass", "Fail"),
    `I2 (%)` = I2,
    `Maximum cohort weight (%)` = max_weight_pct,
    `Direction consistency` = direction_consistency,
    Status = status
  )
write_csv(
  filter(age_table, Analysis != "Age-by-CRC interaction"),
  file.path(output_dir, "Supplementary_Table_Common_Independent_Age_df10_KH.csv"),
  na = "NA"
)
write_csv(
  filter(age_table, Analysis == "Age-by-CRC interaction"),
  file.path(output_dir, "Supplementary_Table_Interaction_df10_KH.csv"),
  na = "NA"
)
