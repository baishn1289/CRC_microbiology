#humann_heatmap
library(patchwork)
library(tidyverse)
library(ggpubr)
library(metafor)

bad_keywords <- paste("plant",  "animal",  "yeast",  "fungi",
                      "mammal",  "eukar",  "mitochond",  "chlorophyll",  "photosynth",
                      "rubisco",  "engineered",  "invertebrate",  "mitochondria",
                      "archaea",  sep = "|")

metacyc_res_EO <-   process_meta_with_rescue(
    meta_res_obj = metacyc_meta.crc_EO,
    group_name = "metacyc_CRC_EO",
    i2_tier1 = 50,      
    i2_tier2 = 80,       
    weight_limit = 40,
    base_consistency = 0.6,  
    rescue_consistency  = 0.75 
  )  %>%  filter(!str_detect(Pathway, regex(bad_keywords, ignore_case = TRUE))) %>% 
  select(-qval.fdr)

metacyc_res_EO_q_val <- metacyc_res_EO %>%   
  filter(str_starts(qc_status,'Pass')) %>%
  filter(!str_detect(Pathway, regex(bad_keywords, ignore_case = TRUE))) %>% 
  mutate( qval.fdr = p.adjust(pval, method = "fdr") ) %>% 
  select(Pathway,qval.fdr)

metacyc_res_EO <- left_join(metacyc_res_EO,metacyc_res_EO_q_val,
                            by = 'Pathway') %>%
  mutate(Group = "EO") 

metacyc_res_LO <- process_meta_with_rescue(
  meta_res_obj = metacyc_meta.crc_LO,
  group_name = "metacyc_CRC_LO",
  i2_tier1 = 50,       
  i2_tier2 = 80,      
  weight_limit = 40,
  base_consistency = 0.6,  
  rescue_consistency  = 0.75 
) %>%  filter(!str_detect(Pathway, regex(bad_keywords, ignore_case = TRUE))) %>% 
  select(-qval.fdr)

metacyc_res_LO_q_val <- metacyc_res_LO %>%   
  filter(str_starts(qc_status,'Pass')) %>%
  filter(!str_detect(Pathway, regex(bad_keywords, ignore_case = TRUE))) %>% 
  mutate( qval.fdr = p.adjust(pval, method = "fdr") ) %>% 
  select(Pathway,qval.fdr)

metacyc_res_LO <- left_join(metacyc_res_LO,metacyc_res_LO_q_val,
                            by = 'Pathway') %>%
  mutate(Group = "LO") 

metacyc_filtered_plot_data <- bind_rows(metacyc_res_EO, metacyc_res_LO) %>%
  mutate(Pathway_Short = str_extract(Pathway, "^[A-Z0-9+]+(-[A-Z0-9+]+)*")) %>% 
  select(Pathway, Pathway_Short, Group, exposure,coef, pval, qval.fdr, I2,pass_qc, 
         stderr,k,max_weight,consistency_score,qc_status) %>%
  group_by(Pathway) %>%
  filter(n_distinct(Group) == 2) %>%
  filter(any(qval.fdr < 0.1  & pval <0.05, na.rm = TRUE)) %>%
  ungroup() %>% 
  mutate( CI_lower = coef - (1.96 * stderr),
          CI_upper = coef + (1.96 * stderr),
          Report = sprintf("β = %.2f; 95%% CI, %.2f–%.2f; adjusted p = %.1e", 
                           coef, CI_lower, CI_upper, qval.fdr)) 

metacyc_plot_data_explortary <- metacyc_filtered_plot_data %>%
  group_by(Pathway) %>%
  mutate(
    coef_EO_val = coef[Group == "EO"],
    coef_LO_val = coef[Group == "LO"],

    coef_abs_ratio = abs(coef_EO_val) / abs(coef_LO_val),

    coef_discussion = case_when(
      coef_abs_ratio < 2/3 ~ "yes", 
      coef_abs_ratio > 3/2 ~ "yes", 
      sign(coef_EO_val) != sign(coef_LO_val) ~ "yes",
      TRUE ~ "no"
    ),
    
    coef_EO_higher = if_else(coef_abs_ratio > 3/2, "yes", "no"),
    coef_opposite  = if_else(sign(coef_EO_val) != sign(coef_LO_val), "yes", "no"),
    coef_LO_higher = if_else(coef_abs_ratio < 2/3 & sign(coef_EO_val) == sign(coef_LO_val), "yes", "no")
  ) %>%
  ungroup()

df_raw_EO <- extract_cohort_stats(metacyc_meta.crc_EO$maaslin_fits, "EO")
df_raw_LO <- extract_cohort_stats(metacyc_meta.crc_LO$maaslin_fits, "LO")

metacyc_combined_cohort_data <- bind_rows(df_raw_EO, df_raw_LO) %>%
  mutate(Pathway_ID = str_trim(str_extract(feature, "^[^:]+")))

target_pathways_interaction <- unique(metacyc_combined_cohort_data$Pathway_ID)

interaction_results <-  map_dfr(
    target_pathways_interaction,
    ~ calc_interaction_p_robust(metacyc_combined_cohort_data, target_id = .x)
  )
  
interaction_results_final <- interaction_results %>%
  filter(!is.na(P_interaction))  %>% 
  mutate(
    Nominal_Significant = ifelse(P_interaction < 0.05, "Yes (*)", "No"),
    `95%_CI_Diff` = sprintf("[%.3f, %.3f]", CI_lower, CI_upper)
  )%>% left_join(metacyc_combined_cohort_data %>% 
                   select(feature,Pathway_ID) %>% unique(),
                 by = 'Pathway_ID')%>% 
  select(Pathway_Short = Pathway_ID, feature,Coef_diff_EO_vs_LO, `95%_CI_Diff`, 
         P_interaction, Nominal_Significant#, everything()
  )  %>%
  filter(feature %in% intersect(metacyc_res_EO[metacyc_res_EO$pass_qc == 'TRUE',]$Pathway, 
                                metacyc_res_LO[metacyc_res_LO$pass_qc == 'TRUE',]$Pathway)) %>%
  mutate( q.fdr.interaction = p.adjust(P_interaction, method = "fdr")) %>% 
  select(-feature) 

metacyc_plot_data_explortary_csv <- metacyc_plot_data_explortary %>% 
  select(Pathway,Age_class=Group,exposure, Report,  I2,pass_qc,qc_status,
         max_weight,consistency_score,
         coef_EO_val,coef_LO_val,coef_abs_ratio,coef_discussion,
         coef_EO_higher,coef_opposite,coef_LO_higher) %>% 
  pivot_wider(
    id_cols   = c(Pathway, exposure),
    names_from = Age_class,
    values_from = c(Report, I2, max_weight,consistency_score,qc_status),
    names_glue = "{Age_class}_{.value}"
  ) %>% mutate(Pathway_Short  = str_extract(Pathway, "^[^:]+")) %>% 
  left_join(interaction_results_final,by = 'Pathway_Short') %>% 
  select(-Pathway_Short) %>% 
  arrange(P_interaction,Pathway)

writexl::write_xlsx(x = metacyc_plot_data_explortary_csv,
          path = paste0(base_symbol,"downstream_analysis/output/table/TableS7.metacyc_pathway_results.xlsx")
          )

metacyc_results_EO_sensi <- process_meta_with_rescue(
  meta_res_obj = metacyc_meta.crc_EO_sensi,
  group_name = "metacyc_CRC_EO_sensi",
  i2_tier1 = 50,       
  i2_tier2 = 80,  
  weight_limit = 40,
  base_consistency = 0.6,  
  rescue_consistency  = 0.75 
) %>%   filter(!str_detect(Pathway, regex(bad_keywords, ignore_case = TRUE))) 

metacyc_results_EO_sensi_qval <- metacyc_results_EO_sensi %>%
  filter(str_starts(qc_status,'Pass')) %>%
  mutate( qval.fdr = p.adjust(pval, method = "fdr") ) %>% 
  select(Pathway,qval.fdr)

metacyc_results_EO_sensi <- metacyc_results_EO_sensi %>% 
  select(-qval.fdr) %>% 
  left_join(metacyc_results_EO_sensi_qval,by='Pathway') %>% 
  select(Pathway,coef,qval.fdr,qc_status)

metacyc_results_EO_csv_sensi <-   left_join(metacyc_res_EO %>% select(Pathway,coef,qval.fdr,qc_status),
                                            metacyc_results_EO_sensi,
            by = 'Pathway',suffix = c('_main','_sensitive')) %>%
  mutate(coef_change_direction = ifelse(coef_sensitive * coef_main <0 ,'different','same')) %>% 
  filter(!is.na(qval.fdr_main)) %>% 
  arrange(coef_change_direction,qval.fdr_main)

write.csv(metacyc_results_EO_csv_sensi,row.names = F,
          file = paste0(base_symbol,
                        "downstream_analysis/output/table/TableS8.Early_onset_samples_metacyc_feature_sensitivity_analysis.csv"))

