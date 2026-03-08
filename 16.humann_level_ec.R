##level ec database
library(lme4)
library(lmerTest)
library(ggpubr)
library(RColorBrewer)
library(KEGGREST)

source("distant_sever_code/19.get_ec_info_tidy_function_use.R")

ec_list <- lapply(dir_list, function(i) {
  file_path <- file.path(read_dir, i, "all_regrouped_level4ec_unstratified_relab.tsv")
  
  if(file.exists(file_path)){
    dt <- fread(file_path, header = TRUE)
    names(dt) <- str_remove(names(dt), '_Abundance-RPKs')
    names(dt)[1] <- "Gene_Family"
    if (i %in% c('A18', 'A17')){
      colnames(dt)[-1] <- paste(i,colnames(dt)[-1],sep = '_')
    }
    inter_samples <- intersect(metadata_allcohorts$Sample_ID,
                               colnames(dt))
    dt <- dt %>% select(any_of(c("Gene_Family",inter_samples)))
    return(dt)
  } else {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
})

ec_world_res <- run_consensus_pathway_EO_LO(
  meta_list = ec_list,   
  meta_EO = world_meta_EO,
  meta_LO = world_meta_LO,
  min_prevalence = 0.1,    
  min_cohorts = 3,         
  abund_threshold = 1e-6,  
  clean_names = TRUE       
)

ec_merged_df <- ec_world_res$merged
dim(ec_merged_df)
ec_mat_raw <- ec_merged_df %>% column_to_rownames("Gene_Family")

setdiff(names(ec_mat_raw),rownames(meta_adj_batch))
setdiff(rownames(meta_adj_batch),names(ec_mat_raw))
ec_mat_raw <- ec_mat_raw[, rownames(meta_adj_batch)]


ec_mat_raw_adjust <- ec_mat_raw[, rownames(meta_adj_batch_adjust)]

fit_ec_adj_adjust <- MMUPHin::adjust_batch(
  feature_abd = ec_mat_raw_adjust,
  batch = "batch",
  covariates =c("Group", "Age","BMI","Sex"),
  data = meta_adj_batch_adjust,
  control = list(verbose = TRUE)
)

ec_mat_adj_adjust <- fit_ec_adj_adjust$feature_abd_adj %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'Gene_Family')

###EOCRC
ec_feat_EO <- ec_world_res$EO %>% 
  column_to_rownames("Gene_Family") %>% 
  as.matrix() %>% 
  apply_pathway_masking( meta_data = cl_metainput_EO,
                         batch_col = "Cohort")

table(cl_metainput_EO$Cohort,cl_metainput_EO$Group) 

ec_metainput_EO = ec_feat_EO %>% 
  as.data.frame() %>% 
  select(any_of(c(rownames(cl_metainput_EO))))

setdiff(names(metacyc_metainput_EO),rownames(cl_metainput_EO))
setdiff(rownames(cl_metainput_EO),names(metacyc_metainput_EO))


ec_meta.crc_EO <- MMUPHin::lm_meta(feature_abd  = ec_metainput_EO,
                                        exposure = 'Group',
                                        batch = 'Cohort',
                                        covariates = c("Age","Sex",'BMI'),
                                        control = list(rma_method="HS", transform="LOG", normalization = 'NONE',
                                                       min_prevalence = 0,min_abundance = 1e-6,
                                                       output= paste0(base_symbol,
                                                       "downstream_analysis/output/humann_meta/EOCRC/ec/")
                                        ),
                                        data = cl_metainput_EO )

###LOCRC
ec_feat_LO <- ec_world_res$LO %>% 
  column_to_rownames("Gene_Family") %>% 
  as.matrix() %>% 
  apply_pathway_masking( meta_data = cl_metainput_LO,
                         batch_col = "Cohort")

table(cl_metainput_LO$Cohort,cl_metainput_LO$Group) 

ec_metainput_LO = ec_feat_LO %>% 
  as.data.frame() %>% 
  select(any_of(c(rownames(cl_metainput_LO))))

setdiff(names(ec_metainput_LO),rownames(cl_metainput_LO))
setdiff(rownames(cl_metainput_LO),names(ec_metainput_LO))

dim(ec_metainput_LO)
dim(cl_metainput_LO)

ec_meta.crc_LO <- MMUPHin::lm_meta(feature_abd = ec_metainput_LO,
                                        exposure = 'Group',
                                        batch = 'Cohort',
                                        covariates = c("Age","Sex",'BMI'),
                                        control = list(rma_method="HS", transform="LOG", normalization = 'NONE',
                                                       min_prevalence = 0,min_abundance = 1e-6,
                                                       output=  paste0(base_symbol,
                                                       "downstream_analysis/output/humann_meta/LOCRC/ec/")
                                        ),
                                        data = cl_metainput_LO )

###敏感性分析
ec_world_EO_LO_combine_sensi <- run_consensus_pathway_EO_LO(
  meta_list = ec_list, 
  meta_EO = world_meta_EO_sensi,
  meta_LO = world_meta_LO_sensi,
  min_prevalence = 0.1,    
  min_cohorts = 3,         
  abund_threshold = 1e-6,  
  clean_names = TRUE     
)

ec_EO_sensi <- ec_world_EO_LO_combine_sensi$EO %>% 
  column_to_rownames("Gene_Family") %>% 
  as.matrix() %>% 
  apply_pathway_masking( meta_data = cl_metainput_EO_sensi,
                         batch_col = "Cohort") 

ec_metainput_EO_sensi = ec_EO_sensi %>% 
  as.data.frame() %>% 
  select(any_of(c(rownames(cl_metainput_EO)))) 

setdiff(names(ec_metainput_EO_sensi),rownames(cl_metainput_EO_sensi))
setdiff(rownames(cl_metainput_EO_sensi),names(ec_metainput_EO_sensi))

ec_meta.crc_EO_sensi <- MMUPHin::lm_meta(feature_abd = ec_metainput_EO_sensi,
                                              exposure = 'Group',
                                              batch = 'Cohort',
                                              covariates = c("Age","Sex",'BMI'),
                                              control = list(rma_method="HS", transform="LOG", normalization = 'NONE',
                                                             min_prevalence = 0,min_abundance = 0,
                                                             output=  paste0(base_symbol,
                                                             "downstream_analysis/output/humann_meta/EOCRC/ec_sensi/"
                                              )),
                                              data = cl_metainput_EO_sensi )


ec_res_EO <- process_meta_with_rescue(
  meta_res_obj = ec_meta.crc_EO,
  group_name = "ec_CRC_EO",
  i2_tier1 = 50,      
  i2_tier2 = 80,      
  weight_limit = 40,
  base_consistency = 0.6,  
  rescue_consistency  = 0.75 
)  %>%  select(-qval.fdr)

ec_res_EO_q_val <- ec_res_EO %>%   
  filter(str_starts(qc_status,'Pass')) %>%
  filter(!str_detect(Pathway, regex(bad_keywords, ignore_case = TRUE))) %>% 
  mutate( qval.fdr = p.adjust(pval, method = "fdr") ) %>% 
  select(Pathway,qval.fdr)

ec_res_EO <- left_join(ec_res_EO,ec_res_EO_q_val,
                            by = 'Pathway') %>%
  mutate(Group = "EO") 

ec_res_LO <- process_meta_with_rescue(
  meta_res_obj = ec_meta.crc_LO,
  group_name = "ec_CRC_LO",
  i2_tier1 = 50,      
  i2_tier2 = 80,   
  weight_limit = 40,
  base_consistency = 0.6,  
  rescue_consistency  = 0.75 
)  %>%  select(-qval.fdr)

ec_res_LO_q_val <- ec_res_LO %>%   
  filter(str_starts(qc_status,'Pass')) %>%
  filter(!str_detect(Pathway, regex(bad_keywords, ignore_case = TRUE))) %>% 
  mutate( qval.fdr = p.adjust(pval, method = "fdr") ) %>% 
  select(Pathway,qval.fdr)

ec_res_LO <- left_join(ec_res_LO,ec_res_LO_q_val,
                            by = 'Pathway') %>%
  mutate(Group = "LO") 

ec_filtered_plot_data <- bind_rows(ec_res_EO, 
                                   ec_res_LO) %>%
  mutate(Pathway_Short = Pathway) %>% 
  select(Pathway, Pathway_Short, Group, exposure,coef, pval, qval.fdr, I2,pass_qc, 
         stderr,k,max_weight,consistency_score,qc_status) %>%
  group_by(Pathway) %>%
  filter(n_distinct(Group) == 2) %>%
  filter(any(qval.fdr < 0.1 & pval <0.05 , na.rm = TRUE)) %>%
  ungroup() %>% 
  mutate( CI_lower = coef - (1.96 * stderr),
          CI_upper = coef + (1.96 * stderr),
          Report = sprintf("β = %.2f; 95%% CI, %.2f–%.2f; adjusted p = %.1e", 
                           coef, CI_lower, CI_upper, qval.fdr)) 

ec_plot_data_explortary <- ec_filtered_plot_data %>%
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

ec_df_raw_EO <- extract_cohort_stats(ec_meta.crc_EO$maaslin_fits, "EO")
ec_df_raw_LO <- extract_cohort_stats(ec_meta.crc_LO$maaslin_fits, "LO")

ec_combined_cohort_data <- bind_rows(ec_df_raw_EO, ec_df_raw_LO) %>%
  mutate(Pathway_ID = str_trim(str_extract(feature, "^[^:]+")))

ec_target_pathways_interaction <- unique(ec_plot_data_explortary$Pathway)

ec_interaction_results <- map_dfr(
  ec_target_pathways_interaction,
  ~ calc_interaction_p_robust(ec_combined_cohort_data, target_id = .x)
)

ec_interaction_results_final <- ec_interaction_results %>%
  filter(!is.na(P_interaction)) %>%
  mutate(
    Nominal_Significant = ifelse(P_interaction < 0.05, "Yes (*)", "No"),
    `95%_CI_Diff` = sprintf("[%.3f, %.3f]", CI_lower, CI_upper)
  ) %>% 
  select(Pathway = Pathway_ID, Coef_diff_EO_vs_LO, `95%_CI_Diff`, 
         P_interaction, Nominal_Significant#, everything()
         ) %>%
  filter(Pathway %in% intersect(ec_res_EO[ec_res_EO$pass_qc == 'TRUE',]$Pathway, 
                                ec_res_LO[ec_res_LO$pass_qc == 'TRUE',]$Pathway)) %>%
  mutate( q.fdr.interaction = p.adjust(P_interaction, method = "fdr"))

ec_plot_data_explortary_csv <- ec_plot_data_explortary %>% 
  select(Pathway,Age_class=Group,exposure, Report,  I2,pass_qc,qc_status,
         max_weight,consistency_score,
         coef_EO_val,coef_LO_val,coef_abs_ratio,coef_discussion,
         coef_EO_higher,coef_opposite,coef_LO_higher) %>% 
  pivot_wider(
    id_cols   = c(Pathway, exposure),
    names_from = Age_class,
    values_from = c(Report, I2, max_weight,consistency_score,qc_status),
    names_glue = "{Age_class}_{.value}"
  ) %>% 
  left_join(ec_interaction_results_final,by = 'Pathway') %>% 
  arrange(P_interaction,Pathway)

writexl::write_xlsx(x = ec_plot_data_explortary_csv,
          path = paste0(base_symbol,
                 "downstream_analysis/output/table/TableS9.ec_enzyme_results.xlsx")
)

ec_results_EO_sensi <- process_meta_with_rescue(
  meta_res_obj = ec_meta.crc_EO_sensi,
  group_name = "ec_CRC_EO_sensi",
  i2_tier1 = 50,       
  i2_tier2 = 80,     
  weight_limit = 40,
  base_consistency = 0.6,  
  rescue_consistency  = 0.75 
) %>%   filter(!str_detect(Pathway, regex(bad_keywords, ignore_case = TRUE))) 

ec_results_EO_sensi_qval <- ec_results_EO_sensi %>%
  filter(str_starts(qc_status,'Pass')) %>%
  mutate( qval.fdr = p.adjust(pval, method = "fdr") ) %>% 
  select(Pathway,qval.fdr)

ec_results_EO_sensi <- ec_results_EO_sensi %>% 
  select(-qval.fdr) %>% 
  left_join(ec_results_EO_sensi_qval,by='Pathway') %>% 
  select(Pathway,coef,qval.fdr,qc_status)

ec_results_EO_csv_sensi <-   left_join(ec_res_EO %>% select(Pathway,coef,qval.fdr,qc_status),
                                            ec_results_EO_sensi,
                                            by = 'Pathway',suffix = c('_main','_sensitive')) %>%
  mutate(coef_change_direction = ifelse(coef_sensitive * coef_main <0 ,'different','same')) %>% 
  filter(!is.na(qval.fdr_main)) %>% 
  arrange(coef_change_direction,qval.fdr_main)

write.csv(ec_results_EO_csv_sensi,row.names = F,
          file = paste0(base_symbol,
                        "downstream_analysis/output/table/TableS10.Early_onset_samples_ec_feature_sensitivity_analysis.csv"))


ec_list_stratified <- lapply(dir_list, function(i) {
  file_path <- file.path(read_dir, i, "all_regrouped_level4ec_stratified_relab.tsv")
  
  if(file.exists(file_path)){
    dt <- fread(file_path, header = TRUE)
    names(dt) <- str_remove(names(dt), '_Abundance-RPKs')
    names(dt)[1] <- "Gene_Family"
    if (i %in% c('A18', 'A17')){
      colnames(dt)[-1] <- paste(i,colnames(dt)[-1],sep = '_')
    }
    inter_samples <- intersect(metadata_allcohorts$Sample_ID,
                               colnames(dt))
    dt <- dt %>% select(any_of(c("Gene_Family",inter_samples)))
    return(dt)
  } else {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
})

names(ec_list_stratified) <- dir_list
