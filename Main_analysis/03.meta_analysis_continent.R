library(data.table)
library(tidyverse)
library(Maaslin2)
library(MMUPHin)
library(purrr)

process_meta_qval <- function(fit_meta,weight_limit_number = 40){
  
  qval_meta = process_meta_with_rescue(
    meta_res_obj = fit_meta,
    group_name = "NA",
    i2_tier1 = 50,       
    i2_tier2 = 80,      
    weight_limit = weight_limit_number,
    base_consistency = 0.6,  
    rescue_consistency  = 0.75 
  ) %>% 
    filter(str_starts(qc_status,'Pass')) %>%
    mutate( qval.fdr = p.adjust(pval, method = "fdr") ) %>% 
    select(Pathway,qval.fdr)
  
  meta_table = process_meta_with_rescue(
    meta_res_obj = fit_meta,
    group_name = "NA",
    i2_tier1 = 50,       
    i2_tier2 = 80,      
    weight_limit = weight_limit_number,
    base_consistency = 0.6,  
    rescue_consistency  = 0.75 
  ) %>%  select(-qval.fdr) %>% 
    left_join(qval_meta,by = 'Pathway')
  
  return(meta_table)
}

###Europe
metadata_Europe <- metadata_allcohorts %>% filter(Continent  == 'Europe' )
table(metadata_Europe$Cohort)

metadata_Europe_EO <- metadata_Europe %>% filter(Age <50)
table(metadata_Europe_EO$Cohort,metadata_Europe_EO$Group)

metadata_Europe_LO <-  metadata_Europe %>% filter(Age >49)
table(metadata_Europe_LO$Cohort,metadata_Europe_LO$Group)


abundance_EO_LO_combine_Europe <- run_consensus_EO_LO_combine(
  meta_list = meta_list_metadata_intersect_samples_Europe,
  meta_EO = metadata_Europe_EO,
  meta_LO = metadata_Europe_LO,min_prevalence = 0.05,abund_threshold = 0.001,min_cohorts = 2
)

feat_EO_Europe <- abundance_EO_LO_combine_Europe$EO 
feat_LO_Europe <- abundance_EO_LO_combine_Europe$LO 

#EOCRC vs EOCTR 
cl_metainput_EO_Europe = metadata_Europe_EO %>% 
  column_to_rownames('Sample_ID') %>% mutate(Cohort = as.factor(Cohort))

feat_EO_Europe <- abundance_EO_LO_combine_Europe$EO %>% 
  column_to_rownames("clade_name") %>% 
  as.matrix() %>% 
  apply_step2_masking( meta_data = cl_metainput_EO,
                       batch_col = "Cohort")

table(cl_metainput_EO_Europe$Cohort,cl_metainput_EO_Europe$Group) 

abundance_metainput_EO_Europe = feat_EO_Europe %>% 
  as.data.frame() %>% 
  select(any_of(c(rownames(cl_metainput_EO_Europe)))) 


fit_meta.crc_Europe_EO <- MMUPHin::lm_meta(feature_abd = abundance_metainput_EO_Europe,
                                           exposure = 'Group',
                                           batch = 'Cohort',
                                           covariates = c("Age","Sex",'BMI'),
                                           control = list(rma_method="HS", transform="LOG", normalization = 'NONE',
                                                          min_prevalence = 0,min_abundance = 0.001,
                                                          output= paste0(base_symbol,
                                                                        "downstream_analysis/meta_analysis_continent/Europe/EOCRC")
                                           ),
                                           data = cl_metainput_EO_Europe )

meta_results_Europe_EO <- process_meta_qval(fit_meta.crc_Europe_EO,
                                            weight_limit_number = 50) %>% 
  mutate(Continent = 'Europe',Age_class = 'EO') 

summary(meta_results_Europe_EO$I2)

#LOCRC vs LOCTR

cl_metainput_LO_Europe = metadata_Europe_LO  %>% 
  column_to_rownames('Sample_ID') %>% mutate(Cohort = as.factor(Cohort))

feat_LO_Europe <- abundance_EO_LO_combine_Europe$LO %>% 
  column_to_rownames("clade_name") %>% 
  as.matrix() %>% 
  apply_step2_masking( meta_data = cl_metainput_LO,
                       batch_col = "Cohort")

table(cl_metainput_LO_Europe$Cohort,cl_metainput_LO_Europe$Group) 

abundance_metainput_LO_Europe = feat_LO_Europe %>% 
  as.data.frame() %>% 
  select(any_of(c(rownames(cl_metainput_LO_Europe))))

setdiff(rownames(cl_metainput_LO_Europe),colnames(abundance_metainput_LO_Europe))

fit_meta.crc_Europe_LO <- MMUPHin::lm_meta(feature_abd = abundance_metainput_LO_Europe,
                                           exposure = 'Group',
                                           batch = 'Cohort',
                                           covariates = c("Age","Sex",'BMI'),
                                           control = list(rma_method="HS", transform="LOG", normalization = 'NONE',
                                                          min_prevalence = 0,min_abundance = 0.001,
                                                          output= paste0(base_symbol,
                                                                        "downstream_analysis/meta_analysis_continent/Europe/LOCRC")
                                           ),
                                           data = cl_metainput_LO_Europe )

meta_results_Europe_LO <- process_meta_qval(fit_meta.crc_Europe_LO,
                                            weight_limit_number = 50) %>% 
  mutate(Continent = 'Europe',Age_class = 'LO') 

###East Asia
metadata_Asia <- metadata_allcohorts %>% filter(Continent  == 'Asia' )
table(metadata_Asia$Cohort)

metadata_Asia_EO <- metadata_Asia %>% filter(Age <50)
table(metadata_Asia_EO$Cohort,metadata_Asia_EO$Group)


metadata_Asia_LO <-  metadata_Asia %>% filter(Age >49)
table(metadata_Asia_LO$Cohort,metadata_Asia_LO$Group)

abundance_Asia_EO_LO_combine <- run_consensus_EO_LO_combine(
  meta_list = meta_list_metadata_intersect_samples_Asia,
  meta_EO = metadata_Asia_EO,
  meta_LO = metadata_Asia_LO,min_prevalence = 0.05,abund_threshold = 0.001,min_cohorts = 2
)
feat_EO_Asia <- abundance_Asia_EO_LO_combine$EO
feat_LO_Asia <- abundance_Asia_EO_LO_combine$LO

#EOCRC vs EOCTR 
cl_metainput_EO_Asia = metadata_Asia_EO %>% 
  column_to_rownames('Sample_ID') %>% 
  mutate(Cohort = as.factor(Cohort))

feat_EO_Asia <- abundance_Asia_EO_LO_combine$EO %>% 
  column_to_rownames("clade_name") %>% 
  as.matrix() %>% 
  apply_step2_masking( meta_data = cl_metainput_EO,
                       batch_col = "Cohort")

table(cl_metainput_EO_Asia$Cohort,cl_metainput_EO_Asia$Group) 

abundance_metainput_EO_Asia = feat_EO_Asia %>% 
  as.data.frame() %>% 
  select(any_of(c(rownames(cl_metainput_EO_Asia)))) 

fit_meta.crc_Asia_EO <- MMUPHin::lm_meta(feature_abd = abundance_metainput_EO_Asia,
                                         exposure = 'Group',
                                         batch = 'Cohort',
                                         covariates = c("Age","Sex",'BMI'),
                                         control = list(rma_method="HS", transform="LOG", normalization = 'NONE',
                                                        min_prevalence = 0,min_abundance = 0.001,
                                                        output= paste0(base_symbol,
                                                        "downstream_analysis/meta_analysis_continent/East_Asia/EOCRC")
                                         ),
                                         data = cl_metainput_EO_Asia )

meta_results_Asia_EO <- process_meta_qval(fit_meta.crc_Asia_EO,
                    weight_limit_number = 50) %>% 
  mutate(Continent = 'East Asia',Age_class = 'EO') 

summary(meta_results_Asia_EO$I2)

#LOCRC vs LOCTR
cl_metainput_LO_Asia = metadata_Asia_LO %>% 
  column_to_rownames('Sample_ID') %>% mutate(Cohort = as.factor(Cohort))

abundance_metainput_LO_Asia = feat_LO_Asia %>% 
  select(any_of(c('clade_name',rownames(cl_metainput_LO_Asia)))) %>% 
  column_to_rownames('clade_name')

setdiff(rownames(cl_metainput_LO_Asia),colnames(abundance_metainput_LO_Asia))

fit_meta.crc_Asia_LO <- MMUPHin::lm_meta(feature_abd = abundance_metainput_LO_Asia,
                                         exposure = 'Group',
                                         batch = 'Cohort',
                                         covariates = c("Age","Sex",'BMI'),
                                         control = list(rma_method="HS", transform="LOG", normalization = 'NONE',
                                                        min_prevalence = 0,min_abundance = 0.001,
                                                        output= paste0(base_symbol,
                                                        "downstream_analysis/meta_analysis_continent/East_Asia/LOCRC")
                                         ),
                                         data = cl_metainput_LO_Asia )

meta_results_Asia_LO <- process_meta_qval(fit_meta.crc_Asia_LO,
                                          weight_limit_number = 50) %>% 
  mutate(Continent = 'East Asia',Age_class = 'LO') 
