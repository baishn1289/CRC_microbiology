#9.humann_meta_analysis

#world_level_metacyc
metacyc_world_res <- run_consensus_pathway_EO_LO(
  meta_list = metacyc_list,  
  meta_EO = world_meta_EO,
  meta_LO = world_meta_LO,
  min_prevalence = 0.1,    
  min_cohorts = 3,         
  abund_threshold = 1e-6,  
  clean_names = TRUE     
)

metacyc_merged_df <- metacyc_world_res$merged
dim(metacyc_merged_df)

###EOCRC
metacyc_feat_EO <- metacyc_world_res$EO %>% 
  column_to_rownames("Feature") %>% 
  as.matrix() %>% 
  apply_pathway_masking( meta_data = cl_metainput_EO,
                       batch_col = "Cohort")

table(cl_metainput_EO$Cohort,cl_metainput_EO$Group) 

metacyc_metainput_EO = metacyc_feat_EO %>% 
  as.data.frame() %>% 
  select(any_of(c(rownames(cl_metainput_EO))))

setdiff(names(metacyc_metainput_EO),rownames(cl_metainput_EO))
setdiff(rownames(cl_metainput_EO),names(metacyc_metainput_EO))

dim(metacyc_metainput_EO)
dim(cl_metainput_EO)

metacyc_meta.crc_EO <- MMUPHin::lm_meta(feature_abd = metacyc_metainput_EO,
                                    exposure = 'Group',
                                    batch = 'Cohort',
                                    covariates = c("Age","Sex",'BMI'),
                                    control = list(rma_method="HS", transform="LOG", normalization = 'NONE',
                                                   min_prevalence = 0,min_abundance = 1e-6,
                                                   output= paste0(base_symbol,"/downstream_analysis/output/humann_meta/EOCRC/metacyc/")
                                    ),
                                    data = cl_metainput_EO )

###LOCRC
metacyc_feat_LO <- metacyc_world_res$LO %>% 
  column_to_rownames("Feature") %>% 
  as.matrix() %>% 
  apply_pathway_masking( meta_data = cl_metainput_LO,
                         batch_col = "Cohort")

table(cl_metainput_LO$Cohort,cl_metainput_LO$Group) 

metacyc_metainput_LO = metacyc_feat_LO %>% 
  as.data.frame() %>% 
  select(any_of(c(rownames(cl_metainput_LO))))

setdiff(names(metacyc_metainput_LO),rownames(cl_metainput_LO))
setdiff(rownames(cl_metainput_LO),names(metacyc_metainput_LO))

dim(metacyc_metainput_LO)
dim(cl_metainput_LO)


metacyc_meta.crc_LO <- MMUPHin::lm_meta(feature_abd = metacyc_metainput_LO,
                                        exposure = 'Group',
                                        batch = 'Cohort',
                                        covariates = c("Age","Sex",'BMI'),
                                        control = list(rma_method="HS", transform="LOG", normalization = 'NONE',
                                                       min_prevalence = 0,min_abundance = 1e-6,
                                                       output= paste0(base_symbol,"downstream_analysis/output/humann_meta/LOCRC/metacyc/")
                                        ),
                                        data = cl_metainput_LO )


#####敏感性分析
world_meta_EO_sensi <- world_meta_EO %>% 
  filter(!Cohort %in% c("Feng_2015",'ThomasAM_2018b')) 

unique(world_meta_EO_sensi$Cohort)

world_meta_LO_sensi <- world_meta_LO

metacyc_world_EO_LO_combine_sensi <- run_consensus_pathway_EO_LO(
    meta_list = metacyc_list, 
    meta_EO = world_meta_EO_sensi,
    meta_LO = world_meta_LO_sensi,
    min_prevalence = 0.1,   
    min_cohorts = 3,         
    abund_threshold = 1e-6,  
    clean_names = TRUE   
  )

metacyc_EO_sensi <- metacyc_world_EO_LO_combine_sensi$EO %>% 
  column_to_rownames("Feature") %>% 
  as.matrix() %>% 
  apply_pathway_masking( meta_data = cl_metainput_EO_sensi,
                       batch_col = "Cohort") 

metacyc_metainput_EO_sensi = metacyc_EO_sensi %>% 
  as.data.frame() %>% 
  select(any_of(c(rownames(cl_metainput_EO)))) 

setdiff(names(metacyc_metainput_EO_sensi),rownames(cl_metainput_EO_sensi))
setdiff(rownames(cl_metainput_EO_sensi),names(metacyc_metainput_EO_sensi))

metacyc_meta.crc_EO_sensi <- MMUPHin::lm_meta(feature_abd = metacyc_metainput_EO_sensi,
                                          exposure = 'Group',
                                          batch = 'Cohort',
                                          covariates = c("Age","Sex",'BMI'),
                                          control = list(rma_method="HS", transform="LOG", normalization = 'NONE',
                                                         min_prevalence = 0,min_abundance = 0,
                                                         output= paste0(base_symbol,"downstream_analysis/output/humann_meta/EOCRC/metacyc_sensi/")
                                          ),
                                          data = cl_metainput_EO_sensi )
