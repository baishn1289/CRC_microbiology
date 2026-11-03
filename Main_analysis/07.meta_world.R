#####heatmap
library(ComplexHeatmap)
library(circlize)

rename_cohorts <- function(x) {
  case_when(
    x == 'A10' ~ 'YachidaS_2019',
    x == 'A12' ~ 'Zeller_2014',
    x == 'A13' ~ 'Feng_2015',
    x == 'A14' ~ 'Vogtmann_2016',
    x == 'A17' ~ 'Hannigan_2018',
    x == 'A18' ~ 'Gupta_2019',
    x == 'Public_study__YuJ_2015' ~ 'Yu_2017',
    x == 'A28' ~ 'Liu_2022',
    x == 'A30_FUSCC' ~ 'FUSCC−SHSD',
    x == 'A33' ~ 'Yang_2020',
    x == 'A37' ~ 'Piccinno_2025_IIGM_TU',
    x == 'A9' ~ 'Wirbel_2019',
    x == 'Public_study__ThomasAM_2018b' ~ 'ThomasAM_2018b',
    x == 'This_study_cohort3__ONCOBIOME_IIGM_CZ' ~ 'Piccinno_2025_IIGM_CZ',
    x == 'This_study_cohort4__ONCOBIOME_IIGM_IT' ~ 'Piccinno_2025_IIGM_IT',
    TRUE ~ x
  )
}

cohort_to_country <- c(
  'YachidaS_2019' = 'Japan', 'Zeller_2014' = 'France', 
  'Feng_2015' = 'Austria',
  'Vogtmann_2016' = 'US&Canada', 'Hannigan_2018' = 'US&Canada', 
  'Gupta_2019' = 'India',
  'Yu_2017' = 'China', 'Liu_2022' = 'China', 
  'FUSCC-SHSD' = 'China', 'Yang_2020' = 'China',
  'Piccinno_2025_IIGM_TU' = 'Turkey', 
  'Wirbel_2019' = 'Germany', 'ThomasAM_2018b' = 'Italy',
  'Piccinno_2025_IIGM_CZ' = 'Cezch', 
  'Piccinno_2025_IIGM_IT' = 'Italy',
  'Overall' = 'Pooled Analysis' 
)

####world level

meta_results_EO_qval <- process_meta_with_rescue(
  meta_res_obj = fit_meta.crc_EO,
  group_name = "CRC_EO",
  i2_tier1 = 50,      
  i2_tier2 = 80,      
  weight_limit = 40,
  base_consistency = 0.6,  
  rescue_consistency  = 0.75 
) %>% 
  filter(str_starts(qc_status,'Pass')) %>%
  mutate( qval.fdr = p.adjust(pval, method = "fdr") ) 


meta_results_LO_qval <- process_meta_with_rescue(
  meta_res_obj = fit_meta.crc_LO,
  group_name = "CRC_LO",
  i2_tier1 = 50,       
  i2_tier2 = 80,      
  weight_limit = 40,
  base_consistency = 0.6,  
  rescue_consistency  = 0.75 
) %>% 
  filter(str_starts(qc_status,'Pass')) %>%
  mutate( qval.fdr = p.adjust(pval, method = "fdr") )

species_interaction_results_qc <- species_interaction_results_final %>% 
  filter(Pathway %in% intersect(meta_results_LO_qval$Pathway,
                                      meta_results_EO_qval$Pathway)) %>%
  mutate( q.fdr.interaction = p.adjust(P_interaction, method = "fdr")) %>% 
  left_join(meta_results_LO_qval %>% select(Pathway,coef,qval.fdr,qc_status),
            by = 'Pathway') %>% 
  left_join(meta_results_EO_qval %>% select(Pathway,coef,qval.fdr,qc_status),
            by = 'Pathway',suffix = c('_LO','_EO')) %>% 
  select(Species,P_interaction,Nominal_Significant,Coef_diff_EO_vs_LO,
         `95%_CI_Diff`,q.fdr.interaction,coef_EO, coef_LO,
         qval.fdr_EO, qval.fdr_LO, qc_status_EO,qc_status_LO)

write.csv(species_interaction_results_qc,row.names = F,
          file = paste0(base_symbol,
          "downstream_analysis/output/table/TableS4.Species_P_interaction.csv"))


meta_results_EO_csv <- process_meta_with_rescue(
  meta_res_obj = fit_meta.crc_EO,
  group_name = "CRC_EO",
  i2_tier1 = 50,      
  i2_tier2 = 80,     
  weight_limit = 40,
  base_consistency = 0.6,  
  rescue_consistency  = 0.75 
) %>% left_join(meta_results_EO_qval %>% 
                  select(feature,qval.fdr),by = 'feature',suffix = c('_original','')) %>% 
  mutate(Species = str_remove(str_extract(feature, 's__[A-Za-z0-9_]+'), '^s__'),
         CI_lower = coef - (1.96 * stderr),
         CI_upper = coef + (1.96 * stderr),
         Report = sprintf("β = %.2f; 95%% CI, %.2f–%.2f; adjusted p = %.1e", 
                          coef, CI_lower, CI_upper, qval.fdr))  %>% 
  arrange(qval.fdr) %>% 
  select(Species, exposure,pval, Report,I2,max_weight,
         stderr,k,pval.tau2,consistency_score,
         pass_qc,qc_status)  

writexl::write_xlsx(x = meta_results_EO_csv,
                    path = paste0(base_symbol,
                                  "downstream_analysis/output/table/TableS2.Species_EOCRC_vs_EOcontorl.xlsx")
)

meta_results_EO <- meta_results_EO_qval %>% 
  filter(qval.fdr <0.1) 

meta_results_EO_heatmap <- meta_results_EO %>% 
  arrange(desc(abs(coef)))%>%
  slice_head(n = 30)


meta_results_LO_csv <- process_meta_with_rescue(
  meta_res_obj = fit_meta.crc_LO,
  group_name = "CRC_LO",
  i2_tier1 = 50,     
  i2_tier2 = 80,   
  weight_limit = 40,
  base_consistency = 0.6,  
  rescue_consistency  = 0.75 
) %>% left_join(meta_results_LO_qval %>% 
                  select(feature,qval.fdr),by = 'feature',suffix = c('_original','')) %>% 
  mutate(Species = str_remove(str_extract(feature, 's__[A-Za-z0-9_]+'), '^s__'),
         CI_lower = coef - (1.96 * stderr),
         CI_upper = coef + (1.96 * stderr),
         Report = sprintf("β = %.2f; 95%% CI, %.2f–%.2f; adjusted p = %.1e", 
                                   coef, CI_lower, CI_upper, qval.fdr)) %>% 
  arrange(qval.fdr) %>% 
  select(Species, exposure,pval, Report,I2,max_weight,
         stderr,k,pval.tau2,consistency_score,
         pass_qc,qc_status)  


writexl::write_xlsx( x = meta_results_LO_csv,
          path = paste0(base_symbol,
                        "downstream_analysis/output/table/TableS3.Species_LOCRC_vs_LOcontorl.xlsx")
)

meta_results_LO <- meta_results_LO_qval %>% 
  filter(qval.fdr <0.1)

meta_results_LO_heatmap <- meta_results_LO %>% 
  arrange(desc(abs(coef))) %>%
    slice_head(n = 30)

table(meta_results_LO$coef < 0)
