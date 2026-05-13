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
  'Overall' = 'Pooled Analysis' # 为 Overall 列单独定义分类
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

###敏感性分析结果比较
meta_results_EO_sensi <- process_meta_with_rescue(
  meta_res_obj = fit_meta.crc_EO_sensi,
  group_name = "meta_CRC_EO_sensi",
  i2_tier1 = 50,      
  i2_tier2 = 80,      
  weight_limit = 40,
  base_consistency = 0.6,  
  rescue_consistency  = 0.75 
) 

meta_results_EO_sensi_qval <- meta_results_EO_sensi %>%
  filter(str_starts(qc_status,'Pass')) %>%
  mutate( qval.fdr = p.adjust(pval, method = "fdr") ) %>% 
  select(Pathway,qval.fdr)

meta_results_EO_sensi <- meta_results_EO_sensi %>% 
  select(-qval.fdr) %>% 
  left_join(meta_results_EO_sensi_qval,by='Pathway') %>% 
  select(Pathway,coef,qval.fdr,qc_status) %>% 
  mutate(Species = str_remove(str_extract(Pathway, 's__[A-Za-z0-9_]+'), '^s__'))

meta_results_EO_csv_sensi <-   left_join(meta_results_EO_csv %>% select(Species,coef,qval.fdr,qc_status),
                                       meta_results_EO_sensi,
                                       by = 'Species',suffix = c('_main','_sensitive')) %>%
  mutate(coef_change_direction = ifelse(coef_sensitive * coef_main <0 ,'different','same')) %>% 
  filter(!is.na(qval.fdr_main)) %>% 
  select(-Pathway) %>% 
  arrange(coef_change_direction,qval.fdr_main)

write.csv(meta_results_EO_csv_sensi,row.names = F,
          file = paste0(base_symbol,
                        "downstream_analysis/output/table/TableS5.Early_onset_samples_Bacteria_feature_sensitivity_analysis.csv")
)

cohort_results <- imap_dfr(fit_meta.crc_EO$maaslin_fits, function(df, name) {
  df %>% 
    filter(metadata == "Group") %>%  
    select(feature, coef, stderr, pval) %>%
    mutate(Cohort = name, Type = "Individual")
})

meta_summary <- fit_meta.crc_EO$meta_fits %>%
  select(feature, coef, stderr, pval = pval) %>% 
  mutate(Cohort = "Overall", Type = "Summary")

plot_data <- bind_rows(cohort_results, meta_summary) %>%
  filter(feature == "k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Oscillospiraceae|g__Flavonifractor|s__Flavonifractor_plautii") 

plot_data <- plot_data %>%
  mutate(lower = coef - 1.96 * stderr,
         upper = coef + 1.96 * stderr)

combined_data <- bind_rows(cohort_results, meta_summary) %>%
  filter(feature %in% meta_results_EO_heatmap$feature) %>%
  mutate(
    # 提取物种名
    feature_short = str_remove(
      str_extract(feature, 's__[A-Za-z0-9_]+'),
      '^s__'
    ),
    # 重命名 Cohort
    Cohort_clean = rename_cohorts(Cohort)
  )


individual_cohorts <- sort(setdiff(unique(combined_data$Cohort_clean), "Overall"))
final_col_order <- c(individual_cohorts, "Overall")

heatmap_df <- combined_data %>%
  select(feature_short, Cohort_clean, coef) %>%
  pivot_wider(names_from = Cohort_clean, values_from = coef) %>%
  column_to_rownames("feature_short") %>%
  select(all_of(final_col_order)) 

pval_df <- combined_data %>%
  select(feature_short, Cohort_clean, pval) %>%
  pivot_wider(names_from = Cohort_clean, values_from = pval) %>%
  column_to_rownames("feature_short") %>%
  select(all_of(final_col_order))

heatmap_mat <- as.matrix(heatmap_df)
pval_mat <- as.matrix(pval_df)

dim(heatmap_df)

country_colors_EO <- c(
  "Japan" = '#313695', "France" = "#fdae61", "Austria" = "#abd9e9", 
  "US&Canada" = "#74add1", "India" = '#e0f3f8', 
  "China" = '#fee090', "Turkey" = "#F781BF", "Germany" = '#4575b4', 
  "Italy" = '#f46d43', "Cezch" = '#d73027', "Pooled Analysis" = "#a50026"
)

current_countries <- cohort_to_country[colnames(heatmap_mat)]
col_ann <- HeatmapAnnotation(
  Country = current_countries,
  col = list(Country = country_colors_EO),
  show_annotation_name = TRUE, 
  annotation_name_side = "left"
)

col_names <- colnames(heatmap_mat)

col_split <- ifelse(col_names == "Overall", "Summary", "Cohorts")
col_split <- factor(col_split, levels = c("Cohorts", "Summary"))

col_fun = colorRamp2(c(-4, 0, 4), c("#313695", "white", "#a50026"))

sig_legend = Legend(
  labels = c("*   p < 0.05", "**  p < 0.01", "*** p < 0.001"),
  title = "Significance",
  type = "points", 
  pch = NA,             
  background = NA,      
  labels_gp = gpar(fontsize = 10, hjust = 0), 
  
  title_gp = gpar(fontsize = 10, fontface = "bold"),
)


target_species <- "Flavonifractor_plautii"
row_fontfaces <- ifelse(rownames(heatmap_mat) == target_species, "italic", "italic")
row_colors <- ifelse(rownames(heatmap_mat) == target_species,"black", "black")
row_fontsizes <- ifelse(rownames(heatmap_mat) == target_species, 8, 8)
ht_opt$HEATMAP_LEGEND_PADDING = unit(1.1, "cm")
ht_opt$ANNOTATION_LEGEND_PADDING = unit(1.1, "cm")

ht_EO <- Heatmap(heatmap_mat, 
              name = "Log Fold Change", 
              col = col_fun,
              top_annotation = col_ann,          
              column_split = col_split,          
              cluster_columns = TRUE,             
              column_gap = unit(2, "mm"),  
              cluster_column_slices = FALSE,      
              row_names_gp = gpar(
                fontface = row_fontfaces, 
                col = row_colors, 
                fontsize = row_fontsizes
              ),
              column_title = "Metagenomic Pooled analysis of EOCRC vs EOControl",
              column_names_rot = 45,
              column_names_gp = gpar(fontsize = 8,fontface ='bold'),
              cell_fun = function(j, i, x, y, width, height, fill) {
                p_val <- pval_mat[i, j]
                if(!is.na(p_val)) {
                  if(p_val < 0.001) {
                    grid.text("***", x, y, gp = gpar(fontsize = 10, fontface = "bold"))
                  } else if(p_val < 0.01) {
                    grid.text("**", x, y, gp = gpar(fontsize = 10, fontface = "bold"))
                  } else if(p_val < 0.05) {
                    grid.text("*", x, y, gp = gpar(fontsize = 10, fontface = "bold"))
                  }
                }
              },
              na_col = "#DDDDDD")

ht_EO
pdf(
  file = paste0(base_symbol,
                "downstream_analysis/output/picture/5.0_heatmap_EOCRCvsEOControl.pdf"),
  width = 13,   
  height = 6, 
  useDingbats = FALSE
)

draw(ht_EO, 
    annotation_legend_list = list(sig_legend), 
    merge_legend = TRUE,             
    heatmap_legend_side = "right",  
    annotation_legend_side = "right" 
)

dev.off()



cohort_results <- imap_dfr(fit_meta.crc_LO$maaslin_fits, function(df, name) {
  df %>% 
    filter(metadata == "Group") %>%  
    select(feature, coef, stderr, pval) %>%
    mutate(Cohort = name, Type = "Individual")
})

meta_summary <- fit_meta.crc_LO$meta_fits %>%
  select(feature, coef, stderr, pval = pval) %>% 
  mutate(Cohort = "Overall", Type = "Summary")

plot_data <- bind_rows(cohort_results, meta_summary) %>%
  filter(feature == "k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Oscillospiraceae|g__Flavonifractor|s__Flavonifractor_plautii")

plot_data <- plot_data %>%
  mutate(lower = coef - 1.96 * stderr,
         upper = coef + 1.96 * stderr)

combined_data <- bind_rows(cohort_results, meta_summary) %>%
  filter(feature %in% meta_results_LO_heatmap$feature) %>%
  mutate(

    feature_short = str_remove(
      str_extract(feature, 's__[A-Za-z0-9_]+'),
      '^s__'
    ),

    Cohort_clean = rename_cohorts(Cohort)
  )


individual_cohorts <- sort(setdiff(unique(combined_data$Cohort_clean), "Overall"))
final_col_order <- c(individual_cohorts, "Overall")

heatmap_df <- combined_data %>%
  select(feature_short, Cohort_clean, coef) %>%
  pivot_wider(names_from = Cohort_clean, values_from = coef) %>%
  column_to_rownames("feature_short") %>%
  select(all_of(final_col_order)) 


pval_df <- combined_data %>%
  select(feature_short, Cohort_clean, pval) %>%
  pivot_wider(names_from = Cohort_clean, values_from = pval) %>%
  column_to_rownames("feature_short") %>%
  select(all_of(final_col_order))

heatmap_mat <- as.matrix(heatmap_df)
pval_mat <- as.matrix(pval_df)

dim(heatmap_df)


country_colors_LO <- country_colors_EO

current_countries <- cohort_to_country[colnames(heatmap_mat)]
col_ann <- HeatmapAnnotation(
  Country = current_countries,
  col = list(Country = country_colors_LO),
  show_annotation_name = TRUE, 
  annotation_name_side = "left"
)

col_names <- colnames(heatmap_mat)

col_split <- ifelse(col_names == "Overall", "Summary", "Cohorts")
col_split <- factor(col_split, levels = c("Cohorts", "Summary"))

col_fun = colorRamp2(c(-4, 0, 4), c("#313695", "white", "#a50026"))

sig_legend = Legend(
  labels = c("*   p < 0.05", "**  p < 0.01", "*** p < 0.001"),
  title = "Significance",
  type = "points", 
  pch = NA,             
  background = NA,   
  labels_gp = gpar(fontsize = 10, hjust = 0), 
  title_gp = gpar(fontsize = 10, fontface = "bold")
)


target_species <- "Fusobacterium_nucleatum"
row_fontfaces <- ifelse(rownames(heatmap_mat) == target_species, "italic", "italic")
row_colors <- ifelse(rownames(heatmap_mat) == target_species, "black", "black")
row_fontsizes <- ifelse(rownames(heatmap_mat) == target_species, 8, 8)
ht_opt$HEATMAP_LEGEND_PADDING = unit(1.1, "cm")
ht_opt$ANNOTATION_LEGEND_PADDING = unit(1.1, "cm")

ht_LO <- Heatmap(heatmap_mat, 
              name = "Log Fold Change", 
              col = col_fun,
              top_annotation = col_ann,          
              column_split = col_split,            
              cluster_columns = TRUE,            
              column_gap = unit(2, "mm"),     
              cluster_column_slices = FALSE,      
              row_names_gp = gpar(
                fontface = row_fontfaces, 
                col = row_colors, 
                fontsize = row_fontsizes
              ),
              column_title = "Metagenomic Pooled Analysis of LOCRC vs LOControl",
              column_names_rot = 45,
              column_names_gp = gpar(fontsize = 8,fontface ='bold'),
              cell_fun = function(j, i, x, y, width, height, fill) {
                p_val <- pval_mat[i, j]
                if(!is.na(p_val)) {
                  if(p_val < 0.001) {
                    grid.text("***", x, y, gp = gpar(fontsize = 10, fontface = "bold"))
                  } else if(p_val < 0.01) {
                    grid.text("**", x, y, gp = gpar(fontsize = 10, fontface = "bold"))
                  } else if(p_val < 0.05) {
                    grid.text("*", x, y, gp = gpar(fontsize = 10, fontface = "bold"))
                  }
                }
              },
              na_col = "#DDDDDD")

ht_LO
dev.off()
pdf(
  file = paste0(base_symbol,
                "downstream_analysis/output/picture/5.0_heatmap_LOCRCvsLOControl.pdf"),
  width = 13,  
  height = 6, 
  useDingbats = FALSE
)

draw(ht_LO, 
     annotation_legend_list = list(sig_legend), 
     merge_legend = TRUE,           
     heatmap_legend_side = "right",  
     annotation_legend_side = "right" 
)

dev.off()
