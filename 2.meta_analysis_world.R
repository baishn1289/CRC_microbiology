library(MMUPHin)
library(metafor)
library(readxl)
library(tidyverse)
library(data.table)
library(patchwork)
library(Maaslin2)
library(ggrastr)

cohorts_heatmap <- metadata_allcohorts %>% 
  filter(Group %in% c( 'CRC', 'Control')) %>% 
  group_by(Cohort,Group,Age_class) %>% 
  reframe(Count = n()) 

total_col <- cohorts_heatmap %>%
  group_by(Group) %>%
  summarise(Count = sum(Count), .groups = "drop") %>%
  mutate(Cohort = "Total")

cohort_levels <- c(
  "YachidaS_2019", "Zeller_2014", "Feng_2015", "Vogtmann_2016",
  "Hannigan_2018",
  "Gupta_2019","Yu_2017","Liu_2022",
  "FUSCC-SHSD","Yang_2020",
  "Piccinno_2025_IIGM_TU","Wirbel_2019",
  "ThomasAM_2018b",
  "Piccinno_2025_IIGM_CZ",
  "Piccinno_2025_IIGM_IT",
  "Total"
)

cohorts_heatmap <- bind_rows(cohorts_heatmap, total_col) %>% 
  mutate(Cohort =  factor(Cohort, levels = cohort_levels),
         Group = factor(Group,levels = c("EOCRC", "EOControl", "LOCRC", "LOControl")))

logic_df <- cohorts_heatmap %>%
  filter(Cohort != "Total") %>% 
  select(Cohort, Group, Count) %>%
  pivot_wider(names_from = Group, values_from = Count,
              values_fill = 0  ) %>%
  mutate(
    EO_Status = ifelse(EOControl >= 3 & `EOCRC` >= 3, "Yes", "No"),
    LO_Status = ifelse(LOControl >= 3 & LOCRC >= 3, "Yes", "No"),
    EO_LO_Status = ifelse(EOCRC >= 3 & LOCRC >= 3, "Yes", "No"),
  ) %>%
  select(Cohort, EO_Status, LO_Status,EO_LO_Status
         ) %>%
  pivot_longer(cols = c("EO_Status", "LO_Status"#,'EO_LO_Status'
                        ), 
               names_to = "Meta_Type", 
               values_to = "Is_Included") %>% 
  mutate(Meta_Type = factor(Meta_Type, 
                             levels = c("LO_Status", "EO_Status"#,'EO_LO_Status'
                             ), 
                             labels = c("LO pooled-analysis", "EO pooled-analysis"#,'CRC_pooled-analysis'
                             )))


logic_df <- bind_rows(
  logic_df,
  expand_grid(
    Cohort   = "Total",
    Meta_Type = levels(logic_df$Meta_Type)
  ) %>%
    mutate(Is_Included = NA)
)

logic_df <- logic_df %>%
  mutate(Cohort =  factor(Cohort,levels = cohort_levels)) %>% 
  mutate(Is_Included = ifelse(Cohort == "Total", NA, Is_Included))

total_x <- which(cohort_levels == "Total")

p_main <- ggplot(cohorts_heatmap,
                 aes(x = Cohort, y = Group, fill = Count)) +
  geom_tile(color = "white", linewidth = 0.5) +  
  geom_text(aes(label = Count), color = "black", size = 4) + 
  geom_vline(xintercept = total_x - 0.5,
             colour = "white", linewidth = 1.4) + 
  scale_x_discrete(drop = FALSE)+
  scale_fill_gradient(low = "#E1F5FE", high = "#01579B", trans = "log10", 
                      name = "Sample Size\n(Log Scale)") + 
  theme_minimal() +
  theme(
    axis.title.x = element_blank(), 
    axis.text.x = element_blank(),  
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),  
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    legend.position = "right"
  ) +
  labs(title = "Cohort Sample Sizes & pooled-analysis Eligibility")

p_bottom <- ggplot(logic_df, aes(x = Cohort, y = Meta_Type, fill = Is_Included)) +
  geom_tile(color = "white", size = 0.5) +
  geom_vline(xintercept = total_x - 0.5,
             colour = "white", linewidth = 1.4) + 
  scale_x_discrete(drop = FALSE)+
  scale_fill_manual(
    values = c("Yes" = "#337C38", "No" = "#E0E0E0"),
    na.value = "white",
    breaks = c("Yes", "No"),
    drop = FALSE,
    name = "Included in\npooled-analysis",
    guide = guide_legend(na.translate = FALSE)
  )+
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1), # 显示队列名
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    legend.position = "right"
  )

heatmap_plot <- p_main / p_bottom + 
  plot_layout(heights = c(4, 1.8), guides = "collect") # guides='collect'合并图例(可选)

print(heatmap_plot)

ggsave(plot = heatmap_plot,width = 12,height = 5,
       filename = paste0(base_symbol,
                         "downstream_analysis/output/picture/1.cohorts_pooled-analysis_survey.pdf")
)

process_meta_with_rescue <- function(meta_res_obj, 
                                     group_name = "CRC", 
                                     i2_tier1 = 50,     
                                     i2_tier2 = 75,   
                                     weight_limit = 40,  
                                     base_consistency = 0.6,
                                     rescue_consistency = 0.75 
) {

  df_meta <- meta_res_obj$meta_fits %>%
    filter(!is.na(coef)) %>% 
    rownames_to_column("Pathway")
  
  df_individual <- bind_rows(meta_res_obj$maaslin_fits) %>%
    select(feature, coef, batch) %>%
    rename(Pathway = feature, cohort_coef = coef, cohort_id = batch) %>%
    filter(!is.na(cohort_coef)) 

  consistency_stats <- df_meta %>% select(Pathway, global_coef = coef) %>%
    left_join(df_individual, by = "Pathway") %>%
    group_by(Pathway) %>%
    filter(cohort_coef != 0 ) %>% 
    summarise(
      n_total_cohorts = n_distinct(cohort_id),
      n_same_dir = sum(sign(cohort_coef) == sign(global_coef), 
                       na.rm = TRUE),
      consistency_score = n_same_dir / n_total_cohorts,
      
      n_opp_dir = n_total_cohorts - n_same_dir,
      pass_min_cohort = n_total_cohorts >= 3,
      .groups = "drop"
    )

  df_final <- df_meta %>%
    left_join(consistency_stats, by = "Pathway") %>%
    rowwise() %>%
    mutate(
      max_weight = ifelse(all(is.na(c_across(starts_with("weight_")))), NA,
                          max(c_across(starts_with("weight_")), na.rm = TRUE)),
      pass_weight = !is.na(max_weight) & max_weight < weight_limit,
      
      is_tier1 = I2 < i2_tier1,
      is_in_rescue_zone = (I2 >= i2_tier1) & (I2 < i2_tier2),

      pass_base_consistency = !is.na(consistency_score) & (consistency_score >= base_consistency),
      pass_rescue_consistency = !is.na(consistency_score) & (consistency_score >= rescue_consistency),

      pass_qc = case_when(
        !pass_weight ~ FALSE,
        !pass_min_cohort ~ FALSE,

        is_tier1 & pass_base_consistency ~ TRUE,
        is_tier1 & !pass_base_consistency ~ FALSE, 
        is_in_rescue_zone & pass_rescue_consistency ~ TRUE,
        
        TRUE ~ FALSE
      ),
      
      qc_status = case_when(
        !pass_weight ~ "Fail_High_Weight",
        !pass_min_cohort ~ "Fail_Too_Few_Cohorts",
        is_tier1 & pass_base_consistency ~ "Pass_Tier1_Robust",
        is_tier1 & !pass_base_consistency ~ "Fail_Tier1_Inconsistent_Direction",
        is_in_rescue_zone & pass_rescue_consistency ~ "Pass_Tier2_Rescued",
        is_in_rescue_zone & !pass_rescue_consistency ~ "Fail_Tier2_Inconsistent",
        I2 >= i2_tier2 ~ "Fail_High_Heterogeneity",
        TRUE ~ "Fail_Other"
      )
    ) %>%
    ungroup()
  
  return(df_final)
}

apply_step2_masking <- function(feature_abd, meta_data, batch_col = "Cohort") {

  MIN_MODEL_COUNT <- 5 
  
  feature_abd_masked <- feature_abd
  cohorts <- unique(meta_data[[batch_col]])
  
  total_masked <- 0
  
  for (cohort in cohorts) {

    s_ids <- rownames(meta_data)[meta_data[[batch_col]] == cohort]
    s_ids <- intersect(s_ids, colnames(feature_abd_masked))
    
    if(length(s_ids) == 0) next

    sub_dat <- feature_abd_masked[, s_ids, drop=FALSE]

    counts <- rowSums(sub_dat > 0)

    risky_taxa <- names(counts)[counts > 0 & counts < MIN_MODEL_COUNT]
    
    if(length(risky_taxa) > 0) {

      feature_abd_masked[risky_taxa, s_ids] <- 0
      total_masked <- total_masked + length(risky_taxa)
    }
  }
  
  return(feature_abd_masked)
}


table(metadata_allcohorts_EO$Cohort,metadata_allcohorts_EO$Group)  

world_meta_EO <- metadata_allcohorts_EO %>% 
  filter(Group %in% c('CRC' ,'Control' ),
         Cohort %in% as.character(logic_df[logic_df$Is_Included =='Yes' & 
                                             logic_df$Meta_Type == 'EO pooled-analysis',]$Cohort)) 
unique(world_meta_EO$Cohort)

world_meta_LO <- metadata_allcohorts_LO  %>% 
  filter(Group %in% c('CRC' ,'Control' ),
         Cohort %in% as.character(logic_df[logic_df$Is_Included =='Yes' & 
                                             logic_df$Meta_Type == 'LO pooled-analysis',]$Cohort)) 

table(world_meta_EO$Cohort,world_meta_EO$Group,is.na(world_meta_EO$BMI))
table(world_meta_EO$Cohort,world_meta_EO$Group,is.na(world_meta_EO$Sex))
table(world_meta_LO$Cohort,world_meta_LO$Group,is.na(world_meta_LO$Sex))
table(world_meta_LO$Cohort,world_meta_LO$Group,is.na(world_meta_LO$BMI))

abundance_world_EO_LO_combine <- run_consensus_EO_LO_combine(
  meta_list = meta_list_metadata_intersect_samples,   
  meta_EO = world_meta_EO,  
  meta_LO = world_meta_LO,min_cohorts = 3,min_prevalence =0.05,
  abund_threshold = 0.001
)

merged_df <- abundance_world_EO_LO_combine$merged   
dim(merged_df)


metadata_allcohorts_EO = metadata_allcohorts %>% filter(Age_class == 'EO')
metadata_allcohorts_LO = metadata_allcohorts %>% filter(Age_class == 'LO')


#EOCRC vs EOCTR 
cl_metainput_EO = world_meta_EO %>%  
  column_to_rownames('Sample_ID') %>% mutate(Cohort = as.factor(Cohort),
                                             Age = as.numeric(Age))

feat_EO <- abundance_world_EO_LO_combine$EO %>% 
  column_to_rownames("clade_name") %>% 
  as.matrix() %>% 
  apply_step2_masking( meta_data = cl_metainput_EO,
                       batch_col = "Cohort")

table(cl_metainput_EO$Cohort,cl_metainput_EO$Group) 

abundance_metainput_EO = feat_EO %>% 
  as.data.frame() %>% 
  select(any_of(c(rownames(cl_metainput_EO)))) 

setdiff(names(abundance_metainput_EO),rownames(cl_metainput_EO))
setdiff(rownames(cl_metainput_EO),names(abundance_metainput_EO))

dim(abundance_metainput_EO)
dim(cl_metainput_EO)


fit_meta.crc_EO <- MMUPHin::lm_meta(feature_abd = abundance_metainput_EO,
                                    exposure = 'Group',
                                    batch = 'Cohort',
                                    covariates = c("Age","Sex",'BMI'),
                                    control = list(rma_method="HS", transform="LOG", normalization = 'NONE',
                                                   min_prevalence = 0,min_abundance = 0.001,
                                                   output= paste0(base_symbol,
                                                                  "downstream_analysis/meta_analysis_world/EOCRC/adjust")
                                    ),
                                    data = cl_metainput_EO )


#LOCRC vs LOCTR 
table(metadata_allcohorts_LO$Cohort,metadata_allcohorts_LO$Group)  

cl_metainput_LO = world_meta_LO %>% 
  column_to_rownames('Sample_ID') %>% mutate(Cohort = as.factor(Cohort),
                                             Age = as.numeric(Age))

feat_LO <- abundance_world_EO_LO_combine$LO %>% 
  column_to_rownames("clade_name") %>% 
  as.matrix() %>% 
  apply_step2_masking( meta_data = cl_metainput_LO,
                       batch_col = "Cohort") 

abundance_metainput_LO = feat_LO %>% 
  as.data.frame() %>% 
  select(any_of(c(rownames(cl_metainput_LO)))) 

setdiff(names(abundance_metainput_LO),rownames(cl_metainput_LO))
setdiff(rownames(cl_metainput_LO),names(abundance_metainput_LO))

fit_meta.crc_LO <- MMUPHin::lm_meta(feature_abd = abundance_metainput_LO,
                                    exposure = 'Group',
                                    batch = 'Cohort',
                                    covariates = c("Age","Sex",'BMI'),
                                    control = list(rma_method="HS", transform="LOG", normalization = 'NONE',
                                                   min_prevalence = 0,min_abundance = 0.001,
                                                   output= paste0(base_symbol,
                                                                  "downstream_analysis/meta_analysis_world/LOCRC/adjust")
                                    ),
                                    data = cl_metainput_LO )

#sensitivity analysis
world_meta_EO_sensi <- world_meta_EO %>% 
  filter(!Cohort %in% c("Feng_2015",'ThomasAM_2018b')) 

unique(world_meta_EO_sensi$Cohort)

world_meta_LO_sensi <- world_meta_LO

abundance_world_EO_LO_combine_sensi <- run_consensus_EO_LO_combine(
  meta_list = meta_list_metadata_intersect_samples,   
  meta_EO = world_meta_EO_sensi,  
  meta_LO = world_meta_LO_sensi,
  min_cohorts = 3,min_prevalence =0.05,
  abund_threshold = 0.001
)

cl_metainput_EO_sensi = world_meta_EO_sensi %>% 
  column_to_rownames('Sample_ID') %>% mutate(Cohort = as.factor(Cohort),
                                             Age = as.numeric(Age))

feat_EO_sensi <- abundance_world_EO_LO_combine_sensi$EO %>% 
  column_to_rownames("clade_name") %>% 
  as.matrix() %>% 
  apply_step2_masking( meta_data = cl_metainput_EO_sensi,
                       batch_col = "Cohort") 

abundance_metainput_EO_sensi = feat_EO_sensi %>% 
  as.data.frame() %>% 
  select(any_of(c(rownames(cl_metainput_EO)))) 

setdiff(names(abundance_metainput_EO_sensi),rownames(cl_metainput_EO_sensi))
setdiff(rownames(cl_metainput_EO_sensi),names(abundance_metainput_EO_sensi))

fit_meta.crc_EO_sensi <- MMUPHin::lm_meta(feature_abd = abundance_metainput_EO_sensi,
                                    exposure = 'Group',
                                    batch = 'Cohort',
                                    covariates = c("Age","Sex",'BMI'),
                                    control = list(rma_method="HS", transform="LOG", normalization = 'NONE',
                                                   min_prevalence = 0,min_abundance = 0,
                                                   output= paste0(base_symbol,
                                                                  "downstream_analysis/meta_analysis_world/EOCRC/adjust_sensi")
                                    ),
                                    data = cl_metainput_EO_sensi )


fit_meta.crc_LO_sensi <- fit_meta.crc_LO

