library(dplyr)
library(purrr)
library(data.table)

meta_list_species <- lapply(meta_list_metadata_intersect_samples, function(df) {
  df_species <- df %>%
    filter(grepl("s__", clade_name) & !grepl("t__", clade_name)) %>%
    mutate(clade_name = as.character(clade_name))
  
  return(df_species)
})

abundance_combined <- meta_list_species %>%
  purrr::reduce(full_join, by = "clade_name") %>%
  mutate(across(everything(), ~ replace_na(.x, 0)))

ml_input_data_raw <- as.data.frame(abundance_combined) 
head(colSums(ml_input_data_raw %>% column_to_rownames('clade_name')))

write.csv( ml_input_data_raw,
  file = paste0(model_dir,'0.ml_input_data_raw.csv'),row.names = F
)

clade_names <- abundance_combined$clade_name
abundance_matrix <- as.matrix(abundance_combined[,-1])
rownames(abundance_matrix) <- clade_names

abundance_transformed <- asin(sqrt(abundance_matrix / 100))

ml_input_data_asin <- as.data.frame(abundance_transformed) %>% 
  rownames_to_column('clade_name')
  
head(colSums(ml_input_data_asin %>% column_to_rownames('clade_name')))

write.csv( ml_input_data_asin,
           file = paste0(model_dir,'0.ml_input_data_asin.csv'),row.names = F
)

abundance_matrix_prop <- abundance_matrix / 100
min(abundance_matrix_prop[abundance_matrix_prop > 0], na.rm = TRUE)
pseudo <- 1e-8

abundance_log <- log10(abundance_matrix_prop + pseudo)
min(abundance_log, na.rm = TRUE)

ml_input_data_log <- as.data.frame(abundance_log) %>%
  rownames_to_column('clade_name')

summary(colSums(ml_input_data_log %>% column_to_rownames('clade_name')))

write.csv(
  ml_input_data_log,
  file = paste0(model_dir,'0.ml_input_data_log.csv'),
  row.names = FALSE
)

run_lodo_nested_meta <- function(meta_list, 
                                 metadata, 
                                 out_dir_features = "LODO_features", 
                                 out_dir_meta = "LODO_meta_results",
                                 min_prevalence = 0.05,
                                 min_train_cohorts = 3,
                                 abund_threshold = 0.001) {

  if(!dir.exists(out_dir_features)) dir.create(out_dir_features, recursive = TRUE)
  if(!dir.exists(out_dir_meta)) dir.create(out_dir_meta, recursive = TRUE)

  all_cohorts <- unique(metadata$Cohort)
  n_cohorts <- length(all_cohorts)
  
  lodo_meta_objects <- list()
  
  for (i in seq_along(all_cohorts)) {
    test_cohort <- all_cohorts[i]
    train_cohorts <- setdiff(all_cohorts, test_cohort)
    
    train_meta_EO <- world_meta_EO %>% filter(!Cohort %in% test_cohort)
    train_meta_LO <- world_meta_LO %>% filter(!Cohort %in% test_cohort)
    
    message(sprintf("\n=== LODO Iteration %d/%d: Test Cohort = [%s] ===", i, n_cohorts, test_cohort))

    train_metadata <- metadata %>% filter(Cohort %in% train_cohorts)
    train_samples <- train_metadata$Sample_ID

    train_abd_list <- lapply(train_cohorts, function(c_name) {
      df <- meta_list[[c_name]]
      target_cols <- intersect(colnames(df), train_samples)
      if(length(target_cols) > 0) {
        df %>% select(clade_name, any_of(target_cols))
      } else {
        NULL
      }
    })

    train_abd_merged <- train_abd_list %>% 
      compact() %>% 
      reduce(function(x, y) full_join(x, y, by = "clade_name")) %>%
      column_to_rownames("clade_name") %>%
      mutate(across(everything(), ~replace_na(as.numeric(.), 0)))

    consensus_res <- run_consensus_EO_LO_combine(
      meta_list =  meta_list_metadata_intersect_samples, 
      meta_EO = train_meta_EO,  
      meta_LO = train_meta_LO,
      min_cohorts = min_train_cohorts,
      min_prevalence = min_prevalence
    )
    
    passed_features <- consensus_res$merged$clade_name

    feature_file <- file.path(out_dir_features, paste0("Features_Test_", test_cohort, ".csv"))
    write.csv(data.frame(Feature = passed_features 
                         ), 
              file = feature_file, row.names = FALSE)
    

    final_train_abd <- train_abd_merged[passed_features, , drop = FALSE]

    cl_train_input <- train_metadata %>% 
      column_to_rownames('Sample_ID') %>% 
      mutate(Cohort = as.factor(Cohort), Age = as.numeric(Age))

    final_train_abd <- final_train_abd[, rownames(cl_train_input), drop = FALSE]
    
    final_train_abd <- apply_step2_masking(
      feature_abd = final_train_abd, 
      meta_data = cl_train_input, 
      batch_col = "Cohort"
    )
    

    fit_meta <- tryCatch({
      MMUPHin::lm_meta(
        feature_abd = final_train_abd,
        exposure = 'Group',
        batch = 'Cohort',
        covariates = c("Age", "Sex", 'BMI'),
        control = list(rma_method = "HS", transform = "LOG", normalization = 'NONE', 
                       min_prevalence = 0, min_abundance = 0),
        data = cl_train_input
      )
    }, error = function(e) {

      return(NULL)
    })

    if(!is.null(fit_meta)) {
      meta_res_df <- fit_meta$meta_fits %>% 
        rownames_to_column("Feature") %>%
        filter(!is.na(pval))

      sig_features <- meta_res_df %>% filter(pval < 0.05)

      res_file <- file.path(out_dir_meta, paste0("MetaRes_Test_", test_cohort, ".csv"))
      write.csv(sig_features, file = res_file, row.names = FALSE)
      
      meta_res_df <-  process_meta_with_rescue(
        meta_res_obj = fit_meta,
        group_name = "CRC_EO",
        i2_tier1 = 50,      
        i2_tier2 = 80,      
        weight_limit = 40,
        base_consistency = 0.6,  
        rescue_consistency  = 0.75 
      ) %>% 
        filter(str_starts(qc_status,'Pass')) %>%
        mutate( qval.fdr = p.adjust(pval, method = "fdr") ) 

      sig_features <- meta_res_df %>% filter(qval.fdr < 0.1)

      res_file <- file.path(out_dir_meta, paste0("Meta_qval_", test_cohort, ".csv"))
      write.csv(sig_features, file = res_file, row.names = FALSE)

      lodo_meta_objects[[test_cohort]] <- fit_meta
    } else {
      lodo_meta_objects[[test_cohort]] <- NA
    }
  }

  return(lodo_meta_objects)
}

current_list_names <- names(meta_list_metadata_intersect_samples)

new_list_names <- case_when(
  current_list_names == 'A10' ~ 'YachidaS_2019',
  current_list_names == 'A12' ~ 'Zeller_2014',
  current_list_names == 'A13' ~ 'Feng_2015',
  current_list_names == 'A14' ~ 'Vogtmann_2016',
  current_list_names == 'A17' ~ 'Hannigan_2018',
  current_list_names == 'A18' ~ 'Gupta_2019',

  current_list_names %in% c('Public_study__YuJ_2015', 'A27') ~ 'Yu_2017', 
  
  current_list_names == 'A28' ~ 'Liu_2022',

  current_list_names %in% c('A30_FUSCC', 'A30') ~ 'FUSCC-SHSD',
  
  current_list_names == 'A33' ~ 'Yang_2020',
  current_list_names == 'A37' ~ 'Piccinno_2025_IIGM_TU',
  current_list_names == 'A9' ~ 'Wirbel_2019',
  current_list_names == 'A15' ~ 'ThomasAM_2018b',
  current_list_names == 'A40' ~ 'Piccinno_2025_IIGM_CZ',
  current_list_names == 'A41' ~ 'Piccinno_2025_IIGM_IT',

  TRUE ~ current_list_names 
)

names(meta_list_metadata_intersect_samples) <- new_list_names

meta_list_metadata_intersect_samples_LODO <- lapply(meta_list_metadata_intersect_samples,function(x){
  x <- x %>% filter(str_detect(clade_name, "\\|s__[^|]+$"))
})
  
lodo_results_EO <- run_lodo_nested_meta(
  meta_list = meta_list_metadata_intersect_samples_LODO,
  metadata = world_meta_EO,
  out_dir_features = paste0(base_symbol,"downstream_analysis/model/LODO_Input_Features_EO"),
  out_dir_meta = paste0(base_symbol,"downstream_analysis/model/LODO_Meta_Sig_EO"),
  min_prevalence = 0.05,
  min_train_cohorts = 3
)


lodo_results_LO <- run_lodo_nested_meta(
  meta_list = meta_list_metadata_intersect_samples_LODO,
  metadata = world_meta_LO,
  out_dir_features = paste0(base_symbol,"downstream_analysis/model/LODO_Input_Features_LO"),
  out_dir_meta = paste0(base_symbol,"downstream_analysis/model/LODO_Meta_Sig_LO"),
  min_prevalence = 0.05,
  min_train_cohorts = 3
)
