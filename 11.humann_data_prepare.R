## humann_results
#metacyc

metacyc_list <- lapply(dir_list, function(i) {
  file_path <- file.path(read_dir, i, "all_samples_pathabundance_unstratified_relab.tsv")
  
  if(file.exists(file_path)){
    dt <- fread(file_path, header = TRUE)
    names(dt) <- str_remove(names(dt), '_Abundance')
    names(dt)[1] <- "Feature"
    if (i %in% c('A18', 'A17')){
      colnames(dt)[-1] <- paste(i,colnames(dt)[-1],sep = '_')
    }
    inter_samples <- intersect(metadata_allcohorts$Sample_ID,
                               colnames(dt))
    dt <- dt %>% select(any_of(c("Feature",inter_samples)))
    return(dt)
  } else {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
})

names(metacyc_list) <- dir_list


run_consensus_pathway_EO_LO <- function(meta_list, 
                                        meta_EO, 
                                        meta_LO, 
                                        min_prevalence = 0.1,    
                                        min_cohorts = 3,   
                                        MIN_STEP1_COUNT =3,
                                        abund_threshold = 1e-6,  
                                        clean_names = TRUE) {    
  
  get_passed_features <- function(cohort_name, df, target_metadata, group_tag) {
    
    sample_cols <- df %>% 
      column_to_rownames(var = colnames(df)[1]) %>% 
      select(any_of(target_metadata$Sample_ID))
    
    n_samples <- ncol(sample_cols)

    if (n_samples == 0) {
      message(sprintf("   [%s] cohort [%s]: no samples", group_tag, cohort_name))
      return(character(0))
    }
    
    if (clean_names) {
      bad_rows <- c("UNMAPPED", "UNINTEGRATED")
      sample_cols <- sample_cols[!rownames(sample_cols) %in% bad_rows, , drop = FALSE]
    }
    
    min_detected <- MIN_STEP1_COUNT
    
    n_detected <- rowSums(sample_cols > abund_threshold)

    min_detected <- MIN_STEP1_COUNT
    
    passed_features <- rownames(sample_cols)[
      (n_detected >= min_detected) & 
        (n_detected / n_samples >= min_prevalence)
    ]

    message(sprintf("   [%s] cohort [%s]: N=%d, all features=%d, satisfied features=%d", 
                    group_tag, 
                    cohort_name, 
                    n_samples,        
                    nrow(sample_cols),  
                    length(passed_features)))
    
    return(passed_features)
  }

  cohort_votes_EO <- lapply(names(meta_list), function(n) {
    get_passed_features(n, meta_list[[n]], meta_EO, "EO")
  })
  
  cohort_votes_LO <- lapply(names(meta_list), function(n) {
    get_passed_features(n, meta_list[[n]], meta_LO, "LO")
  })
  

  all_votes_EO <- unlist(cohort_votes_EO)
  keep_EO <- names(table(all_votes_EO))[table(all_votes_EO) >= min_cohorts]
  
  all_votes_LO <- unlist(cohort_votes_LO)
  keep_LO <- names(table(all_votes_LO))[table(all_votes_LO) >= min_cohorts]
  
  final_keep_features <- unique(c(keep_EO, keep_LO))

  merged_df <- meta_list %>%
    reduce(function(x, y) full_join(x, y, by = colnames(x)[1])) %>% 
    filter(!!sym(colnames(.)[1]) %in% final_keep_features) %>%
    mutate(across(where(is.numeric), ~replace_na(., 0)))
 
  feat_col <- colnames(merged_df)[1]
  
  EO_sub <- merged_df %>% select(all_of(feat_col), any_of(meta_EO$Sample_ID))
  LO_sub <- merged_df %>% select(all_of(feat_col), any_of(meta_LO$Sample_ID))
  
  return(list(merged = merged_df, EO = EO_sub, LO = LO_sub))
}


apply_pathway_masking <- function(feature_abd, meta_data, batch_col = "Cohort", min_count = 5) {

  feature_abd_masked <- feature_abd
  cohorts <- unique(meta_data[[batch_col]])
  total_masked_entries <- 0
  
  for (cohort in cohorts) {
    s_ids <- rownames(meta_data)[meta_data[[batch_col]] == cohort]
    s_ids <- intersect(s_ids, colnames(feature_abd_masked))
    
    if(length(s_ids) == 0) next

    sub_dat <- feature_abd_masked[, s_ids, drop=FALSE]

    counts <- rowSums(sub_dat > 0)

    risky_features <- names(counts)[counts > 0 & counts < min_count]
    
    if(length(risky_features) > 0) {
      feature_abd_masked[risky_features, s_ids] <- 0
      total_masked_entries <- total_masked_entries + length(risky_features)
    }
  }

  return(feature_abd_masked)
}

