library(MMUPHin)
library(metafor)
library(readxl)
library(tidyverse)
library(data.table)

disk_symbol <- "J:/"
base_symbol <- "/home/global_EOCRC_microbiology/"


#metadata_allcohorts
#meta_list_metadata_intersect_samples
read_dir <- paste0(base_symbol,"data")
meta_list <- lapply(dir_list, function(i) {
  fread(file.path(read_dir, i, "/merged_metaphlan_taxonomic_profile.tsv"), header = T,
        skip = "#mpa_vJun23_CHOCOPhlAnSGB_202403")
})
names(meta_list) <- dir_list

colnames(meta_list[["A17"]])[-1] <- paste0('A17_',colnames(meta_list[["A17"]])[-1])
colnames(meta_list[["A18"]])[-1] <- paste0('A18_',colnames(meta_list[["A18"]])[-1])


meta_list_metadata_intersect_samples <- meta_list

{
  make_dirs <- function(paths) {
    for (p in paths) {
      dir.create(p, recursive = TRUE, showWarnings = FALSE)
    }
  }
  

  base_down <- file.path(base_symbol, "downstream_analysis")
  
  dirs_basic <- file.path(
    base_down,
    c("", "generated_data", "output", "output/picture", "output/table")
  )
  

  world_base <- file.path(base_down, "meta_analysis_world")
  
  groups <- c("EOCRC", "LOCRC", "interCRC")
  modes  <- c("adjust", "crude")
  
  dirs_world <- c(
    world_base,
    file.path(world_base, groups),
    as.vector(outer(groups, modes,
                    function(g, m) file.path(world_base, g, m)))
  )

  continent_base <- file.path(base_down, "meta_analysis_continent")
  
  continents <- c("Europe", "East_Asia", "North_America")
  groups2    <- c("EOCRC", "LOCRC")
  
  dirs_continent <- c(
    continent_base,
    file.path(continent_base, continents),
    as.vector(outer(continents, groups2,
                    function(cn, g)
                      file.path(continent_base, cn, g)))
  )
  

  humann_base <- file.path(base_down, "output", "humann_meta")
  
  dirs_humann_metacyc <- c(
    humann_base,
    file.path(humann_base, c("EOCRC", "LOCRC")),
    file.path(humann_base, "EOCRC", "metacyc"),
    file.path(humann_base, "LOCRC", "metacyc")
  )

  dirs_humann_ec <- file.path(
    humann_base,
    c("EOCRC/ec","EOCRC/ec_crude",
      "LOCRC/ec","LOCRC/ec_crude")
  )
  
  lapply(c(dirs_basic,dirs_world,dirs_continent,
           dirs_humann_metacyc,dirs_humann_ec),
         make_dirs)
  
}

run_consensus_EO_LO_combine <- function(meta_list, 
                                        meta_EO,  
                                        meta_LO,  
                                        min_prevalence = 0.05, 
                                        min_cohorts = 2, 
                                        MIN_STEP1_COUNT =3,
                                        abund_threshold = 1e-3, 
                                        renormalize = TRUE) {

  get_passed_taxa <- function(cohort_name, full_df, target_samples) {
    sample_cols <- full_df %>% 
      filter(str_detect(clade_name, "\\|s__[^|]+$")) %>% 
      column_to_rownames('clade_name') %>% 
      select(any_of(target_samples)) %>% 
      mutate(across(everything(), as.numeric))
    
    n_samples <- ncol(sample_cols)
    
    if (n_samples == 0) {
      message(sprintf("   cohort [%s]: no samples", cohort_name))
      return(character(0))
    }

    n_detected <- rowSums(sample_cols > abund_threshold)

    min_detected <- MIN_STEP1_COUNT

    
    passed_taxa <- rownames(sample_cols)[
      (n_detected >= min_detected) & 
        (n_detected / n_samples >= min_prevalence)
    ]
    
    message(sprintf("   cohort [%s] (N=%d): number of species=%d (satisfy >=%d samples & >=%.0f%%)", 
                    cohort_name, n_samples, length(passed_taxa),min_detected, min_prevalence*100))
    
    return(passed_taxa)
  }
  
  cohort_votes_EO <- lapply(names(meta_list), function(cn) {
    get_passed_taxa(cn, meta_list[[cn]], meta_EO$Sample_ID)
  })
  
  cohort_votes_LO <- lapply(names(meta_list), function(cn) {
    get_passed_taxa(cn, meta_list[[cn]], meta_LO$Sample_ID)
  })
  
  # ---(Consensus Counting) ---
  all_votes_EO <- unlist(cohort_votes_EO)
  keep_taxa_EO <- names(table(all_votes_EO))[table(all_votes_EO) >= min_cohorts]
  
  all_votes_LO <- unlist(cohort_votes_LO)
  keep_taxa_LO <- names(table(all_votes_LO))[table(all_votes_LO) >= min_cohorts]
  
  keep_taxa_final <- unique(c(keep_taxa_EO, keep_taxa_LO))
  
  message(sprintf("\n>>> consensus results: retain %d species (EO contribution=%d, LO contribution=%d)", 
                  length(keep_taxa_final), length(keep_taxa_EO), length(keep_taxa_LO)))

  merged_df <- meta_list %>%
    reduce(function(x, y) full_join(x, y, by = "clade_name")) %>%
    filter(clade_name %in% keep_taxa_final) %>%
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  if (renormalize) {
    taxa_names <- merged_df$clade_name
    mat <- merged_df %>% select(-clade_name) %>% as.matrix()
    col_sums <- colSums(mat)
    col_sums[col_sums == 0] <- 1 
    mat_norm <- sweep(mat, 2, col_sums, "/") * 100
    merged_df <- as.data.frame(mat_norm) %>% mutate(clade_name = taxa_names) %>% select(clade_name, everything())
  }
  
  EO_sub <- merged_df %>% select(clade_name, any_of(meta_EO$Sample_ID))
  LO_sub <- merged_df %>% select(clade_name, any_of(meta_LO$Sample_ID))
  
  return(list(merged = merged_df, EO = EO_sub, LO = LO_sub))
}

