library(tidyverse)

plot_data_discussion <- plot_data %>% group_by(Pathway) %>% 
  filter(
    any(P_interaction < 0.05, na.rm = TRUE) #|

  ) %>%
  ungroup() %>% pull(Pathway_Short)  %>% unique() 

merged_stratified_df_pattern <- paste(c(plot_data_discussion), collapse = "|")

metacyc_merged_global <- purrr::reduce(metacyc_list_stratified, full_join, by = "Feature") %>%
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
  column_to_rownames("Feature") 

col_sums_global <- colSums(metacyc_merged_global, na.rm = TRUE)
col_sums_global[col_sums_global == 0] <- 1
metacyc_norm_global <- sweep(metacyc_merged_global, 2, col_sums_global, "/")

target_pattern <- paste(plot_data_discussion, collapse = "|")
metacyc_subset <- metacyc_norm_global %>% 
  rownames_to_column('Feature') %>% 
  filter(str_starts(Feature,merged_stratified_df_pattern)) %>% 
  column_to_rownames('Feature')

common_samples <- intersect(colnames(metacyc_subset), rownames(meta_adj_batch_adjust))
metacyc_subset <- metacyc_subset[, common_samples, drop = FALSE]
meta_ready <- meta_adj_batch_adjust[common_samples, , drop = FALSE]

fit_metacyc_adj <- MMUPHin::adjust_batch(
  feature_abd = metacyc_subset, # 使用全正的 Log 数据
  batch = "batch",
  covariates = c("Group", "Age", "BMI", "Sex"),
  data = meta_ready,
  control = list(verbose = TRUE)
)

# metacyc_adj_log <- fit_metacyc_adj$feature_abd_adj
metacyc_abd_adj_adjust <- fit_metacyc_adj$feature_abd_adj %>%
  as.data.frame()  %>% rownames_to_column('Feature')

plot_data_long_adjust <- metacyc_abd_adj_adjust %>%
  filter(
    str_starts(Feature,merged_stratified_df_pattern)
    # str_starts(Feature,'PWY-7688')
  ) %>% 
  as.data.table() %>% 
  pivot_longer(-Feature, names_to = "Sample_ID", values_to = "Abundance") %>%
  mutate(Species = str_extract(Feature,'(s__[A-Za-z0-9-_|]+)|unclassified'),
         genus = str_extract(Feature,'(g__[A-Za-z0-9-_|]+)|unclassified'),
         Pathway = str_extract(Feature,'[A-Za-z0-9: -_]+'
         )
  ) %>% 
  inner_join(metadata_allcohorts, by = "Sample_ID")


plot_species_contribution <- function(data, target_pwy_id) {

  df_pwy <- data %>% filter(grepl(paste0("^", target_pwy_id), Pathway))
  if(nrow(df_pwy) == 0) return(NULL)

  top_species <- df_pwy %>%
    group_by(Species) %>% summarise(total = sum(Abundance)) %>%
    arrange(desc(total)) %>% slice_head(n=8) %>% pull(Species)

  df_summary <- df_pwy %>%
    mutate(Species_Plot = ifelse(Species %in% top_species, Species, "Others")) %>%
    mutate(Species_Label = str_remove(Species_Plot, ".*s__")) %>%
    group_by(Age_class, Group, Species_Label) %>%
    summarise(Mean_Abundance = mean(Abundance), .groups = "drop") %>%
    mutate(Plot_Group = paste(Age_class, Group, sep = "_")) %>%
    mutate(Plot_Group = factor(Plot_Group, levels = c("EO_Control", "EO_CRC", "LO_Control", "LO_CRC")))

  full_name <- unique(df_pwy$Pathway)[1]
  
  ggplot(df_summary, aes(x = Plot_Group, y = Mean_Abundance, fill = Species_Label)) +
    geom_col(position = "stack", width = 0.7, color = "black", linewidth = 0.1) +
    scale_fill_manual(values = c(brewer.pal(8, "Set2"), "grey80")) +
    theme_classic() +
    labs(title = full_name %>%
           gsub("^([0-9])", "EC \\1", .) %>%
           str_replace(":.*$", ""),
         y = "Mean Relative Abundance (CPM/RelAb)",
         x = NULL,
         fill = "Contributor") +
    
    scale_x_discrete(labels = c("EOControl", "EOCRC", "LOControl", "LOCRC")) +
    
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold",size =14),
      axis.text.x = element_text(size = 11, face = "bold"),
      legend.position = "right"
    )
  
  
}


out_dir_stack_adjust <- "D:/R_workspace/global_EOCRC_microbiology/downstream_analysis/output/picture/21.metacyc_species_contribution_adjust/"
if(!dir.exists(out_dir_stack_adjust)) dir.create(out_dir_stack_adjust, recursive = T)

for (pwy_id in plot_data_discussion) {
  p <- plot_species_contribution(plot_data_long_adjust, pwy_id)
  
  if (!is.null(p)) {
    safe_name <- str_replace_all(pwy_id, "[^A-Za-z0-9]", "_")
    ggsave(paste0(out_dir_stack_adjust, "Stack_", safe_name, ".pdf"), p, 
           width = 8, height = 5)
    print(paste("Saved:", pwy_id))
  }
}


