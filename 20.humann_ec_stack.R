
ec_merged_global <- purrr::reduce(ec_list_stratified, full_join, by = "Gene_Family") %>%
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>%
  column_to_rownames("Gene_Family") # rows = Gene_Familys, cols = samples

col_sums_global <- colSums(ec_merged_global, na.rm = TRUE)
col_sums_global[col_sums_global == 0] <- 1
ec_norm_global <- sweep(ec_merged_global, 2, col_sums_global, "/")

target_pattern <- paste(plot_data_discussion, collapse = "|")
ec_subset <- ec_norm_global %>% 
  rownames_to_column('Gene_Family') %>% 
  filter(str_starts(Gene_Family,merged_stratified_df_pattern)) %>% 
  column_to_rownames('Gene_Family')

common_samples <- intersect(colnames(ec_subset), rownames(meta_adj_batch_adjust))
ec_subset <- ec_subset[, common_samples, drop = FALSE]
meta_ready <- meta_adj_batch_adjust[common_samples, , drop = FALSE]

fit_ec_adj <- MMUPHin::adjust_batch(
  feature_abd = ec_subset, 
  batch = "batch",
  covariates = c("Group", "Age", "BMI", "Sex"),
  data = meta_ready,
  control = list(verbose = TRUE)
)

ec_abd_adj_adjust <- fit_ec_adj$feature_abd_adj %>%
  as.data.frame()  %>% rownames_to_column('Gene_Family')

ec_plot_data_long <- ec_abd_adj_adjust %>%
  as.data.table() %>% 
  pivot_longer(-Gene_Family, names_to = "Sample_ID", values_to = "Abundance") %>%
  mutate(Species = str_extract(Gene_Family,'(s__[A-Za-z0-9- _|]+)|unclassified'),
         genus = str_extract(Gene_Family,'(g__[A-Za-z0-9-_|]+)|unclassified'),
         Pathway = str_extract(Gene_Family,'[A-Za-z0-9: -_]+'
         )
  ) %>% 
  inner_join(metadata_allcohorts, by = "Sample_ID")


ec_plot_data_split <- ec_plot_data_long %>%
  split(.$Pathway)  

out_dir_stack_ec <-  paste0(base_symbol,
                            "downstream_analysis/output/picture/22.ec_species_contribution_adjust/")
if(!dir.exists(out_dir_stack_ec)) dir.create(out_dir_stack_ec, recursive = T)

for (pwy_id in plot_data_discussion) {
  p <- plot_species_contribution(ec_plot_data_long, pwy_id)
  
  if (!is.null(p)) {
    safe_name <- str_replace_all(pwy_id, "[^A-Za-z0-9]", "_")
    ggsave(paste0(out_dir_stack_ec, "Stack_", safe_name, ".pdf"), p, width = 8, height = 5)
    print(paste("Saved:", pwy_id))
  }
}
