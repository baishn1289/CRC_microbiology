
dim(metadata_allcohorts)
dim(merged_df)

meta_adj_batch <- metadata_allcohorts %>%
  mutate(
    batch = factor(Cohort),        
    Group = factor(Group),         
    Age_class = factor(Age_class), 
    Sex = factor(Sex),            
    BMI = as.numeric(BMI),
    Age = as.numeric(Age)       
  ) %>% 
  column_to_rownames("Sample_ID")

print(sum(is.na(meta_adj_batch$BMI)))
print(sum(is.na(meta_adj_batch$Sex)))
table(meta_adj_batch$Sex)
print(table(meta_adj_batch$Cohort, is.na(meta_adj_batch$BMI)))

prop.table(table(metadata_allcohorts$Cohort, is.na(metadata_allcohorts$BMI)), 1)
table_missing_crc <- table(is.na(metadata_allcohorts$BMI), metadata_allcohorts$Group)
chisq.test(table_missing_crc)


merged_df_batch <- merged_df %>%     
  column_to_rownames("clade_name") %>%
  as.matrix()

merged_df_batch <- merged_df_batch/100   
# 确保列名与 meta 行名对齐
merged_df_batch <- merged_df_batch[, rownames(meta_adj_batch)]
dim(merged_df_batch)
dim(meta_adj_batch)
if (!is.numeric(merged_df_batch)) {
  mode(merged_df_batch) <- "numeric"
}

table(colSums(merged_df_batch))

merged_df_batch <- as.data.frame(merged_df_batch)

setdiff(rownames(meta_adj_batch),names(merged_df_batch))
setdiff(names(merged_df_batch),rownames(meta_adj_batch))

dim(merged_df_batch)


meta_adj_batch_adjust <- meta_adj_batch %>% 
  filter(!is.na(BMI) & !is.na(Age))
merged_df_batch_adjust <- merged_df_batch %>% 
  dplyr::select(any_of(rownames(meta_adj_batch_adjust))
  )
dim(merged_df_batch_adjust)
dim(meta_adj_batch_adjust)            
setdiff(rownames(meta_adj_batch_adjust),names(merged_df_batch_adjust))
setdiff(names(merged_df_batch_adjust),rownames(meta_adj_batch_adjust))

fit_tax_adj_adjust <- MMUPHin::adjust_batch(
  feature_abd = merged_df_batch_adjust,
  batch = "batch",
  covariates = c("Group", "BMI", "Sex", "Age"), 
  
  data = meta_adj_batch_adjust,
  control = list(verbose = TRUE)
)

tax_mat_adj_adjust <- fit_tax_adj_adjust$feature_abd_adj %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'clade_name')
