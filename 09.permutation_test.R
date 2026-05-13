library(showtext)
library(furrr)      
library(progress)  
library(scales)
library(glue)

target_25_species <- df_compare %>% 
  filter(Status =='Shared (Consistent)') %>% 
  pull(Species)

abd_27 <- merged_df %>%
  filter(str_detect(clade_name,paste(target_25_species, collapse = "|"))) %>%
  column_to_rownames("clade_name") %>%
  as.matrix()

missing_samples <- setdiff(metadata_allcohorts$Sample_ID, colnames(abd_27))
if (length(missing_samples) > 0) {
  warning("The following samples can not be found in abd_27：", paste(head(missing_samples, 10), collapse = ", "),
          if(length(missing_samples) > 10) paste0(" ... ( ", length(missing_samples)-10, " )"))
}

meta_all <- metadata_allcohorts %>% filter(Sample_ID %in% colnames(abd_27))

if(!"Age_class" %in% colnames(meta_all)) stop("metadata doesn't have Age_class column")
meta_all <- meta_all %>%
  mutate(Age = as.numeric(Age),
         Sex = as.factor(Sex),
         BMI = as.numeric(BMI),
         Age_class = as.character(Age_class)) %>%
  column_to_rownames(var = "Sample_ID")

compute_meta_coefs <- function(meta_df_rownames, abd_matrix,Age_Group) {
  meta_df_rownames <- meta_df_rownames %>% 
    mutate(Group = factor(Group, levels = c("Control","CRC")),
           Cohort = as.factor(Cohort))

  sample_ids <- rownames(meta_df_rownames)
  missing <- setdiff(sample_ids, colnames(abd_matrix))
  if(length(missing) > 0) {
    stop("compute_meta_coefs(): The following sample can not be found in abundance matrix: ", paste(missing, collapse = ", "))
  }
  abd_sub <- abd_matrix[, sample_ids, drop = FALSE]
  
  meta_dir <- file.path(meta_dir_base, Age_Group)
  
  control_list = list(rma_method="HS", transform="LOG", 
                      normalization = 'NONE', min_prevalence = 0, 
                      min_abundance = 0.001,output =  meta_dir)
  
  fit <- tryCatch({
    suppressMessages({
      MMUPHin::lm_meta(feature_abd = abd_sub, exposure = 'Group', batch = 'Cohort',
                       covariates = c("Age","Sex","BMI"), data = meta_df_rownames,
                       control = control_list)
    })
  }, error = function(e) {
    message("lm_meta failed: ", e$message)
    return(NULL)
  })
  
  if(is.null(fit)) {
    return(tibble(feature = rownames(abd_sub), coef = NA_real_))
  }
  
  if(!"meta_fits" %in% names(fit)) {
    stop("lm_meta output does not include meta_fits; please inspect the object structure.")
  }
  
  res <- fit$meta_fits %>%
    as_tibble() %>%
    select(feature, coef) %>%
    mutate(coef = as.numeric(coef))
  
  missing_feats <- setdiff(rownames(abd_sub), res$feature)
  if(length(missing_feats) > 0) {
    res <- bind_rows(res, tibble(feature = missing_feats, coef = NA_real_))
  }

  res <- res %>% slice(match(rownames(abd_sub), feature))
  return(res %>% select(feature, coef))
}


meta_EO <- meta_all %>%
  filter(Age_class == "EO",Cohort %in% world_meta_EO$Cohort)
meta_LO <- meta_all %>%
  filter(Age_class == "LO",Cohort %in% world_meta_LO$Cohort)

meta_dir_base <- paste0(base_symbol,"downstream_analysis/meta_permutation")
dir.create(meta_dir_base, recursive = TRUE, showWarnings = FALSE)

dir.create(file.path(meta_dir_base,"EO"), showWarnings = FALSE)
dir.create(file.path(meta_dir_base,"LO"), showWarnings = FALSE)

res_EO_obs <- compute_meta_coefs(meta_df_rownames = meta_EO, 
                                 abd_matrix = abd_27, 
                                 Age_Group ='EO')
res_LO_obs <- compute_meta_coefs(meta_df_rownames = meta_LO, 
                                 abd_matrix = abd_27, 
                                 Age_Group ='LO')


res_obs <- inner_join(res_EO_obs, res_LO_obs, by = "feature", suffix = c("_EO","_LO")) %>%
  mutate(delta_abs = abs(coef_EO) - abs(coef_LO))

observed_delta <- median(res_obs$delta_abs, na.rm = TRUE)
message("Observed median delta_abs for target features: ", observed_delta)

run_single_permutation <- function(seed_val, meta_df_with_rownames, abd_matrix) {
  # meta_df_with_rownames: metadata with rownames = Sample_ID
  set.seed(seed_val)
  meta_df <- meta_df_with_rownames %>% 
    tibble::rownames_to_column(var = "Sample_ID")

  meta_perm <- meta_df %>%
    group_by(Cohort, Group) %>%
    mutate(Age_class_perm = sample(Age_class)) %>%
    ungroup()
  
  meta_perm2 <- meta_perm %>%
    select(-Sample_ID) %>%
    as.data.frame()
  
  rownames(meta_perm2) <- meta_perm$Sample_ID

  meta_EO_perm <- meta_perm2[meta_perm2$Age_class_perm == "EO", , drop = FALSE] %>% 
    filter(Cohort %in% world_meta_EO$Cohort)
  meta_LO_perm <- meta_perm2[meta_perm2$Age_class_perm == "LO", , drop = FALSE] %>% 
    filter(Cohort %in% world_meta_LO$Cohort)
  
  if(nrow(meta_EO_perm) < 2 || nrow(meta_LO_perm) < 2) {
    return(NA_real_)
  }
  
  res_EO <- tryCatch({
    compute_meta_coefs(meta_df_rownames = meta_EO_perm, 
                       abd_matrix = abd_matrix,
                       Age_Group ='EO')
  }, error = function(e) {
    message("compute_meta_coefs EO error: ", e$message)
    return(tibble(feature = rownames(abd_matrix), coef = NA_real_))
  })
  
  res_LO <- tryCatch({
    compute_meta_coefs(meta_df_rownames = meta_LO_perm, 
                       abd_matrix = abd_matrix,
                       Age_Group ='LO')
  }, error = function(e) {
    message("compute_meta_coefs LO error: ", e$message)
    return(tibble(feature = rownames(abd_matrix), coef = NA_real_))
  })
  
  joined <- inner_join(res_EO, res_LO, by = "feature", suffix = c("_EO","_LO")) %>%
    mutate(delta_abs = abs(coef_EO) - abs(coef_LO))
  
  stat <- median(joined$delta_abs, na.rm = TRUE)
  return(stat)
}

plan(multisession, workers = 12)  

n_permutations <- 1000
perm_null_distribution <- future_map_dbl(1:n_permutations, ~run_single_permutation(.x + 1000, 
                                                                                   meta_df_with_rownames= meta_all,
                                                                                   abd_matrix= abd_27),
                                         .options = furrr_options(seed = TRUE),
                                         .progress = TRUE)
future::plan(sequential)

perm_null_distribution_clean <- perm_null_distribution[!is.na(perm_null_distribution)]
n_failed <- sum(is.na(perm_null_distribution))


empirical_p <- (sum(perm_null_distribution_clean >= observed_delta) + 1) / (length(perm_null_distribution_clean) + 1)
message("Empirical permutation p-value = ", empirical_p)

df_null <- tibble(Delta = perm_null_distribution_clean)

p_pub <- df_null %>%
  ggplot(aes(x = Delta)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "#d9d9d9", 
                 color = "#666666", size = 0.2) +
  geom_density(aes(y = after_stat(density)), color = "#2b8cbe", linewidth = 0.9, alpha = 0.3) +
  geom_vline(xintercept = observed_delta, color = "#de2d26", 
             linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = median(perm_null_distribution_clean, na.rm = TRUE), 
             color = "#252525", linetype = "solid", linewidth = 0.6) +
  annotate("text", x = observed_delta, y = max(density(perm_null_distribution_clean)$y)*0.95,
           label = glue("Observed = {round(observed_delta,3)}\nempirical p = {formatC(empirical_p, digits=3)}"),
           hjust = -0.05, vjust = 1, color = "#de2d26", fontface = "bold", size = 3.6) +
  annotate("text", x = Inf, y = Inf, label = glue("n_perm = {length(perm_null_distribution_clean)}"),
           hjust = 1.05, vjust = 2, size = 3.5) +
  labs(title = "Cohort-stratified permutation null distribution",
       subtitle = glue("Global statistic: median(|coef_EO| - |coef_LO|)"),
       x = expression("Median(" * abs(beta[EO]) - abs(beta[LO]) * ")"),
       y = "Density") +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    axis.title = element_text(size = 12),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.4)
  )
p_pub
showtext_auto()
ggsave(paste0(base_symbol,
              "downstream_analysis/output/picture/5.4.permutation_null_publication.pdf"), 
       p_pub, width = 11, height = 5)


ecdf_plot <- tibble(stat = perm_null_distribution_clean) %>%
  ggplot(aes(x = stat)) +
  stat_ecdf(geom = "step", linewidth = 0.8) +
  geom_vline(xintercept = observed_delta, color = "#de2d26", linetype = "dashed") +
  annotate("text", x = observed_delta, y = 0.05, 
           label = paste0("obs=",round(observed_delta,3)), hjust = -0.05) +
  labs(x = "Permuted statistic", y = "ECDF", title = "Empirical CDF of permuted stats")

ggsave(paste0(base_symbol,
              "downstream_analysis/output/picture/5.4.ecdf_plot.pdf"), 
       ecdf_plot, width = 11, height = 5)
