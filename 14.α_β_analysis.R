library(patchwork)
library(vegan)
library(tidyverse)
library(MMUPHin)
library(ape)
library(broom)
library(metafor)

calc_alpha_from_list <- function(data_list) {
  alpha_res_list <- lapply(names(data_list), function(cohort_name) {
    df <- data_list[[cohort_name]]

    if ("clade_name" %in% colnames(df)) {
      df <- df %>% column_to_rownames("clade_name")
    } else {
      df <- df %>% column_to_rownames(colnames(df)[1])
    }

    df_t <- t(df)

    shannon <- diversity(df_t, index = "shannon")
    simpson <- diversity(df_t, index = "simpson")
    richness <- specnumber(df_t) 
    
    data.frame(
      Sample_ID = names(shannon),
      Shannon = shannon,
      Simpson = simpson,
      Richness = richness,
      Cohort = cohort_name
    )
  })
  bind_rows(alpha_res_list)
}


alpha_all <- calc_alpha_from_list(meta_list_metadata_intersect_samples)


alpha_data_EO <- alpha_all %>%
  inner_join(world_meta_EO, by = "Sample_ID") %>% 
  select(Sample_ID, Shannon, Simpson, Richness, Group, Cohort = Cohort.y, Age, Sex, BMI)

alpha_data_LO <- alpha_all %>%
  inner_join(world_meta_LO, by = "Sample_ID") %>% 
  select(Sample_ID, Shannon, Simpson, Richness, Group, Cohort = Cohort.y, Age, Sex, BMI)

run_alpha_meta_standard <- function(alpha_df, metric = "Shannon") {
  
  alpha_df$Age <- as.numeric(as.character(alpha_df$Age))
  
  cohorts <- unique(as.character(alpha_df$Cohort))

  cohort_stats_list <- list()
  
  for (coh in cohorts) {

    sub_df <- alpha_df %>% 
      filter(Cohort == coh) %>%
      select(Sample_ID, Group, Age, Sex, BMI, all_of(metric)) %>%
      drop_na() 
    
    counts <- table(sub_df$Group)
    if(sum(counts >= 3) < 2) next 

    formula_str <- paste(metric, "~ Group + Age + Sex + BMI")
    fit <- lm(as.formula(formula_str), data = sub_df)

    target_level <- levels(sub_df$Group)[2]  
    target_term  <- paste0("Group", target_level)

    res <- tryCatch({
      broom::tidy(fit) %>% 
        filter(term == target_term) %>% 
        select(estimate, std.error, p.value)
    }, error = function(e) NULL)
    
    if(!is.null(res) && nrow(res) > 0) {
      cohort_stats_list[[coh]] <- data.frame(
        Cohort = coh,
        feature = metric,
        coef = res$estimate,
        stderr = res$std.error,
        pval = res$p.value,
        Type = "Individual"
      )
    }
  }

  cohort_res_df <- bind_rows(cohort_stats_list)

  if(nrow(cohort_res_df) < 2) {
    warning("Not enough cohorts for meta-analysis (N < 2). Returning NULL.")
    return(NULL)
  }

  meta_fit <- rma(yi = cohort_res_df$coef, 
                  sei = cohort_res_df$stderr, 
                  method = "REML")

  meta_sum <- summary(meta_fit)
  
  summary_res <- data.frame(
    Cohort = "Overall",
    feature = metric,
    coef = as.numeric(meta_sum$beta),  
    stderr = as.numeric(meta_sum$se), 
    pval = as.numeric(meta_sum$pval), 
    Type = "Summary",
    I2 = meta_fit$I2            
  )
  
  return(list(
    cohort_res = cohort_res_df, 
    meta_res = summary_res,     
    model = meta_fit            
  ))
}

dir.create(paste0(base_symbol,"downstream_analysis/output/EOCRC_Alpha"), recursive = T)
fit_alpha_EO <- run_alpha_meta_standard(alpha_df = alpha_data_EO, 
                                        metric = "Shannon")
dir.create(paste0(base_symbol,"downstream_analysis/output/LOCRC_Alpha"), recursive = T)
fit_alpha_LO <- run_alpha_meta_standard(alpha_df = alpha_data_LO, 
                                        metric = "Shannon")


purrr::reduce(list(fit_alpha_EO$cohort_res %>% mutate(Age_class = 'EO'),
                                    fit_alpha_LO$cohort_res %>% mutate(Age_class = 'LO'),
                                    fit_alpha_EO$meta_res %>% mutate(Age_class = 'EO') ,
                                    fit_alpha_LO$meta_res %>% mutate(Age_class = 'LO')), 
                               bind_rows) %>% 
  write.csv(.,row.names = F,
            file = paste0(base_symbol,
                          "downstream_analysis/output/table/Table.alpha_analysis_results.csv"))



get_alpha_data_v2 <- function(alpha_index, fit_obj) {
  
  cohort_res <- fit_obj$cohort_res
  meta_res   <- fit_obj$meta_res
  use_rename <- exists("rename_cohorts")
  use_country_map <- exists("cohort_to_country")
  
  res_df <- bind_rows(cohort_res, meta_res) %>%
    mutate(
      lower = coef - 1.96 * stderr,
      upper = coef + 1.96 * stderr,
      weight = 1 / (stderr^2),
      size_scaled = if_else(Type == "Summary", 0, sqrt(weight / max(weight, na.rm=T)) * 3.5 + 1),
      
      label_val = sprintf("%.2f", coef),
      label_ci  = sprintf("(%.2f, %.2f)", lower, upper),
      label_pval = case_when(
        pval < 0.001 ~ "< 0.001",
        pval < 0.01 ~ sprintf("%.3f", pval),
        TRUE ~ sprintf("%.2f", pval)
      ),

      Cohort_clean = if(use_rename) rename_cohorts(Cohort) else Cohort
    ) 

  if (use_country_map) {
    mapped_countries <- cohort_to_country[res_df$Cohort_clean]
    res_df$Country <- ifelse(is.na(mapped_countries), "Unknown", mapped_countries)
    res_df$Country[res_df$Cohort == "Overall"] <- "Global"
  } else {
    res_df$Country <- ifelse(res_df$Cohort == "Overall", "Global", "Unknown")
  }
  
  res_df <- res_df %>%
    mutate(feature = if_else(str_detect(feature, "s__"), 
                             str_remove(str_extract(feature, 's__[A-Za-z0-9_]+'), '^s__'), 
                             feature))
  
  return(res_df)
}

forest_alpha_fun <- function(alpha_index, Age_class) {

  obj_name <- paste0("fit_alpha_", Age_class)
  if (!exists(obj_name)) {
    warning(paste("Object", obj_name, "not found."))
    return(NULL)
  }
  fit_obj <- get(obj_name)

  if(is.null(fit_obj)) {
    warning(paste("Fit object", obj_name, "is NULL (likely insufficient cohorts). Skipping plot."))
    return(NULL)
  }
  
  plot_df_base <- get_alpha_data_v2(alpha_index, fit_obj) %>% 
    mutate(feature = str_remove(str_extract(feature, 's__[A-Za-z0-9_]+'), '^s__'))
  
  if (nrow(plot_df_base) == 0) {
    warning(paste("No data found for species:", alpha_index, "in", Age_class))
    return(NULL)
  }

  df_text <- plot_df_base %>%
    mutate(is_bold = (pval < 0.05)) %>%
    select(Cohort_clean, Country, Type, coef, lower, upper, size_scaled, label_ci, label_pval, is_bold)

  header_row <- data.frame(
    Cohort_clean = "Cohort", Type = "Header",
    label_ci = "Effect Size (95% CI)", label_pval = "P-value",
    Country = NA, coef = NA, lower = NA, upper = NA, size_scaled = NA, is_bold = FALSE
  )
  

  ind_order <- df_text %>% 
    filter(Type == "Individual") %>% 
    arrange(Country, Cohort_clean) %>%
    pull(Cohort_clean) %>% 
    unique() 

  final_levels <- c("Overall", rev(ind_order), "Cohort")
  
  plot_df_final <- bind_rows(df_text, header_row) %>%
    mutate(Cohort_clean = factor(Cohort_clean, levels = final_levels))

  p_left <- plot_df_final %>%
    ggplot(aes(y = Cohort_clean)) +
    geom_text(aes(x = 0, label = Cohort_clean), hjust = 0, size = 3.8, color = "black",
              fontface = ifelse(plot_df_final$Type == "Header", "bold", "plain")) +
    theme_void() +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + # 修改 expand 防止截断
    theme(plot.margin = margin(t = 0, r = 5, b = 0, l = 0))

  x_min <- min(plot_df_base$lower, na.rm = T)
  x_max <- max(plot_df_base$upper, na.rm = T)
  x_padding <- (x_max - x_min) * 0.2 
  xlims <- c(x_min - x_padding, x_max + x_padding)
  x_label <- paste0("Log Fold Change (", Age_class, "CRC vs ", Age_class, "Control)")
  
  p_mid <- plot_df_final %>%
    filter(Type != "Header") %>% 
    ggplot(aes(y = Cohort_clean, x = coef)) +
    theme_classic() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbar(data = . %>% filter(Type == "Individual"), aes(xmin = lower, xmax = upper), 
                  width = 0.2, linewidth = 0.6, color = "black") + 
    geom_point(data = . %>% filter(Type == "Individual"), aes(size = size_scaled), fill = "black", shape = 22, color = "black") +
    geom_errorbar(data = . %>% filter(Type == "Summary"), aes(xmin = lower, xmax = upper), 
                  width = 0.2, linewidth = 0.6, color = "black") + 
    geom_point(data = . %>% filter(Type == "Summary"), shape = 23, size = 5, 
               fill = "#D62728", color = "black", stroke = 1) +
    scale_size_identity() +
    coord_cartesian(xlim = xlims) + 
    labs(x = x_label, y = NULL) +
    theme(
      axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
      legend.position = "none", panel.grid.major.y = element_blank(),
      axis.title.x = element_text(size = 11, margin = margin(t = 10)),
      plot.margin = margin(t = 0, l = 5, r = 5, b = 0) 
    ) +
    scale_y_discrete(drop = FALSE) # 关键：防止空行被丢弃

  p_right <- plot_df_final %>%
    ggplot(aes(y = Cohort_clean)) +
    geom_text(aes(x = 0, label = label_ci), hjust = 0.5, size = 3.8, color = "black",
              fontface = ifelse(plot_df_final$Type == "Header", "bold", "plain")) +
    geom_text(aes(x = 1.2, label = label_pval), hjust = 0.5, size = 3.8, color = "black",
              fontface = ifelse(plot_df_final$Type == "Header" | plot_df_final$is_bold, "bold", "plain")) +
    theme_void() +
    scale_x_continuous(limits = c(-0.5, 1.7)) + 
    theme(plot.margin = margin(l = 5)) 

  layout_widths <- c(3, 5, 4) 
  alpha_index_short <- paste0(alpha_index,' index')
  
  final_plot <- p_left + p_mid + p_right +
    plot_layout(widths = layout_widths) + 
    plot_annotation(
      title = paste0("Meta-analysis (", Age_class, "): ", alpha_index_short),
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 20)))
    )

  save_dir <- paste0(base_symbol,"downstream_analysis/output/picture/")
  if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  file_name <- paste0(save_dir, "13.forestPlot_", Age_class, "_", alpha_index_short, ".pdf")
  ggsave(filename = file_name, plot = final_plot, width = 11, height = 6)
  
  return(final_plot)
}

p_alpha_forest_eo <- forest_alpha_fun('Shannon', "EO")
p_alpha_forest_lo <- forest_alpha_fun('Shannon', "LO")


run_beta_meta_final <- function(feature_abd, meta_data, output_dir
                                ) {
  
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  common_samples <- intersect(colnames(feature_abd), meta_data$Sample_ID)
  
  meta_base <- meta_data %>% 
    filter(Sample_ID %in% common_samples) %>%
    column_to_rownames("Sample_ID") %>%
    mutate(Cohort = as.factor(Cohort),
           Age = as.numeric(Age)) 

  feat_base <- feature_abd[, rownames(meta_base), drop=FALSE]

  stopifnot(all(colnames(feat_base) == rownames(meta_base)))

  sample_sums <- colSums(feat_base)
  if (any(sample_sums == 0)) {
    zero_samples <- names(sample_sums[sample_sums == 0])
    feat_base <- feat_base[, !colnames(feat_base) %in% zero_samples]
    meta_base <- meta_base[!rownames(meta_base) %in% zero_samples, ]
  }
  
  cohorts <- unique(as.character(meta_base$Cohort))
  res_list <- list()
  
  message(">>> Starting PERMANOVA & PERMDISP Analysis...")

  for (coh in cohorts) {
    sub_meta_raw <- meta_base %>% filter(Cohort == coh)

    vars_covariates <- c() 
    if (sum(!is.na(sub_meta_raw$Age)) > 5 && var(sub_meta_raw$Age, na.rm=T) > 0) vars_covariates <- c(vars_covariates, "Age")
    if (sum(!is.na(sub_meta_raw$Sex)) > 5 && length(unique(na.omit(sub_meta_raw$Sex))) > 1) vars_covariates <- c(vars_covariates, "Sex")
    if (sum(!is.na(sub_meta_raw$BMI)) > 5 && var(sub_meta_raw$BMI, na.rm=T) > 0) vars_covariates <- c(vars_covariates, "BMI")
    
    vars_to_use <- c(vars_covariates, "Group")
    
    sub_meta <- sub_meta_raw %>% 
      mutate(Group = factor(Group, levels = c("Control", "CRC"))) %>%
      select(all_of(vars_to_use)) %>% 
      drop_na() %>% droplevels()
    
    if (nrow(sub_meta) < 5 || length(unique(sub_meta$Group)) < 2) next
    
    sub_feat <- feat_base[, rownames(sub_meta), drop=FALSE]
    sub_feat <- sub_feat[rowSums(sub_feat) > 0, , drop=FALSE]
    
    if(nrow(sub_feat) == 0) next
    if(ncol(sub_feat) == 0 || nrow(sub_meta) != ncol(sub_feat)) next
    
    dist_mat <- vegan::vegdist(t(sub_feat), method = "bray")
    if(any(is.na(dist_mat))) next
    
    formula_str <- paste("dist_mat ~", paste(vars_to_use, collapse = " + "))
    
    tryCatch({
      set.seed(123)
      ad_res <- vegan::adonis2(as.formula(formula_str), data = sub_meta, permutations = 999)

      if("Group" %in% rownames(ad_res)) {
        r2_val <- ad_res["Group", "R2"]
        p_val  <- ad_res["Group", "Pr(>F)"]
      } else {

        r2_val <- ad_res[1, "R2"] 
        p_val  <- ad_res[1, "Pr(>F)"]
      }

    
      res_list[[coh]] <- data.frame(
        Cohort = coh, feature = "Beta_Diversity",
        coef = r2_val, pval = p_val,
        Type = "Individual"
      )
    }, error = function(e) { message(paste("Error:", coh, e$message)) })
  }
  
  message(">>> Running Overall Stratified Analysis...")

  meta_overall <- meta_base %>%
    mutate(Group = factor(Group, levels = c("Control", "CRC")),
           Cohort = as.factor(Cohort)) %>% 
    drop_na(Group, Age, Sex, BMI) %>% droplevels()
  
  feat_overall <- feat_base[, rownames(meta_overall)]
  feat_overall_norm <- sweep(feat_overall, 2, colSums(feat_overall), "/")
  
  fit_adjust <- MMUPHin::adjust_batch(
    feature_abd = feat_overall_norm,
    batch = "Cohort",
    covariates = c("Group", "Age", "Sex", "BMI"),
    data = meta_overall,
    control = list(verbose = FALSE)
  )
  feat_adj <- fit_adjust$feature_abd_adj
  feat_adj[feat_adj < 0] <- 0
  
  dist_global <- vegan::vegdist(t(feat_adj), method = "bray")
  
  set.seed(123)
  ad_global <- vegan::adonis2(dist_global ~ Age + Sex + BMI + Group, 
                              data = meta_overall, strata = meta_overall$Cohort, permutations = 999)
  
  r2_global <- if("Group" %in% rownames(ad_global)) ad_global["Group", "R2"] else ad_global[1, "R2"]
  p_global  <- if("Group" %in% rownames(ad_global)) ad_global["Group", "Pr(>F)"] else ad_global[1, "Pr(>F)"]
  
  res_global <- data.frame(
    Cohort = "Overall", feature = "Beta_Diversity",
    coef = r2_global, pval = p_global,
    Type = "Summary"
  )
  
  final_df <- bind_rows(res_list) %>% bind_rows(res_global)
  write.csv(final_df, file.path(output_dir, "beta_meta_stats_with_permdisp.csv"), row.names = F)
  return(final_df)
}

fit_beta_EO_df <- run_beta_meta_final(feat_EO, world_meta_EO, 
                                      paste0(base_symbol,"downstream_analysis/output/EOCRC_Beta"))
fit_beta_EO_df
fit_beta_LO_df <- run_beta_meta_final(feat_LO, world_meta_LO, 
                                      paste0(base_symbol,"downstream_analysis/output/LOCRC_Beta"))
fit_beta_LO_df

bind_rows(fit_beta_EO_df %>% mutate(Age_class ='EO'),
          fit_beta_LO_df %>% mutate(Age_class ='LO')) %>% 
  write.csv(.,row.names =F,
            file = paste0(base_symbol,
                          "downstream_analysis/output/table/Table.beta_analysis_results.csv"))


plot_beta_lollipop <- function(Age_class) {

  obj_name <- paste0("fit_beta_", Age_class, "_df")
  if (!exists(obj_name)) stop(paste("Error: Data object", obj_name, "not found."))
  beta_stats <- get(obj_name)

  use_rename <- exists("rename_cohorts")
  
  plot_df <- beta_stats %>%
    mutate(
      stars = case_when(
        pval < 0.001 ~ "  ***",
        pval < 0.01  ~ "  **",
        pval < 0.05  ~ "  *",
        TRUE         ~ "" 
      ),
      Cohort_clean = case_when(
        Cohort == "Overall" ~ "Meta-Analysis",
        use_rename ~ rename_cohorts(Cohort),
        TRUE ~ as.character(Cohort)
      ),
      Point_Type = if_else(Cohort == "Overall", "Summary", "Individual")
    )

  ind_order <- plot_df %>% 
    filter(Cohort != "Overall") %>% 
    arrange(coef) %>% 
    pull(Cohort_clean) %>% 
    unique() 

  if(length(ind_order) == 0) ind_order <- character(0)

  final_levels <- unique(c("Meta-Analysis", ind_order))
  
  plot_df <- plot_df %>%
    mutate(Cohort_clean = factor(Cohort_clean, levels = final_levels))

  x_max <- max(plot_df$coef, na.rm = T) * 1.15 
  if(x_max == 0) x_max <- 0.1
  
  p <- ggplot(plot_df, aes(x = coef, y = Cohort_clean)) +
    theme_classic() +
    
    geom_segment(aes(x = 0, xend = coef, y = Cohort_clean, yend = Cohort_clean), 
                 color = "grey70", size = 0.8) +
    geom_point(aes(shape = Point_Type, color = Point_Type, size = Point_Type, fill = Point_Type), 
               stroke = 1) +
    geom_text(aes(label = stars), hjust = -0.3, vjust = 0.75, size = 5, color = "black") +
    
    scale_shape_manual(values = c("Individual" = 21, "Summary" = 23)) + 
    scale_color_manual(values = c("Individual" = "black", "Summary" = "black")) +
    scale_fill_manual(values = c("Individual" = "white", "Summary" = "#D62728")) + 
    scale_size_manual(values = c("Individual" = 3.5, "Summary" = 5)) +
    
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    
    labs(
      title = paste0("Beta Diversity Effect Size (", Age_class, "CRC vs ", Age_class, 'Control)'),
      subtitle = "PERMANOVA R2 based on Bray-Curtis dissimilarity",
      x = expression(paste("Effect Size (", italic("R")^2, ")")),
      y = NULL
    ) +
    
    theme(
      legend.position = "none", 
      axis.line.y = element_blank(), 
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(size = 11, color = "black", face = "bold"), 
      axis.text.x = element_text(size = 10, color = "black"),
      plot.margin = margin(20, 20, 20, 20)
    )

  meta_idx <- which(levels(plot_df$Cohort_clean) == "Meta-Analysis")
  if (length(meta_idx) == 1) {
    sep_pos <- meta_idx + 0.5
    p <- p + geom_hline(yintercept = sep_pos, linetype = "dashed", color = "grey80")
  }

  save_dir <- paste0(base_symbol,"downstream_analysis/output/picture/")
  if(!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  file_name <- paste0(save_dir, "13.Lollipop_Beta_", Age_class, ".pdf")
  ggsave(filename = file_name, plot = p, width = 8, height = 6) 
  
  return(p)
}

p_lollipop_eo <- plot_beta_lollipop("EO")
p_lollipop_lo <- plot_beta_lollipop("LO")

