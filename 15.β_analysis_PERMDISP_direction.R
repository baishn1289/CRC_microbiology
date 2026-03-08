library(vegan)
library(boot)     
library(effsize)

check_dispersion_direction_v3 <- function(feat, meta, group_col = "Group", title_prefix = "") {

  common <- intersect(rownames(feat), rownames(meta))
  feat <- feat[common, ]
  meta <- meta[common, ]

  if(any(feat < 0)) {
    feat[feat < 0] <- 0
  }
  
  groups <- factor(meta[[group_col]], levels = c("Control", "CRC"))
  
  dist_mat <- vegan::vegdist(feat, method = "bray")
  bd <- vegan::betadisper(dist_mat, groups)

  set.seed(123)
  perm_res <- vegan::permutest(bd, permutations = 999)
  perm_p <- perm_res$tab["Groups", "Pr(>F)"]

  df_dist <- data.frame(
    Sample_ID = rownames(meta),
    Group = groups,
    Distance = bd$distances
  )
  
  vec_crc <- df_dist$Distance[df_dist$Group == "CRC"]
  vec_ctrl <- df_dist$Distance[df_dist$Group == "Control"]
  
  mean_crc <- mean(vec_crc)
  mean_ctrl <- mean(vec_ctrl)
  mean_diff <- mean_crc - mean_ctrl  

  d_res <- effsize::cohen.d(vec_crc, vec_ctrl)

  diff_fun <- function(data, indices) {
    d <- data[indices, ] 
    m_ctrl <- mean(d$Distance[d$Group == "Control"])
    m_crc <- mean(d$Distance[d$Group == "CRC"])
    return(m_crc - m_ctrl) 
  }
  
  set.seed(123)
  boot_res <- boot::boot(data = df_dist, statistic = diff_fun, R = 1000)
  boot_ci <- boot::boot.ci(boot_res, type = "perc")
  ci_lower <- boot_ci$percent[4]
  ci_upper <- boot_ci$percent[5]
  

  p <- ggplot(df_dist, aes(x = Group, y = Distance, fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.5) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1.5, color = "grey30") +
    scale_fill_manual(values = c("Control" = "#3B4992", "CRC" = "#EE0000")) +
    
    annotate("text", x = 1.5, y = max(df_dist$Distance) * 1.08, 
             label = paste0("PERMDISP P = ", sprintf("%.3f", perm_p), "\n",
                            "Mean Diff (CRC-Ctrl) = ", sprintf("%.3f", mean_diff), "\n",
                            "Cohen's d = ", sprintf("%.2f", d_res$estimate)),
             size = 5, lineheight = 1) +
    
    labs(title = paste(title_prefix, "- Distance to Centroid"),
         y = "Distance to Centroid (Bray-Curtis)",
         x = NULL) +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          axis.text.x  = element_text(size = 16),
          axis.text.y  = element_text(size = 16),
          axis.title.y = element_text(size = 18))
  
  print(p)
  
  return(list(
    mean_diff = mean_diff,
    cohen_d = d_res$estimate,
    p_val = perm_p,
    plot = p,
    ci_lower = ci_lower,  
    ci_upper = ci_upper
  ))
}


keep_cohorts_fun <- function(meta, min_n = 3) {
  cohort_counts <- meta %>%
    count(Cohort, Group) %>%
    group_by(Cohort) %>%
    summarise(keep = all(c("CRC", "Control") %in% Group) && all(n >= min_n)) %>%
    filter(keep) %>%
    pull(Cohort) %>%
    as.character()
  return(cohort_counts)
}


meta_adj_batch_adjust_EO <- meta_adj_batch_adjust %>% filter(Age_class == "EO")
valid_cohorts_EO <- keep_cohorts_fun(meta_adj_batch_adjust_EO, min_n = 3)
valid_cohorts_EO
meta_clean_final_EO <- meta_adj_batch_adjust_EO %>% filter(Cohort %in% valid_cohorts_EO)
feat_clean_final_EO <- tax_mat_adj_adjust %>% 
  column_to_rownames("clade_name") %>% 
  t() %>% 
  .[rownames(meta_clean_final_EO), ]

setdiff(valid_cohorts_EO,unique(world_meta_EO$Cohort))
setdiff(unique(world_meta_EO$Cohort),valid_cohorts_EO)

meta_adj_batch_adjust_LO <- meta_adj_batch_adjust %>% filter(Age_class == "LO")
valid_cohorts_LO <- keep_cohorts_fun(meta_adj_batch_adjust_LO, min_n = 3)
valid_cohorts_LO
meta_clean_final_LO <- meta_adj_batch_adjust_LO %>% filter(Cohort %in% valid_cohorts_LO)
feat_clean_final_LO <- tax_mat_adj_adjust %>% 
  column_to_rownames("clade_name") %>% 
  t() %>% 
  .[rownames(meta_clean_final_LO), ]

setdiff(valid_cohorts_LO,unique(world_meta_LO$Cohort))
setdiff(unique(world_meta_LO$Cohort),valid_cohorts_LO)

res_eo_v3 <- check_dispersion_direction_v3(feat_clean_final_EO, 
                                           meta_clean_final_EO, 
                                           title_prefix = "Early-Onset (EO)")

ggsave(plot = res_eo_v3$plot,height = 8,width = 6,
       filename = paste0(base_symbol,
                          "downstream_analysis/output/picture/14.EO_Distance_to_Centroid.pdf"))

# LO
res_lo_v3 <- check_dispersion_direction_v3(feat_clean_final_LO, 
                                           meta_clean_final_LO, 
                                           title_prefix = "Late-Onset (LO)")
ggsave(plot = res_lo_v3$plot,height = 8,width = 6,
       filename =paste0(base_symbol,
                        "downstream_analysis/output/picture/14.LO_Distance_to_Centroid.pdf"))



dist_eo <- res_eo_v3$plot$data %>% mutate(Age_class = "EO")
dist_lo <- res_lo_v3$plot$data %>% mutate(Age_class = "LO")


dist_combined <- bind_rows(dist_eo, dist_lo)

dist_combined$Group <- factor(dist_combined$Group, levels = c("Control", "CRC"))
dist_combined$Age_class <- factor(dist_combined$Age_class, levels = c("LO", "EO"))


fit_disp_interaction <- lm(Distance ~ Group * Age_class, data = dist_combined)

summary(fit_disp_interaction)
