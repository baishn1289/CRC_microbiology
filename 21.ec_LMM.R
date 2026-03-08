library(lme4)
library(broom.mixed)
library(emmeans)
library(ggsignif) 
library(ggbreak)
library(broom.mixed)
library(emmeans)
library(patchwork)
library(scales)
library(lmerTest)

ec_LMM_input <- ec_subset %>%
  rownames_to_column('Gene_Family') %>% 
  as.data.table() %>% 
  pivot_longer(-Gene_Family, names_to = "Sample_ID", values_to = "Abundance") %>%
  mutate(Species = str_extract(Gene_Family,'(s__[A-Za-z0-9- _|]+)|unclassified'),
         genus = str_extract(Gene_Family,'(g__[A-Za-z0-9-_|]+)|unclassified'),
         Pathway = str_extract(Gene_Family,'[A-Za-z0-9: -_]+'
         )
  ) %>% 
  inner_join(metadata_allcohorts, by = "Sample_ID")


run_lmm_for_pathway_species_v2_mod <- function(pathway,
                                               species,
                                               df,
                                               out_dir = ".",
                                               pseudocount_default = 1e-10,
                                               models_store_name = "models_list",
                                               results_store_name = "models_results_df",
                                               save_plot = TRUE,
                                               plot_width = 4,
                                               plot_height = 9.5) {
  
  safe_name <- function(x) {
    x %>% as.character() %>%
      str_replace_all("[^A-Za-z0-9]+", "_") %>%
      str_replace_all("^_+|_+$", "")
  }
  elem_name <- paste0(safe_name(pathway), "__", safe_name(species))
  file_base <- paste0(safe_name(pathway), "__", safe_name(species))
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  df_sspa <- df %>% filter(Pathway == pathway, Species == species)
  if (nrow(df_sspa) == 0) stop("No rows found for given pathway and species.")
  all_nonzero_values <- df_sspa$Abundance[!is.na(df_sspa$Abundance) & df_sspa$Abundance > 0]
  min_nonzero_val <- if (length(all_nonzero_values) > 0) min(all_nonzero_values, na.rm = TRUE) else pseudocount_default

  df_lmm_ready <- df_sspa %>%
    ungroup() %>%
    mutate(
      BMI = as.numeric(BMI),
      Sex = as.factor(Sex),
      Cohort = as.factor(Cohort)
    ) %>%
    mutate(Compare_Group = paste(Age_class, Group, sep = "")) %>%
    mutate(Compare_Group = factor(Compare_Group, levels = c("LOCRC", "EOCRC", "LOControl", "EOControl"))) %>%
    mutate(
      Abundance_pc = Abundance + (min_nonzero_val / 2),
      Log_SSPA = log10(Abundance_pc)
    )

  model_all <- tryCatch({
    lmer(Log_SSPA ~ Compare_Group + Sex + BMI + (1 | Cohort), data = df_lmm_ready)
  }, error = function(e) stop("Model fitting failed: ", e$message))
  
  res_all <- broom.mixed::tidy(model_all, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate, std.error, p.value)
  
  emm <- tryCatch(emmeans(model_all, ~ Compare_Group), error = function(e) NULL)
  pairs_res <- if (!is.null(emm)) tryCatch(as.data.frame(pairs(emm, adjust = "none")), error = function(e) NULL) else NULL
  
  extract_pair_p <- function(df_pairs, a, b) {
    if (is.null(df_pairs)) return(NA_real_)
    found <- df_pairs %>% filter(str_detect(contrast, fixed(a)) & str_detect(contrast, fixed(b)))
    if (nrow(found) == 0) found <- df_pairs %>% filter(str_detect(contrast, fixed(b)) & str_detect(contrast, fixed(a)))
    if (nrow(found) == 0) return(NA_real_)
    found$p.value[1]
  }
  p_main <- extract_pair_p(pairs_res, "EOCRC", "LOCRC")
  p_eo   <- extract_pair_p(pairs_res, "EOCRC", "EOControl")
  p_lo   <- extract_pair_p(pairs_res, "LOCRC", "LOControl")
  
  format_p_text <- function(p) {
    if (is.na(p)) return("NA")
    if (p < 0.001) return("P < 0.001")
    sprintf("P = %.3g", p)
  }

  df_plot <- df_lmm_ready %>%
    mutate(Compare_Group = factor(Compare_Group, levels = c("EOControl", "EOCRC", "LOControl", "LOCRC")),
           Group = ifelse(str_detect(as.character(Compare_Group), "CRC"), "CRC", "Control"))

  df_prev <- df_plot %>%
    group_by(Compare_Group) %>%
    summarise(prevalence = mean(Abundance > 0, na.rm = TRUE), n = n(), .groups = "drop") %>% 
    mutate(Group = ifelse(str_detect(as.character(Compare_Group), "CRC"),
                          "CRC", "Control"))

  df_mean_ab <- df_plot %>%
    group_by(Compare_Group, Group) %>%
    summarise(mean_ab = mean(Abundance_pc, na.rm = TRUE), .groups = "drop") %>%
    mutate(Compare_Group = factor(Compare_Group, levels = c("EOControl", "EOCRC", "LOControl", "LOCRC")))

  y_max <- max(df_plot$Log_SSPA, na.rm = TRUE)
  y_min <- min(df_plot$Log_SSPA, na.rm = TRUE)
  y_range <- y_max - y_min
  tick_height <- 0.03 * y_range 
  
  annot_df <- tibble::tibble(
    comparison = c("EOControl vs EOCRC",
                   "LOControl vs LOCRC",
                   "EOCRC vs LOCRC"),
    x_start = c(1, 3, 2),
    x_end   = c(2, 4, 4),
    y       = y_max + c(0.2, 0.4, 0.6) * y_range,
    label   = c(format_p_text(p_eo),
                format_p_text(p_lo),
                format_p_text(p_main))
  )

  p1 <- ggplot(df_prev, aes(x = Compare_Group, y = prevalence, fill = Group)) +
    geom_col(width = 0.6, color = "black", alpha = 0.8) +
    geom_text(
      aes(label = scales::percent(prevalence, accuracy = 0.1)),
      vjust = -0.4, size = 3
    ) +
    scale_fill_manual(values = c("Control" = "#3B4992",
                                 "CRC" = "#EE0000")) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, min(1, max(df_prev$prevalence, na.rm = TRUE) * 1.3))
    ) +
    labs(title = "Prevalence (Abundance > 0)",
         x = NULL, y = "Prevalence") +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 10, face = "bold"),
      legend.position = "none"  
    )
  

  p2 <- ggplot(df_plot, aes(x = Compare_Group, y = Log_SSPA, fill = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
    geom_jitter(width = 0.18, size = 1.5, alpha = 0.5, color = "grey30") +
    scale_fill_manual(values = c("Control" = "#3B4992", "CRC" = "#EE0000")) +
    labs(title = "Log10(Abundance + pseudocount)",
         subtitle = NULL,
         x = NULL, y = expression(Log[10] ~ "(Abundance + pseudocount)")) +
    theme_classic() +
    theme(
      axis.text.x = element_text(size = 10, face = "bold"),
      legend.position = "none"  
    ) +
    geom_segment(  
      data = annot_df,
      aes(x = x_start, xend = x_end, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.6
    )+   geom_segment(
      data = annot_df,
      aes(x = x_start, xend = x_start,
          y = y, yend = y - tick_height),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +

    geom_segment(
      data = annot_df,
      aes(x = x_end, xend = x_end,
          y = y, yend = y - tick_height),
      inherit.aes = FALSE,
      linewidth = 0.6
    ) +
    geom_text( 
      data = annot_df,
      aes(x = (x_start + x_end) / 2, y = y + 0.08 * y_range, label = label),
      inherit.aes = FALSE,
      size = 3.2
    )
  

  y_max_mean <- max(df_mean_ab$mean_ab, na.rm = TRUE)
  p3 <- ggplot(df_mean_ab, aes(x = Compare_Group, y = mean_ab, fill = Group)) +
    geom_col(width = 0.6, color = "black", alpha = 0.9) +
    geom_text(aes(label = formatC(mean_ab, format = "e", digits = 2)), vjust = -0.4, size = 3) +
    scale_y_continuous(
      limits = c(0, y_max_mean * 1.2), 
      expand = expansion(mult = c(0, 0))
    ) +
    scale_fill_manual(values = c("Control" = "#3B4992", "CRC" = "#EE0000")) +
    labs(title = "Mean Abundance (per-group)",
         x = NULL, y = "Mean Abundance (with pseudocount)") +
    theme_classic() +
    theme(axis.text.x = element_text(size = 10, face = "bold"),
          legend.position = "none")
  p_combined <- p1  / p3 / p2 + plot_layout(heights = c(1,1.3, 1.8)) +
    plot_annotation(title = paste0(str_remove(species,'s__')),
                    subtitle = pathway %>%
                      gsub("^([0-9])", "EC \\1", .) %>%
                      str_replace(":.*$", ""),
                    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                                  plot.subtitle = element_text(hjust = 0.5,size = 10),
                                  plot.caption = element_text(hjust = 0.5) )
    )

  plot_filepath <- NA_character_
  if (save_plot) {
    plot_filepath <- file.path(out_dir, paste0(file_base, ".pdf"))
    tryCatch({
      ggsave(plot_filepath, p_combined, width = plot_width, height = plot_height)
    }, error = function(e) message("Failed to save plot: ", e$message))
  }

  res_row <- tibble::tibble(
    Pathway = pathway,
    Species = species,
    p_EO_vs_LO_CRC = p_main,
    p_EOCRC_vs_EOControl = p_eo,
    p_LOCRC_vs_LOControl = p_lo,
    n_samples = nrow(df_lmm_ready),
    model_term_count = nrow(res_all)
  )
  if (!exists(models_store_name, envir = .GlobalEnv)) assign(models_store_name, list(), envir = .GlobalEnv)
  models_list <- get(models_store_name, envir = .GlobalEnv)
  models_list[[elem_name]] <- model_all
  assign(models_store_name, models_list, envir = .GlobalEnv)
  
  if (!exists(results_store_name, envir = .GlobalEnv)) {
    assign(results_store_name, res_row, envir = .GlobalEnv)
  } else {
    existing <- get(results_store_name, envir = .GlobalEnv)
    existing <- bind_rows(existing, res_row)
    assign(results_store_name, existing, envir = .GlobalEnv)
  }
  
  invisible(list(
    model = model_all,
    res_all = res_all,
    emmeans_pairs = pairs_res,
    contrast_p = tibble::tibble(
      Pathway = pathway,
      Species = species,
      p_EO_vs_LO_CRC = p_main,
      p_EOCRC_vs_EOControl = p_eo,
      p_LOCRC_vs_LOControl = p_lo
    ),
    plots = list(p1 = p1, p2 = p2, p3 = p3, combined = p_combined),
    plot_file = plot_filepath
  ))
}

out_dir_LMM <- paste0(base_symbol,"downstream_analysis/output/picture/23_LMM/")
dir.create(out_dir_LMM)

pairs_to_run <- list(
  list(p = "2.5.1.30", s = "s__Coprococcus_eutactus"),
  list(p = "2.5.1.30", s = "s__Enterococcus_faecalis"),
  list(p = "2.5.1.30", s = "s__Enterococcus_faecium"),
  list(p = "2.5.1.30", s = "s__Flavonifractor_plautii"),
  list(p = "2.5.1.30", s = "s__Mitsuokella_multacida"),
  list(p = "2.5.1.30", s = "s__Porphyromonas_asaccharolytica"),
  list(p = "2.5.1.30", s = "s__Prevotella_intermedia"),
  list(p = "1.2.7.3", s = "s__Alistipes_shahii"),
  list(p = "1.2.7.3", s = "s__Odoribacter_splanchnicus"),
  list(p = "1.2.7.3", s = "s__Prevotella_buccae"),
  list(p = "1.2.7.3", s = "s__Acidaminococcus_intestini"),
  list(p = "1.2.7.3", s = "s__Porphyromonas_asaccharolytica"),
  list(p = "1.2.7.3", s = "s__Porphyromonas_uenonis"),
  list(p = "1.2.7.3", s = "s__Prevotella_intermedia"),
  list(p = "6.3.2.49", s = "s__Bacteroides_finegoldii"),
  list(p = "6.3.2.49", s = "s__Bacteroides_uniformis"),
  list(p = "6.3.2.49", s = "s__Bacteroides_xylanisolvens"),
  list(p = "6.3.2.49", s = "s__Eubacterium_hallii"),
  list(p = "6.3.2.49", s = "s__Hungatella_hathewayi"),
  list(p = "6.3.2.49", s = "s__Roseburia_faecis"),
  list(p = "6.3.2.49", s = "s__Ruminococcus_torques")
  )

rm(models_results_df)

for (it in pairs_to_run) {
  run_lmm_for_pathway_species_v2_mod(it$p, it$s, df = ec_LMM_input,
                                     out_dir = out_dir_LMM)
}

ec_LMM_resultes_df <- unique(models_results_df )


rm(models_results_df)
pairs_to_run_metacyc <- list(
  list (p = 'PWY-7688: dTDP-&alpha;-D-ravidosamine and dTDP-4-acetyl-&alpha;-D-ravidosamine biosynthesis',
        s = 's__Bacteroides_cellulosilyticus'),
  list (p= 'PWY-7688: dTDP-&alpha;-D-ravidosamine and dTDP-4-acetyl-&alpha;-D-ravidosamine biosynthesis',
        s = 's__Bacteroides_fragilis'),
  list (p= 'PWY-7688: dTDP-&alpha;-D-ravidosamine and dTDP-4-acetyl-&alpha;-D-ravidosamine biosynthesis',
        s ='s__Bacteroides_intestinalis')
)

for (it in pairs_to_run_metacyc) {
  run_lmm_for_pathway_species_v2_mod(it$p, it$s, df = metacyc_LMM_input,
                                     out_dir = out_dir_LMM)
}

metacyc_LMM_resultes_df <- unique(models_results_df )

LMM_results_ouptput <- ec_LMM_resultes_df %>% 
  select(Pathway,Species,p_EO_vs_LO_CRC,
         p_EOCRC_vs_EOControl,p_LOCRC_vs_LOControl)

write.csv(LMM_results_ouptput,row.names = F,
          file = paste0(base_symbol,
          "downstream_analysis/output/table/TableS11.ec_species_LMM.csv"))
