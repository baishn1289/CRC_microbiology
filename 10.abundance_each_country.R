# country_abundance
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(ggsci)
library(gghalves)
library(ggh4x)

target_features_EO <- meta_results_EO$feature
target_features_LO <- meta_results_LO$feature

plot_species_abundance <- function(data,
                                   target_species_full_name,
                                   Age_class,
                                   output_dir,adjust_mode,
                                   save_table = TRUE) {
  
  sub_df <- data %>%
    filter(clade_name == target_species_full_name)
  
  if (nrow(sub_df) == 0) return(NULL)
  
  short_name <- unique(sub_df$Species_Short)

  valid_countries <- sub_df %>%
    group_by(Country) %>%
    summarise(
      n_crc  = sum(Group == "CRC", na.rm = TRUE),
      n_ctrl = sum(Group == "Control", na.rm = TRUE),
      n_ctrl_pos = sum(Group == "Control" & LogAbundance>-10, na.rm = TRUE),
      n_crc_pos = sum(Group == "CRC" & LogAbundance>-10, na.rm = TRUE),
      # all_zero = all(LogAbundance == -10, na.rm = TRUE), 
      .groups = "drop"
    ) %>%
    filter( n_crc >= 3,n_ctrl >= 3) %>%
    pull(Country)
  
  if (length(valid_countries) == 0) {
    message(paste("Skipping", short_name,
                  "- No country has >= 3 samples in both groups."))
    return(NULL)
  }
  
  sub_df_filtered <- sub_df %>%
    filter(Country %in% valid_countries)
  
  if (save_table) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    file_name <- paste0(
      gsub(" ", "_", short_name),
      "_", adjust_mode,"_",Age_class,
      "_by_country_boxplot_data.csv"
    )
    
    output_path <- file.path(output_dir, file_name)

    sub_df_filtered %>%
      select(
        Sample_ID,
        Species = clade_name,
        Species_Short,
        Age_class,
        Country,
        Group,
        LogAbundance
      ) %>%
      arrange(Country, Group) %>%
      write.csv(output_path, row.names = FALSE)
  }
  
  Cohorts_name <- paste0("(", Age_class, " Cohorts)")
  
  p <- ggplot(sub_df_filtered,
              aes(x = Group, y = LogAbundance, fill = Group)) +
    
    # geom_violin(trim = F, alpha = 0.3, color = NA) +
    geom_boxplot(width = 0.4, alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 1, alpha = 0.4,
                color = "grey30") +
    facet_wrap(~Country, scales = "free_y", nrow = 1) +
    
    stat_compare_means(
      comparisons = list(c("Control", "CRC")),
      method = "wilcox.test",
      label = "p.format",
      tip.length = 0.01,
      size = 3.5
    ) +
    
    scale_fill_npg() +
    theme_bw() +
    labs(
      title = bquote(italic(.(short_name)) ~ .(Cohorts_name)),
      y = "Log10(Relative Abundance)",
      x = NULL
    ) +
    
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text = element_text(size = 10, color = "black"),
      strip.background = element_rect(fill = "#E5E5E5"),
      strip.text = element_text(face = "bold", size = 11),
      legend.position = "none"
    )
  
  return(p)
}

plot_species_abundance_ALL <- function(data, adjust_mode,top_species_short) {
  
  data_plot <- data %>%
    filter(!(Age_class == "EO" & Country == "Turkey"))

  country_levels_LO <- data_plot %>%
    filter(Age_class == "LO") %>%
    distinct(Country) %>%
    arrange(Country) %>%
    pull(Country)
  
  if ("Turkey" %in% country_levels_LO) {
    country_levels_LO <- c(
      setdiff(country_levels_LO, "Turkey"),
      "Turkey"
    )
  }

  data_plot <- data_plot %>%
    mutate(
      Country = factor(Country, levels = country_levels_LO),
      Age_class = factor(Age_class, levels = c( "LO","EO"))
    )

  p <- ggplot(
    data_plot,
    aes(x = Group, y = LogAbundance, fill = Group)
  ) +
    # geom_violin(trim = FALSE, alpha = 0.3, color = NA) +
    geom_jitter(width = 0.2, height = 0,size = 1, alpha = 0.4,
                color = "grey30") +
    geom_boxplot(
      width = 0.4,
      alpha = 0.8,
      outlier.shape = NA
    ) +
    
    facet_grid(
      rows = vars(Age_class),
      cols = vars(Country),
      scales = "free_y",
      drop = FALSE, switch = "y"  
    ) +
    
    stat_compare_means(
      comparisons = list(c("Control", "CRC")),
      method = "wilcox.test",
      label = "p.format",
      size = 3
    )+  scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.1))
    )+
    scale_fill_manual(
      values = c(
        "Control" = "#386cb0",
        "CRC"     = "#e31a1c"
      )
    )+
    theme_bw() +
    
    labs(
      title = paste0("Abundance by country (", adjust_mode, "ed by MMUPHin)"),
      y = paste0("Log10(Relative Abundance)\n(",top_species_short,")"),
      x = NULL
    ) +
    
    theme(
      legend.position = "none",
      panel.spacing = unit(4, "mm"),
      
      strip.background = element_rect(fill = "#E5E5E5"),
      strip.text = element_text(size = 11, face = "bold"),
      
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      )
    )+theme(
      strip.placement = "outside",   
      
      strip.text.y.left = element_text(     
        angle = 0,
        size = 12,
        face = "bold"),
      strip.background.y =element_blank(),
      panel.spacing = unit(4, "mm")
    )
  
  return(p)
}

country_abundance_boxplot <- function(mat_data_input,adjust_mode,plot_species_abundance_ALL_mode,
                                      top_species_EO_input,top_species_LO_input){
  
  plot_data_ready <- mat_data_input %>%
    # filter(clade_name %in% target_features) %>%
    pivot_longer(
      cols = -clade_name, 
      names_to = "Sample_ID", 
      values_to = "Abundance"
    ) %>%
    
    inner_join(metadata_allcohorts, by = "Sample_ID") %>%
    mutate(
      Species_Short = str_remove(str_extract(clade_name, 's__[A-Za-z0-9_]+'), '^s__'),
      LogAbundance = log10(Abundance + 1e-7),
      Group = factor(Group, levels = c("Control", "CRC"))
    ) 
  
  all_countries <- plot_data_ready %>%
    distinct(Country) %>%
    arrange(Country) %>%
    pull(Country)
  
  plot_data_ready_ALL <- plot_data_ready %>%
    filter(
      (Age_class == "EO" & clade_name %in% top_species_EO_input) |
        (Age_class == "LO" & clade_name %in% top_species_LO_input)
    ) %>%
    mutate(
      Country = factor(Country, levels = all_countries),
      Age_class = factor(Age_class, levels = c("EO", "LO"))
    )
  
  plot_data_ready_EO <- plot_data_ready %>% 
    filter(Age_class == "EO",
           clade_name %in% top_species_EO_input ) 
  
  plot_data_ready_LO <- plot_data_ready %>% 
    filter(Age_class == "LO",
           clade_name %in% top_species_LO_input ) 
  
  for (bug in top_species_EO_input) {
    p <- plot_species_abundance(data = plot_data_ready_EO, 
                                target_species_full_name = bug,
                                Age_class = 'EO',adjust_mode = adjust_mode,
                                output_dir = output_dir_EO)
    short_name <- str_remove(str_extract(bug, 's__[A-Za-z0-9_]+'), '^s__')
    file_name <- paste0(output_dir_EO, "EO_Boxplot_", adjust_mode, '_',short_name, ".pdf")
    ggsave(file_name, p, width = 12, height = 6)
    
    print(paste("Saved:", short_name))
  }
  
  
  for (bug in top_species_LO_input) {
    p <- plot_species_abundance(data = plot_data_ready_LO, 
                                target_species_full_name = bug,
                                Age_class = "LO",adjust_mode = adjust_mode,
                                output_dir = output_dir_LO)
    short_name <- str_remove(str_extract(bug, 's__[A-Za-z0-9_]+'), '^s__')
    file_name <- paste0(output_dir_LO, "LO_Boxplot_",  adjust_mode, '_',short_name, ".pdf")
    ggsave(file_name, p, width = 16, height = 2.5)
    
    print(paste("Saved:", short_name))
  }
  
  if (plot_species_abundance_ALL_mode == FALSE) {
    next
  }
  
  top_species_short <- str_remove(str_extract(top_species_EO_input, 's__[A-Za-z0-9_]+'), '^s__')
  
  p_all <- plot_species_abundance_ALL(
    data = plot_data_ready_ALL,
    adjust_mode = adjust_mode,
    top_species_short = top_species_short
  )
  
  file_all <- paste0(base_symbol,"downstream_analysis/output/picture/",
    "EO_LO_ALL_Boxplots_", top_species_short,'_',adjust_mode, ".pdf"
  )
  
  ggsave(
    file_all,
    p_all,width = 16,height = 4
  )
  
}

output_dir_EO <-  paste0(base_symbol,"downstream_analysis/output/picture/Boxplots_EO/")
if(!dir.exists(output_dir_EO)) dir.create(output_dir_EO, recursive = TRUE)

top_species_EO <- c(target_fp,target_fn,
                    target_features_EO[
                      str_detect(
                        target_features_EO,
                        "Akkermansia_muciniphila|Parabacteroides_distasonis|Faecalibacterium_prausnitzii|Roseburia_intestinalis")
                    ])

output_dir_LO <- paste0(base_symbol,"downstream_analysis/output/picture/Boxplots_LO/")
if(!dir.exists(output_dir_LO)) dir.create(output_dir_LO, recursive = TRUE)

top_species_LO <- c(target_fn,target_fp,
                    target_features_LO[
                      str_detect(
                        target_features_LO,
                        "Roseburia_hominis|Faecalibacterium_prausnitzii|Roseburia_intestinalis")
                    ])


country_abundance_boxplot(mat_data_input = tax_mat_adj_adjust,adjust_mode = 'adjust',
                          plot_species_abundance_ALL_mode = TRUE,
                          top_species_EO_input = target_fn,
                          top_species_LO_input = target_fn )

country_abundance_boxplot(mat_data_input = tax_mat_adj_adjust,adjust_mode = 'adjust',
                          plot_species_abundance_ALL_mode = TRUE,
                          top_species_EO_input = 'k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Roseburia|s__Roseburia_inulinivorans',
                          top_species_LO_input = 'k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Roseburia|s__Roseburia_inulinivorans')

country_abundance_boxplot(mat_data_input = tax_mat_adj_adjust,adjust_mode = 'adjust',
                          plot_species_abundance_ALL_mode = TRUE,
                          top_species_EO_input = 'k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Roseburia|s__Roseburia_faecis',
                          top_species_LO_input = 'k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Roseburia|s__Roseburia_faecis'
)


country_abundance_boxplot(mat_data_input = tax_mat_adj_adjust,adjust_mode = 'adjust',
                          plot_species_abundance_ALL_mode = TRUE,
                          top_species_EO_input = 'k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Enterocloster|s__Enterocloster_bolteae',
                          top_species_LO_input = 'k__Bacteria|p__Firmicutes|c__Clostridia|o__Eubacteriales|f__Lachnospiraceae|g__Enterocloster|s__Enterocloster_bolteae'
)

