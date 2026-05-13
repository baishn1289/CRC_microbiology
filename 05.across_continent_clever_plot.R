library(ggrepel)
library(ggprism)
library(patchwork)

df_all <- bind_rows(list( meta_results_Europe_EO,meta_results_Europe_LO,
                          meta_results_Asia_EO,meta_results_Asia_LO)) %>%
  mutate(Continent = ifelse(Continent =='East Asia','East_Asia',Continent))%>%
  mutate(Continent = factor(Continent, levels = c("Europe", "East_Asia")),
         feature = str_remove(str_extract(feature, 's__[A-Za-z0-9_]+'), '^s__'),
         Age_Group = Age_class,
         CI_lower = coef - (1.96 * stderr),
         CI_upper = coef + (1.96 * stderr),
         Report = sprintf("β = %.2f; 95%% CI, %.2f–%.2f; adjusted p = %.1e", 
                          coef, CI_lower, CI_upper, qval.fdr))%>% 
  arrange(Continent,Age_Group,qval.fdr) %>%
  select(feature, exposure, Continent, Age_Group,
         pval, Report, I2,consistency_score,max_weight,qc_status) 

writexl::write_xlsx(x = df_all,
       path = paste0(base_symbol,
                        "downstream_analysis/output/table/TableS6.continent_all_species.xlsx")
          )

plot_dot_design <- function(data_subset, title_text) {

  order_levels <- data_subset %>%
    group_by(feature) %>%
    summarise(mean_coef = mean(coef)) %>%
    arrange(mean_coef) %>%
    pull(feature)
  
  plot_data <- data_subset %>%
    mutate(feature = factor(feature, levels = order_levels))
  
  ggplot(plot_data, aes(x = coef, y = feature)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_segment(aes(x = 0, xend = coef, y = feature, yend = feature), color = "grey90", size = 0.5) +
    geom_point(aes(color = Continent, shape = Continent), size = 3, alpha = 0.9) +
    scale_color_manual(values = c("Europe" = "#E69F00", "East_Asia" = "#56B4E9", "North_America" = "#009E73")) +
    coord_cartesian(xlim = c(-4, 4)) +
    theme_bw() +
    theme(
      axis.text.y = element_text(face = "italic", size = 8,colour = 'black'), # 菌名斜体
      axis.title.y = element_blank(),
      panel.grid.major.y = element_line(color = "grey95"), # 辅助横线
      legend.position = "bottom",
      legend.text = element_text(size = 12,colour = 'black'),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(x = "Meta-analysis Coefficient (Log scale)", title = title_text)
}

df_filtered <- df_all %>%
  filter(qval.fdr < 0.1) %>%
  group_by(feature, Age_Group) %>%
  filter(n_distinct(Continent) >= 2 | max(abs(coef)) > 2) %>%
  ungroup() 

p_eo_dot <- plot_dot_design(filter(df_filtered, Age_Group == "EO"), "Early-Onset")
p_lo_dot <- plot_dot_design(filter(df_filtered, Age_Group == "LO"), "Late-Onset")

final_plot_B <- p_eo_dot + p_lo_dot + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(legend.position = "bottom")

final_plot_B
ggsave(final_plot_B,width = 13,height = 8,
       filename =paste0(base_symbol,"downstream_analysis/output/picture/4.1.Cleveland_Dot_Plot.pdf"))
